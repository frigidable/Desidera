"""SC1: native-FITS WCS scene editor, explicitly unadopted.

This workbench uses the already-tested immutable ICRS SkyMaskDraft and
read-only FITS Explorer raster/coordinate conventions. It does not
estimate a final isophote, sky, PSF wings, galaxy flux or a science mask.

Original cutouts load on a worker; Qt graphics and edits stay on the
GUI thread. The preview is bounded and never downsampled for mask QA.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import math
from pathlib import Path
from uuid import uuid4

import numpy as np

from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtGui import QColor, QPainterPath, QPen, QShortcut, QKeySequence
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDoubleSpinBox, QGraphicsPathItem,
    QGraphicsPixmapItem, QGraphicsScene, QGraphicsView, QHBoxLayout, QLabel,
    QLineEdit, QListWidget, QPushButton, QSplitter, QVBoxLayout, QWidget,
    QGroupBox, QScrollArea, QSlider, QToolButton, QFileDialog,
)

from lsbpipe.gui.fits_explorer import FITSGraphicsView, FITSExplorer
from lsbpipe.gui.fits_raster import fits_pixel_from_scene, rgb_for_layer
from lsbpipe.photometry.mask_drafts import (
    read_active_mask_draft, save_mask_draft,
)
from lsbpipe.photometry.sky_masks import (
    SkyMaskDraft, SkyMaskEdit, SkyMaskShape, rasterize_sky_shape,
)
from lsbpipe.photometry.ellipse_seed import (
    propose_local_isophote, seed_as_sky_geometry,
    assess_local_isophote_sensitivity, sga_informed_radial_seed,
)
from lsbpipe.photometry.scene_contaminants import prepare_scene_contamination
from lsbpipe.photometry.sga_priors import (
    read_local_sga_priors, fetch_viewer_sga2025_priors,
)
from lsbpipe.survey import legacy_validity
from lsbpipe.surveys.legacy import SkyField, BrickHit
from lsbpipe.surveys.object_projects import read_object_project
from lsbpipe.gui.fits_raster import prepare_fits_explorer


PURPOSE_COLORS = {
    "protect_target": "#a7ee87",
    "exclude": "#f17d77",
    "restore_auto": "#6dc9e5",
    "sky_train": "#e9c96e",
    "sky_holdout": "#b99ded",
}
PURPOSE_LABELS = (
    ("Protect galaxy light", "protect_target"),
    ("Exclude contaminant", "exclude"),
    ("Restore editable auto mask", "restore_auto"),
    ("Sky training region (draft)", "sky_train"),
    ("Sky holdout region (draft)", "sky_holdout"),
)
MAX_RENDER_PIXELS = 6_000_000


def load_scene_original_cutout(project_root: str, band: str,
                               width_arcmin: float) -> dict:
    """Read one original selected-brick quartet; no network or writes.

    Uses the same strict no-padding, source-bound and WCS-aligned reader
    as FITS Explorer. Width is explicitly limited for a full-pixel
    inspection; this is NOT coverage of a larger galaxy outside view.
    """
    project = read_object_project(project_root)
    if project.field is None or project.selected_dataset is None:
        raise ValueError("Project must contain target coordinates and native source")
    if project.data_root is None:
        raise ValueError("Original FITS data root is not configured")
    release, region, brick, _selected_band = project.selected_dataset
    if band not in ("g", "r", "i", "z") or (region == "north" and band == "i"):
        raise ValueError("Selected optical band is unavailable")
    if not 0.5 <= width_arcmin <= min(9., project.field.width_arcmin):
        raise ValueError("Cutout width must remain inside the saved project field")
    field = SkyField(project.field.ra_deg, project.field.dec_deg, width_arcmin)
    hit = BrickHit(release, region, brick, 0., 0., 0., 0., {band: 1})
    original = prepare_fits_explorer(
        directory=project.data_root, field=field, hit=hit, band=band)
    data = original["data"]
    shape = original["pixel_shape"]
    if np.prod(shape) > MAX_RENDER_PIXELS:
        raise ValueError("Cutout too large for native-pixel editor; use a smaller view")
    validity, _ = legacy_validity(
        data["image"], data["invvar"], data["maskbits"], band)
    hard = (~validity | ~np.isfinite(data["nexp"]) |
            (data["nexp"] <= 0))
    rgb = rgb_for_layer(data["image"], "image", band, contrast="faint")
    rgba = np.zeros((*shape, 4), dtype=np.uint8)
    rgba[hard] = (238, 89, 104, 160)
    rgba.setflags(write=False)
    return {
        "status": "NATIVE_SCENE_CUTOUT_NOT_SCIENCE_ADOPTED",
        "data": data, "wcs": original["wcs"],
        "pixel_shape": shape, "rgb": rgb, "hard_rgba": rgba,
        "native_hard": hard, "band": band,
        "source": original["source"],
        "display_width_arcmin": width_arcmin,
        "project_field_width_arcmin": project.field.width_arcmin,
        "project_id": project.project_id,
        "hdu_indices": original["hdu_indices"],
    }


def scene_shape_from_world_points(kind: str, points: tuple, *,
                                  axis_ratio: float = 1.) -> SkyMaskShape:
    """Translate user-selected ICRS points to the existing sky geometry.

    Two points define a circle or an ellipse centre/major-axis endpoint.
    Polygon vertices remain sky coordinates; angle is north through east.
    """
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    if kind == "polygon":
        return SkyMaskShape("polygon", vertices=tuple(points))
    if kind not in ("circle", "ellipse") or len(points) != 2:
        raise ValueError("Circle and ellipse need two ICRS positions")
    center = SkyCoord(*points[0], unit="deg", frame="icrs")
    endpoint = SkyCoord(*points[1], unit="deg", frame="icrs")
    radius = float(center.separation(endpoint).to_value(u.arcsec))
    if not math.isfinite(radius) or not 0 < radius <= 1800.:
        raise ValueError("Selected radius outside safe 0–1800 arcsec range")
    if kind == "circle":
        return SkyMaskShape("circle", center=tuple(points[0]),
                           radius_arcsec=radius)
    east, north = center.spherical_offsets_to(endpoint)
    e, n = east.to_value(u.arcsec), north.to_value(u.arcsec)
    return SkyMaskShape(
        "ellipse", center=tuple(points[0]), a_arcsec=radius,
        q=float(axis_ratio), pa_deg=float(math.degrees(math.atan2(e, n)) % 180.))


def scene_shape_outline(shape: SkyMaskShape, wcs, ny: int):
    """Native WCS -> Qt image path, including rotation and RA wrap."""
    if shape.kind == "polygon":
        ra = np.array([p[0] for p in shape.vertices])
        dec = np.array([p[1] for p in shape.vertices])
    else:
        angle = np.linspace(0, 2*np.pi, 97)[:-1]
        if shape.kind == "circle":
            east = shape.radius_arcsec*np.sin(angle)
            north = shape.radius_arcsec*np.cos(angle)
        else:
            pa = math.radians(shape.pa_deg)
            a = shape.a_arcsec
            b = a*shape.q
            east = a*np.sin(pa)*np.cos(angle)+b*np.cos(pa)*np.sin(angle)
            north = a*np.cos(pa)*np.cos(angle)-b*np.sin(pa)*np.sin(angle)
        # Great-circle destination from east/north angular offsets. Using
        # a direct spherical formula avoids WCS NaNs at the RA=0 meridian
        # seen in the Windows/Astropy SC1 GUI regression. No pixel
        # resampling, and PA remains measured from north through east.
        bearing = np.arctan2(east, north)
        distance = np.deg2rad(np.hypot(east, north) / 3600.)
        ra0 = math.radians(shape.center[0])
        dec0 = math.radians(shape.center[1])
        dec_rad = np.arcsin(np.clip(
            math.sin(dec0)*np.cos(distance)
            + math.cos(dec0)*np.sin(distance)*np.cos(bearing),
            -1., 1.))
        ra_rad = ra0 + np.arctan2(
            np.sin(bearing)*np.sin(distance)*math.cos(dec0),
            np.cos(distance)-math.sin(dec0)*np.sin(dec_rad))
        ra = np.mod(np.rad2deg(ra_rad), 360.)
        dec = np.rad2deg(dec_rad)
    xx, yy = wcs.all_world2pix(ra, dec, 0)
    if not np.isfinite(xx).all() or not np.isfinite(yy).all():
        raise ValueError("Native WCS cannot display the selected sky shape")
    path = QPainterPath()
    for i, (x, y) in enumerate(zip(xx, yy)):
        px, py = float(x)+.5, ny-float(y)-.5
        if i == 0:
            path.moveTo(px, py)
        else:
            path.lineTo(px, py)
    path.closeSubpath()
    return path


def prepare_scene_neighbour_masks(dataset, project_root, draft, center_icrs,
                                  sga_prior=None):
    """Prepare native Tractor masks, optionally protecting an explicit SGA core."""
    project = read_object_project(project_root)
    return prepare_scene_contamination(
        dataset, project, draft, center_icrs, sga_prior=sga_prior)


def propose_scene_isophote(dataset, center_icrs, threshold, smoothing,
                           contamination, bridge_guard_pix=0,
                           sga_prior=None):
    """Require the prepared source-bound mask before extracting an isophote."""
    if (contamination is None or
            contamination.get("status") !=
            "CATALOGUE_CONTAMINATION_MASK_NOT_ADOPTED" or
            contamination["project_id"] != dataset["project_id"] or
            contamination["band"] != dataset["band"] or
            contamination["dataset_source"] != dataset["source"] or
            tuple(contamination["anchor_icrs"]) != tuple(center_icrs)):
        raise ValueError("Current-band Tractor neighbour mask required")
    wcs = dataset["wcs"]
    x, y = wcs.all_world2pix(
        float(center_icrs[0]), float(center_icrs[1]), 0)
    if not (math.isfinite(float(x)) and math.isfinite(float(y))):
        raise ValueError("Selected ICRS galaxy center cannot be projected")
    if np.shape(contamination["mask"]) != tuple(dataset["pixel_shape"]):
        raise ValueError("Neighbour mask not on native pixel grid")
    valid = ~(dataset["native_hard"] | contamination["mask"])
    model = propose_local_isophote(
        dataset["data"]["image"], valid,
        (float(x), float(y)),
        threshold_contrast=float(threshold), smooth_sigma_pix=float(smoothing),
        bridge_guard_pix=int(bridge_guard_pix))
    sensitivity = assess_local_isophote_sensitivity(
        dataset["data"]["image"], valid, (float(x), float(y)), model,
        threshold_contrast=float(threshold),
        smooth_sigma_pix=float(smoothing),
        bridge_guard_pix=int(bridge_guard_pix))
    def sky_shapes(candidate):
        geometry = seed_as_sky_geometry(candidate.seed, wcs)
        ellipse = SkyMaskShape(
            "ellipse",
            center=(geometry["ra_deg"], geometry["dec_deg"]),
            a_arcsec=geometry["sma_arcsec"],
            q=geometry["q_sky"],
            pa_deg=geometry["pa_north_through_east_deg"])
        coords = np.asarray(candidate.contour_xy)
        ra, dec = wcs.all_pix2world(coords[:, 0], coords[:, 1], 0)
        if not (np.isfinite(ra).all() and np.isfinite(dec).all()):
            raise ValueError("WCS cannot project trial isophotal contour")
        contour = SkyMaskShape(
            "polygon", vertices=tuple((float(a)%360., float(b))
                                      for a,b in zip(ra,dec)))
        return {"ellipse":ellipse, "contour":contour,
                "model":candidate}

    alternatives = tuple(
        (sky_shapes(candidate) if candidate is not None else None)
        for candidate in sensitivity["trial_models"])
    ellipse, contour = alternatives[0]["ellipse"], alternatives[0]["contour"]
    sga_result = None
    sga_seed_shape = None
    sga_warning = None
    sga_comparison = None
    if sga_prior is not None:
        if (sga_prior.core_shape() is None or
                contamination["project_id"] != dataset["project_id"]):
            raise ValueError("Selected SGA prior has no valid geometry")
        from astropy.coordinates import SkyCoord
        from astropy import units as u
        prior_shape = sga_prior.core_shape()
        from_centre = SkyCoord(*prior_shape.center, unit="deg")
        anchor_centre = SkyCoord(*center_icrs, unit="deg")
        center_sep = float(from_centre.separation(anchor_centre).arcsec)
        if center_sep > 60.:
            raise ValueError("SGA prior centre too far from target anchor")
        # SGA/Tractor are not statistically independent. The radial
        # diagnostic itself uses the SAME original pixels after masking.
        try:
            sga_result = sga_informed_radial_seed(
                dataset["data"]["image"], valid, wcs, prior_shape,
                local_sky=model.local_sky,
                smoothed_rms=model.smoothed_sky_rms,
                smooth_sigma_pix=float(smoothing))
            if sga_result.detected_a_arcsec is not None:
                sga_seed_shape = SkyMaskShape(
                    "ellipse", center=sga_result.centre_icrs,
                    a_arcsec=sga_result.detected_a_arcsec,
                    q=sga_result.q_sky,
                    pa_deg=sga_result.pa_deg)
                # Explicitly INVENTORY masked source centers, not
                # automatic 'shred' assignments or a new mask.
                footprint = rasterize_sky_shape(
                    sga_seed_shape,wcs,dataset["pixel_shape"])
                candidates = []
                for source in contamination["stats"].get(
                        "masked_catalog_centres", ()):
                    ix, iy = int(round(source["x"])),int(round(source["y"]))
                    if (0 <= ix < footprint.shape[1] and
                            0 <= iy < footprint.shape[0] and
                            footprint[iy,ix] > 0):
                        candidates.append(source)
                sga_comparison = {
                    "masked_catalog_centres_in_initial_ellipse": len(candidates),
                    "review_candidates": tuple(candidates[:25]),
                    "catalogue_inventory_truncated":
                        contamination["stats"].get(
                            "masked_centres_truncated", 0),
                }
        except (ValueError, TypeError, ArithmeticError) as exc:
            # Independent connected-component results stay usable when
            # SGA-informed analysis cannot support an outer extent.
            sga_warning = str(exc)[:250]
        model_sky = seed_as_sky_geometry(model.seed,wcs)
        pa_delta = abs(model_sky["pa_north_through_east_deg"]
                       -prior_shape.pa_deg)%180.
        sga_comparison = {
            **(sga_comparison or {}),
            "reference_a_arcsec":prior_shape.a_arcsec,
            "reference_q":prior_shape.q,
            "reference_pa_deg":prior_shape.pa_deg,
            "independent_q":model_sky["q_sky"],
            "independent_pa_deg":
                model_sky["pa_north_through_east_deg"],
            "delta_q":abs(model_sky["q_sky"]-prior_shape.q),
            "delta_pa_deg":min(pa_delta,180.-pa_delta),
            "center_offset_arcsec":center_sep,
            "sga_ref":sga_prior.sga_id,
            "sga_source_sha256":sga_prior.source_sha256,
        }
    return {"status":model.state, "ellipse":ellipse, "contour":contour,
            "proposal":model, "alternatives":alternatives,
            "sga_radial":sga_result, "sga_seed_shape":sga_seed_shape,
            "sga_warning":sga_warning, "sga_comparison":sga_comparison,
            "band":dataset["band"],
            "project_id":dataset["project_id"],
            "source":dataset["source"],
            "catalog_sha256":contamination["catalog_sha256"],
            "draft_digest":contamination["draft_digest"],
            "contamination":contamination["stats"],
            "sensitivity":sensitivity,
            "bridge_guard_pix":int(bridge_guard_pix)}


class SceneDrawingView(FITSGraphicsView):
    pointChosen = Signal(float, float)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.drawing = False
        self.setObjectName("scene_native_fits_view")

    def set_drawing(self, drawing: bool):
        self.drawing = bool(drawing)
        self.setDragMode(
            QGraphicsView.NoDrag if self.drawing
            else QGraphicsView.ScrollHandDrag)

    def mousePressEvent(self, event):
        if self.drawing and event.button() == Qt.LeftButton:
            point = self.mapToScene(event.position().toPoint())
            if self.scene() is not None and self.scene().sceneRect().contains(point):
                self.pointChosen.emit(point.x(), point.y())
            event.accept()
            return
        super().mousePressEvent(event)


class SceneEditor(QWidget):
    """Real editable scene DRAFT in the same Photometry Workspace.

    Only explicit save creates a new immutable project draft snapshot.
    Sky coordinate geometry is shared when requested, but the original
    native optical pixel grids and hard masks are independent.
    """
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("scene_editor_sc1")
        self._project_root = None
        self._project_id = None
        self._source_key = None
        self._draft = None
        self._saved_digest = None
        self._undo = []
        self._redo = []
        self._pending_world = []
        self._pending_scene = []
        self._preview_path = None
        self._replace_index = None
        self._loaded = None
        self._isophote_proposal = None
        self._isophote_preview = None
        self._contamination = None
        self._contamination_overlay = None
        self._conflict_paths = []
        self._sga_shred_paths = []
        self._isophote_center_icrs = None
        self._picking_isophote_center = False
        self._project_target = None
        self._sga_lookup = None
        self._sga_candidate = None
        self._sga_item = None
        self._future = None
        self._work_kind = None
        self._load_generation = 0
        self._executor = ThreadPoolExecutor(
            max_workers=1, thread_name_prefix="desidera-scene-fits")
        self._poll_timer = QTimer(self)
        self._poll_timer.setInterval(45)
        self._poll_timer.timeout.connect(self._poll_load)
        # Human-drawn drafts are automatically checkpointed, never
        # adopted or promoted to an approved segmentation.
        self._autosave_timer = QTimer(self)
        self._autosave_timer.setSingleShot(True)
        self._autosave_timer.setInterval(650)
        self._autosave_timer.timeout.connect(self.flush_autosave)
        self._paths = []
        self._available = True
        self._band = "g"
        self._image = None
        self._native_overlay = None
        self._work_cache = {}

        root = QVBoxLayout(self)
        root.setContentsMargins(9, 8, 9, 8)
        root.setSpacing(7)
        self.heading = QLabel("Scene and Smart Masks · SC1")
        self.heading.setStyleSheet("font-size: 17px; font-weight: 600;")
        root.addWidget(self.heading)
        self.warning = QLabel(
            "Original native FITS; editable ICRS region drafts only. "
            "Isophotal galaxy extent, stellar wings, residual sky and "
            "photometric measurements are NOT validated.")
        self.warning.setWordWrap(True)
        root.addWidget(self.warning)

        toolbar = QHBoxLayout()
        self.band = QComboBox()
        self.band.setObjectName("scene_band")
        for b in "griz":
            self.band.addItem(b.upper(), b)
        self.band.currentIndexChanged.connect(self._band_changed)
        toolbar.addWidget(QLabel("Band"))
        toolbar.addWidget(self.band)
        self.width = QDoubleSpinBox()
        self.width.setObjectName("scene_cutout_width")
        self.width.setRange(.5, 9.)
        self.width.setDecimals(1)
        self.width.setSingleStep(.5)
        self.width.setValue(6.)
        self.width.setSuffix(" arcmin")
        self.width.setToolTip(
            "Only an original native-pixel cutout is rendered. "
            "A larger target may extend beyond this explicit view.")
        toolbar.addWidget(QLabel("View width"))
        toolbar.addWidget(self.width)
        self.load_btn = QPushButton("Load original FITS")
        self.load_btn.setObjectName("scene_load_original")
        self.load_btn.clicked.connect(self.load_original)
        toolbar.addWidget(self.load_btn)
        self.display_mode = QComboBox()
        self.display_mode.setObjectName("scene_image_stretch")
        self.display_mode.addItem("Asinh · faint", ("faint", "asinh"))
        self.display_mode.addItem("Asinh · standard", ("standard", "asinh"))
        self.display_mode.addItem("Linear · faint", ("faint", "linear"))
        self.display_mode.currentIndexChanged.connect(self._rerender_image)
        self.display_mode.setToolTip(
            "Display-only stretch on original IMAGE values; no science data changed.")
        self.fit_btn = QPushButton("Fit")
        self.fit_btn.clicked.connect(lambda: self.view.fit_image())
        toolbar.addWidget(self.fit_btn)
        toolbar.addStretch(1)
        root.addLayout(toolbar)
        # Keep the narrow-screen minimum width independent of display options.
        display_row = QHBoxLayout()
        display_row.addWidget(QLabel("Display stretch"))
        display_row.addWidget(self.display_mode)
        display_row.addStretch(1)
        root.addLayout(display_row)

        splitter = QSplitter(Qt.Horizontal, self)
        self.graphics = QGraphicsScene(self)
        self.view = SceneDrawingView()
        self.view.setScene(self.graphics)
        self.view.pointChosen.connect(self._on_click)
        self.view.pointHovered.connect(self._inspect)
        splitter.addWidget(self.view)

        controls = QWidget()
        control_layout = QVBoxLayout(controls)
        control_layout.setContentsMargins(7, 4, 7, 4)
        self.mode = QComboBox()
        self.mode.setObjectName("scene_draw_mode")
        self.mode.addItem("Pan / inspect", "pan")
        self.mode.addItem("Circle · two clicks", "circle")
        self.mode.addItem("Ellipse · two clicks", "ellipse")
        self.mode.addItem("Polygon · vertices", "polygon")
        self.mode.currentIndexChanged.connect(self._mode_changed)
        control_layout.addWidget(QLabel("Tool"))
        control_layout.addWidget(self.mode)
        self.purpose = QComboBox()
        self.purpose.setObjectName("scene_mask_purpose")
        for label, purpose in PURPOSE_LABELS:
            self.purpose.addItem(label, purpose)
        control_layout.addWidget(QLabel("Region meaning"))
        control_layout.addWidget(self.purpose)
        self.scope = QComboBox()
        self.scope.addItem("All optical bands (ICRS geometry)", "all")
        self.scope.addItem("This band only", "band")
        self.scope.setObjectName("scene_mask_scope")
        control_layout.addWidget(self.scope)
        self.axis_ratio = QDoubleSpinBox()
        self.axis_ratio.setObjectName("scene_ellipse_q")
        self.axis_ratio.setRange(.05, 1.)
        self.axis_ratio.setDecimals(2)
        self.axis_ratio.setSingleStep(.05)
        self.axis_ratio.setValue(.7)
        control_layout.addWidget(QLabel("Ellipse minor/major axis ratio"))
        control_layout.addWidget(self.axis_ratio)
        self.comment = QLineEdit()
        self.comment.setObjectName("scene_edit_comment")
        self.comment.setPlaceholderText("Optional reason / source identity")
        control_layout.addWidget(self.comment)

        # A SGA ellipse is an externally fitted central-scale PRIOR.
        # It never replaces a native FITS result or an outer LSB contour.
        sga_box = QGroupBox("SGA · initial geometry prior (optional)")
        sga_box.setObjectName("scene_sga_prior_group")
        sga_layout = QVBoxLayout(sga_box)
        self.sga_find = QPushButton("Find SGA-2025 · DESI Viewer")
        self.sga_find.setObjectName("scene_sga_find_desi")
        self.sga_find.setToolTip(
            "Query the official DESI Legacy Surveys SGA-2025 overlay "
            "in a 3-arcmin sky region: one small JSON request, no 1.4-GB "
            "catalogue download. Viewer ellipses are provisional geometry "
            "and may differ from the merged SGA2025 scientific catalogue.")
        self.sga_find.clicked.connect(self._find_desi_sga2025)
        sga_layout.addWidget(self.sga_find)
        self.sga_parent_find = QPushButton(
            "Try SGA-2025 Parent · initial geometry")
        self.sga_parent_find.setObjectName("scene_sga_find_parent")
        self.sga_parent_find.setToolTip(
            "Separate DESI Viewer SGA2025 parent overlay. May provide "
            "an initial catalogue diameter even when fitted SGA2025 "
            "ellipse is missing. NOT a measured isophote or fitted SGAID.")
        self.sga_parent_find.clicked.connect(
            self._find_desi_sga2025_parent)
        sga_layout.addWidget(self.sga_parent_find)
        self.sga_load = QPushButton("Or choose local SGA FITS…")
        self.sga_load.setObjectName("scene_sga_open_local")
        self.sga_load.setToolTip(
            "Explicit local SGA-2020 or SGA-2025 catalogue. "
            "Read-only, no network or file copy; large file hashed in worker.")
        self.sga_load.clicked.connect(self._open_sga_file)
        sga_layout.addWidget(self.sga_load)
        self.sga_candidates = QComboBox()
        self.sga_candidates.setObjectName("scene_sga_candidates")
        self.sga_candidates.setToolTip(
            "Choose the right galaxy among positional matches; matching "
            "catalogues and Tractor sources are NOT independent evidence.")
        self.sga_candidates.currentIndexChanged.connect(
            self._sga_candidate_changed)
        sga_layout.addWidget(self.sga_candidates)
        self.sga_use_prior = QCheckBox(
            "Protect selected SGA core against Tractor shredding")
        self.sga_use_prior.setObjectName("scene_sga_use_protection")
        self.sga_use_prior.setChecked(False)
        self.sga_use_prior.toggled.connect(self._sga_protection_changed)
        self.sga_use_prior.setToolTip(
            "Opt-in: protects the selected SGA central ellipse during "
            "preliminary source screening only. DESI Viewer geometry is "
            "an overlay prior, not necessarily catalogue D26 and never "
            "the full LSB radius.")
        sga_layout.addWidget(self.sga_use_prior)
        self.sga_use_seed = QCheckBox(
            "Use SGA centre/q/PA for initial geometry")
        self.sga_use_seed.setObjectName("scene_sga_use_geometry_seed")
        self.sga_use_seed.setToolTip(
            "Use selected SGA central geometry as an orientation prior. "
            "Estimate radial extent from masked native DESI pixels, "
            "with no arbitrary diameter multiplier or science adoption.")
        self.sga_use_seed.setChecked(False)
        self.sga_use_seed.toggled.connect(self._clear_isophote_preview)
        sga_layout.addWidget(self.sga_use_seed)
        self.sga_status = QLabel(
            "Use Find SGA-2025 to query the DESI Viewer by current target "
            "coordinates. It returns a provisional central overlay ellipse, "
            "not a calibrated outer LSB boundary. Source selection is manual.")
        self.sga_status.setObjectName("scene_sga_status")
        self.sga_status.setWordWrap(True)
        sga_layout.addWidget(self.sga_status)
        control_layout.addWidget(sga_box)

        # A threshold contour and a moment-ellipse are offered separately.
        # The researcher must inspect and explicitly accept either one.
        isogroup = QGroupBox("Suggest galaxy isophote · preliminary")
        isogroup.setObjectName("scene_isophote_group")
        iso_layout = QVBoxLayout(isogroup)
        settings_row = QHBoxLayout()
        settings_row.addWidget(QLabel("Contrast"))
        self.iso_threshold = QDoubleSpinBox()
        self.iso_threshold.setObjectName("scene_isophote_threshold")
        self.iso_threshold.setRange(1., 8.)
        self.iso_threshold.setSingleStep(.3)
        self.iso_threshold.setDecimals(1)
        self.iso_threshold.setValue(2.2)
        self.iso_threshold.setToolTip(
            "Local smoothed-background contrast, NOT calibrated significance.")
        settings_row.addWidget(self.iso_threshold)
        settings_row.addWidget(QLabel("Smooth px"))
        self.iso_smooth = QDoubleSpinBox()
        self.iso_smooth.setObjectName("scene_isophote_smoothing")
        self.iso_smooth.setRange(1.,12.)
        self.iso_smooth.setSingleStep(1.)
        self.iso_smooth.setDecimals(1)
        self.iso_smooth.setValue(4.)
        settings_row.addWidget(self.iso_smooth)
        iso_layout.addLayout(settings_row)
        self.iso_bridge_guard = QCheckBox(
            "Conservative contour · suppress thin bridges (2 px)")
        self.iso_bridge_guard.setObjectName("scene_isophote_bridge_guard")
        self.iso_bridge_guard.setChecked(True)
        self.iso_bridge_guard.setToolTip(
            "Binary morphological opening; can remove real clumpy or "
            "tidal features. Compare with guard OFF before interpretation.")
        iso_layout.addWidget(self.iso_bridge_guard)
        center_row = QHBoxLayout()
        self.iso_pick_center = QPushButton("Pick galaxy centre")
        self.iso_pick_center.setObjectName("scene_isophote_pick_center")
        self.iso_pick_center.clicked.connect(self._pick_isophote_center)
        center_row.addWidget(self.iso_pick_center)
        self.iso_reset_center = QPushButton("Use target ICRS")
        self.iso_reset_center.clicked.connect(self._reset_isophote_center)
        center_row.addWidget(self.iso_reset_center)
        iso_layout.addLayout(center_row)
        self.iso_prepare = QPushButton("1 · Prepare neighbour masks")
        self.iso_prepare.setObjectName("scene_prepare_neighbour_masks")
        self.iso_prepare.setToolTip(
            "Read same-brick local DR11 South Tractor catalogue; "
            "make provisional source envelopes and inspect MASKBITS. "
            "No model subtraction, PSF wing validation or downloads.")
        self.iso_prepare.clicked.connect(self._prepare_neighbour_masks)
        iso_layout.addWidget(self.iso_prepare)
        self.iso_mask_toggle = QCheckBox("Show proposed contaminant mask")
        self.iso_mask_toggle.setObjectName("scene_show_neighbor_mask")
        self.iso_mask_toggle.setChecked(True)
        self.iso_mask_toggle.toggled.connect(self._update_contamination_overlay)
        iso_layout.addWidget(self.iso_mask_toggle)
        opacity_row = QHBoxLayout()
        opacity_row.addWidget(QLabel("Mask opacity"))
        self.iso_opacity = QSlider(Qt.Horizontal)
        self.iso_opacity.setObjectName("scene_neighbour_mask_opacity")
        self.iso_opacity.setRange(0, 100)
        self.iso_opacity.setValue(45)
        self.iso_opacity.valueChanged.connect(self._update_contamination_overlay)
        opacity_row.addWidget(self.iso_opacity)
        iso_layout.addLayout(opacity_row)
        self.iso_show_conflicts = QCheckBox(
            "Show 20 nearest unresolved neighbours")
        self.iso_show_conflicts.setChecked(False)
        self.iso_show_conflicts.setToolTip(
            "Hide clutter by default. Enable to inspect the 20 closest "
            "catalogued sources left unmasked near the galaxy; counts "
            "still include all unresolved sources.")
        self.iso_show_conflicts.toggled.connect(
            self._update_contamination_overlay)
        iso_layout.addWidget(self.iso_show_conflicts)
        self.iso_show_sga_shreds = QCheckBox(
            "Show masked candidates inside SGA initial ellipse")
        self.iso_show_sga_shreds.setObjectName(
            "scene_sga_initial_shred_candidates")
        self.iso_show_sga_shreds.setChecked(False)
        self.iso_show_sga_shreds.setToolTip(
            "Blue review rings mark previously masked Tractor centres "
            "inside the SGA-supported initial ellipse. These are "
            "NOT proven galaxy shreds or automatically restored.")
        self.iso_show_sga_shreds.toggled.connect(
            self._update_contamination_overlay)
        iso_layout.addWidget(self.iso_show_sga_shreds)
        self.iso_mask_status = QLabel(
            "Neighbour mask not prepared. Local DR11 South Tractor "
            "catalogue required; obtain it through PSF-0 if absent.")
        self.iso_mask_status.setObjectName("scene_neighbor_mask_status")
        self.iso_mask_status.setWordWrap(True)
        iso_layout.addWidget(self.iso_mask_status)
        self.iso_mask_details_toggle = QToolButton()
        self.iso_mask_details_toggle.setText("Mask provenance and cautions")
        self.iso_mask_details_toggle.setCheckable(True)
        self.iso_mask_details_toggle.setChecked(False)
        self.iso_mask_details_toggle.setObjectName("scene_mask_details_toggle")
        iso_layout.addWidget(self.iso_mask_details_toggle)
        self.iso_mask_details = QLabel("")
        self.iso_mask_details.setWordWrap(True)
        self.iso_mask_details.setVisible(False)
        self.iso_mask_details_toggle.toggled.connect(
            self.iso_mask_details.setVisible)
        iso_layout.addWidget(self.iso_mask_details)
        self.iso_generate = QPushButton("2 · Find masked isophote")
        self.iso_generate.setObjectName("scene_isophote_generate")
        self.iso_generate.clicked.connect(self._generate_isophote)
        iso_layout.addWidget(self.iso_generate)
        iso_layout.addWidget(QLabel("Compare threshold / smoothing"))
        self.iso_candidate = QComboBox()
        self.iso_candidate.setObjectName("scene_isophote_trial_selector")
        self.iso_candidate.setToolTip(
            "Switch between four actual ICRS contours from nearby detection "
            "thresholds/scales. A stable area alone does not guarantee "
            "spatial agreement. No trial is science-adopted.")
        self.iso_candidate.currentIndexChanged.connect(
            self._show_isophote_preview)
        iso_layout.addWidget(self.iso_candidate)
        self.iso_shape = QComboBox()
        self.iso_shape.addItem("Native isophote contour", "contour")
        self.iso_shape.addItem("Equivalent ICRS ellipse", "ellipse")
        self.iso_shape.currentIndexChanged.connect(self._show_isophote_preview)
        iso_layout.addWidget(self.iso_shape)
        self.iso_share_bands = QCheckBox(
            "Share across g/r/i/z (explicit)")
        self.iso_share_bands.setObjectName("scene_isophote_share_bands")
        self.iso_share_bands.setChecked(False)
        self.iso_share_bands.setToolTip(
            "Default is native-band only. A common ICRS boundary needs "
            "deliberate confirmation and independent per-band QA.")
        iso_layout.addWidget(self.iso_share_bands)
        self.iso_low_conf_override = QCheckBox(
            "Allow low-confidence draft after review")
        self.iso_low_conf_override.setObjectName(
            "scene_isophote_review_override")
        self.iso_low_conf_override.setChecked(False)
        self.iso_low_conf_override.toggled.connect(self._refresh_controls)
        iso_layout.addWidget(self.iso_low_conf_override)
        proposal_actions = QHBoxLayout()
        self.iso_accept = QPushButton("Add proposed region")
        self.iso_accept.setObjectName("scene_isophote_accept")
        self.iso_accept.clicked.connect(self._accept_isophote)
        proposal_actions.addWidget(self.iso_accept)
        self.iso_discard = QPushButton("Discard")
        self.iso_discard.clicked.connect(self._clear_isophote_preview)
        proposal_actions.addWidget(self.iso_discard)
        iso_layout.addLayout(proposal_actions)
        self.iso_summary = QLabel(
            "No contour yet. Run source screening, then a multi-scale proposal.")
        self.iso_summary.setObjectName("scene_isophote_summary")
        self.iso_summary.setWordWrap(True)
        iso_layout.addWidget(self.iso_summary)
        self.iso_sga_summary = QLabel(
            "SGA-informed geometry not yet calculated.")
        self.iso_sga_summary.setObjectName(
            "scene_sga_informed_initial_geometry_qa")
        self.iso_sga_summary.setWordWrap(True)
        iso_layout.addWidget(self.iso_sga_summary)
        self.iso_details_toggle = QToolButton()
        self.iso_details_toggle.setText("Isophote QA and all sensitivity trials")
        self.iso_details_toggle.setCheckable(True)
        self.iso_details_toggle.setChecked(False)
        self.iso_details_toggle.setObjectName("scene_isophote_details_toggle")
        iso_layout.addWidget(self.iso_details_toggle)
        self.iso_details = QLabel(
            "Uses local provisional sky and native FITS. No PSF/sky/stellar "
            "correction or science adoption. Review before accepting.")
        self.iso_details.setWordWrap(True)
        self.iso_details.setObjectName("scene_isophote_evidence")
        self.iso_details.setVisible(False)
        self.iso_details_toggle.toggled.connect(self.iso_details.setVisible)
        iso_layout.addWidget(self.iso_details)
        self.iso_threshold.valueChanged.connect(self._clear_isophote_preview)
        self.iso_smooth.valueChanged.connect(self._clear_isophote_preview)
        self.iso_bridge_guard.toggled.connect(self._clear_isophote_preview)
        control_layout.addWidget(isogroup)

        buttons = QHBoxLayout()
        self.finish_btn = QPushButton("Finish polygon")
        self.finish_btn.setObjectName("scene_finish_polygon")
        self.finish_btn.clicked.connect(self.finish_polygon)
        buttons.addWidget(self.finish_btn)
        self.cancel_btn = QPushButton("Cancel shape")
        self.cancel_btn.clicked.connect(self._cancel_pending)
        buttons.addWidget(self.cancel_btn)
        control_layout.addLayout(buttons)
        control_layout.addWidget(QLabel("Unadopted edits (click one to remove)"))
        self.edits = QListWidget()
        self.edits.setObjectName("scene_edits_list")
        self.edits.setMinimumHeight(95)
        self.edits.currentRowChanged.connect(self._refresh_controls)
        control_layout.addWidget(self.edits, 1)
        history = QHBoxLayout()
        self.undo_btn = QPushButton("Undo")
        self.undo_btn.setObjectName("scene_undo")
        self.undo_btn.clicked.connect(self.undo)
        history.addWidget(self.undo_btn)
        self.redo_btn = QPushButton("Redo")
        self.redo_btn.setObjectName("scene_redo")
        self.redo_btn.clicked.connect(self.redo)
        history.addWidget(self.redo_btn)
        self.edit_btn = QPushButton("Redraw selected")
        self.edit_btn.setObjectName("scene_edit_selected")
        self.edit_btn.setToolTip(
            "Replace the selected outline by drawing it again on the "
            "original image; keeps source/edit identity and undo history.")
        self.edit_btn.clicked.connect(self.edit_selected)
        history.addWidget(self.edit_btn)
        self.remove_btn = QPushButton("Remove selected")
        self.remove_btn.setObjectName("scene_remove")
        self.remove_btn.clicked.connect(self.remove_selected)
        history.addWidget(self.remove_btn)
        control_layout.addLayout(history)
        self.show_hard = QCheckBox("Show native unusable pixels")
        self.show_hard.setChecked(True)
        self.show_hard.stateChanged.connect(self._update_hard_overlay)
        control_layout.addWidget(self.show_hard)
        self.save_btn = QPushButton("Save scene draft")
        self.save_btn.setObjectName("scene_save_draft")
        self.save_btn.clicked.connect(self.save_draft)
        self.save_btn.setToolTip(
            "Explicitly save one new immutable project mask-draft revision. "
            "This is NOT scientific approval of segmentation or photometry.")
        control_layout.addWidget(self.save_btn)
        self.status = QLabel("Open an object project to begin.")
        self.status.setWordWrap(True)
        self.status.setObjectName("scene_status")
        control_layout.addWidget(self.status)
        controls_scroll = QScrollArea(self)
        controls_scroll.setObjectName("scene_controls_scroll")
        controls_scroll.setWidgetResizable(True)
        controls_scroll.setWidget(controls)
        controls_scroll.setMinimumWidth(220)
        self.view.setMinimumWidth(0)
        splitter.addWidget(controls_scroll)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 0)
        splitter.setSizes([920, 345])
        root.addWidget(splitter, 1)
        self.pixel_info = QLabel("Hover original FITS to inspect native pixel and ICRS.")
        self.pixel_info.setObjectName("scene_pixel_info")
        self.pixel_info.setWordWrap(True)
        root.addWidget(self.pixel_info)
        QShortcut(QKeySequence("Ctrl+Z"), self).activated.connect(self.undo)
        QShortcut(QKeySequence("Ctrl+Shift+Z"), self).activated.connect(self.redo)
        QShortcut(QKeySequence("Escape"), self).activated.connect(
            self._cancel_pending)
        self._mode_changed()
        self._refresh_controls()

    @property
    def draft(self):
        return self._draft

    def bind_project(self, project_root):
        """Metadata and saved draft only: never open FITS or download on entry."""
        if project_root is None:
            self.flush_autosave()
            self._switch_source(None, None)
            self._project_root = None
            self._project_id = None
            self._draft = None
            self._saved_digest = None
            self.status.setText("Open an object project to prepare a scene.")
            self._refresh_controls()
            return
        project = read_object_project(project_root)
        self._project_target = (
            (project.field.ra_deg, project.field.dec_deg)
            if project.field is not None else None)
        key = (str(Path(project.root).resolve()), project.project_id,
               project.selected_dataset, str(project.data_root))
        if key == self._source_key:
            return
        self.flush_autosave()
        self._sga_lookup = None
        self._sga_candidate = None
        self.sga_candidates.blockSignals(True)
        self.sga_candidates.clear()
        self.sga_candidates.blockSignals(False)
        self.sga_use_prior.setChecked(False)
        self.sga_use_seed.setChecked(False)
        self.sga_status.setText(
            "No SGA prior selected for this galaxy; local FITS only. "
            "SGA D26 is a central scale, never an outer LSB boundary.")
        # Keep unsaved revisions in memory if researcher navigates to a
        # different project and later comes back in the same GUI session.
        if self._draft is not None and self._project_root is not None:
            self._work_cache[self._source_key] = (
                self._draft, self._saved_digest,
                list(self._undo), list(self._redo))
        self._switch_source(key, None)
        self._project_root = Path(project.root)
        self._project_id = project.project_id
        if project.field is not None:
            max_width = min(9., max(.5, project.field.width_arcmin))
            self.width.setMaximum(max_width)
            self.width.setValue(min(6., max_width))
        if key in self._work_cache:
            self._draft, self._saved_digest, self._undo, self._redo = (
                self._work_cache[key])
        else:
            try:
                self._draft, self._saved_digest = read_active_mask_draft(
                    project.root)
                self._available = True
            except (OSError, ValueError, KeyError, TypeError) as exc:
                self._draft = None
                self._saved_digest = None
                self._available = False
                self.status.setText(
                    f"Cannot read saved mask history: {exc}. "
                    "Editing blocked; original project not overwritten.")
                self._refresh_controls()
                return
            self._undo, self._redo = [], []
        self._available = True
        self.status.setText(
            f"Scene draft ready for {project.display_name}: "
            f"{len(self._draft.edits)} edits. Load native FITS explicitly. "
            "No science result is adopted.")
        self._refresh_controls()
        self._redraw_outlines()

    def _switch_source(self, key, dataset):
        self._load_generation += 1
        self._source_key = key
        self._loaded = dataset
        self._pending_world = []
        self._pending_scene = []
        self._preview_path = None
        self._replace_index = None
        self._isophote_proposal = None
        self._isophote_preview = None
        self._contamination = None
        self._contamination_overlay = None
        self._conflict_paths = []
        self._sga_shred_paths = []
        self._sga_item = None
        self._isophote_center_icrs = None
        self._picking_isophote_center = False
        self.graphics.clear()
        self._image = None
        self._native_overlay = None
        self._paths.clear()
        # Preserve a running stale worker until it finishes, then drop
        # its result by generation. Otherwise a band/project switch could
        # stop polling forever and leave the Load button blocked.
        if self._future is None or self._future.done():
            self._poll_timer.stop()

    def set_band(self, band: str):
        if band not in "griz" or len(band) != 1:
            raise ValueError("Unsupported scene band")
        index = self.band.findData(band)
        if self.band.currentIndex() != index:
            self.band.setCurrentIndex(index)

    def _band_changed(self, *_):
        new_band = self.band.currentData()
        if new_band is None or new_band == self._band:
            return
        self._band = new_band
        self._switch_source(self._source_key, None)
        self.status.setText(
            f"Switched to {new_band.upper()} band. Load its own native "
            "FITS/WCS to display the same ICRS draft; no mask pixels copied.")
        self._refresh_controls()

    def load_original(self):
        if self._project_root is None or not self._available:
            self.status.setText("Open a valid galaxy project before loading FITS.")
            return
        if self._future is not None and not self._future.done():
            self.status.setText("Original FITS load in progress.")
            return
        self._switch_source(self._source_key, None)
        width = float(self.width.value())
        self.status.setText(
            f"Reading original {self._band.upper()} quartet, "
            f"{width:g} arcmin cutout… No source files will be written.")
        self.load_btn.setEnabled(False)
        self._future = self._executor.submit(
            load_scene_original_cutout, str(self._project_root),
            self._band, width)
        self._work_kind = "load"
        self._active_generation = self._load_generation
        self._poll_timer.start()

    def _poll_load(self):
        future = self._future
        if future is None or not future.done():
            return
        self._poll_timer.stop()
        self._future = None
        kind = self._work_kind
        self._work_kind = None
        if self._active_generation != self._load_generation:
            self._refresh_controls()
            return
        try:
            result = future.result()
            if kind == "isophote":
                self._accept_isophote_result(result)
            elif kind == "sga":
                self._accept_sga_lookup(result)
            elif kind == "contamination":
                self._accept_contamination_result(result)
            else:
                self._accept_dataset(result)
        except (OSError, ValueError, KeyError, TypeError, RuntimeError,
                MemoryError) as exc:
            if kind == "isophote":
                self.iso_details.setText(
                    "No proposal: " + str(exc)
                    + ". Try adjusting contrast/smoothing or pick the "
                      "diffuse component centre. Nothing was saved.")
            elif kind == "sga":
                self.sga_status.setText(
                    "SGA lookup failed: " + str(exc)
                    + ". No SGA region adopted; original FITS untouched.")
            elif kind == "contamination":
                self.iso_mask_status.setText(
                    "Mask preparation FAILED: " + str(exc)
                    + ". No unmasked isophote will be offered. "
                      "For DR11 South obtain the local Tractor catalogue "
                      "in PSF-0; no automatic download occurred.")
            else:
                self.status.setText(
                    f"Original FITS scene could not be displayed: {exc}")
        self._refresh_controls()

    def _accept_dataset(self, data):
        """Main-thread Qt creation from a native worker's immutable arrays."""
        if data["status"] != "NATIVE_SCENE_CUTOUT_NOT_SCIENCE_ADOPTED":
            raise ValueError("Scene workbench refuses science-adopted input")
        if data["band"] != self._band or data["project_id"] != self._project_id:
            raise ValueError("Scene image identity does not match active project")
        ny, nx = map(int, data["pixel_shape"])
        if ny*nx > MAX_RENDER_PIXELS or ny <= 0 or nx <= 0:
            raise ValueError("Scene native pixel count out of bounds")
        if any(np.shape(data[key]) != (ny, nx) for key in ("native_hard",)):
            raise ValueError("Scene invalidity grid mismatch")
        self.graphics.clear()
        self._paths.clear()
        self.graphics.setSceneRect(0., 0., float(nx), float(ny))
        self._image = QGraphicsPixmapItem(FITSExplorer._pixmap(data["rgb"]))
        self.graphics.addItem(self._image)
        self._native_overlay = QGraphicsPixmapItem(
            FITSExplorer._pixmap(data["hard_rgba"], alpha=True))
        self._native_overlay.setZValue(2.)
        self.graphics.addItem(self._native_overlay)
        self._loaded = data
        self._redraw_sga_prior()
        self._update_hard_overlay()
        self._redraw_outlines()
        self.view.fit_image()
        self.status.setText(
            f"{data['source']} · original {nx}×{ny} px cutout, "
            f"{data['display_width_arcmin']:g}′ view of "
            f"{data['project_field_width_arcmin']:g}′ project field. "
            "Native bad pixels shown; all masks are DRAFT only.")
        self._refresh_controls()

    def _rerender_image(self, *_):
        if self._loaded is None or self._image is None:
            return
        contrast, scaling = self.display_mode.currentData()
        raster = rgb_for_layer(
            self._loaded["data"]["image"], "image", self._band,
            contrast=contrast, scaling=scaling)
        self._image.setPixmap(FITSExplorer._pixmap(raster))

    def _update_hard_overlay(self, *_):
        if self._native_overlay is not None:
            self._native_overlay.setVisible(self.show_hard.isChecked())

    def _clear_isophote_preview(self, *_):
        """Discard display-only candidate, without changing saved scene edits."""
        if self._isophote_preview is not None:
            if self._isophote_preview.scene() is self.graphics:
                self.graphics.removeItem(self._isophote_preview)
            self._isophote_preview = None
        self._isophote_proposal = None
        for item in self._sga_shred_paths:
            if item.scene() is self.graphics:
                self.graphics.removeItem(item)
        self._sga_shred_paths.clear()
        self.iso_shape.blockSignals(True)
        seed_index = self.iso_shape.findData("sga_seed")
        if seed_index >= 0:
            self.iso_shape.removeItem(seed_index)
        self.iso_shape.blockSignals(False)
        self.iso_sga_summary.setText(
            "No SGA-informed radial proposal on this band.")
        self.iso_candidate.blockSignals(True)
        self.iso_candidate.clear()
        self.iso_candidate.blockSignals(False)
        self.iso_summary.setText(
            "No current proposal. Run masked detection to compare contours.")
        self.iso_low_conf_override.setChecked(False)
        self._refresh_controls()

    def _discard_prepared_masks(self):
        """Invalidate candidate masks whenever target, source or edits change."""
        self._clear_isophote_preview()
        if self._contamination_overlay is not None:
            if self._contamination_overlay.scene() is self.graphics:
                self.graphics.removeItem(self._contamination_overlay)
            self._contamination_overlay = None
        for item in self._conflict_paths:
            if item.scene() is self.graphics:
                self.graphics.removeItem(item)
        self._conflict_paths.clear()
        self._contamination = None
        self.iso_mask_details.setText("")
        self.iso_mask_status.setText(
            "Neighbour masks need preparation for the current "
            "source, band, target centre and scene-draft revision.")
        self._refresh_controls()

    def _find_desi_sga2025(self):
        """On-demand official SGA2025 fitted-ellipse Viewer overlay."""
        self._request_sga2025_viewer_layer("ellipse")

    def _find_desi_sga2025_parent(self):
        """Explicit SGA2025 parent initial geometry, not a fitted isophote."""
        self._request_sga2025_viewer_layer("parent")

    def _request_sga2025_viewer_layer(self, layer):
        if (self._project_root is None or self._project_target is None
                or self._future is not None):
            return
        project = read_object_project(self._project_root)
        if (project.selected_dataset is None or
                project.selected_dataset[0] != "dr11" or
                project.selected_dataset[1] not in ("north", "south")):
            self.sga_status.setText(
                "Select a DR11 North/South source in the object project "
                "before querying the matching SGA2025 Viewer region.")
            return
        region = project.selected_dataset[1]
        self._sga_lookup = None
        self._sga_candidate = None
        self.sga_candidates.blockSignals(True)
        self.sga_candidates.clear()
        self.sga_candidates.blockSignals(False)
        self.sga_use_prior.setChecked(False)
        self.sga_use_seed.setChecked(False)
        self._redraw_sga_prior()
        self._discard_prepared_masks()
        options = {"region":region, "radius_arcmin":3.}
        if layer == "parent":
            options["layer"] = "parent"
        self._future = self._executor.submit(
            fetch_viewer_sga2025_priors, self._project_target, **options)
        self._work_kind = "sga"
        self._active_generation = self._load_generation
        kind = ("parent initial geometry" if layer == "parent"
                else "fitted-ellipse overlay")
        self.sga_status.setText(
            f"Querying official DESI Viewer SGA2025 {kind} by "
            f"{region} sky coordinates. Version-mutable small JSON; "
            "no full FITS download, and nothing science-adopted.")
        self._poll_timer.start()
        self._refresh_controls()

    def _open_sga_file(self):
        """Select an existing SGA catalog: no downloads, no extra folders."""
        if (self._project_root is None or self._project_target is None
                or self._future is not None):
            return
        path, _ = QFileDialog.getOpenFileName(
            self, "Select local SGA-2020 or SGA-2025 FITS",
            "", "FITS catalog (*.fits *.fit *.fits.gz);;All files (*)")
        if not path:
            return
        self._sga_lookup = None
        self._sga_candidate = None
        self.sga_candidates.blockSignals(True)
        self.sga_candidates.clear()
        self.sga_candidates.blockSignals(False)
        self.sga_use_prior.setChecked(False)
        self.sga_use_seed.setChecked(False)
        self._redraw_sga_prior()
        self._discard_prepared_masks()
        self._future = self._executor.submit(
            read_local_sga_priors, path, self._project_target)
        self._work_kind = "sga"
        self._active_generation = self._load_generation
        self.sga_status.setText(
            "Reading and hashing selected local SGA FITS in background. "
            "No galaxy automatically matched or adopted.")
        self._poll_timer.start()
        self._refresh_controls()

    def _accept_sga_lookup(self, lookup):
        if self._project_id is None or self._project_target is None:
            raise ValueError("SGA result belongs to an inactive project")
        if lookup.status != "SGA_LOCAL_MATCHES_REQUIRE_RESEARCHER_SELECTION":
            raise ValueError("SGA geometry result has unexpected status")
        self._sga_lookup = lookup
        self._sga_candidate = None
        self.sga_candidates.blockSignals(True)
        self.sga_candidates.clear()
        self.sga_candidates.addItem(
            f"Select matched galaxy ({len(lookup.matches)} candidates)", None)
        for row in lookup.matches:
            identifier = (
                f"Parent PGC {row.sga_id}"
                if row.d26_ref == "VIEWER_PARENT_DIAMETER" else
                f"Viewer ref {row.sga_id}"
                if row.d26_ref == "VIEWER_RADIUS_ARCSEC" else
                f"SGAID {row.sga_id}")
            label = (
                f"{row.name} · {row.release} · {identifier} · "
                f"{row.offset_arcsec:.1f} arcsec")
            self.sga_candidates.addItem(label, row)
        self.sga_candidates.setCurrentIndex(0)
        self.sga_candidates.blockSignals(False)
        unusable = getattr(lookup, "unusable_count", 0)
        invalid_notes = getattr(lookup, "unusable_examples", ())
        if lookup.matches:
            origin = ("DESI Viewer overlay" if
                      lookup.source_path.startswith("https://")
                      else "Local FITS")
            self.sga_status.setText(
                f"{origin} {lookup.release}: "
                f"{len(lookup.matches)} usable geometric candidate(s). "
                + (f"{unusable} nearby Viewer entry/entries were rejected "
                   "because no valid ellipse exists. " if unusable else "")
                + "Select the correct galaxy explicitly; the overlay "
                  "is not independently validated outer LSB geometry.")
        elif unusable:
            self.sga_status.setText(
                f"SGA2025 Viewer returned {unusable} nearby record(s), "
                "but ZERO usable initial ellipses. "
                + ("Examples: " + "; ".join(invalid_notes) + ". "
                   if invalid_notes else "")
                + "Do not use these records as galaxy geometry or "
                  "enable SGA-based Tractor protection. "
                  "Try the separate SGA2025 Parent initial geometry "
                  "layer, or continue without SGA.")
        else:
            self.sga_status.setText(
                "No SGA match with usable geometry within 3 arcmin of the saved "
                "object centre. Continue without SGA; native Tractor "
                "screening still works.")
        self._redraw_sga_prior()
        self._refresh_controls()

    def _sga_candidate_changed(self, *_):
        prior = self.sga_candidates.currentData()
        if (prior is not None and (
                self._sga_lookup is None or
                prior.source_sha256 != self._sga_lookup.source_sha256
                or prior not in self._sga_lookup.matches)):
            prior = None
        self._sga_candidate = prior
        # Explicit row selection permits a read-only geometry preview;
        # mask protection remains a separate, unchecked opt-in.
        self.sga_use_seed.setChecked(
            prior is not None and prior.core_shape() is not None)
        self.sga_use_prior.setChecked(False)
        self._discard_prepared_masks()
        self._redraw_sga_prior()
        if prior is None:
            return
        d26 = (f"{prior.d26_arcmin:.3f} arcmin"
               if prior.d26_arcmin is not None else "unavailable")
        radius_label = (
            "Parent initial diameter"
            if prior.d26_ref == "VIEWER_PARENT_DIAMETER" else
            "Viewer overlay diameter"
            if prior.d26_ref == "VIEWER_RADIUS_ARCSEC" else "D26")
        identifier = (
            "Parent PGC"
            if prior.d26_ref == "VIEWER_PARENT_DIAMETER" else
            "Viewer ref"
            if prior.d26_ref == "VIEWER_RADIUS_ARCSEC" else "SGAID")
        pa = (f"{prior.pa_deg:.1f} deg"
              if prior.pa_deg is not None else "unavailable")
        q = (f"{prior.ba:.3f}"
             if prior.ba is not None else "unavailable")
        self.sga_status.setText(
            f"{prior.release} {identifier} {prior.sga_id} · "
            f"{radius_label}={d26} "
            f"(ref={prior.d26_ref or 'unknown'}); "
            f"q={q}, PA={pa}; "
            f"offset={prior.offset_arcsec:.1f} arcsec. "
            + ("Gold dashed ellipse = CENTRAL PRIOR only. "
               if prior.core_shape() is not None
               else "NO USABLE SGA ELLIPSE: nothing drawn or protected. ")
            + "May miss LSB outskirts. SGA/Tractor not independent. "
            + ("Flags: " + ", ".join(prior.flags)
               if prior.flags else "Published geometry still needs QA."))
        self._refresh_controls()

    def _sga_protection_changed(self, *_):
        self._discard_prepared_masks()

    def _redraw_sga_prior(self):
        if self._sga_item is not None:
            if self._sga_item.scene() is self.graphics:
                self.graphics.removeItem(self._sga_item)
            self._sga_item = None
        prior = self._sga_candidate
        if prior is None or self._loaded is None:
            return
        shape = prior.core_shape()
        if shape is None:
            return
        try:
            path = scene_shape_outline(
                shape, self._loaded["wcs"], self._loaded["pixel_shape"][0])
        except (TypeError, ValueError) as exc:
            self.sga_status.setText(
                "SGA prior cannot be projected on this WCS: " + str(exc))
            return
        item = QGraphicsPathItem(path)
        pen = QPen(QColor("#e2bf62"), 2.25, Qt.DashLine)
        pen.setCosmetic(True)
        item.setPen(pen)
        item.setZValue(3.3)
        self.graphics.addItem(item)
        self._sga_item = item

    def _prepare_neighbour_masks(self):
        if (self._loaded is None or self._future is not None or
                self._project_root is None or self._draft is None or
                self._project_target is None):
            return
        self._discard_prepared_masks()
        centre = (self._isophote_center_icrs
                  if self._isophote_center_icrs is not None
                  else self._project_target)
        prior = (self._sga_candidate
                 if self.sga_use_prior.isChecked() else None)
        arguments = (
            self._loaded, str(self._project_root), self._draft, centre)
        # Keep existing non-SGA work flow and its test contracts intact.
        self._future = self._executor.submit(
            prepare_scene_neighbour_masks, *arguments,
            *((prior,) if prior is not None else ()))
        self._work_kind = "contamination"
        self._active_generation = self._load_generation
        self.iso_mask_status.setText(
            "Reading local Tractor source identities and preparing "
            "provisional, non-destructive native masks in background…")
        self._poll_timer.start()
        self._refresh_controls()

    def _accept_contamination_result(self, result):
        if (self._loaded is None or self._draft is None or
                result.get("status") !=
                "CATALOGUE_CONTAMINATION_MASK_NOT_ADOPTED" or
                result["project_id"] != self._project_id or
                result["dataset_source"] != self._loaded["source"] or
                result["band"] != self._band or
                result["draft_digest"] != self._draft.digest or
                result.get("sga_prior_id") != (
                    self._sga_candidate.sga_id
                    if self.sga_use_prior.isChecked()
                    and self._sga_candidate is not None else None) or
                result["anchor_icrs"] !=
                (self._isophote_center_icrs
                 if self._isophote_center_icrs is not None
                 else self._project_target)):
            raise ValueError("Prepared masks do not match current scene")
        if np.shape(result["mask"]) != tuple(self._loaded["pixel_shape"]):
            raise ValueError("Contaminant mask grid mismatch")
        self._contamination = result
        ny, nx = self._loaded["pixel_shape"]
        layer = np.zeros((ny, nx, 4), dtype=np.uint8)
        layer[result["auto_mask"]] = (175, 106, 238, 180)
        self._contamination_overlay = QGraphicsPixmapItem(
            FITSExplorer._pixmap(layer, alpha=True))
        self._contamination_overlay.setZValue(2.5)
        self.graphics.addItem(self._contamination_overlay)
        # Amber review rings are source identities whose provisional
        # footprint touches protected target light. They are NOT masked.
        self._conflict_paths.clear()
        anchor_x, anchor_y = self._loaded["wcs"].all_world2pix(
            *result["anchor_icrs"], 0)
        nearest = sorted(
            result["stats"].get("near_target_sources", ()),
            key=lambda item:
                (item["x"]-anchor_x)**2+(item["y"]-anchor_y)**2)
        for candidate in nearest[:20]:
            radius = min(12., max(4., candidate["radius_pix"]*.30))
            cx, cy = candidate["x"]+.5, ny-candidate["y"]-.5
            ring = QPainterPath()
            ring.addEllipse(cx-radius, cy-radius, 2.*radius, 2.*radius)
            item = QGraphicsPathItem(ring)
            pen = QPen(QColor("#ffba65"), 2., Qt.DashLine)
            pen.setCosmetic(True)
            item.setPen(pen)
            item.setZValue(3.5)
            self.graphics.addItem(item)
            self._conflict_paths.append(item)
        self._update_contamination_overlay()
        counts = result["stats"]
        self.iso_mask_status.setText(
            f"PRELIMINARY {self._band.upper()} · "
            f"{counts['masked_catalog_sources']} masked neighbours "
            f"({counts.get('elliptical_catalog_sources', 0)} morphology-based), "
            f"{counts['skipped_target_overlap']} unresolved overlaps"
            + (f"; {counts.get('possible_sga_shreds', 0)} inside "
               "chosen SGA central ellipse" if result.get("sga_prior_id")
               is not None else "")
            + ". Purple = candidate excludes; enable amber rings "
            "to inspect the nearest 20.")
        self.iso_mask_details.setText(
            f"{counts['excluded_pixels']} proposed excluded pixels; "
            f"{counts['native_hard_pixels']} hard-invalid original pixels. "
            f"Catalogue SHA-256 {result['catalog_sha256']}; "
            f"SGA {result.get('sga_prior_id') or 'none'} "
            f"SHA-256 {result.get('sga_prior_sha256') or 'none'}; "
            f"{counts.get('missing_shape_fallbacks', 0)} shape fallbacks. "
            "Original FITS unchanged. NOT an exact Tractor pixel model, "
            "PSF-wing subtraction or a calibrated science mask. Flags: "
            + ", ".join(result["flags"]))
        self._refresh_controls()

    def _update_contamination_overlay(self, *_):
        if self._contamination_overlay is not None:
            self._contamination_overlay.setVisible(
                self.iso_mask_toggle.isChecked())
            self._contamination_overlay.setOpacity(
                self.iso_opacity.value()/100.)
        for item in self._conflict_paths:
            item.setVisible(self.iso_show_conflicts.isChecked())
        for item in self._sga_shred_paths:
            item.setVisible(self.iso_show_sga_shreds.isChecked())

    def _pick_isophote_center(self):
        if self._loaded is None:
            return
        self._discard_prepared_masks()
        self._cancel_pending()
        self._picking_isophote_center = True
        self.view.set_drawing(True)
        self.iso_details.setText(
            "Click the diffuse galaxy centre once on the ORIGINAL FITS; "
            "this is an anchor for connected emission, not a fitted centre.")

    def _reset_isophote_center(self):
        self._picking_isophote_center = False
        self._isophote_center_icrs = None
        self.view.set_drawing(self.mode.currentData() != "pan")
        self._discard_prepared_masks()
        self.iso_details.setText(
            "Detection anchor reset to saved galaxy ICRS coordinates. "
            "No image pixels or masks changed.")

    def _generate_isophote(self):
        if (self._loaded is None or self._future is not None or
                self._project_target is None or self._draft is None or
                self._contamination is None):
            self.iso_details.setText(
                "Prepare and inspect the current Tractor catalogue "
                "neighbour mask before isophotal detection. "
                "Raw-image detection is disabled.")
            return
        self._picking_isophote_center = False
        self.view.set_drawing(self.mode.currentData() != "pan")
        self._clear_isophote_preview()
        center = (self._isophote_center_icrs
                  if self._isophote_center_icrs is not None
                  else self._project_target)
        sga = (self._sga_candidate if self.sga_use_seed.isChecked()
               else None)
        options = {"sga_prior":sga} if sga is not None else {}
        self._future = self._executor.submit(
            propose_scene_isophote, self._loaded, center,
            float(self.iso_threshold.value()),
            float(self.iso_smooth.value()), self._contamination,
            2 if self.iso_bridge_guard.isChecked() else 0, **options)
        self._work_kind = "isophote"
        self._active_generation = self._load_generation
        self.iso_summary.setText("Analysing four nearby thresholds/scales…")
        self.iso_details.setText(
            "Finding target-anchored emission and comparing three nearby "
            "threshold/smoothing variants on the same masked native FITS. "
            "All results are preliminary.")
        self._poll_timer.start()
        self._refresh_controls()

    def _accept_isophote_result(self, result):
        if (self._loaded is None or
                result["status"] != "LOCAL_ISOPHOTE_PROPOSAL_NOT_ADOPTED" or
                result["project_id"] != self._project_id or
                result["band"] != self._band or
                result["source"] != self._loaded["source"]):
            raise ValueError("Proposed isophote belongs to an old source/band")
        for kind in ("contour", "ellipse"):
            scene_shape_outline(
                result[kind], self._loaded["wcs"],
                self._loaded["pixel_shape"][0])
        if (self._contamination is None or
                self._contamination["catalog_sha256"] != result["catalog_sha256"]
                or self._contamination["draft_digest"] != result["draft_digest"]):
            raise ValueError("Neighbour masks changed during detection")
        comparison = result.get("sga_comparison")
        if comparison is not None and (
                not self.sga_use_seed.isChecked() or
                self._sga_candidate is None or
                self._sga_candidate.sga_id != comparison["sga_ref"] or
                self._sga_candidate.source_sha256 !=
                comparison["sga_source_sha256"]):
            raise ValueError("SGA geometry changed during radial analysis")
        self._isophote_proposal = result
        for item in self._sga_shred_paths:
            if item.scene() is self.graphics:
                self.graphics.removeItem(item)
        self._sga_shred_paths.clear()
        candidates = ((comparison or {}).get("review_candidates") or ())
        ny = self._loaded["pixel_shape"][0]
        for entry in candidates:
            x = float(entry["x"])+.5
            y = ny-float(entry["y"])-.5
            radius = 6.
            ring = QPainterPath()
            ring.addEllipse(x-radius,y-radius,2*radius,2*radius)
            item = QGraphicsPathItem(ring)
            pen = QPen(QColor("#66c9ff"),2.,Qt.DashLine)
            pen.setCosmetic(True)
            item.setPen(pen)
            item.setZValue(3.7)
            self.graphics.addItem(item)
            self._sga_shred_paths.append(item)
        self._update_contamination_overlay()
        self.iso_low_conf_override.setChecked(False)
        self.iso_shape.blockSignals(True)
        if result.get("sga_seed_shape") is not None:
            if self.iso_shape.findData("sga_seed") < 0:
                self.iso_shape.addItem(
                    "SGA-informed initial ellipse · supported radius",
                    "sga_seed")
        else:
            idx = self.iso_shape.findData("sga_seed")
            if idx >= 0:
                self.iso_shape.removeItem(idx)
        self.iso_shape.blockSignals(False)
        if comparison is not None:
            radial = result.get("sga_radial")
            outer = (f"{radial.detected_a_arcsec:.1f} arcsec"
                     if radial is not None
                     and radial.detected_a_arcsec is not None
                     else "NOT detected")
            masked = comparison.get(
                "masked_catalog_centres_in_initial_ellipse", 0)
            self.iso_sga_summary.setText(
                f"SGA INITIAL: a={comparison['reference_a_arcsec']:.1f}\" "
                f"→ data-supported a={outer}; "
                f"SGA q={comparison['reference_q']:.3f}, "
                f"PA={comparison['reference_pa_deg']:.1f}°. "
                f"Independent moments: q={comparison['independent_q']:.3f}, "
                f"PA={comparison['independent_pa_deg']:.1f}° "
                f"(Δq={comparison['delta_q']:.3f}, "
                f"ΔPA={comparison['delta_pa_deg']:.1f}°). "
                f"SGA–target centre offset="
                f"{comparison['center_offset_arcsec']:.1f}\". "
                f"{masked} masked Tractor centres inside proposal "
                "need source-by-source review; none restored."
                + ("" if radial is not None else
                   " No supported extension: " +
                   str(result.get('sga_warning') or 'insufficient evidence')))
        else:
            self.iso_sga_summary.setText(
                "No SGA geometry used. Select an SGA prior to compare "
                "a data-supported initial ellipse with independent moments.")
        self.iso_candidate.blockSignals(True)
        self.iso_candidate.clear()
        for index, trial in enumerate(result["sensitivity"]["trials"]):
            variant = (result.get("alternatives") or ())[index] if (
                result.get("alternatives")) else (
                {"model":result["proposal"]} if index == 0 else None)
            if variant is None:
                continue
            label = (
                f"{'Reference' if index == 0 else 'Alternative '+str(index)}"
                f" · {trial['threshold']:.2f} contrast / "
                f"{trial['smooth_pix']:.1f} px / "
                f"{trial['connected_pixels']:,} pixels")
            self.iso_candidate.addItem(label, index)
        self.iso_candidate.setCurrentIndex(0)
        self.iso_candidate.blockSignals(False)
        model = result["proposal"]
        sensitivity = result["sensitivity"]
        low = sensitivity["confidence"].startswith("LOW_CONFIDENCE")
        label = ("LOW CONFIDENCE · manual review required" if low
                 else "PRELIMINARY · stable enough to inspect")
        self.iso_summary.setText(
            f"{label}. {self._band.upper()} connected area "
            f"{model.seed.component_pixels:,} px; q={model.seed.q:.2f}; "
            f"centroid offset {model.seed.centroid_offset_pix:.1f} px. "
            f"Sensitivity area ratio {sensitivity['area_ratio']:.2f}; "
            f"min spatial IoU "
            f"{sensitivity['minimum_spatial_iou'] if sensitivity['minimum_spatial_iou'] is not None else float('nan'):.2f} "
            "across four trials. "
            + ("Review flags; acceptance blocked until checkbox confirmed."
               if low else "Review before adding to draft."))
        trial_text = "; ".join(
            f"{trial['threshold']:.2g} / {trial['smooth_pix']:.2g}px: "
            + (f"{trial['connected_pixels']}px, "
               f"q={trial['q']:.2f}, "
               f"drift={trial['centroid_offset_pix']:.1f}px"
               if trial["state"]=="DETECTED"
               else "NO COMPONENT")
            for trial in sensitivity["trials"])
        self.iso_details.setText(
            f"PROVISIONAL {self._band.upper()} band-specific: "
            f"local sky {model.local_sky:.4g} nMgy/pixel; "
            f"smoothed RMS {model.smoothed_sky_rms:.4g} nMgy/pixel; "
            f"masked neighbours "
            f"{result['contamination']['masked_catalog_sources']}. "
            f"Bridge opening {result['bridge_guard_pix']} native px. "
            f"Sensitivity trials: {trial_text}. "
            f"QA warnings: {', '.join(sensitivity['flags']) or 'none'}; "
            f"{', '.join(model.flags)}. NOT a calibrated outer isophote; "
            "PSF wings, bright neighbors, variable sky, "
            "incomplete contamination and correlated-noise significance "
            "must be independently reviewed.")
        self._show_isophote_preview()

    def _selected_isophote(self):
        """Return active trial in memory, without saving or refitting FITS."""
        proposal = self._isophote_proposal
        if proposal is None:
            return None
        alternatives = proposal.get("alternatives")
        if not alternatives:
            # Legacy test payloads and initial standalone SC1 proposals.
            return {"model": proposal["proposal"],
                    "contour":proposal["contour"],
                    "ellipse":proposal["ellipse"]}
        index = self.iso_candidate.currentData()
        if not isinstance(index, int) or not 0 <= index < len(alternatives):
            return None
        return alternatives[index]

    def _show_isophote_preview(self, *_):
        if self._isophote_preview is not None:
            if self._isophote_preview.scene() is self.graphics:
                self.graphics.removeItem(self._isophote_preview)
            self._isophote_preview = None
        if self._isophote_proposal is None or self._loaded is None:
            self._refresh_controls()
            return
        candidate = self._selected_isophote()
        if candidate is None:
            self._refresh_controls()
            return
        kind = self.iso_shape.currentData()
        shape = (self._isophote_proposal.get("sga_seed_shape")
                 if kind == "sga_seed" else candidate[kind])
        if shape is None:
            self._refresh_controls()
            return
        try:
            path = scene_shape_outline(
                shape, self._loaded["wcs"], self._loaded["pixel_shape"][0])
        except (ValueError, TypeError, IndexError) as exc:
            self.iso_details.setText("Cannot display proposal: " + str(exc))
            self._isophote_proposal = None
            self._refresh_controls()
            return
        item = QGraphicsPathItem(path)
        low = self._isophote_proposal["sensitivity"][
            "confidence"].startswith("LOW_CONFIDENCE")
        pen = QPen(QColor("#ffa767" if low else "#63e0ef"),
                   2.5, Qt.DashLine)
        pen.setCosmetic(True)
        item.setPen(pen)
        item.setZValue(4.)
        self.graphics.addItem(item)
        self._isophote_preview = item
        if self.iso_candidate.count():
            model = candidate["model"]
            summary = self._isophote_proposal["sensitivity"]
            label = ("LOW CONFIDENCE" if low else "PRELIMINARY")
            self.iso_summary.setText(
                ("SGA-derived PA/q with empirical radius · "
                 if kind == "sga_seed" else "")
                + f"{label} · trial {self.iso_candidate.currentData()} · "
                f"{model.seed.component_pixels:,} pixels; "
                f"q={model.seed.q:.2f}; "
                f"centroid offset={model.seed.centroid_offset_pix:.1f} px. "
                f"Across four trials: area ratio "
                f"{summary['area_ratio']:.2f}, min IoU "
                f"{summary.get('minimum_spatial_iou') if summary.get('minimum_spatial_iou') is not None else float('nan'):.2f}. "
                + ("Manual-review override required." if low
                   else "No science adoption."))
            original_qa = self._isophote_proposal["sensitivity"]
            trial_text = "; ".join(
                f"{trial['threshold']:.2f}/{trial['smooth_pix']:.1f}px:"
                + (f"{trial['connected_pixels']}px, "
                   f"IoU={trial.get('overlap_to_reference')}"
                   if trial["state"] == "DETECTED"
                   else "FAILED")
                for trial in original_qa["trials"])
            self.iso_details.setText(
                f"Selected trial {self.iso_candidate.currentData()} "
                f"in {self._band.upper()}: contrast "
                f"{model.threshold_contrast:.2f}, "
                f"Gaussian smoothing {model.smoothing_pix:.1f}px, "
                f"local sky {model.local_sky:.4g} nMgy/pixel, "
                f"smoothed RMS {model.smoothed_sky_rms:.4g}. "
                f"All trials: {trial_text}. "
                f"QA flags: {', '.join(original_qa['flags']) or 'none'}. "
                + ("SGA radial QA: "
                   + ", ".join(self._isophote_proposal["sga_radial"].flags)
                   + ". " if self._isophote_proposal.get("sga_radial")
                   is not None else "")
                + "NOT a calibrated outer isophote: nearby stars, "
                "PSF wings, spatial sky structure and unmeasured "
                "correlated noise may still bias the boundary.")
        self._refresh_controls()

    def _accept_isophote(self):
        """One explicit user action adds a DRAFT, never adopted photometry."""
        if (self._isophote_proposal is None or
                self._loaded is None or self._draft is None or
                self._isophote_preview is None):
            return
        proposed = self._isophote_proposal
        if (proposed["band"] != self._band or
                proposed["project_id"] != self._project_id or
                proposed["source"] != self._loaded["source"] or
                self._contamination is None or
                proposed["catalog_sha256"] !=
                self._contamination["catalog_sha256"] or
                proposed["draft_digest"] != self._draft.digest):
            self._clear_isophote_preview()
            self.iso_details.setText("Stale isophote: recalculate on current FITS.")
            return
        sensitivity = proposed["sensitivity"]
        kind = self.iso_shape.currentData()
        if ((sensitivity["confidence"].startswith("LOW_CONFIDENCE")
             or kind == "sga_seed")
                and not self.iso_low_conf_override.isChecked()):
            self.iso_summary.setText(
                "LOW CONFIDENCE: contour cannot be added until "
                "explicit manual-review override is checked.")
            return
        candidate = self._selected_isophote()
        if candidate is None:
            return
        model = candidate["model"]
        shape_for_edit = (
            proposed.get("sga_seed_shape")
            if kind == "sga_seed" else candidate[kind])
        if shape_for_edit is None:
            return
        scope = None if self.iso_share_bands.isChecked() else self._band
        provenance = (
            f"AUTO_ISOPHOTE_PRELIMINARY from_band={self._band} "
            f"trial={self.iso_candidate.currentData()} "
            f"scope={'all_optical_EXPLICIT' if scope is None else self._band} "
            f"confidence={sensitivity['confidence']} "
            f"manual_override={int(self.iso_low_conf_override.isChecked())} "
            f"kind={kind} threshold={model.threshold_contrast:g} "
            f"smooth_pix={model.smoothing_pix:g} "
            f"bridge_pix={proposed['bridge_guard_pix']} "
            f"min_iou={sensitivity.get('minimum_spatial_iou')} "
            + (f"SGAref={proposed['sga_comparison']['sga_ref']} "
               f"SGAsha={proposed['sga_comparison']['sga_source_sha256'][:16]} "
               f"SGAinitial_a={shape_for_edit.a_arcsec:.2f} "
               f"SGAinitial_q={shape_for_edit.q:.3f} "
               f"SGAinitial_PA={shape_for_edit.pa_deg:.1f} "
               if kind == "sga_seed" else "")
            + f"sky={model.local_sky:.7g} "
            f"rms={model.smoothed_sky_rms:.7g} "
            f"TractorSHA256={proposed['catalog_sha256']} "
            f"source={proposed['source']}")
        if kind == "sga_seed":
            cmp = proposed["sga_comparison"]
            # Compact origin-first metadata must survive the 512-character
            # draft-comment bound. Keep FULL, not abbreviated, catalogue
            # hashes for independently checkable provenance.
            provenance = (
                f"SGA_INITIAL_GEOMETRY_DRAFT kind=sga_seed "
                f"from_band={self._band} "
                f"scope={'all_EXPLICIT' if scope is None else self._band} "
                f"SGAref={cmp['sga_ref']} "
                f"SGAsha256={cmp['sga_source_sha256']} "
                f"TractorSHA256={proposed['catalog_sha256']} "
                f"a_arcsec={shape_for_edit.a_arcsec:.2f} "
                f"q={shape_for_edit.q:.4f} "
                f"PA_deg={shape_for_edit.pa_deg:.2f} "
                f"masked_candidates={cmp.get('masked_catalog_centres_in_initial_ellipse',0)} "
                f"threshold={model.threshold_contrast:g} "
                f"smooth_pix={model.smoothing_pix:g} "
                "MANUAL_REVIEW_NO_SCIENCE_ADOPTION")
            comment = (provenance +
                       (" | " + self.comment.text().strip()
                        if self.comment.text().strip() else ""))[:512]
        else:
            comment = (
                self.comment.text().strip() + " | " + provenance
            ).strip(" |")[:512]
        try:
            edit = SkyMaskEdit(
                "sc1_iso_" + uuid4().hex,
                "protect_target",
                shape_for_edit,
                band=scope,
                comment=comment)
            self._push_edits(self._draft.edits + (edit,))
        except (ValueError, TypeError) as exc:
            self.iso_details.setText("Proposed region rejected: " + str(exc))
            return
        self._clear_isophote_preview()
        self.iso_summary.setText(
            "Added provisional " + kind + " as "
            + ("explicit shared ICRS geometry" if scope is None
               else f"{self._band.upper()}-only geometry")
            + ". No science is adopted.")
        self.iso_details.setText(
            "Provisional " + kind + " added to the scene DRAFT. "
            "Undo/Remove or Redraw selected can revise it; it will be "
            "autosaved as an unadopted user-reviewed mask. No science adopted.")

    def _mode_changed(self, *_):
        primitive = self.mode.currentData()
        self.view.set_drawing(primitive != "pan")
        self.finish_btn.setEnabled(primitive == "polygon")
        self.axis_ratio.setEnabled(primitive == "ellipse")
        self._cancel_pending()

    def _cancel_pending(self):
        self._pending_world.clear()
        self._pending_scene.clear()
        self._replace_index = None
        if self._preview_path is not None:
            self.graphics.removeItem(self._preview_path)
            self._preview_path = None

    def _update_pending_path(self):
        if self._preview_path is not None:
            self.graphics.removeItem(self._preview_path)
            self._preview_path = None
        if not self._pending_scene:
            return
        path = QPainterPath()
        for j, (x, y) in enumerate(self._pending_scene):
            if j == 0:
                path.moveTo(x, y)
            else:
                path.lineTo(x, y)
        self._preview_path = QGraphicsPathItem(path)
        pen = QPen(QColor("#f6dc78"), 2.1, Qt.DashLine)
        pen.setCosmetic(True)
        self._preview_path.setPen(pen)
        self._preview_path.setZValue(5.)
        self.graphics.addItem(self._preview_path)

    def _on_click(self, x_scene, y_scene):
        if self._loaded is None or self._draft is None:
            return
        pixel = fits_pixel_from_scene(
            x_scene, y_scene, self._loaded["pixel_shape"])
        if pixel is None:
            return
        primitive = self.mode.currentData()
        if primitive == "pan" and not self._picking_isophote_center:
            return
        ny = self._loaded["pixel_shape"][0]
        xpix, ypix = float(x_scene)-.5, float(ny)-float(y_scene)-.5
        ra, dec = self._loaded["wcs"].all_pix2world(xpix, ypix, 0)
        if not np.isfinite(ra) or not np.isfinite(dec) or abs(dec) > 90:
            self.status.setText("Invalid WCS position; shape not added.")
            return
        if self._picking_isophote_center:
            self._picking_isophote_center = False
            self._isophote_center_icrs = (float(ra)%360., float(dec))
            self.view.set_drawing(self.mode.currentData() != "pan")
            self._discard_prepared_masks()
            self.iso_details.setText(
                f"Galaxy candidate centre: RA {ra%360:.7f}°, "
                f"Dec {dec:+.7f}°. Prepare masks, then find isophote.")
            return
        self._pending_world.append((float(ra)%360., float(dec)))
        self._pending_scene.append((float(x_scene),float(y_scene)))
        self._update_pending_path()
        if primitive != "polygon" and len(self._pending_world) == 2:
            self._finish_shape(primitive)
        elif primitive == "polygon":
            self.status.setText(
                f"Polygon: {len(self._pending_world)} vertices selected. "
                "Use Finish polygon after at least three vertices.")
        else:
            self.status.setText("Choose the region's radius or major-axis endpoint.")

    def finish_polygon(self):
        if self.mode.currentData() != "polygon":
            return
        if len(self._pending_world) < 3:
            self.status.setText("Polygon needs at least three native-image clicks.")
            return
        self._finish_shape("polygon")

    def _finish_shape(self, kind):
        points = tuple(self._pending_world)
        try:
            shape = scene_shape_from_world_points(
                kind, points, axis_ratio=float(self.axis_ratio.value()))
            if self._loaded is not None:
                scene_shape_outline(
                    shape, self._loaded["wcs"], self._loaded["pixel_shape"][0])
            index = self._replace_index
            edits = list(self._draft.edits)
            edit = SkyMaskEdit(
                edit_id=(edits[index].edit_id if index is not None
                         else "sc1_" + uuid4().hex),
                purpose=self.purpose.currentData(),
                shape=shape,
                band=None if self.scope.currentData() == "all" else self._band,
                comment=self.comment.text()[:512])
            if index is None:
                edits.append(edit)
            else:
                edits[index] = edit
            self._push_edits(edits)
        except (ValueError, TypeError, IndexError) as exc:
            self.status.setText(f"Region rejected: {exc}")
        finally:
            self._cancel_pending()

    def _push_edits(self, edits):
        if self._draft is None:
            raise RuntimeError("No object-linked scene draft")
        self._undo.append(self._draft.edits)
        self._redo.clear()
        self._draft = SkyMaskDraft(
            self._draft.project_id, tuple(edits), self._draft.revision + 1)
        self._discard_prepared_masks()
        self.status.setText(
            f"Changed scene draft: {len(edits)} regions; unsaved, "
            "not a scientifically adopted source mask.")
        self._redraw_outlines()
        self._refresh_controls()
        self._autosave_timer.start()

    def undo(self):
        if self._draft is None or not self._undo:
            return
        self._cancel_pending()
        previous = self._undo.pop()
        self._redo.append(self._draft.edits)
        self._draft = SkyMaskDraft(
            self._draft.project_id, previous, self._draft.revision + 1)
        self._discard_prepared_masks()
        self._redraw_outlines()
        self._refresh_controls()
        self._autosave_timer.start()

    def redo(self):
        if self._draft is None or not self._redo:
            return
        self._cancel_pending()
        upcoming = self._redo.pop()
        self._undo.append(self._draft.edits)
        self._draft = SkyMaskDraft(
            self._draft.project_id, upcoming, self._draft.revision + 1)
        self._discard_prepared_masks()
        self._redraw_outlines()
        self._refresh_controls()
        self._autosave_timer.start()

    def edit_selected(self):
        if self._draft is None or self._loaded is None:
            return
        index = self.edits.currentRow()
        if not 0 <= index < len(self._draft.edits):
            return
        edit = self._draft.edits[index]
        if not edit.applies_to(self._band):
            self.status.setText(
                "Selected region is local to another band. "
                "Switch to that band to redraw its native geometry.")
            return
        self.mode.setCurrentIndex(self.mode.findData(edit.shape.kind))
        self._cancel_pending()
        self._replace_index = index
        self.purpose.setCurrentIndex(self.purpose.findData(edit.purpose))
        self.scope.setCurrentIndex(
            self.scope.findData("all" if edit.band is None else "band"))
        self.comment.setText(edit.comment)
        if edit.shape.kind == "ellipse":
            self.axis_ratio.setValue(edit.shape.q)
        self.status.setText(
            f"Redraw selected {edit.shape.kind} on original FITS "
            "using two clicks or polygon vertices. Edit identity "
            "is preserved; Cancel shape aborts without changing draft.")

    def remove_selected(self):
        if self._draft is None:
            return
        index = self.edits.currentRow()
        if index < 0 or index >= len(self._draft.edits):
            return
        self._push_edits(
            self._draft.edits[:index] + self._draft.edits[index+1:])

    def save_draft(self):
        self._autosave_timer.stop()
        if self._project_root is None or self._draft is None:
            return
        try:
            digest = save_mask_draft(
                self._project_root, self._draft,
                expected_digest=self._saved_digest)
        except (OSError, ValueError, KeyError, TypeError) as exc:
            self.status.setText(
                f"Draft NOT saved: {exc}. Refresh/reconcile elsewhere; "
                "no source FITS or previous mask revision overwritten.")
            return
        self._saved_digest = digest
        self.status.setText(
            f"Saved immutable scene DRAFT {digest[:12]} · "
            f"revision {self._draft.revision}. Not science adopted.")
        self._refresh_controls()

    def flush_autosave(self):
        """Persist changed ICRS drafts only; refuse stale concurrent HEAD."""
        self._autosave_timer.stop()
        if (self._project_root is None or self._draft is None or
                not self._available or
                self._draft.digest == self._saved_digest):
            return
        try:
            digest = save_mask_draft(
                self._project_root, self._draft,
                expected_digest=self._saved_digest)
        except (OSError, ValueError, KeyError, TypeError) as exc:
            self.status.setText(
                "SCENE DRAFT AUTOSAVE FAILED: " + str(exc)
                + ". Your unsaved edits remain in this session; "
                  "reconcile the external revision before closing.")
            return
        self._saved_digest = digest
        self.status.setText(
            f"Scene draft autosaved · {digest[:12]} · "
            f"{len(self._draft.edits)} regions · NOT SCIENCE ADOPTED.")
        self._refresh_controls()

    def _redraw_outlines(self):
        for item in self._paths:
            self.graphics.removeItem(item)
        self._paths.clear()
        self.edits.clear()
        if self._draft is None:
            return
        for edit in self._draft.edits:
            self.edits.addItem(
                f"{edit.purpose} · {edit.shape.kind} · "
                f"{'all bands' if edit.band is None else edit.band.upper()}")
            if self._loaded is None or not edit.applies_to(self._band):
                continue
            try:
                path = scene_shape_outline(
                    edit.shape, self._loaded["wcs"],
                    self._loaded["pixel_shape"][0])
            except (ValueError, IndexError, TypeError):
                continue
            item = QGraphicsPathItem(path)
            pen = QPen(QColor(PURPOSE_COLORS[edit.purpose]), 1.9)
            pen.setCosmetic(True)
            item.setPen(pen)
            item.setZValue(3.)
            self.graphics.addItem(item)
            self._paths.append(item)

    def _inspect(self, x_scene, y_scene):
        if self._loaded is None:
            return
        pixel = fits_pixel_from_scene(
            x_scene, y_scene, self._loaded["pixel_shape"])
        if pixel is None:
            self.pixel_info.setText("Outside original FITS cutout.")
            return
        x, y = pixel
        ra, dec = self._loaded["wcs"].all_pix2world(float(x), float(y), 0)
        self.pixel_info.setText(
            f"Pixel ({x},{y}) · ICRS {ra%360:.7f}°, {dec:+.7f}° · "
            f"IMAGE {float(self._loaded['data']['image'][y,x]):.6g} nMgy/pixel · "
            f"{'UNUSABLE' if self._loaded['native_hard'][y,x] else 'native valid'}")

    def _refresh_controls(self):
        ready = self._draft is not None and self._available
        active = ready and self._loaded is not None
        self.load_btn.setEnabled(ready and (
            self._future is None or self._future.done()))
        self.mode.setEnabled(active)
        self.purpose.setEnabled(active)
        self.scope.setEnabled(active)
        self.finish_btn.setEnabled(active and self.mode.currentData() == "polygon")
        self.comment.setEnabled(active)
        self.save_btn.setEnabled(ready and (
            self._saved_digest is None or self._draft.digest != self._saved_digest))
        self.undo_btn.setEnabled(ready and bool(self._undo))
        self.redo_btn.setEnabled(ready and bool(self._redo))
        self.remove_btn.setEnabled(
            ready and 0 <= self.edits.currentRow() < len(self._draft.edits))
        self.edit_btn.setEnabled(
            active and 0 <= self.edits.currentRow() < len(self._draft.edits))
        worker_idle = self._future is None or self._future.done()
        self.sga_find.setEnabled(
            ready and worker_idle and self._project_target is not None)
        self.sga_parent_find.setEnabled(
            ready and worker_idle and self._project_target is not None)
        self.sga_load.setEnabled(ready and worker_idle
                                 and self._project_target is not None)
        self.sga_candidates.setEnabled(
            ready and worker_idle and self.sga_candidates.count() > 0)
        self.sga_use_prior.setEnabled(
            ready and worker_idle and self._sga_candidate is not None
            and self._sga_candidate.core_shape() is not None)
        self.sga_use_seed.setEnabled(
            ready and worker_idle and self._sga_candidate is not None
            and self._sga_candidate.core_shape() is not None)
        self.iso_prepare.setEnabled(
            active and worker_idle and self._project_target is not None)
        self.iso_mask_toggle.setEnabled(self._contamination is not None)
        self.iso_show_sga_shreds.setEnabled(
            self._isophote_proposal is not None
            and bool(self._sga_shred_paths))
        self.iso_generate.setEnabled(
            active and worker_idle and self._project_target is not None
            and self._contamination is not None
            and self._draft.digest == self._contamination["draft_digest"])
        self.iso_pick_center.setEnabled(active and worker_idle)
        self.iso_reset_center.setEnabled(active and worker_idle)
        confidence_ok = (
            self._isophote_proposal is not None
            and (self.iso_low_conf_override.isChecked()
                 or (self.iso_shape.currentData() != "sga_seed"
                     and not self._isophote_proposal["sensitivity"][
                         "confidence"].startswith("LOW_CONFIDENCE"))))
        self.iso_accept.setEnabled(
            active and worker_idle and confidence_ok
            and self._isophote_preview is not None)
        self.iso_discard.setEnabled(self._isophote_proposal is not None)
        self.iso_shape.setEnabled(self._isophote_proposal is not None)
        self.iso_candidate.setEnabled(
            active and self._isophote_proposal is not None
            and self.iso_candidate.count() > 1)
        self._update_hard_overlay()

    def shutdown(self):
        self.flush_autosave()
        self._load_generation += 1
        self._poll_timer.stop()
        self._executor.shutdown(wait=False, cancel_futures=True)

    def closeEvent(self, event):
        self.shutdown()
        super().closeEvent(event)
