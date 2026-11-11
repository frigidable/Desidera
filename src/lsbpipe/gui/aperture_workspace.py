"""Stage 04.2 read-only project-local finite-aperture diagnostics.

Reuses *currently loaded* Stage 02 native pixels and prepared Scene draft.
No FITS reads/downloads on opening this tab, no saved galaxy photometry,
no sky adoption, no inferred total magnitude or interpolated masked flux.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import json
import math

import numpy as np
from PySide6.QtCore import Qt, QTimer, QRectF
from PySide6.QtGui import QColor, QPainter, QPen
from PySide6.QtWidgets import (
    QApplication, QGroupBox, QHBoxLayout, QLabel, QPlainTextEdit,
    QPushButton, QSizePolicy, QVBoxLayout, QWidget,
)

from lsbpipe.photometry.aperture_flux import (
    ApertureProvenance, ApertureSkyScenario,
)
from lsbpipe.photometry.aperture_growth import (
    geometry_hypotheses, predeclared_growth_radii,
    measure_growth_diagnostic,
)
from lsbpipe.photometry.contracts import NativeBandFrame
from lsbpipe.photometry.sky_masks import evaluate_mask_draft


def build_project_aperture_growth(
        dataset, draft, contamination, *, source_triplet,
        target_icrs, sga_q=None, sga_pa_deg=None):
    """Explicit Stage 04 worker on the one already loaded Scene source.

    The source triplet comes from the bound project manifest, not from
    parsing its visual label. Does not acquire another brick or read FITS.
    """
    if (dataset.get("status") != "NATIVE_SCENE_CUTOUT_NOT_SCIENCE_ADOPTED"
            or contamination.get("status") !=
            "CATALOGUE_CONTAMINATION_MASK_NOT_ADOPTED"
            or dataset["project_id"] != draft.project_id
            or contamination.get("project_id") != dataset["project_id"]
            or contamination.get("band") != dataset["band"]
            or contamination.get("dataset_source") != dataset["source"]
            or contamination.get("draft_digest") != draft.digest):
        raise ValueError("Stage 04 refuses stale project/Scene evidence")
    if len(source_triplet) != 3:
        raise ValueError("Source identity must be release/region/brick")
    release,region,brick = source_triplet
    if dataset["band"] not in "griz":
        raise ValueError("Unknown native optical band")
    data=dataset["data"]
    f=NativeBandFrame(
        image_nMgy=data["image"],
        invvar_nMgy2=data["invvar"],
        nexp=data["nexp"],
        maskbits=data["maskbits"],
        native_wcs=dataset["wcs"],release=release,region=region,
        brick=brick,band=dataset["band"],
        source_evidence={
            "scene_dataset_source":dataset["source"],
            "project_id":dataset["project_id"],
            "data_status":"SCENE_ORIGINAL_NATIVE_CUTOUT_NOT_FILE_HASH_VERIFIED",
            "local_hashes":{},
        })
    if f.shape != tuple(dataset["pixel_shape"]):
        raise ValueError("Original Scene native image grid changed")
    hard=np.asarray(dataset["native_hard"])
    if not np.array_equal(~f.valid,hard):
        raise ValueError(
            "Native quartet validity disagrees with previously loaded Scene")
    preview=evaluate_mask_draft(
        draft,dataset["wcs"],f.shape,dataset["band"],
        native_hard=hard,
        automatic_contaminants=contamination["auto_mask"])
    # Crucial: sky target-protection and provisional isophote layers
    # MUST NOT enter the photometric invalidity masks.
    source=np.array(preview["effective_auto"],copy=True)
    manual=np.array(preview["manual_exclude"],copy=True)
    if preview["draft_sha256"] != draft.digest:
        raise ValueError("Scene mask revision changed during photometry")
    provenance=ApertureProvenance(
        project_id=dataset["project_id"],
        scene_draft_sha256=draft.digest,
        catalogue_sha256=contamination["catalog_sha256"],
        source_policy_id="PROVISIONAL_SCENE_EFFECTIVE_AUTO_PLUS_EXCLUDE_V1")
    protected=int(np.count_nonzero(contamination["protected"]))
    radii,radial_qa=predeclared_growth_radii(
        f,target_icrs,
        target_protection_native_pixels=(
            protected if protected >= 25 else None))
    hypotheses=geometry_hypotheses(
        sga_q=sga_q,sga_pa_deg=sga_pa_deg)
    curves=[]
    for hypo in hypotheses:
        curves.append(measure_growth_diagnostic(
            f,anchor_icrs=target_icrs,provenance=provenance,
            effective_source_mask=source,manual_exclude=manual,
            radii_arcsec=radii,hypothesis=hypo,
            sky=ApertureSkyScenario()))
    result={
        "schema":"desidera.photometry.project-aperture-growth/1",
        "status":"DIRECT_FINITE_APERTURE_DIAGNOSTICS_NOT_ADOPTED",
        "project_id":dataset["project_id"],
        "dataset_source":dataset["source"],
        "dataset_id":f.dataset_id,
        "grid_fingerprint":curves[0]["grid_fingerprint"],
        "scene_draft_sha256":draft.digest,
        "catalog_sha256":contamination["catalog_sha256"],
        "source_policy_id":provenance.source_policy_id,
        "band":dataset["band"],
        "target_anchor_icrs":list(map(float,target_icrs)),
        "selected_sga_hypothesis_available":(
            sga_q is not None and sga_pa_deg is not None),
        "radial_policy":radial_qa,
        "n_provisional_catalogue_source_pixels":
            int(np.count_nonzero(contamination["auto_mask"])),
        "n_effective_source_exclusion_pixels":
            int(np.count_nonzero(source)),
        "n_explicit_manual_exclusion_pixels":
            int(np.count_nonzero(manual)),
        "n_sky_target_protection_pixels_not_masked_for_photometry":
            protected,
        "sky_method":"zero_additional_correction",
        "sky_stage_science_adopted":False,
        "uses_only_existing_loaded_native_scene":True,
        "source_fits_hashes_verified":False,
        "original_fits_modified":False,
        "background_subtracted_fits_created":False,
        "masked_light_reconstructed":False,
        "total_galaxy_magnitude_inferred":False,
        "noise_covariance_validated":False,
        "source_membership_validated":False,
        "psf_science_validated":False,
        "adopted_for_science":False,
        "curves":curves,
        "warnings":[
            "FIRST_STAGE04_DIAGNOSTIC_CURVE_NOT_FINAL_GALAXY_PHOTOMETRY",
            "PROJECT_ANCHOR_FIXED_NOT_A_FITTED_PHOTOMETRIC_CENTRE",
            "PROVISIONAL_SOURCE_MEMBERSHIP_CAN_REMOVE_TARGET_KNOTS",
            "NATIVE_SCENE_CUTOUT_SOURCE_FILE_HASHES_UNVERIFIED",
            "RESIDUAL_SKY_SYSTEMATICS_NOT_CALIBRATED",
            "COADD_NOISE_COVARIANCE_NOT_CALIBRATED",
            "EMPIRICAL_PSF_WINGS_NOT_SCIENTIFICALLY_VALIDATED",
        ],
    }
    json.dumps(result,ensure_ascii=False,allow_nan=False)
    return result


class ApertureCurvePlot(QWidget):
    """Simple read-only native-curve preview; no fitting or smoothing."""

    def __init__(self,parent=None, *, mode="cumulative"):
        super().__init__(parent)
        if mode not in ("cumulative","annular"):
            raise ValueError("Unknown signed aperture plot mode")
        self._mode=mode
        self.setObjectName("aperture_stage_"+mode+"_plot")
        self.setMinimumHeight(175)
        self._curves=()
        self.setSizePolicy(QSizePolicy.Ignored,QSizePolicy.Fixed)

    def set_curves(self,curves):
        self._curves=tuple(curves)
        self.update()

    def _entries(self,curve):
        if self._mode=="annular":
            data=zip(
                curve["radii_arcsec"][1:],
                curve["independently_integrated_annuli"])
            key="signed_mean_intensity_nMgy_per_arcsec2"
        else:
            data=zip(curve["radii_arcsec"],curve["cumulative"])
            key="signed_valid_flux_nMgy"
        return [
            (float(rad),float(row[key]))
            for rad,row in data if row[key] is not None
        ]

    def paintEvent(self,event):
        super().paintEvent(event)
        p=QPainter(self)
        p.setRenderHint(QPainter.Antialiasing,True)
        p.fillRect(self.rect(),self.palette().base())
        w,h=self.width(),self.height()
        left,right,top,bottom=66.,15.,16.,35.
        if w<110 or h<110:
            return
        bounds=QRectF(left,top,w-left-right,h-top-bottom)
        p.setPen(QPen(self.palette().mid().color(),1.))
        p.drawRect(bounds)
        points=[]
        for curve in self._curves:
            data=self._entries(curve)
            points.extend(data)
        if not points:
            p.setPen(self.palette().text().color())
            p.drawText(bounds,Qt.AlignCenter,
                       "No direct finite-aperture curve yet")
            return
        xmax=max(x for x,y in points)
        ymin=min(min(y for x,y in points),0.)
        ymax=max(max(y for x,y in points),0.)
        if ymax<=ymin:
            ymax=ymin+1.
        gap=(ymax-ymin)*.08
        ymin-=gap
        ymax+=gap
        colors=(QColor("#527cdd"),QColor("#cf8344"),QColor("#48a68c"))
        def pos(x,y):
            return (left+bounds.width()*x/xmax,
                    top+bounds.height()*(1.-(y-ymin)/(ymax-ymin)))
        if ymin<0.<ymax:
            x0,y0=pos(0.,0.)
            p.setPen(QPen(self.palette().mid().color(),1.,Qt.DashLine))
            p.drawLine(int(left),int(y0),int(left+bounds.width()),int(y0))
        for index,curve in enumerate(self._curves):
            entries=self._entries(curve)
            col=colors[index%len(colors)]
            p.setPen(QPen(col,2.))
            for i,(rad,flux) in enumerate(entries):
                x,y=pos(rad,flux)
                p.drawEllipse(QRectF(x-2.5,y-2.5,5.,5.))
                if i:
                    x0,y0=pos(*entries[i-1])
                    p.drawLine(int(x0),int(y0),int(x),int(y))
        p.setPen(self.palette().text().color())
        p.drawText(5,18,(
            "nMgy" if self._mode=="cumulative" else "nMgy/arcsec²"))
        p.drawText(3,int(top+12),f"{ymax:.3g}")
        p.drawText(3,int(top+bounds.height()),f"{ymin:.3g}")
        p.drawText(int(left+bounds.width()-72),h-8,
                   f"a ≤ {xmax:.1f}″")


class ApertureWorkspace(QWidget):
    """Explicit Stage 04 diagnostic UI connected to current Scene source."""

    def __init__(self,scene_editor,sky_stage,parent=None):
        super().__init__(parent)
        self.setObjectName("photometry_aperture_stage")
        self.setMinimumWidth(0)
        self.setSizePolicy(QSizePolicy.Ignored,QSizePolicy.Expanding)
        self._scene=scene_editor
        self._sky=sky_stage
        self._executor=ThreadPoolExecutor(
            max_workers=1,thread_name_prefix="desidera-aperture")
        self._future=None
        self._active_identity=None
        self._result_identity=None
        self._result=None
        self._timer=QTimer(self)
        self._timer.setInterval(100)
        self._timer.timeout.connect(self._poll_worker)
        root=QVBoxLayout(self)
        root.setContentsMargins(15,13,15,13)
        root.setSpacing(8)
        title=QLabel("04 · Object aperture photometry")
        title.setObjectName("aperture_stage_title")
        title.setStyleSheet("font-size: 19px; font-weight: 600;")
        root.addWidget(title)
        note=QLabel(
            "Direct finite-aperture flux only. No total galaxy magnitude, "
            "source filling, sky adoption or calibrated noise covariance. "
            "Original DR11 sky-subtracted pixels are retained.")
        note.setWordWrap(True)
        root.addWidget(note)
        setup=QGroupBox("Current native Scene → automatic fixed ellipses")
        box=QVBoxLayout(setup)
        self.context_label=QLabel(
            "In Stage 02 load the original FITS and prepare neighbour masks.")
        self.context_label.setObjectName("aperture_stage_context")
        self.context_label.setWordWrap(True)
        box.addWidget(self.context_label)
        root.addWidget(setup)
        buttons=QHBoxLayout()
        self.run_btn=QPushButton("Build aperture curves")
        self.run_btn.setObjectName("aperture_stage_run")
        self.run_btn.setToolTip(
            "Use project-anchored circles and explicit provisional ellipses "
            "on the current native FITS scene. Worker never reads a new "
            "brick or writes an image. No manual aperture drawing required.")
        self.run_btn.clicked.connect(self.run_diagnostic)
        buttons.addWidget(self.run_btn)
        self.copy_btn=QPushButton("Copy aperture QA")
        self.copy_btn.setObjectName("aperture_stage_copy_qa")
        self.copy_btn.setEnabled(False)
        self.copy_btn.clicked.connect(self.copy_qa)
        buttons.addWidget(self.copy_btn)
        buttons.addStretch(1)
        root.addLayout(buttons)
        self.status=QLabel("Not measured · science HOLD.")
        self.status.setObjectName("aperture_stage_status")
        self.status.setWordWrap(True)
        root.addWidget(self.status)
        plots=QHBoxLayout()
        self.plot=ApertureCurvePlot(self,mode="cumulative")
        plots.addWidget(self.plot,1)
        self.annulus_plot=ApertureCurvePlot(self,mode="annular")
        plots.addWidget(self.annulus_plot,1)
        root.addLayout(plots)
        self.readout=QPlainTextEdit(self)
        self.readout.setObjectName("aperture_stage_readout")
        self.readout.setMinimumWidth(0)
        self.readout.setSizePolicy(QSizePolicy.Ignored,QSizePolicy.Expanding)
        self.readout.setReadOnly(True)
        self.readout.setPlainText(
            "No aperture measurement. Open Stage 02, prepare native "
            "source masks, then explicitly build diagnostic curves.")
        root.addWidget(self.readout,1)
        self.refresh_context()

    def _snapshot(self):
        # Source and drafts are owned by Stage 02. Stage 03's project
        # evidence guard is reused, but its sky fit is never adopted.
        snap=self._sky._snapshot()
        scene=self._scene
        if (snap is None or scene._project_target is None
                or scene._source_key is None or
                len(scene._source_key)<3 or
                scene._source_key[2] is None):
            return None
        source_tuple=tuple(scene._source_key[2][:3])
        if len(source_tuple)!=3:
            return None
        prior=scene._sga_candidate if (
            scene.sga_use_prior.isChecked()) else None
        q=getattr(prior,"ba",None) if prior is not None else None
        pa=getattr(prior,"pa_deg",None) if prior is not None else None
        if q is None or pa is None:
            q,pa=None,None
        identity=snap[0]+(
            tuple(source_tuple),tuple(scene._project_target),
            (q,pa))
        return (identity,snap[1],snap[2],snap[3],
                source_tuple,tuple(scene._project_target),q,pa)

    def refresh_context(self):
        snap=self._snapshot()
        if snap is None:
            self.context_label.setText(
                "No current source-bound Scene. Load original native FITS "
                "and prepare provisional neighbour masks in Stage 02.")
        else:
            self.context_label.setText(
                f"{snap[1]['source']} · project {snap[1]['project_id']} · "
                f"draft {snap[2].digest[:12]} · "
                "same project anchor; no sky subtraction.")
        if self._result is not None and (
                snap is None or self._result_identity!=snap[0]):
            self._result=None
            self._result_identity=None
            self.copy_btn.setEnabled(False)
            self.plot.set_curves(())
            self.annulus_plot.set_curves(())
            self.readout.setPlainText(
                "Previous aperture curves invalidated by project, band, "
                "Scene draft, source or geometric hypothesis change.")
            self.status.setText("Stale aperture result · run again.")
        self.run_btn.setEnabled(snap is not None and self._future is None)
        return snap is not None

    def run_diagnostic(self):
        if self._future is not None:
            return
        snap=self._snapshot()
        if snap is None:
            self.refresh_context()
            self.status.setText("Stage 02 original FITS/masks required.")
            return
        identity,dataset,draft,contamination,triplet,anchor,q,pa=snap
        self._active_identity=identity
        self._result=None
        self._result_identity=None
        self.copy_btn.setEnabled(False)
        self.plot.set_curves(())
        self.annulus_plot.set_curves(())
        self._future=self._executor.submit(
            build_project_aperture_growth,
            dataset,draft,contamination,
            source_triplet=triplet,target_icrs=anchor,
            sga_q=q,sga_pa_deg=pa)
        self.run_btn.setEnabled(False)
        self.status.setText(
            "Integrating exact observed native-pixel ellipses and "
            "independent rings in background…")
        self._timer.start()

    def _poll_worker(self):
        future=self._future
        if future is None or not future.done():
            return
        self._timer.stop()
        self._future=None
        identity=self._active_identity
        self._active_identity=None
        snap=self._snapshot()
        if snap is None or snap[0]!=identity:
            self.status.setText(
                "Aperture result discarded: Scene evidence changed "
                "during calculation.")
            self.refresh_context()
            return
        try:
            result=future.result()
        except (OSError,ValueError,KeyError,TypeError,ArithmeticError,
                RuntimeError,ImportError) as exc:
            self.status.setText(f"Aperture diagnostic refused: {exc}")
            self.readout.setPlainText(
                "No adopted photometry. Review original Scene, "
                "native FITS/WCS and exact Photutils dependency.")
            self.refresh_context()
            return
        if (result.get("project_id")!=identity[0]
                or result.get("dataset_source")!=identity[1]
                or result.get("band")!=identity[2]
                or result.get("scene_draft_sha256")!=identity[3]
                or result.get("catalog_sha256")!=identity[4]):
            self.status.setText(
                "Aperture report provenance mismatch — discarded.")
            self.refresh_context()
            return
        self._result=result
        self._result_identity=identity
        self.copy_btn.setEnabled(True)
        self.plot.set_curves(result["curves"])
        self.annulus_plot.set_curves(result["curves"])
        self.readout.setPlainText(self._format_readout(result))
        self.status.setText(
            f"{len(result['curves'])} fixed-geometry curves on "
            f"{result['band'].upper()} · direct valid-area flux · "
            "UNADOPTED. Original FITS unchanged.")
        self.refresh_context()

    @staticmethod
    def _format_readout(result):
        lines=[
            "04 · Native signed aperture photometry — UNADOPTED",
            f"Source: {result['dataset_source']}",
            f"Project anchor (ICRS): {result['target_anchor_icrs']}",
            f"Radial policy: {result['radial_policy']['strategy']}",
            f"Radius basis: {result['radial_policy']['basis']}",
            "Sky: ZERO additional correction to original DR11 coadd",
            "All masked galaxy light remains UNRECOVERED.",
            "Diagonal INVVAR errors omit coadd covariance.",
            "",
        ]
        for curve in result["curves"]:
            h=curve["geometry_hypothesis"]
            lines.extend([
                f"Geometry: {h['name']} · q={h['q']:.3g}, "
                f"PA={h['pa_deg']:.1f}° · {h['source']}",
                "a [arcsec] | valid flux [nMgy] | coverage | "
                "formal σ [nMgy]",
            ])
            for radius,row in zip(curve["radii_arcsec"],curve["cumulative"]):
                fl=row["signed_valid_flux_nMgy"]
                sn=row["formal_diagonal_sigma_nMgy"]
                cov=row["area_native_pixel2"]["valid_fraction_of_ideal"]
                lines.append(
                    f"{radius:8.2f} | "
                    f"{fl if fl is not None else '—':>15} | "
                    f"{cov:.3f} | {sn if sn is not None else '—'}")
            lines += [
                "Direct annuli: a_in–a_out [arcsec] | signed intensity "
                "[nMgy/arcsec²] | AB μ if positive",
            ])
            for ring in curve["independently_integrated_annuli"]:
                inner=ring["a_inner_arcsec"]
                outer=ring["a_outer_arcsec"]
                intensity=ring["signed_mean_intensity_nMgy_per_arcsec2"]
                mu=ring["ab_surface_brightness_mag_arcsec2"]
                lines.append(
                    f"{inner:6.2f}–{outer:6.2f} | "
                    f"{intensity if intensity is not None else '—'} | "
                    f"{mu if mu is not None else 'undefined'}")
            lines.extend([
                "Independent annulus vs CoG closure [nMgy]: "
                + str(curve["max_abs_annulus_closure_nMgy"]),
                "Outermost valid-area AB magnitude is conditional "
                "on an observed *partial* aperture, NOT total light.",
                "",
            ]
        lines += [
            "SCIENCE HOLD: no missing-area interpolation, no PSF-corrected "
            "total, no calibrated sky/noise covariance.",
            "Copy aperture QA contains complete per-ring signed intensities "
            "and mask provenance.",
        ]
        return "\n".join(lines)

    def copy_qa(self):
        snap=self._snapshot()
        if (self._result is None or snap is None or
                snap[0]!=self._result_identity):
            self.refresh_context()
            self.status.setText("No current aperture QA to copy.")
            return
        payload={
            **self._result,
            "export_generated_utc":datetime.now(timezone.utc).isoformat(),
        }
        QApplication.clipboard().setText(json.dumps(
            payload,ensure_ascii=False,indent=2,allow_nan=False))
        self.status.setText(
            "Copied current Stage 04 source-bound QA. "
            "Original image and Scene masks unchanged.")

    def shutdown(self):
        self._timer.stop()
        self._executor.shutdown(wait=False,cancel_futures=True)

    def closeEvent(self,event):
        self.shutdown()
        super().closeEvent(event)
