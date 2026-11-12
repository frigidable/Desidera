"""Conservative, read-only positional linking of independent band detections.

An image detection is not a physical object. In particular, one source may
produce multiple Tractor entries and a blend may contain multiple objects.
This module does not classify sources, change masks, or infer photometry.

Stage SA2a: downstream callers must provide detections extracted from the
*original* band images with explicit provenance. Duplicating rows of the
single griz Tractor catalogue is not independent multiband evidence.
"""
from __future__ import annotations

import hashlib
import json
import math
import re

import numpy as np
from scipy.spatial import cKDTree


ASSOCIATION_SCHEMA = "desidera.scene-positional-association/2"
_ALLOWED_BANDS = frozenset(("g", "r", "i", "z"))


def _xyz(ra_deg, dec_deg):
    ra, dec = math.radians(ra_deg), math.radians(dec_deg)
    return (math.cos(dec)*math.cos(ra),
            math.cos(dec)*math.sin(ra), math.sin(dec))


def _offset_arcsec(a, b):
    """Great-circle angular distance, robust near RA wrap and the poles."""
    x, y = np.asarray(a["xyz"]), np.asarray(b["xyz"])
    chord = min(2., float(np.linalg.norm(x-y)))
    return math.degrees(2.*math.asin(chord/2.))*3600.


def associate_band_detections(band_data, *, match_radius_arcsec,
                              max_detections=12000, max_edges=100000):
    """Group positional *candidates*, never declare astrophysical identity.

    band_data: {band: {dataset_id, image_sha256, segmentation_id,
                       detections: [{detection_id, ra_deg, dec_deg}, ...]}}.

    Every detection is retained. Cross-band pairs inside the explicit
    angular radius form an undirected candidate graph. A connected component
    containing repeated bands or non-transitive pairwise separations is
    AMBIGUOUS and must never be silently merged into a physical source.
    Independent Tractor model fluxes/identities are *not* constructed here.
    """
    if not isinstance(band_data, dict) or not band_data:
        raise ValueError("Nonempty band-to-native-detections mapping required")
    if (isinstance(match_radius_arcsec, bool) or
            not isinstance(match_radius_arcsec, (int, float)) or
            not math.isfinite(match_radius_arcsec) or
            not .05 <= match_radius_arcsec <= 10.):
        raise ValueError("Explicit angular match radius must be 0.05..10 arcsec")
    if (type(max_detections) is not int or max_detections < 1 or
            type(max_edges) is not int or max_edges < 1):
        raise ValueError("Positive association safety limits required")
    records = []
    provenance = {}
    seen_native_images = {}
    for band in sorted(band_data):
        if band not in _ALLOWED_BANDS:
            raise ValueError(f"Unsupported optical band: {band!r}")
        block = band_data[band]
        if not isinstance(block, dict):
            raise ValueError("Band evidence must be a mapping")
        meta = {}
        for key in ("dataset_id", "image_sha256", "segmentation_id"):
            value = block.get(key)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"Missing {key} provenance for band {band}")
            meta[key] = value
        if not meta["dataset_id"].endswith("/" + band):
            raise ValueError(
                "Band provenance dataset_id disagrees with detection band")
        image_sha = meta["image_sha256"]
        if re.fullmatch(r"[0-9a-f]{64}", image_sha) is None:
            raise ValueError(
                "image_sha256 must be a complete lowercase SHA-256")
        if image_sha in seen_native_images:
            raise ValueError(
                "Identical original image SHA across bands cannot constitute "
                "independent image-detection evidence")
        seen_native_images[image_sha] = band
        provenance[band] = meta
        detections = block.get("detections")
        if not isinstance(detections, (list, tuple)):
            raise ValueError("Native band detections must be a sequence")
        seen = set()
        for detection in detections:
            if not isinstance(detection, dict):
                raise ValueError("Each image detection must be a mapping")
            name = detection.get("detection_id")
            if not isinstance(name, str) or not name or len(name) > 256:
                raise ValueError("Stable band-local detection_id required")
            if name in seen:
                raise ValueError("Duplicate detection_id within one band")
            seen.add(name)
            try:
                ra, dec = float(detection["ra_deg"]), float(detection["dec_deg"])
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError("Detection needs finite ICRS RA/Dec") from exc
            if not (math.isfinite(ra) and math.isfinite(dec)
                    and -90. <= dec <= 90.):
                raise ValueError("Invalid ICRS position")
            ra = ra % 360.
            records.append({
                "band":band, "detection_id":name,
                "ra_deg":ra, "dec_deg":dec,
                "xyz":_xyz(ra,dec),
            })
            if len(records) > max_detections:
                raise ValueError("Image-detection count exceeds safety cap")
    # Stable ordering, independent of detection/GUI iteration order.
    records.sort(key=lambda r:(r["band"],r["detection_id"]))
    n = len(records)
    adjacency = [set() for _ in range(n)]
    edges = []
    # Query *only* cross-band neighbours. The previous one-tree scan
    # visited O(N²) near-identical same-band detections even though
    # every same-band pair was then discarded. Its edge cap could not
    # constrain that work. Build trees for the larger side of each
    # nonempty band pair; query only the smaller side, fail before
    # storing excessive actual cross-band edges.
    indexed_bands = {
        band: [i for i, r in enumerate(records) if r["band"] == band]
        for band in sorted(band_data)
    }
    if sum(bool(rows) for rows in indexed_bands.values()) >= 2:
        vectors = np.asarray([r["xyz"] for r in records], dtype=float)
        radius = 2.*math.sin(math.radians(match_radius_arcsec/3600.)/2.)
        bands = tuple(sorted(indexed_bands))
        for position, left_band in enumerate(bands):
            for right_band in bands[position+1:]:
                left, right = indexed_bands[left_band], indexed_bands[right_band]
                if not left or not right:
                    continue
                query_indices, tree_indices = (
                    (left, right) if len(left) <= len(right)
                    else (right, left)
                )
                tree = cKDTree(vectors[tree_indices])
                for index in query_indices:
                    for match in sorted(tree.query_ball_point(
                            vectors[index], r=radius)):
                        other = tree_indices[match]
                        a, b = sorted((index, other))
                        if len(edges) >= max_edges:
                            raise ValueError("Too many positional candidate edges")
                        adjacency[a].add(b)
                        adjacency[b].add(a)
                        edges.append((a, b))
    visited = set()
    groups = []
    for seed in range(n):
        if seed in visited:
            continue
        stack = [seed]
        visited.add(seed)
        indices = []
        while stack:
            node = stack.pop()
            indices.append(node)
            for other in sorted(adjacency[node]):
                if other not in visited:
                    visited.add(other)
                    stack.append(other)
        indices.sort()
        members = [records[i] for i in indices]
        bands = [r["band"] for r in members]
        if len(set(bands)) != len(bands):
            # A one-to-many component may contain thousands of entries;
            # do not compute an unbounded pairwise distance matrix.
            max_sep = None
            status = "AMBIGUOUS_MULTIPLE_DETECTIONS_IN_BAND"
        else:
            max_sep = max((_offset_arcsec(records[a],records[b])
                           for j,a in enumerate(indices)
                           for b in indices[j+1:]), default=0.)
            if len(indices) == 1:
                status = "SINGLE_BAND_DETECTION"
            elif max_sep > match_radius_arcsec + 1e-8:
                status = "AMBIGUOUS_NONTRANSITIVE_OFFSETS"
            else:
                status = "POSITIONAL_CANDIDATE_ONLY"
        # Stable GROUP identity must bind the source bytes, detector
        # provenance, sky positions and the matching rule. Reusing
        # detection_id after a new segmentation/image cannot reattach
        # an old group's downstream notes or science decisions.
        identity = {
            "schema": ASSOCIATION_SCHEMA,
            "match_radius_arcsec": float(match_radius_arcsec),
            "members": [
                {
                    "band": r["band"],
                    "detection_id": r["detection_id"],
                    "ra_deg": r["ra_deg"],
                    "dec_deg": r["dec_deg"],
                    **provenance[r["band"]],
                } for r in members
            ],
        }
        digest = hashlib.sha256(json.dumps(
            identity, sort_keys=True, separators=(",", ":"),
            allow_nan=False).encode("utf-8")).hexdigest()[:24]
        groups.append({
            "group_id":f"pos-{digest}",
            "status":status,
            "members":tuple({key:r[key] for key in (
                "band", "detection_id", "ra_deg", "dec_deg")}
                for r in members),
            "bands":tuple(sorted(set(bands))),
            "max_pairwise_offset_arcsec":max_sep,
            "physical_membership":"UNDETERMINED",
            "photometry_action":"NO_NEW_ACTION_READ_ONLY",
        })
    return {
        "schema":ASSOCIATION_SCHEMA,
        "match_radius_arcsec":float(match_radius_arcsec),
        "band_provenance":provenance,
        "input_detection_count":n,
        "candidate_edge_count":len(edges),
        "candidate_neighbour_search":"CROSS_BAND_ONLY_KDTREE_SMALLER_SIDE_V2",
        "same_band_neighbours_never_queried":True,
        "group_id_policy":"INPUT_IMAGE_SEGMENTATION_POSITION_AND_RADIUS_SHA256_V2",
        "group_id_is_not_a_physical_source_identity":True,
        "groups":tuple(groups),
        "status":"UNVALIDATED_POSITIONAL_REVIEW_ONLY",
        "cautions":(
            "IMAGE_DETECTION_NOT_A_PHYSICAL_OBJECT",
            "POSITIONAL_MATCH_NOT_MEMBERSHIP_OR_ASSOCIATION_PROOF",
            "NO_TRACTOR_GRIZ_DUPLICATE_AS_INDEPENDENT_EVIDENCE",
            "IMAGE_SHA256_DECLARED_BY_CALLER_NOT_RECOMPUTED_IN_ASSOCIATION",
            "NO_AUTOMATIC_MASK_OR_PHOTOMETRY_ACTION",
            "REQUIRES_INDEPENDENT_PSF_ASTROMETRY_AND_REAL_DATA_VALIDATION",
        ),
    }
