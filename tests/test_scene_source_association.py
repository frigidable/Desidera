"""Stage SA2a: conservative, provenance-bound linking only."""
from copy import deepcopy

import pytest

from lsbpipe.photometry.scene_source_association import (
    associate_band_detections, ASSOCIATION_SCHEMA,
)


def block(band, *detections):
    return {"dataset_id":f"dr11/south/1554p180/{band}",
            "image_sha256":{"g":"a","r":"b","i":"c","z":"d"}[band]*64,
            "segmentation_id":f"native-test-{band}",
            "detections":[{"detection_id":str(label),
                           "ra_deg":ra,"dec_deg":dec}
                          for label, ra, dec in detections]}


def test_ra_wrap_one_to_one_never_infers_membership_or_changes_input():
    data={"g":block("g",("G",359.99995,0.)),
          "r":block("r",("R",.00005,0.))}
    original=deepcopy(data)
    out=associate_band_detections(data,match_radius_arcsec=.5)
    assert data==original
    assert out["schema"]==ASSOCIATION_SCHEMA
    assert out["candidate_edge_count"]==1
    assert len(out["groups"])==1
    g=out["groups"][0]
    assert g["status"]=="POSITIONAL_CANDIDATE_ONLY"
    assert g["physical_membership"]=="UNDETERMINED"
    assert g["photometry_action"]=="NO_NEW_ACTION_READ_ONLY"
    assert g["max_pairwise_offset_arcsec"]==pytest.approx(.36,rel=1e-5)


def test_blend_one_r_to_two_g_sources_is_explicitly_ambiguous():
    data={"g":block("g",("g1",10.,1.),("g2",10.+.4/3600.,1.)),
          "r":block("r",("r1",10.+.2/3600.,1.))}
    out=associate_band_detections(data,match_radius_arcsec=.3)
    g=out["groups"][0]
    assert g["status"]=="AMBIGUOUS_MULTIPLE_DETECTIONS_IN_BAND"
    assert len(g["members"])==3


def test_chain_of_cross_band_offsets_does_not_implicitly_merge():
    data={"g":block("g",("g1",10.,0.)),
          "r":block("r",("r1",10.+.4/3600.,0.)),
          "i":block("i",("i1",10.+.8/3600.,0.))}
    out=associate_band_detections(data,match_radius_arcsec=.5)
    assert len(out["groups"])==1
    assert out["groups"][0]["status"]=="AMBIGUOUS_NONTRANSITIVE_OFFSETS"


def test_single_band_not_invented_multiband_counterpart_and_determinism():
    data={"g":block("g",("a",11.,0.),("b",15.,0.)),
          "z":block("z",("c",11.,0.))}
    out=associate_band_detections(data,match_radius_arcsec=1.)
    reversed_data={"z":data["z"],"g":dict(data["g"],
                                detections=list(reversed(data["g"]["detections"]))) }
    again=associate_band_detections(reversed_data,match_radius_arcsec=1.)
    assert out==again
    assert len(out["groups"])==2
    assert sorted(g["status"] for g in out["groups"])==[
        "POSITIONAL_CANDIDATE_ONLY","SINGLE_BAND_DETECTION"]


def test_provenance_duplicate_identity_and_resource_limits_fail_closed():
    valid={"g":block("g",("g1",10.,0.))}
    with pytest.raises(ValueError, match="safety limits"):
        associate_band_detections(valid,match_radius_arcsec=.3,
                                  max_detections=0)
    with pytest.raises(ValueError,match="Missing image_sha256"):
        associate_band_detections({"g":dict(valid["g"],image_sha256="")},
                                  match_radius_arcsec=.5)
    duplicate={"g":block("g",("a",10.,0.),("a",12.,0.))}
    with pytest.raises(ValueError,match="Duplicate"):
        associate_band_detections(duplicate,match_radius_arcsec=.5)
    close={"g":block("g",("a",10.,0.)),
           "r":block("r",("b",10.,0.))}
    close["i"]=block("i",("c",10.,0.))
    with pytest.raises(ValueError,match="candidate edges"):
        associate_band_detections(close,match_radius_arcsec=.5,max_edges=1)
    with pytest.raises(ValueError,match="Invalid ICRS"):
        associate_band_detections({"g":block("g",("a",10.,91.))},
                                  match_radius_arcsec=.5)


def test_association_rejects_unverifiable_or_duplicate_image_evidence():
    valid = {"g":block("g",("a",10.,0.)),
             "r":block("r",("b",10.,0.))}
    assert associate_band_detections(
        valid,match_radius_arcsec=.3)["candidate_edge_count"] == 1
    for invalid in ("g"*64,"A"*64,"a"*63,"anything","a"*65):
        with pytest.raises(ValueError,match="lowercase SHA-256"):
            associate_band_detections(
                {"g":dict(valid["g"],image_sha256=invalid)},
                match_radius_arcsec=.3)
    with pytest.raises(ValueError,match="Identical original image SHA"):
        associate_band_detections(
            {"g":valid["g"],
             "r":dict(valid["r"],image_sha256=valid["g"]["image_sha256"])},
            match_radius_arcsec=.3)
    with pytest.raises(ValueError,match="dataset_id disagrees"):
        associate_band_detections(
            {"g":dict(valid["g"],dataset_id="dr11/south/1554p180/z")},
            match_radius_arcsec=.3)
    with pytest.raises(ValueError,match="Explicit angular match"):
        associate_band_detections(valid,match_radius_arcsec=True)
    with pytest.raises(ValueError,match="safety limits"):
        associate_band_detections(valid,match_radius_arcsec=.3,
                                  max_detections=True)


def test_group_ids_bind_dataset_image_segmentation_position_and_radius():
    """No old group ID may survive a changed native observation/segmentation."""
    original = {"g":block("g",("a",10.,0.)),
                "r":block("r",("b",10.+.2/3600.,0.))}
    def gid(data, radius=.4):
        report = associate_band_detections(
            data,match_radius_arcsec=radius)
        assert report["candidate_edge_count"] == 1
        assert report["groups"][0]["status"] == "POSITIONAL_CANDIDATE_ONLY"
        return report["groups"][0]["group_id"]
    baseline = gid(original)
    assert baseline.startswith("pos-") and len(baseline) == 28
    assert gid(original) == baseline
    changed_image = deepcopy(original)
    changed_image["g"]["image_sha256"] = "f"*64
    assert gid(changed_image) != baseline
    changed_segmentation = deepcopy(original)
    changed_segmentation["r"]["segmentation_id"] += "-new"
    assert gid(changed_segmentation) != baseline
    changed_position = deepcopy(original)
    changed_position["r"]["detections"][0]["ra_deg"] += .02/3600.
    assert gid(changed_position) != baseline
    assert gid(original, radius=.5) != baseline
    changed_dataset = deepcopy(original)
    changed_dataset["g"]["dataset_id"] = "dr11/south/1555p180/g"
    assert gid(changed_dataset) != baseline
