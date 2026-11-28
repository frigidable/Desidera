import csv
import hashlib
import io
import json
import zipfile

import pytest

from lsbpipe.photometry.joint_archive_review import (
    read_joint_archive, save_review_receipt, reopen_review_receipt,
)


CONTEXT=dict(project_id="offline-project",source_triplet=("dr11","north","0965p780"),
             target_icrs=(96.,78.),target_radius_arcsec=60.)


def make_archive(path, *, anchor=(96.,78.), tamper=False, source=None):
    """Explicitly synthetic archive-format fixture, no observational evidence."""
    rows=io.StringIO()
    w=csv.writer(rows)
    w.writerow(["band","a_inner_arcsec","a_outer_arcsec","observed_nMgy",
                "model_same_support_nMgy","conditional_sigma_nMgy","eligible","train","test"])
    for band in ("g","r"):
        w.writerow([band,0,5,2,2.1,.2,True,True,False])
        w.writerow([band,5,10,-.2,.1,.3,True,False,True])
    import numpy as np
    arr=io.BytesIO()
    np.savez_compressed(arr,fixture=np.arange(3))
    fit={"status":"DIAGNOSTIC_NOT_ADOPTED","success":True,"boundary":False,
         "band_order":["g","r"],"record_sha256":["a"*64,"b"*64],
         "covariance_sha256":"c"*64,"covariance_calibrated":False,
         "train_rings":[True,False],"test_rings":[False,True],
         "chi2":2.3,"dof":1,"covariance_source":"SYNTHETIC FIXTURE ONLY"}
    files={"rings.csv":rows.getvalue().encode(),"arrays.npz":arr.getvalue(),
           "fit.json":json.dumps(fit).encode()}
    triplet=source or CONTEXT["source_triplet"]
    manifest={"status":"DIAGNOSTIC_NOT_ADOPTED","schema":"desidera.annular-diagnostic/1",
              "band_order":fit["band_order"],"record_sha256":fit["record_sha256"],
              "covariance_sha256":fit["covariance_sha256"],"covariance_calibrated":False,
              "full_geometry_containment_verified":False,
              "inputs":[{"dataset_id":"/".join((*triplet,b)),"native_input_sha256":"d"*64,
                         "grid_id":"wcs-sha256:"+"e"*64,"anchor_icrs":list(anchor)}
                        for b in ("g","r")],
              "files_sha256":{k:hashlib.sha256(v).hexdigest() for k,v in files.items()}}
    if tamper:
        files["rings.csv"]+=b"tampered"
    with zipfile.ZipFile(path,"w",compression=zipfile.ZIP_DEFLATED) as z:
        for k,v in files.items():
            z.writestr(k,v)
        z.writestr("manifest.json",json.dumps(manifest))
    return path


def test_archive_signed_flux_association_and_receipt(tmp_path):
    path=make_archive(tmp_path/"synthetic.zip")
    review=read_joint_archive(path,**CONTEXT)
    assert review["rows"][1]["observed_nMgy"]==-.2
    assert review["status"]=="DIAGNOSTIC_NOT_ADOPTED"
    assert not review["original_fits_reverified"]
    receipt=tmp_path/"review.json"
    with pytest.raises(ValueError):
        save_review_receipt(receipt,review,confirmed=False)
    assert not receipt.exists()
    save_review_receipt(receipt,review,confirmed=True)
    assert reopen_review_receipt(receipt,**CONTEXT)["archive_sha256"]==review["archive_sha256"]
    with pytest.raises(FileExistsError):
        save_review_receipt(receipt,review,confirmed=True)
    with pytest.raises(ValueError,match="current project"):
        reopen_review_receipt(receipt,**dict(CONTEXT,project_id="different"))
    with path.open("ab") as f:
        f.write(b"modified")
    with pytest.raises(ValueError,match="changed"):
        reopen_review_receipt(receipt,**CONTEXT)


def test_reject_mismatched_source_anchor_and_corruption(tmp_path):
    with pytest.raises(ValueError,match="selected native source"):
        read_joint_archive(make_archive(tmp_path/"wrong.zip",source=("dr9","north","0965p780")),**CONTEXT)
    with pytest.raises(ValueError,match="target radius"):
        read_joint_archive(make_archive(tmp_path/"far.zip",anchor=(96.,79.)),**CONTEXT)
    with pytest.raises(ValueError,match="checksum mismatch"):
        read_joint_archive(make_archive(tmp_path/"tamper.zip",tamper=True),**CONTEXT)
    path=make_archive(tmp_path/"unexpected.zip")
    with zipfile.ZipFile(path,"a") as z:
        z.writestr("../unexpected",b"not extracted")
    with pytest.raises(ValueError,match="Unexpected archive"):
        read_joint_archive(path,**CONTEXT)
    assert not (tmp_path.parent/"unexpected").exists()


def test_joint_review_dialog_offline_states(tmp_path,monkeypatch):
    pytest.importorskip("PySide6")
    monkeypatch.setenv("QT_QPA_PLATFORM","offscreen")
    from PySide6.QtWidgets import QApplication
    from lsbpipe.gui.joint_archive_dialog import JointArchiveDialog
    app=QApplication.instance() or QApplication([])
    dialog=JointArchiveDialog(project_name="Synthetic fixture — not science",
                              project_root=tmp_path,context=CONTEXT)
    try:
        assert not dialog.save_btn.isEnabled()
        assert dialog.rings.rowCount()==0
        path=make_archive(tmp_path/"synthetic.zip")
        dialog.start_load(path)
        assert "Verifying" in dialog.status.text()
        assert not dialog.open_btn.isEnabled()
        dialog._future.result(timeout=10)
        dialog._poll()
        assert dialog.rings.rowCount()==4
        assert dialog.rings.item(1,2).text()=="-0.2"
        assert "Original FITS reverified: no" in dialog.summary.toPlainText()
        assert not dialog.save_btn.isEnabled()
        dialog.confirm.setChecked(True)
        assert dialog.save_btn.isEnabled()
        dialog.resize(640,480)
        dialog.show()
        app.processEvents()
        assert dialog.open_btn.isVisible()
        assert dialog.save_btn.isVisible()
        dialog.start_load(tmp_path/"missing.zip")
        assert dialog.rings.rowCount()==0
        assert not dialog.save_btn.isEnabled()
        with pytest.raises(FileNotFoundError):
            dialog._future.result(timeout=10)
        dialog._poll()
        assert "Review failed" in dialog.status.text()
        assert dialog.open_btn.isEnabled()
        assert not dialog.confirm.isEnabled()
    finally:
        dialog.reject()
        dialog.deleteLater()
        app.processEvents()
