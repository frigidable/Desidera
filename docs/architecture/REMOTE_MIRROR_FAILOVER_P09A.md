# DESIDERA P0.9a — controlled astronomy-service mirror failover

**Scope:** SIMBAD mirror failover for manually initiated object identification. FITS archive substitution is *not* implemented or implied.

## Verified authoritative service endpoints

SIMBAD is served by two explicitly supported `astroquery.simbad` mirrors:

- Strasbourg: `simbad.cds.unistra.fr`
- Harvard: `simbad.harvard.edu`

References:

- https://docs.astropy.org/projects/astroquery/en/latest/api/astroquery.simbad.Conf.html
- https://simbad.cfa.harvard.edu/

The code uses the maintained `Simbad.server` API, never alters process-global `astroquery` configuration, and never permits arbitrary hostnames or HTTP transport.

## UI and failure semantics

Project Hub → **Galaxy Knowledge · SIMBAD…** → **SIMBAD server**:

- **Automatic · Strasbourg → Harvard** (default).
- **Automatic · Harvard → Strasbourg** (reverse priority).
- **Strasbourg only**.
- **Harvard only**.

A search is manual, runs in a worker, and attempts at most two approved endpoints. A network timeout, connection failure, HTTP 403/404/408/425/429/5xx or a recognised service transport error triggers trying the other approved endpoint. A valid empty result is a final result, not a reason to query a second site. Invalid user input, impossible coordinates, a malformed query or an unexpected software error does **not** silently produce a fallback result. If both servers fail, the GUI displays a clear failure and no new candidate snapshot is saved. Older immutable snapshots and manually confirmed identities remain intact.

Each successful immutable snapshot records `endpoint_evidence` with the actually used hostname, policy, attempted hosts, error *types* (not full potentially sensitive exception messages), and whether failover occurred. The dialogue displays this provenance, including from a saved offline snapshot. It continues to require explicit human confirmation of catalogue identity; no automatic adoption of coordinates or other science quantities.

## Distinction from archive data access

DESI Legacy original coadd products, CCD catalogues and calibration evidence must never silently switch to a different data product. The official NERSC archive is documented at:

- https://www.legacysurvey.org/dr9/files/
- https://www.legacysurvey.org/dr11/files/

NOIRLab Data Lab, Legacy Survey viewer cutout services and raw-image archive provide **other access methods and/or different products**, not independently verified byte-equivalent mirrors of the `legacysurvey-<brick>-image-<band>.fits.fz` original path. They must remain separate, explicitly labelled sources until an original-format, complete, release-specific mirror is verified and integrity checked.

A future archive failover may use only an explicit source registry with release/region/brick/band/product-kind identity, exact relative paths, FITS metadata checks, persisted download URL and SHA-256 of each downloaded file, atomically checked original data, and per-source failure receipts. HTTP 404 for missing survey coverage must not be turned into a fabricated image from a different release or a viewer crop. Likewise NED's new API is a distinct data provider, not a SIMBAD mirror.

## Windows acceptance

From a Git worktree based on this feature branch and Python 3.12:

```powershell
& .\.venv\Scripts\python.exe -m pip install -e ".[exact,gui,knowledge,test]" pytest-qt
& .\.venv\Scripts\python.exe -m pytest -q tests/test_remote_failover.py tests/test_galaxy_knowledge.py
& .\.venv\Scripts\python.exe -m pytest -q tests/gui/test_simbad_mirror_gui.py tests/gui/test_galaxy_knowledge_gui.py
& .\.venv\Scripts\python.exe -m lsbpipe.gui
```

Open the project's Galaxy Knowledge dialog, start with **Automatic · Strasbourg → Harvard**, then try the one-server modes and compare behavior. Verify the displayed actual host, the ability to re-open the same snapshot offline, and that local FITS and project coordinates are unchanged. A failure reaching both servers is an accepted diagnostic outcome, not a completed scientific identity result.

## Limitations

Network availability and mirror freshness are not demonstrated by offline CI. MIRROR FALLBACK IS NOT A SCIENTIFIC MEASUREMENT; records are still labelled `CANDIDATES_UNVERIFIED` / `NOT_SCIENTIFICALLY_ADOPTED`. A mirror result might differ in update latency. A valid query response from one mirror is not silently reconciled against the other, and a result from one is never claimed to validate the other.

The feature is intentionally isolated from the science dataset manager and existing original-FITS downloading.
