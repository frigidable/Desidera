# DESIDERA P0.7a — SIMBAD Galaxy Knowledge Hub (identity-only)

**Status:** Draft implementation and QA; NOT a scientific measurement or source of authoritative galaxy positions.

## Windows installation and manual GUI acceptance

From the project repository, using PowerShell 7:

```powershell
cd C:\Users\frigi\Documents\PhD\Desidera
git fetch origin
git switch --track origin/feature/p07a-galaxy-knowledge-simbad
python -m pip install -e ".[gui,knowledge,test]"
python -m pytest -q tests/test_galaxy_knowledge.py tests/gui/test_galaxy_knowledge_gui.py
python -m lsbpipe.gui
```

If the working directory is under a different folder, change only the `cd` path. Do not switch branches with uncommitted local changes: save or stash them first.

In **Galaxy Projects**, select a galaxy with a saved name and (preferably) ICRS coordinates. Click **Galaxy Knowledge · SIMBAD…**. On opening, only locally cached results are read; no network request starts. The **Search SIMBAD…** button performs explicit name lookup and, if project coordinates are known, a cone search (default 60 arcsec, configurable to 300 arcsec). Results are presented separately, with source positions and angular separations.

Select the correct SIMBAD row, verify the proposed association, click **Confirm selected identity…**, and confirm in the dialogue. An immutable confirmation receipt is saved under `knowledge/simbad/`. Project UUID, target ICRS coordinates, project metadata, local original FITS, dataset registry, analysis QA, and photometry are unchanged. Reopening the same object works offline, including the last cached list and the confirmed catalogue identity. A missing SIMBAD counterpart is a normal result.

### Technical contracts

- `src/lsbpipe/surveys/galaxy_knowledge.py`: opt-in astroquery SIMBAD adapter, candidate normalization, angular separation, content-addressed query snapshot and user confirmation.
- `src/lsbpipe/gui/galaxy_knowledge_dialog.py`: nonblocking lookup UI, candidate review, explicit confirmation; no remote work during project startup.
- `src/lsbpipe/gui/project_hub.py`: Galaxy Knowledge button and read-only confirmed-association preview.
- `tests/test_galaxy_knowledge.py`: no-network fake SIMBAD service tests, identity ambiguity, saved hash validation, cross-project isolation and original-file preservation.
- `tests/gui/test_galaxy_knowledge_gui.py`: offline dialog, explicit search, selection confirmation and offline reopen.
- `pyproject.toml`: optional `knowledge` extra installs astroquery >= 0.4.8.

### P0.7a scope and limitations

Candidate matching is **not** astrophysical identification or resolution of conflicting catalogue identities. A nearest-position match can be a foreground star or a background galaxy. A user-confirmed alias is merely a catalogue assertion and does not modify the project's saved ICRS target. No astronomical measurement, distance, redshift, velocity, or literature reference is imported as authoritative. There is no automatic position update or background synchronization. Coordinate errors and uncertainty in LSB centroids require separate evaluation.

SIMBAD service availability, network proxies, service rate limits, and package versions can cause failures; the UI must display them without modifying local evidence. Never use selected candidates as an independent validation of the study object's nature.

### Follow-up scope

P0.7b should add NED using the **new** API (legacy NED services no longer contain data ingested after January 2026), followed by explicit cross-source identity disagreement handling and publication-linked measurements. P0.8 should add immutable analysis-configuration snapshots bound to project UUID, dataset identities, parameters, and current verified data fingerprints, without automatically adopting results.

No original FITS or existing scientific results are changed or overwritten by this feature.
