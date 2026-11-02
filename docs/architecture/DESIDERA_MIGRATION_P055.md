# DESIDERA — project identity and compatibility (P0.5.5)

The canonical human-facing name and installable Python distribution is
**DESIDERA**. Use `python -m desidera.gui` for the desktop Project Hub and
`python -m desidera --help` for the CLI. Editable installation:
`python -m pip install -e ".[gui,test]"`.

The legacy Python implementation namespace `lsbpipe` deliberately remains
importable so ongoing research notebooks, original photometric QA manifests,
and older scripts are not broken. Old entry points `python -m lsbpipe.gui`
and `python -m lsbpipe` remain supported. New source code should use
the DESIDERA name in user-visible UI, notices and documentation; avoid
changing published result provenance.

**On-disk migration is non-destructive.** Existing
`lsbpipe-object.json`, `lsbpipe-project.json`,
`lsbpipe.local-project/0.1`, `lsbpipe.object-project/1`, project UUIDs,
cached Legacy Survey metadata, old `LSBpipeData` and `LSBpipeProjects`
directories, and saved original FITS remain under their pre-existing paths.
Only the default directory proposed for *new* projects uses
`DESIDERAProjects`. New archive source FITS must never be re-downloaded
merely because the visible product title changes.

The repository GitHub slug may remain `frigidable/LSBpipe` until an
authorized repository rename is performed in GitHub settings. GitHub
repositories, external links, package upload identities, Python import
paths, and on-disk science schemas are separate namespaces and should not
be silently rewritten in a code-only rename.

The small semi-transparent overdraw seen at the top left of Windows
WebEngine workspaces has not yet been conclusively attributed. For an
A/B test use `DESIDERA_WEBENGINE_SOFTWARE=1` before
`python -m desidera.gui`. This passes Chromium
`--disable-gpu --disable-gpu-compositing` before PySide6 QtWebEngine
initialization. Treat it as diagnostic only, not proof of the cause.
