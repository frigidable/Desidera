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

## Windows top-left overlay isolation (P0.5.6)

The overlapping translucent controls near the workspace's top-left corner
are a reported real-machine rendering defect, not a FITS/science problem.
The original `DESIDERA_WEBENGINE_SOFTWARE=1` switch disables Chromium
GPU/compositing but may provoke GLES context initialization failures and
is **not** a proven fix. The new `DESIDERA_DISABLE_WEBENGINE=1` switch
does not instantiate a `QWebEngineView` at all; the embedded remote tab
becomes an explicit diagnostic placeholder. It never claims that the hosted
Viewer loaded. Legacy catalog discovery, archive FITS selection, Dataset
Manager and local analysis remain enabled; an external browser can still
open the official Viewer.

Compare screenshots of the exact same galaxy workspace, with and without
this switch and with no other graphics switches. If artifacts disappear
without any WebEngine view, embedded Chromium/Qt texture composition is
implicated; if they remain, inspect overlays/external software/Qt QWidget
painting before trying GPU or driver changes. It is not necessary to
reinstall the Python environment or redownload survey FITS.
