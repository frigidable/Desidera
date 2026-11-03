# P0.4 UX: Native FITS Explorer and interactive Four-layer QA

This is a UI refactor: original FITS, pixels, WCS, masks, QA status and
scientific calibration remain unchanged.

## Single-layer FITS Explorer

The image now occupies the main viewport, with scrollable display controls
in a narrow left sidebar. Layer, contrast, stretch, numerical limits,
Fit, smooth-preview and aperture-outline controls remain available.
Individual SATUR, ALLMASK, Gaia BRIGHT/MEDIUM, SGA GALAXY, and invalid
pixel overlays are available under **Mask overlays**. HDU provenance,
cache statistics, numeric legend and range details are disclosed under
**Details & provenance**. The original FITS are not modified. Mouse
wheel zoom, pan and per-pixel ICRS/value inspection continue to work.

## Interactive Four-layer QA

The initial main tabs are **Remote sky** and **FITS Explorer**.
**Compare 4 layers** is a separate button above the tabs and an
action in Science QA. It adds an interactive **Four-layer QA** tab
**after FITS Explorer**, only when requested. If no native quartet is
loaded, the existing frozen-context read-only FITS preparation runs
first. The third tab is removed with **Close comparison**.

The comparison consists of four independent, pannable, zoomable
QGraphicsViews showing IMAGE, INVVAR, NEXP, and optical MASKBITS on the
same native pixel grid. **Link zoom** synchronizes transforms; **Fit all**
resets the views. A shared hover readout reports original pixel
coordinates, ICRS and native values. Source/sky outlines are drawn in
each QGraphicsScene. No PNG is used for rendering.

NumPy conversion runs in one dedicated background worker; QPixmaps
are created on the GUI thread. Changed datasets, removal and workspace
shutdown prevent old displays from publishing. The same full-resolution
6-million-pixel limit applies; no silent downsampling or reprojection.

The former static figure action remains **Export static PNG (optional)**
as a separate output artifact, not as the four-layer on-screen interface.

## Remaining gates

Offline real Qt tests cover the sidebar and disclosures, rendered
interactive panels, linked zoom, native readout, read-only input arrays,
on-demand tab lifecycle and refusal of unverified states. Large
original DR9 FITS and high-DPI/VRAM stress tests are still required on
Windows; the comparison is diagnostic, not adopted science. A future
tiled renderer may handle larger images without excessive memory use.
