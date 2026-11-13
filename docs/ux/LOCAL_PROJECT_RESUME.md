# Reopening previously downloaded DESI Legacy galaxy projects

The standard GUI now has **Open existing FITS project…**, an offline-only
entry point that does not require the user to start with Find data.

## Existing project migration

For example, choose:

\`\`\`text
C:\Users\frigi\Documents\PhD\LSBpipeData\UGC_2162
\`\`\`

The folder can either contain \`dr9/north/coadd/...\` directly or a
\`data/dr9/north/coadd/...\` child. Supported releases are DR9, DR10, DR11
with their existing North/South conventions.

The local project scanner reads original **filenames and lengths** to
recognize per-region, per-brick, per-band IMAGE/INVVAR/NEXP/MASKBITS
quartets. A zero-byte, missing or out-of-place file does not count as
present. This is **not** an archival checksum, WCS QA or photometric
acceptance. The existing strict science FITS loader still checks the
four matching original HDUs/WCS before running a diagnostic.

Original target coordinates are taken, in order, from an existing
\`lsbpipe-project.json\`, the original \`download-report.json\`, or explicit
ICRS values typed in the GUI before opening the folder. Never guess
a galaxy center from the DESI brick name. If the download report is
missing and the target was not explicitly entered, opening fails with an
actionable message instead of silently using a default target.

The local brick list displays **unknown** median NEXP, and does not
draw authoritative survey brick outlines. The archive is still available
through explicit **Find data**, which can re-query official geometry.
Double-click on an offline-restored brick never masquerades a placeholder
as a real brick center.

Any earlier QA run in \`qa-runs/\` with a matching receipt, selected
release, region, brick, band and exact target field can be presented as
**restored history**. Its source SHA-256 is *not* revalidated at reopen;
running the relevant Science QA action generates a new, independently
fingerprinted run. The last action is never silently reexecuted.

## Session storage

At first import, the app atomically creates \`lsbpipe-project.json\`
beside the original data folder (outside \`data/\` if present). Its fields
include original ICRS center and square width, selected release/region/
brick/band, diagnostic target and sky radii, and the last report link.
Science FITS are never overwritten. No absolute path from a foreign
project format is automatically trusted.

The app also writes a user-local \`last-project.json\` under
\`%LOCALAPPDATA%\LSBpipe\` (or \`~/.config/lsbpipe/\` outside Windows).
On the next launch, the latest project is scanned and reopened from the
real filesystem *without an archive search or automatic download*.
The old saved JSON report remains read-only. If files are moved or
missing, import warns or reports incomplete local quartets rather than
inventing them.

During automated tests, \`LSBPIPE_SESSION_DIR\` isolates this pointer
from the actual Windows user's active project.

## File provenance and future work

The files table retains the independent HTTP evidence column. An existing
local product receives an explicit **LOCAL PRESENT (unverified)**
tooltip and disk-size display; a missing file receives **NOT FOUND**.
Download selected is a no-op when every selected file already exists.
When a completed download finishes, the new session is registered for
future startup recovery.

Full science-workflow checkpointing, pixel-level memory paging, resume
of *unfinished* downloads and automatic data-integrity revalidation
are separate later milestones. Currently the saved project restores the
known inputs and last completed QA report, not an interrupted in-memory
calculation or a certified ready-to-publish photometric reduction.


## Object-bound Photometry Workspace recovery (SC1 / WP06, 2026-11-09)

The newer **Galaxy Projects** hub opens each galaxy using its immutable
object-project UUID. Opening a previously saved object project restores
the last **Photometry → PSF-0 / Scene & Masks** tab and active native
`g/r/i/z` band when the remembered release/region/brick still matches.
This pointer is an atomic, project-local
`photometry-workspace-state.json`, not an adopted scientific result.
Explicitly returning to Survey Explorer records that location instead.

Completed, *experimental* PSF-0 models are checkpointed in the existing
immutable `analyses/psf-pilots` archive after both single-band and batch
builds. On Photometry entry a background worker calls
`read_saved_psf0_bands`: it checks matching original IMAGE/INVVAR/NEXP/
MASKBITS SHA-256, saved Gaia/Tractor catalogue hash where applicable,
the actual PSF recipe, saved kernel integrity and the dataset identity.
Verified models and QA are restored in their original bands without
reconstructing stellar cutouts, rerunning ePSF fitting or downloading FITS.
Bad or stale checkpoints remain unusable and produce visible warnings.
A failed checkpoint write leaves its preview accessible in the current
session, but does not claim persistence. No original FITS is changed.

`Scene & Masks` restores the revisioned ICRS mask draft from its
existing store and only loads the original image cutout after a deliberate
user action. PSF/scene restoration does not imply scientific adoption,
accepted segmentation or final photometry. Previously computed kernels
and independent bands do not become invalid simply because the GUI's
default viewing band changes.

The Project Hub offers an **opt-in** `Open last galaxy on startup`
preference stored by Qt's standard user-local `QSettings`. With the
option disabled (default), launch shows Galaxy Projects. With it enabled,
launch enters the most recent accessible object project and follows its
saved workspace position, without remote access or automatic calculations.
The `Galaxy Projects` navigation remains available at all times.

This implements recovery of completed, verified *experimental products*
and scene drafts. Full crash-resumable numerical jobs, cache DAG
invalidation and recovery of in-progress photometry remain WP06 follow-up
work. The legacy standalone last-FITS-project pointer described above is
separate from this object-bound startup preference.
