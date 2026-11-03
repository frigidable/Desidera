"""Immutable GUI science context selection and stale-result guards.

Pure selection logic: no file I/O, no QThread, no publication. The owning
SurveyExplorer retains public compatibility wrappers.
"""
from __future__ import annotations

from lsbpipe.gui.qa_records import ScienceContext


class ScienceContextController:
    def __init__(self, owner):
        self.owner = owner

    def current(self, action):
        w = self.owner
        if not w._brick_result_matches_target():
            return None
        row = w.bricks.currentRow()
        if not (w.bricks.selectedItems() and 0 <= row < len(w.hits)):
            return None
        try:
            field = w._field()
        except ValueError:
            return None
        return ScienceContext.capture(
            action=action, field=field, hit=w.hits[row],
            band=w.science_qa.current_band(),
            directory=w.science_qa.directory.text().strip(),
            north_directory=w.science_qa.north_directory.text().strip(),
            south_directory=w.science_qa.south_directory.text().strip(),
            target_radius_arcsec=w.target_radius.value(),
            sky_inner_arcsec=w.sky_inner.value(),
            sky_outer_arcsec=w.sky_outer.value(),
            project_id=(w._object_project.project_id
                        if w._object_project is not None else None))

    def ticket_current(self, *, context, ticket):
        w = self.owner
        return (ticket == w._qa_generation and
                self.current(context.action) == context)

    def raster_current(self, dataset):
        """Compare displayed science field and explicit brick/band identity."""
        w = self.owner
        try:
            field = w._field()
        except ValueError:
            return False
        saved_field = dataset["field"]
        if any(abs(a - b) > 1e-9 for a, b in (
                (field.ra_deg, saved_field.ra_deg),
                (field.dec_deg, saved_field.dec_deg),
                (field.width_arcmin, saved_field.width_arcmin))):
            return False
        row = w.bricks.currentRow()
        if not (w.bricks.selectedItems() and 0 <= row < len(w.hits)):
            return False
        hit = w.hits[row]
        identity = (
            f"{hit.release.upper()} {hit.region} · {hit.name} · "
            f"{w.science_qa.current_band()}")
        return dataset["source"] == identity
