"""Presentation-only workflow guidance for the three essential galaxy steps.

Does not grant science acceptance, mutate a project, infer a catalog
selection, or start work. It reflects explicit scientific readiness.
"""
from __future__ import annotations
from dataclasses import dataclass


@dataclass(frozen=True)
class WorkflowCue:
    step: int
    title: str
    instruction: str

    def __post_init__(self):
        if self.step not in (1, 2, 3):
            raise ValueError("Invalid workflow stage")


def workflow_cue(*, has_target: bool, has_catalog: bool,
                 selected_brick: bool, has_local_fits: bool,
                 has_interactive_image: bool) -> WorkflowCue:
    """Derive an honest next action from current UI, without inference."""
    if not has_target:
        return WorkflowCue(1, "01 / 03 · Target",
                           "Enter ICRS coordinates, then find survey data.")
    if not has_catalog and not has_local_fits:
        return WorkflowCue(1, "01 / 03 · Target",
                           "Find data to see which survey bricks cover this field.")
    if not selected_brick:
        return WorkflowCue(2, "02 / 03 · Data",
                           "Select an available brick to continue.")
    if not has_local_fits:
        return WorkflowCue(2, "02 / 03 · Data",
                           "Choose the four science FITS products, then download.")
    if not has_interactive_image:
        return WorkflowCue(3, "03 / 03 · Inspect",
                           "Original FITS found. Open the selected band to inspect.")
    return WorkflowCue(3, "03 / 03 · Inspect",
                       "Explore FITS pixels and masks; diagnostics are optional.")
