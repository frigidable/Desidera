# DESIDERA — GitHub development and review workflow

## Operating model
The owner supplies scientific requirements and accepts the final
scientific scope. ChatGPT may choose routine design details without
serial approval: layout, typography, spacing, component grouping,
icons, colors, responsive behavior and implementation refactoring
that does not alter scientific meaning.

Tooling: connected GitHub read/write, PRs, GitHub Actions, Qt UI
regression tests and recorded screenshots. Codex, Figma and paid
design services are **optional**, not required.

This is **session-driven autonomy**. GitHub Actions runs on configured
events automatically. ChatGPT cannot continue working or review new
artifacts after a chat ends unless a separate supported task/agent
is explicitly configured. Report this distinction honestly.

## Work stages
1. READ: fetch current branch SHA, source, docs, open PRs, test logs.
2. AUDIT: classify confirmed defects vs hypotheses; capture source
   evidence and screenshots where possible.
3. DESIGN: independently propose and select a coherent solution that
   fits PRODUCT_BRIEF and DESIGN_SYSTEM; preserve science truth.
4. IMPLEMENT: make small reviewable commits on a branch based on the
   correct scientific parent branch.
5. VALIDATE: execute deterministic tests and obtain GitHub Actions
   status; capture offline screenshots and logs.
6. REVIEW: inspect **actual** screenshot images when tools can retrieve
   them; record what was observed and what remains unverified.
7. REPAIR: correct failures/regressions and repeat checks.
8. HANDOFF: leave a draft PR with objective evidence, unresolved gaps
   and explicit user acceptance boundary. Never merge automatically.

## Branch discipline
Check the latest `main` SHA, open PRs and the specific task's code
dependencies before choosing a base. Routine documentation, CI and
compatibility maintenance should normally branch from `main`.
If a change depends on an unmerged scientific work package, start from
that work package's exact parent branch and target that parent in the PR,
without silently including the rest of the experimental stack.

Use **one short-lived branch per coherent change**; repair failures on the
same branch and in the same PR instead of making successive housekeeping
branches. Do not rebase, merge or retarget independent scientific drafts
without explicit approval. GitHub Actions on the exact PR head and manual
review are required before proposing a merge. ChatGPT must not auto-merge.

## Reproducibility contract
- Offline tests must not query the survey archive or remote Viewer.
- Separate `preview`, `archive candidate`, `HEAD evidence`,
  `downloaded`, `diagnostic` and `science-adopted` statuses.
- Publish screen-resolution, state, screenshot filenames, test logs
  and the exact branch commit with the CI job.
- Document rendering limitations of offscreen QtWebEngine;
  a blank placeholder must never be described as remote imagery.
- Do not claim success from a workflow that has not actually completed.
- If connector cannot expose image bytes, label screenshot visual
  review `unverified` and retain artifact links for later inspection.

## Review deliverable template
Branch and base | changes | test command / run | passed / failed /
not run | screenshot artifacts captured? | screenshots **inspected**? |
scope exclusions | remaining defects | PR link.
