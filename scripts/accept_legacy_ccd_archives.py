"""P0.6c: fail-closed real FITS archive parser acceptance (no science adoption).

The three fixed official CCD-table URLs are downloaded afresh; success requires
all three and independent Astropy row agreement with Desidera's parser.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import time
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from astropy.io import fits
from lsbpipe.surveys.legacy_ccd_catalog import _read_table

BASE_URL = "https://portal.nersc.gov/cfs/cosmo/data/legacysurvey"
SAMPLES = (("dr9", "south", "1554p180"),
           ("dr11", "south", "1554p180"),
           ("dr9", "north", "2100p542"))
REQUIRED = frozenset(("camera", "expnum", "ccdname", "filter",
                      "brick_x0", "brick_x1", "brick_y0", "brick_y1"))
MAX_BYTES = 32 * 1024 * 1024


def archive_url(release: str, region: str, brick: str,
                base_url: str = BASE_URL) -> str:
    if (release not in ("dr9", "dr10", "dr11") or
            region not in ("north", "south") or
            len(brick) != 8 or not brick[:4].isdigit() or
            brick[4] not in "pm" or not brick[5:].isdigit()):
        raise ValueError("Invalid archive sample identity")
    return (f"{base_url.rstrip('/')}/{release}/{region}/coadd/{brick[:3]}/"
            f"{brick}/legacysurvey-{brick}-ccds.fits")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def download(url: str, destination: Path, timeout: int = 75) -> None:
    if not url.startswith("https://"):
        raise ValueError("Only HTTPS archive URLs accepted")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temp = destination.with_suffix(".fits.incomplete")
    try:
        for attempt in range(2):
            try:
                request = Request(url, headers={
                    "User-Agent": "Desidera-P0.6c-research-audit/1.0"})
                with urlopen(request, timeout=timeout) as response, temp.open("wb") as out:
                    if not response.geturl().startswith("https://"):
                        raise RuntimeError("Unsafe archive redirect")
                    if response.headers.get("Content-Type", "").startswith("text/html"):
                        raise RuntimeError("Archive returned HTML rather than FITS")
                    total = 0
                    while chunk := response.read(1024 * 1024):
                        total += len(chunk)
                        if total > MAX_BYTES:
                            raise RuntimeError("CCD FITS larger than 32 MiB")
                        out.write(chunk)
                if total == 0:
                    raise ValueError("Empty archive response")
                temp.replace(destination)
                return
            except (HTTPError, URLError, TimeoutError) as exc:
                if isinstance(exc, HTTPError) and exc.code in (400, 401, 403, 404, 410):
                    raise
                if attempt:
                    raise
                time.sleep(2)
    finally:
        if temp.exists():
            temp.unlink()


def _text(value) -> str:
    if isinstance(value, bytes):
        value = value.decode("ascii")
    return str(value).strip()


def independent_rows(path: Path) -> tuple[dict, list[tuple]]:
    """Independent Astropy-only extraction, with no Desidera parser helpers."""
    with fits.open(path, mode="readonly", memmap=False) as hdul:
        indexes = [i for i, h in enumerate(hdul) if isinstance(h, fits.BinTableHDU)]
        if len(indexes) != 1:
            raise ValueError("Independent reader requires exactly one FITS BINTABLE")
        table = hdul[indexes[0]]
        columns = {name.lower(): name for name in table.columns.names}
        if not REQUIRED.issubset(columns):
            raise ValueError("Missing required CCD columns in independent reader")
        formats = {name: str(table.columns[columns[name]].format)
                   for name in sorted(REQUIRED)}
        rows = [(_text(row[columns["camera"]]).lower(),
                 str(int(row[columns["expnum"]])),
                 _text(row[columns["ccdname"]]).upper(),
                 _text(row[columns["filter"]]).lower(), index)
                for index, row in enumerate(table.data)]
        return {"hdu": indexes[0], "count": len(rows),
                "formats": formats,
                "column_names": list(table.columns.names)}, rows


def audit_file(path: Path, *, release: str, region: str, brick: str,
               source_url: str) -> tuple[dict, list[tuple]]:
    if path.name != f"legacysurvey-{brick}-ccds.fits":
        raise ValueError("Not the exact original per-brick CCD filename")
    independent, expected = independent_rows(path)
    if not expected:
        raise ValueError("Empty original CCD catalogue")
    counts = Counter(r[3] for r in expected)
    by_band = {}
    for band, count in sorted(counts.items()):
        parsed, n, meta = _read_table(path, band)
        actual = Counter((r["instrument"], r["exposure_id"],
                          r["detector_id"].upper(), r["band"],
                          int(r["record_locator"].split(":")[-1]))
                         for r in parsed)
        reference = Counter(r for r in expected if r[3] == band)
        if actual != reference or n != count or meta["catalogue_rows"] != len(expected):
            raise AssertionError(f"{band}: Desidera/independent row disagreement")
        by_band[band] = {"count": n, "fits_integrity": meta["fits_integrity"]}
    physical = [(camera, exp, ccd, band) for camera, exp, ccd, band, _ in expected]
    if len(physical) != len(set(physical)):
        raise AssertionError("Duplicate physical CCD")
    stable = sorted(physical)
    return ({
        "release": release, "region": region, "brick": brick,
        "url": source_url, "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
        "row_count": len(physical),
        "unique_ccd_count": len({x[:3] for x in physical}),
        "unique_exposure_count": len({x[:2] for x in physical}),
        "fits_hdu": independent["hdu"],
        "column_formats": independent["formats"],
        "column_names": independent["column_names"], "bands": by_band,
        "ccd_rows_sha256": hashlib.sha256(
            json.dumps(stable, separators=(",", ":")).encode()).hexdigest(),
        "parser_vs_independent_astropy": "PASS",
        "release_and_region_binding": "URL_PATH_ASSERTION_NOT_INDEPENDENTLY_CERTIFIED",
        "ccd_completeness": "NOT_ASSESSED",
        "pixel_level_overlap": "NOT_ASSESSED",
        "science_status": "NOT_SCIENTIFICALLY_ADOPTED",
    }, stable)


def report_markdown(report: dict) -> str:
    lines = ["# Desidera P0.6c — original archive CCD acceptance", "",
             f"UTC: {report['created_at']}", f"Status: **{report['status']}**", "",
             "| Release | Region | Brick | Rows | CCDs | Exposures | Parser |",
             "|---|---|---|---:|---:|---:|---|"]
    for x in report["samples"]:
        lines.append(f"| {x['release']} | {x['region']} | {x['brick']} | "
                     f"{x['row_count']} | {x['unique_ccd_count']} | "
                     f"{x['unique_exposure_count']} | {x['parser_vs_independent_astropy']} |")
    if report["failures"]:
        lines += ["", "## Failures"] + [f"- {e}" for e in report["failures"]]
    lines += ["", "## DR9/DR11 CCD comparison", "",
              "```json", json.dumps(report["cross_release"], indent=2), "```",
              "", "This is parser/schema acceptance only. Archive path is not",
              "independent release authentication. CCD pixel contributions,",
              "coadd IMAGE hashes, WCS, PSF, calibration, noise, and",
              "statistical independence have not been certified.", ""]
    return "\n".join(lines)


def run(cache: Path, report_dir: Path, base_url: str = BASE_URL,
        *, skip_download: bool = False) -> bool:
    report_dir.mkdir(parents=True, exist_ok=True)
    samples, rows_by_sample, failures = [], {}, []
    for release, region, brick in SAMPLES:
        url = archive_url(release, region, brick, base_url)
        path = cache / release / region / brick / f"legacysurvey-{brick}-ccds.fits"
        try:
            if not skip_download:
                download(url, path)
            if not path.is_file():
                raise FileNotFoundError(f"Archive sample missing: {path}")
            result, rows = audit_file(path, release=release, region=region,
                                      brick=brick, source_url=url)
            samples.append(result)
            rows_by_sample[(release, region, brick)] = rows
            print(f"PASS {release}/{region}/{brick}: {result['row_count']} rows")
        except Exception as exc:
            failure = f"{release}/{region}/{brick}: {type(exc).__name__}: {exc}"
            print("FAIL " + failure, file=sys.stderr)
            failures.append(failure)
    cross = {"status": "NOT_EVALUATED", "independence_certified": False,
             "joint_scientific_combination_allowed": False}
    left = ("dr9", "south", "1554p180")
    right = ("dr11", "south", "1554p180")
    if left in rows_by_sample and right in rows_by_sample:
        a, b = set(rows_by_sample[left]), set(rows_by_sample[right])
        shared = sorted(a & b)
        cross = {
            "status": "COMPARISON_COMPLETE_NO_INDEPENDENCE_CLAIM",
            "dr9_ccd_count": len(a), "dr11_ccd_count": len(b),
            "shared_ccd_count": len(shared), "shared_ccd_keys": shared,
            "dr9_only_count": len(a - b), "dr11_only_count": len(b - a),
            "shared_exposure_count": len({r[:2] for r in a} & {r[:2] for r in b}),
            "independence_certified": False,
            "joint_scientific_combination_allowed": False}
    if len(samples) != len(SAMPLES):
        failures.append("Not all required real archive cases passed")
    report = {
        "schema": "desidera.real-ccd-archive-acceptance/1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": "INCOMPLETE_OR_FAILED" if failures else "PASS_LIMITED_PARSER_ACCEPTANCE",
        "samples": samples, "cross_release": cross, "failures": failures,
        "certifies_real_archive_parser_compatibility": not bool(failures),
        "certifies_scientific_independence": False}
    (report_dir / "legacy_ccd_real_archive_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (report_dir / "legacy_ccd_real_archive_report.md").write_text(
        report_markdown(report), encoding="utf-8")
    return not failures


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, default=Path(".real-ccd-downloads"))
    parser.add_argument("--report-dir", type=Path, default=Path("real-ccd-evidence"))
    parser.add_argument("--base-url", default=BASE_URL)
    parser.add_argument("--skip-download", action="store_true")
    args = parser.parse_args()
    return 0 if run(args.cache_dir, args.report_dir, args.base_url,
                    skip_download=args.skip_download) else 1


if __name__ == "__main__":
    raise SystemExit(main())
