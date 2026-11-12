"""Parse native QC reports without changing their metrics."""
from __future__ import annotations
import csv
import re
import math
from pathlib import Path
from typing import Optional

SUMMARY_PATTERN = re.compile(
    r"C:(?P<complete>[0-9.]+)%"
    r"\[S:(?P<single>[0-9.]+)%,D:(?P<duplicated>[0-9.]+)%\],"
    r"F:(?P<fragmented>[0-9.]+)%,"
    r"M:(?P<missing>[0-9.]+)%,"
    r"n:(?P<n>[0-9]+)"
)
COUNT_PATTERNS = {
    "complete": re.compile(r"(?m)^\s*(\d+)\s+Complete BUSCOs \(C\)\s*$"),
    "single": re.compile(
        r"(?m)^\s*(\d+)\s+Complete and single-copy BUSCOs \(S\)\s*$"
    ),
    "duplicated": re.compile(
        r"(?m)^\s*(\d+)\s+Complete and duplicated BUSCOs \(D\)\s*$"
    ),
    "fragmented": re.compile(r"(?m)^\s*(\d+)\s+Fragmented BUSCOs \(F\)\s*$"),
    "missing": re.compile(r"(?m)^\s*(\d+)\s+Missing BUSCOs \(M\)\s*$"),
    "n": re.compile(r"(?m)^\s*(\d+)\s+Total BUSCO groups searched\s*$"),
}
BUSCO_COLUMNS = (
    "busco_status",
    "busco_summary",
    "busco_complete",
    "busco_single",
    "busco_duplicated",
    "busco_fragmented",
    "busco_missing",
    "busco_n",
)
METRIC_COLUMNS = (
    "complete",
    "single",
    "duplicated",
    "fragmented",
    "missing",
    "n",
)



def quality_row(csv_path: Path) -> tuple[Optional[dict[str, str]], str]:
    try:
        with csv_path.open("rt", encoding="utf-8-sig", newline="") as handle:
            header = handle.readline()
            if not header:
                return None, "empty_csv"
            delimiter = "\t" if header.count("\t") > header.count(",") else ","
            handle.seek(0)
            reader = csv.DictReader(handle, delimiter=delimiter)
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error):
        return None, "unreadable_csv"
    if len(rows) != 1:
        return None, f"csv_rows_{len(rows)}"
    lowered = {column.lower() for column in (reader.fieldnames or [])}
    if not any("complet" in column for column in lowered):
        return None, "not_quality_csv"
    if not any("contam" in column for column in lowered):
        return None, "not_quality_csv"
    return rows[0], "ok"

def parse_summary(path: Path) -> Optional[dict[str, str]]:
    text = path.read_text(encoding="utf-8", errors="replace")
    match = SUMMARY_PATTERN.search(text.replace(" ", ""))
    if match:
        values = match.groupdict()
        metrics = [float(values[key]) for key in METRIC_COLUMNS[:-1]]
        if int(values["n"]) <= 0 or any(
            not math.isfinite(value) or not 0 <= value <= 100 for value in metrics
        ):
            return None
        complete, single, duplicated, fragmented, missing = metrics
        if abs(complete - single - duplicated) > 0.21:
            return None
        if abs(complete + fragmented + missing - 100) > 0.31:
            return None
        return values

    # BUSCO 5/6 also writes a line-oriented count section.  This is the
    # format used by the user's existing Notion workflow and remains useful
    # when the compact C:[S,D],F,M,n summary is absent.
    counts: dict[str, int] = {}
    for name, pattern in COUNT_PATTERNS.items():
        count_match = pattern.search(text)
        if count_match is None:
            return None
        counts[name] = int(count_match.group(1))
    total = counts["n"]
    if (total <= 0
        or counts["complete"] != counts["single"] + counts["duplicated"]
        or counts["complete"] + counts["fragmented"] + counts["missing"] != total):
        return None

    def percentage(count: int) -> str:
        return f"{100.0 * count / total:.6f}".rstrip("0").rstrip(".")

    return {
        "complete": percentage(counts["complete"]),
        "single": percentage(counts["single"]),
        "duplicated": percentage(counts["duplicated"]),
        "fragmented": percentage(counts["fragmented"]),
        "missing": percentage(counts["missing"]),
        "n": str(total),
    }
