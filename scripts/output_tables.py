"""Readable QC and representative tables; no selection decisions are made here."""
from __future__ import annotations

import csv
from pathlib import Path

from merge_eukcc import atomic_write


def read(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    import io
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    atomic_write(path, buffer.getvalue())


def export_tables(out: Path, candidates: Path, eukcc: Path, busco: Path,
                  qc_tools: tuple[str, ...]) -> None:
    columns = ["bin_id", "sample", "binner", "original_bin_id", "contigs", "bp",
               "n50", "longest_contig", "sha256", "bin_path"]
    _, candidates_rows = read(candidates)
    rows = [{key: row[key] for key in columns} for row in candidates_rows]
    for tool, path in (("eukcc", eukcc), ("busco", busco)):
        if tool not in qc_tools:
            continue
        fields, tool_rows = read(path)
        by_id = {row["bin_id"]: row for row in tool_rows}
        # All fields beyond candidate metadata are native QC results/statuses.
        metadata = set(read(candidates)[0])
        metrics = [key for key in fields if key not in metadata]
        names = {key: key if key.startswith(tool + "_") else tool + "_" + key
                 for key in metrics}
        columns.extend(names.values())
        for row in rows:
            values = by_id[row["bin_id"]]
            row.update({name: values.get(key, "") for key, name in names.items()})
    write(out / "2_qc/summary.tsv", columns, rows)

    qc_by_id = {row["bin_id"]: row for row in rows}
    rep_columns, representatives = read(out / "0_log/selection/representatives.tsv")
    _, alternatives = read(out / "0_log/selection/alternatives.tsv")
    _, manifests = read(out / "3_selection/bins/manifest.tsv")
    paths = {row["representative_bin_id"]: row["bin_path"] for row in manifests}
    represented: dict[str, list[str]] = {}
    for row in alternatives:
        represented.setdefault(row["representative_bin_id"], []).append(row["alternative_bin_id"])
    extra = [key for key in columns if key not in rep_columns and key != "bin_id"]
    final_columns = [*rep_columns, *extra, "selection_reason", "represented_bins"]
    for row in representatives:
        bin_id = row["representative_bin_id"]
        row.update({key: qc_by_id[bin_id][key] for key in extra})
        row["bin_path"] = paths[bin_id]
        row["selection_reason"] = "passed_thresholds_and_no_overlap_with_higher_ranked_cross_binner_bin"
        # Each alternative has one primary representative, as chosen by the selector.
        row["represented_bins"] = ";".join([bin_id, *represented.get(bin_id, [])])
    write(out / "3_selection/representatives.tsv", final_columns, representatives)
