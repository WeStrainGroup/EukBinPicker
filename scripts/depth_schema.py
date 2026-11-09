#!/usr/bin/env python3
"""Identify generic and JGI depth-table coverage dimensions.

The public CLI accepts generic coverage matrices, then materializes a canonical
JGI table.  Downstream preparation/filtering consumes that canonical table.
Keeping header interpretation here prevents the binners from silently choosing
different coverage columns.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence


ID_COLUMNS = ("contigName", "contig_id", "contigname", "contig")
JGI_CORE_COLUMNS = ("contigName", "contigLen", "totalAvgDepth")


class DepthSchemaError(ValueError):
    """The depth header is ambiguous or incompatible with the contract."""


@dataclass(frozen=True)
class DepthSchema:
    columns: tuple[str, ...]
    id_index: int
    length_index: Optional[int]
    total_index: Optional[int]
    coverage_names: tuple[str, ...]
    source_coverage_names: tuple[str, ...]
    coverage_indices: tuple[int, ...]
    variance_indices: tuple[Optional[int], ...]
    source_format: str
    total_avg_depth_policy: str

    @property
    def coverage_column_count(self) -> int:
        return len(self.coverage_names)


def _validate_columns(columns: Sequence[str]) -> tuple[str, ...]:
    normalized = tuple(columns)
    if not normalized:
        raise DepthSchemaError("depth table header is empty")
    unsafe = [
        repr(column)
        for column in normalized
        if not column or column != column.strip() or "\x00" in column
    ]
    if unsafe:
        raise DepthSchemaError(
            "depth table has empty or unsafe column name(s): "
            + ", ".join(unsafe)
        )
    duplicates = sorted(
        {column for column in normalized if normalized.count(column) > 1}
    )
    if duplicates:
        raise DepthSchemaError(
            "depth table has duplicate column name(s): "
            + ", ".join(duplicates)
        )
    return normalized


def _validate_normalized_names(coverage_names: Sequence[str]) -> None:
    normalized_columns = list(JGI_CORE_COLUMNS)
    for name in coverage_names:
        normalized_columns.extend((name, f"{name}-var"))
    duplicates = sorted(
        {
            column
            for column in normalized_columns
            if normalized_columns.count(column) > 1
        }
    )
    if duplicates:
        raise DepthSchemaError(
            "coverage names collide after JGI normalization: "
            + ", ".join(duplicates)
        )


def identify_depth_schema(
    columns: Sequence[str],
    sample: str,
    *,
    require_jgi: bool = False,
    single_coverage_name: Optional[str] = None,
    contig_column: Optional[str] = None,
    coverage_columns: Optional[Sequence[str]] = None,
    variance_columns: Optional[Sequence[str]] = None,
) -> DepthSchema:
    """Return deterministic indexes and names for one depth-table header."""

    normalized = _validate_columns(columns)
    column_set = set(normalized)
    if contig_column is not None or coverage_columns is not None or variance_columns is not None:
        if not contig_column or not coverage_columns:
            raise DepthSchemaError("custom tables require --contig-column and --coverage-columns together")
        if variance_columns is not None and len(variance_columns) != len(coverage_columns):
            raise DepthSchemaError("variance columns must pair one-to-one with coverage columns in the same order")
        selected = [contig_column, *coverage_columns, *(variance_columns or ())]
        if len(set(selected)) != len(selected):
            raise DepthSchemaError("contig, coverage and variance columns must be distinct")
        missing = set(selected) - column_set
        if missing:
            raise DepthSchemaError("unknown selected column(s): " + ", ".join(sorted(missing)))
        if "contigLen" in selected:
            raise DepthSchemaError("contigLen is reserved for FASTA length validation")
        names = tuple(coverage_columns)
        if single_coverage_name is not None:
            if len(names) != 1:
                raise DepthSchemaError("coverage aliases require exactly one coverage column")
            names = (single_coverage_name,)
        _validate_normalized_names(names)
        return DepthSchema(
            columns=normalized, id_index=normalized.index(contig_column),
            length_index=normalized.index("contigLen") if "contigLen" in column_set else None,
            total_index=None, coverage_names=names,
            source_coverage_names=tuple(coverage_columns),
            coverage_indices=tuple(normalized.index(x) for x in coverage_columns),
            variance_indices=tuple(normalized.index(x) for x in variance_columns) if variance_columns else tuple(None for _ in names),
            source_format="custom_coverage", total_avg_depth_policy="sum_of_coverage_columns",
        )
    is_jgi = all(column in column_set for column in JGI_CORE_COLUMNS)

    if is_jgi:
        id_index = normalized.index("contigName")
        length_index = normalized.index("contigLen")
        total_index = normalized.index("totalAvgDepth")
        extra_indices = tuple(
            index
            for index, name in enumerate(normalized)
            if name not in JGI_CORE_COLUMNS
        )
        coverage_names: list[str] = []
        coverage_indices: list[int] = []
        variance_indices: list[Optional[int]] = []
        if extra_indices:
            cursor = 0
            while cursor < len(extra_indices):
                depth_index = extra_indices[cursor]
                depth_name = normalized[depth_index]
                if depth_name.endswith("-var"):
                    raise DepthSchemaError(
                        "JGI variance column lacks its preceding depth column: "
                        f"{depth_name}"
                    )
                if cursor + 1 >= len(extra_indices):
                    raise DepthSchemaError(
                        f"JGI depth column lacks variance pair: {depth_name}"
                    )
                variance_index = extra_indices[cursor + 1]
                expected_variance = f"{depth_name}-var"
                if normalized[variance_index] != expected_variance:
                    raise DepthSchemaError(
                        f"JGI depth column {depth_name} must be followed by "
                        f"{expected_variance}"
                    )
                coverage_names.append(depth_name)
                coverage_indices.append(depth_index)
                variance_indices.append(variance_index)
                cursor += 2
            source_coverage_names = tuple(coverage_names)
            source_format = "jgi_depth_pairs"
            if single_coverage_name is not None:
                if len(coverage_names) != 1:
                    raise DepthSchemaError("coverage aliases require one JGI depth pair")
                coverage_names = [single_coverage_name]
        else:
            fallback_name = single_coverage_name or f"{sample}.bam"
            coverage_names = [fallback_name]
            coverage_indices = [total_index]
            variance_indices = [None]
            source_coverage_names = ("totalAvgDepth",)
            source_format = "jgi_total_only"
        _validate_normalized_names(coverage_names)
        return DepthSchema(
            columns=normalized,
            id_index=id_index,
            length_index=length_index,
            total_index=total_index,
            coverage_names=tuple(coverage_names),
            source_coverage_names=source_coverage_names,
            coverage_indices=tuple(coverage_indices),
            variance_indices=tuple(variance_indices),
            source_format=source_format,
            total_avg_depth_policy="provided",
        )

    if require_jgi:
        missing = [
            column for column in JGI_CORE_COLUMNS if column not in column_set
        ]
        raise DepthSchemaError(
            "depth table is missing required JGI column(s): "
            + ", ".join(missing)
        )

    id_column = next(
        (name for name in ID_COLUMNS if name in column_set),
        normalized[0],
    )
    id_index = normalized.index(id_column)
    length_index = (
        normalized.index("contigLen") if "contigLen" in column_set else None
    )
    metadata_indices = {id_index}
    if length_index is not None:
        metadata_indices.add(length_index)
    coverage_indices = tuple(
        index for index in range(len(normalized)) if index not in metadata_indices
    )
    if not coverage_indices:
        raise DepthSchemaError(
            "generic depth table requires at least one coverage column"
        )
    source_names = tuple(normalized[index] for index in coverage_indices)
    if len(coverage_indices) == 1:
        coverage_names = (single_coverage_name or f"{sample}.bam",)
        source_format = "generic_single_coverage"
    else:
        if "totalAvgDepth" in source_names:
            raise DepthSchemaError(
                "generic multi-coverage input must omit totalAvgDepth; it is "
                "derived as the sum of coverage columns"
            )
        coverage_names = source_names
        source_format = "generic_multi_coverage"
    _validate_normalized_names(coverage_names)
    return DepthSchema(
        columns=normalized,
        id_index=id_index,
        length_index=length_index,
        total_index=None,
        coverage_names=coverage_names,
        source_coverage_names=source_names,
        coverage_indices=coverage_indices,
        variance_indices=tuple(None for _ in coverage_indices),
        source_format=source_format,
        total_avg_depth_policy="sum_of_coverage_columns",
    )


def write_metabat_depth(source, destination):
    """Keep observed JGI variance, or write MetaBAT2's native mean-only format."""
    import csv
    from pathlib import Path
    with Path(source).open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        schema = identify_depth_schema(header, "input", require_jgi=True)
        rows = list(reader)
    present = {row[index] != "NA" if index is not None else False
               for row in rows for index in schema.variance_indices}
    if len(present) > 1:
        raise DepthSchemaError("MetaBAT2 cannot mix observed and missing variance; provide all variances or omit metabat2")
    if present == {True}:
        import shutil
        shutil.copyfile(source, destination)
        return "jgi"
    with Path(destination).open("w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["contigName", *schema.coverage_names])
        for row in rows:
            writer.writerow([row[schema.id_index], *(row[i] for i in schema.coverage_indices)])
    return "means"


if __name__ == "__main__":
    import sys
    try:
        print(write_metabat_depth(sys.argv[1], sys.argv[2]))
    except (ValueError, OSError) as error:
        sys.exit(str(error))
