"""Quality-ranked selection using matching contig IDs and normalized sequences."""
from __future__ import annotations
from collections import Counter
from dataclasses import dataclass
import gzip
import hashlib
import math
from pathlib import Path
from typing import TextIO

REVERSE_COMPLEMENT = str.maketrans("ACGTRYMKBDHVNWS", "TGCAYRKMVHDBNWS")
@dataclass
class Candidate:
    row: dict[str, str]
    completeness: float
    contamination: float
    contig_count: int
    n50: int
    quality_score: float
    contig_bp: dict[tuple[str, str], int]
    fasta_bp: int

    @property
    def bin_id(self) -> str:
        return self.row["bin_id"]

    @property
    def binner(self) -> str:
        return self.row["binner"]

    @property
    def path(self) -> Path:
        return Path(self.row["bin_path"])


def open_text(path: Path) -> TextIO:
    with path.open("rb") as probe:
        compressed = probe.read(2) == b"\x1f\x8b"
    if compressed:
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("rt", encoding="utf-8", newline="")


def parse_float(value: str, field: str, bin_id: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{bin_id}: invalid {field}: {value!r}") from error


def parse_int(value: str, field: str, bin_id: str) -> int:
    try:
        return int(str(value).replace(",", "").strip())
    except (TypeError, ValueError) as error:
        raise ValueError(f"{bin_id}: invalid {field}: {value!r}") from error


def canonical_sequence(sequence: str) -> str:
    normalized = sequence.upper().replace("U", "T")
    reverse = normalized.translate(REVERSE_COMPLEMENT)[::-1]
    return min(normalized, reverse)


def fasta_contigs(path: Path) -> dict[tuple[str, str], int]:
    if not path.is_file():
        raise ValueError(f"bin FASTA is missing: {path}")
    contigs: dict[tuple[str, str], int] = {}
    current_header = ""
    sequence_parts: list[str] = []
    observed_ids: set[str] = set()

    def finish() -> None:
        nonlocal current_header, sequence_parts
        if not current_header:
            return
        sequence = "".join(sequence_parts).replace(" ", "").replace("\t", "")
        if not sequence:
            raise ValueError(f"{path}: empty FASTA record {current_header!r}")
        canonical = canonical_sequence(sequence)
        digest = hashlib.sha256(canonical.encode("ascii")).hexdigest()
        identifier = current_header.split()[0]
        if identifier in observed_ids:
            raise ValueError(
                f"{path}: duplicate contig ID at {current_header!r}"
            )
        observed_ids.add(identifier)
        contigs[(identifier, digest)] = len(sequence)

    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish()
                current_header = line[1:].strip()
                if not current_header:
                    raise ValueError(f"{path}: empty FASTA header")
                sequence_parts = []
            else:
                if not current_header:
                    raise ValueError(f"{path}: sequence before first FASTA header")
                sequence_parts.append(line)
    finish()
    if not contigs:
        raise ValueError(f"{path}: no FASTA records")
    return contigs


def quality_sort_key(
    *,
    completeness: float,
    contamination: float,
    contig_count: int,
    n50: int,
    bin_id: str,
    contamination_penalty: float,
) -> tuple[float, float, float, int, int, str]:
    """Return the agreed 1,3,2,5,4,6 deterministic ordering."""
    quality_score = completeness - contamination_penalty * contamination
    return (
        -quality_score,
        -completeness,
        contamination,
        contig_count,
        -n50,
        bin_id,
    )


def candidate_sort_key(
    candidate: Candidate,
    contamination_penalty: float,
) -> tuple[float, float, float, int, int, str]:
    return quality_sort_key(
        completeness=candidate.completeness,
        contamination=candidate.contamination,
        contig_count=candidate.contig_count,
        n50=candidate.n50,
        bin_id=candidate.row.get("sort_id", candidate.bin_id),
        contamination_penalty=contamination_penalty,
    )


def shared_overlap(left: Candidate, right: Candidate) -> tuple[int, float]:
    shared = sum(
        bp for digest, bp in left.contig_bp.items() if digest in right.contig_bp
    )
    denominator = min(left.fasta_bp, right.fasta_bp)
    return shared, (shared / denominator if denominator else 0.0)

def load_candidates(
    rows: list[dict],
    contamination_penalty: float,
    min_completeness: float,
    max_contamination: float,
    expected_rank: int,
    expected_study: str,
    expected_sample: str,
) -> list[Candidate]:
    columns = set(rows[0]) if rows else set()
    required = {
        "bin_id",
        "sample_rank",
        "study",
        "sample",
        "binner",
        "bin_path",
        "decision_qc",
        "contamination_source",
        "qc_pass",
        "completeness",
        "contamination",
        "contigs",
        "n50",
        "bp",
    }
    missing = sorted(required.difference(columns))
    if rows and missing:
        raise ValueError(f"QC table is missing columns: {', '.join(missing)}")

    candidates: list[Candidate] = []
    observed_ids: set[str] = set()
    for row in rows:
        observed_identity = (
            parse_int(row["sample_rank"], "sample_rank", row.get("bin_id", "")),
            row["study"],
            row["sample"],
        )
        expected_identity = (expected_rank, expected_study, expected_sample)
        if observed_identity != expected_identity:
            raise ValueError(
                "QC table contains a row outside the expected sample: "
                f"{observed_identity!r} != {expected_identity!r}"
            )
        if row["qc_pass"].strip().lower() != "yes":
            continue
        bin_id = row["bin_id"]
        if not bin_id:
            raise ValueError("eligible QC row has an empty bin_id")
        if bin_id in observed_ids:
            raise ValueError(f"duplicate eligible bin_id: {bin_id}")
        observed_ids.add(bin_id)
        completeness = parse_float(row["completeness"], "completeness", bin_id)
        contamination = parse_float(row["contamination"], "contamination", bin_id)
        if (
            not math.isfinite(completeness)
            or not 0.0 <= completeness <= 100.0
            or not math.isfinite(contamination)
            or contamination < 0.0
            or completeness < min_completeness
            or contamination > max_contamination
        ):
            raise ValueError(
                f"{bin_id}: qc_pass=yes conflicts with configured thresholds"
            )
        reported_contigs = parse_int(row["contigs"], "contigs", bin_id)
        n50 = parse_int(row["n50"], "n50", bin_id)
        reported_bp = parse_int(row["bp"], "bp", bin_id)
        contig_bp = fasta_contigs(Path(row["bin_path"]))
        observed_bp = sum(contig_bp.values())
        cumulative = 0
        observed_n50 = 0
        for length in sorted(contig_bp.values(), reverse=True):
            cumulative += length
            if cumulative * 2 >= observed_bp:
                observed_n50 = length
                break
        if n50 != observed_n50:
            raise ValueError(f"{bin_id}: FASTA N50 {observed_n50} != QC N50 {n50}")
        if len(contig_bp) != reported_contigs:
            raise ValueError(
                f"{bin_id}: FASTA contigs {len(contig_bp)} != QC contigs "
                f"{reported_contigs}"
            )
        if observed_bp != reported_bp:
            raise ValueError(
                f"{bin_id}: FASTA bp {observed_bp} != QC bp {reported_bp}"
            )
        candidates.append(
            Candidate(
                row=row,
                completeness=completeness,
                contamination=contamination,
                contig_count=reported_contigs,
                n50=n50,
                quality_score=completeness
                - contamination_penalty * contamination,
                contig_bp=contig_bp,
                fasta_bp=observed_bp,
            )
        )
    return candidates

def select_rows(rows, overlap_threshold, contamination_penalty, min_completeness,
                max_contamination, expected_rank, expected_study, expected_sample):
    if not 0.0 < overlap_threshold <= 1.0:
        raise ValueError("overlap threshold must be in (0, 1]")
    if not math.isfinite(contamination_penalty) or contamination_penalty < 0.0:
        raise ValueError("contamination penalty must be non-negative")
    if (
        not math.isfinite(min_completeness)
        or not 0.0 <= min_completeness <= 100.0
    ):
        raise ValueError("minimum completeness must be in [0, 100]")
    if (
        not math.isfinite(max_contamination)
        or not 0.0 <= max_contamination <= 100.0
    ):
        raise ValueError("maximum contamination must be in [0, 100]")
    candidates = load_candidates(rows, contamination_penalty, min_completeness,
                                 max_contamination, expected_rank, expected_study, expected_sample)
    representative_rows: list[dict[str, object]] = []
    alternative_rows: list[dict[str, object]] = []
    representative_paths: list[str] = []
    representative_by_binner: Counter[str] = Counter()
    alternative_by_binner: Counter[str] = Counter()
    ambiguous = 0

    ordered = sorted(
        candidates,
        key=lambda candidate: candidate_sort_key(
            candidate,
            contamination_penalty,
        ),
    )
    selected: list[Candidate] = []
    selection_rank = 0
    for candidate in ordered:
        matches: list[tuple[float, Candidate, int]] = []
        for representative in selected:
            # Deliberately do not inspect bins from the same binner.
            if representative.binner == candidate.binner:
                continue
            shared_bp, overlap = shared_overlap(candidate, representative)
            if overlap >= overlap_threshold:
                matches.append((overlap, representative, shared_bp))
        if not matches:
            selection_rank += 1
            selected.append(candidate)
            representative_by_binner[candidate.binner] += 1
            representative_paths.append(str(candidate.path.resolve()))
            representative_rows.append(
                {
                    "sample_rank": expected_rank,
                    "study": expected_study,
                    "sample": expected_sample,
                    "selection_rank": selection_rank,
                    "representative_bin_id": candidate.bin_id,
                    "selected_binner": candidate.binner,
                    "bin_path": str(candidate.path.resolve()),
                    "decision_qc": candidate.row["decision_qc"],
                    "contamination_source": candidate.row[
                        "contamination_source"
                    ],
                    "qc_pass": "yes",
                    "completeness": f"{candidate.completeness:g}",
                    "contamination": f"{candidate.contamination:g}",
                    "quality_score": f"{candidate.quality_score:g}",
                    "contigs": candidate.contig_count,
                    "n50": candidate.n50,
                    "bp": candidate.fasta_bp,
                }
            )
            continue

        matches.sort(
            key=lambda match: (
                -match[0],
                candidate_sort_key(match[1], contamination_penalty),
            )
        )
        best_overlap, representative, shared_bp = matches[0]
        reason = (
            "ambiguous_merge"
            if len(matches) > 1
            else "cross_binner_overlap_suppressed"
        )
        if len(matches) > 1:
            ambiguous += 1
        alternative_by_binner[candidate.binner] += 1
        alternative_rows.append(
            {
                "sample_rank": expected_rank,
                "study": expected_study,
                "sample": expected_sample,
                "representative_bin_id": representative.bin_id,
                "representative_binner": representative.binner,
                "alternative_bin_id": candidate.bin_id,
                "alternative_binner": candidate.binner,
                "decision_qc": candidate.row["decision_qc"],
                "shared_bp": shared_bp,
                "overlap": f"{best_overlap:.6f}",
                "threshold": f"{overlap_threshold:.6f}",
                "reason": reason,
                "matching_representative_count": len(matches),
                "matching_representatives": ";".join(
                    match[1].bin_id
                    for match in sorted(
                        matches,
                        key=lambda match: (
                            -match[0],
                            candidate_sort_key(
                                match[1],
                                contamination_penalty,
                            ),
                        ),
                    )
                ),
                "representative_quality_score": f"{representative.quality_score:g}",
                "alternative_quality_score": f"{candidate.quality_score:g}",
            }
        )

    # Count only final retained bins, including overlaps within one binner.
    contig_uses = Counter(key for candidate in selected for key in candidate.contig_bp)
    for candidate, row in zip(selected, representative_rows):
        shared_bp = sum(bp for key, bp in candidate.contig_bp.items() if contig_uses[key] > 1)
        row["shared_contig_bp"] = shared_bp
        row["shared_contig_fraction"] = f"{shared_bp / candidate.fasta_bp:.6f}"

    return representative_rows, alternative_rows
