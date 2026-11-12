"""Read existing bin sets without changing their contigs or native bin boundaries."""
from pathlib import Path

from select_within_sample_representatives import fasta_contigs
from seqkit_bin_stats import is_fasta, original_bin_id, safe_component, sha256_file


def discover_bin_sets(directories: list[Path]) -> dict[str, list[Path]]:
    sources = {}
    seen = set()
    for directory in directories:
        directory = directory.expanduser().resolve()
        if not directory.is_dir():
            raise ValueError(f"bin directory does not exist: {directory}")
        if directory in seen:
            raise ValueError(f"duplicate bin directory: {directory}")
        seen.add(directory)
        base = safe_component(directory.name)
        label = base
        suffix = 2
        while label in sources:
            label = f"{base}_{suffix}"
            suffix += 1
        paths = sorted(path for path in directory.iterdir() if is_fasta(path))
        if any(any(character in str(path) for character in "\t\r\n") for path in paths):
            raise ValueError(f"bin paths cannot contain tabs or newlines: {directory}")
        if not paths:
            raise ValueError(f"no bin FASTA files directly inside: {directory}")
        identifiers = [safe_component(original_bin_id(path)) for path in paths]
        if len(set(identifiers)) != len(identifiers):
            raise ValueError(f"bin filenames produce duplicate identifiers: {directory}")
        sources[label] = paths
    return sources


def validate_bin_sets(sources: dict[str, list[Path]], assembly: Path) -> list[dict]:
    reference = fasta_contigs(assembly)
    reference_ids = {key[0] for key in reference}
    inputs = []
    for source, paths in sources.items():
        for path in paths:
            records = fasta_contigs(path)
            for key in records:
                if key[0] not in reference_ids:
                    raise ValueError(f"{path}: contig {key[0]!r} is absent from --contigs")
                if key not in reference:
                    raise ValueError(f"{path}: sequence differs from --contigs for {key[0]!r}")
            inputs.append(dict(source=source, path=str(path), sha256=sha256_file(path),
                               contigs=len(records), bp=sum(records.values())))
    return inputs
