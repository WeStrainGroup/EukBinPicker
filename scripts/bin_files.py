"""FASTA filename handling and checksums shared by import and workflow."""
import hashlib
import re
from pathlib import Path

FASTA_SUFFIXES = ('.fa', '.fna', '.fasta')

def safe_component(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._")
    if not cleaned:
        raise ValueError(f"cannot form a safe identifier from {value!r}")
    return cleaned


def original_bin_id(path: Path) -> str:
    name = path.name
    lower = name.lower()
    if lower.endswith(".gz"):
        name = name[:-3]
        lower = lower[:-3]
    for suffix in FASTA_SUFFIXES:
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    raise ValueError(f"not a recognized FASTA path: {path}")


def is_fasta(path: Path) -> bool:
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    return path.is_file() and any(name.endswith(suffix) for suffix in FASTA_SUFFIXES)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
