#!/usr/bin/env python3
"""Validate indexed BAMs and calculate native JGI contig coverage."""
import hashlib
import os
from pathlib import Path
import subprocess


def tool_command(command, environment=None, conda="conda"):
    if not environment:
        return command
    selector = "-p" if "/" in environment else "-n"
    return [conda, "run", "--no-capture-output", selector, environment, *command]


def checked_output(command):
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode or "[E::" in result.stderr or "reverting to slow method" in result.stderr:
        raise ValueError(f"{command[0]} failed: {result.stderr.strip()}")
    return result.stdout


def file_identity(path):
    state = path.stat()
    return {"path": str(path), "bytes": state.st_size, "mtime_ns": state.st_mtime_ns}


def validate_bams(paths, lengths, fasta_handle, environment=None, conda="conda"):
    paths = [Path(p).expanduser().resolve() for p in paths]
    if len(set(paths)) != len(paths) or len({p.name for p in paths}) != len(paths):
        raise ValueError("BAM paths and basenames must be unique (one BAM per sample)")
    identities = []
    reference = None
    required_m5 = {}
    missing_m5 = False
    for path in paths:
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"BAM does not exist or is empty: {path}")
        if path.suffix != ".bam" or any(c in str(path) for c in "\t\r\n"):
            raise ValueError(f"expected a .bam path without tabs/newlines: {path}")
        indices = list(dict.fromkeys(p for p in [Path(str(path)+".bai"), path.with_suffix(".bai"),
                                                 Path(str(path)+".csi"), path.with_suffix(".csi")]
                                    if p.is_file() and p.stat().st_size))
        if not indices:
            raise ValueError(f"BAM requires a nonempty BAI or CSI index: {path}")
        def samtools(*args):
            return checked_output(tool_command(["samtools", *args, str(path)], environment, conda))
        samtools("quickcheck", "-v")
        header = samtools("view", "-H")
        dictionary = []
        md5 = {}
        sorted_by_coordinate = False
        for line in header.splitlines():
            fields = line.split("\t")
            tags = dict(field.split(":", 1) for field in fields[1:] if ":" in field)
            if fields[0] == "@HD":
                sorted_by_coordinate = tags.get("SO") == "coordinate"
            if fields[0] == "@SQ":
                if "SN" not in tags or "LN" not in tags:
                    raise ValueError(f"invalid BAM reference header: {path}")
                dictionary.append((tags["SN"], int(tags["LN"])))
                if "M5" in tags:
                    md5[tags["SN"]] = tags["M5"].lower()
        if not sorted_by_coordinate:
            raise ValueError(f"BAM must declare coordinate sorting (@HD SO:coordinate): {path}")
        if not dictionary or len(dict(dictionary)) != len(dictionary):
            raise ValueError(f"BAM has missing or duplicate reference IDs: {path}")
        if reference is not None and dictionary != reference:
            raise ValueError(f"BAM reference dictionaries differ: {path}")
        reference = dictionary
        refs = dict(dictionary)
        for name, length in lengths.items():
            if refs.get(name) != length:
                raise ValueError(f"FASTA/BAM reference ID or length mismatch: {name} in {path}")
            if name in md5:
                if name in required_m5 and required_m5[name] != md5[name]:
                    raise ValueError(f"BAM reference M5 conflict: {name}")
                required_m5[name] = md5[name]
            else:
                missing_m5 = True
        # Exercise the index; a stale or mismatching dictionary must not be accepted.
        indexed = samtools("idxstats")
        indexed_refs = [(f[0], int(f[1])) for row in indexed.splitlines()
                        if (f := row.split("\t"))[0] != "*"]
        if indexed_refs != dictionary:
            raise ValueError(f"BAM index reference dictionary differs: {path}")
        identities.append({**file_identity(path), "header_sha256": hashlib.sha256(header.encode()).hexdigest(),
                           "indices": [{**file_identity(i), "sha256": hashlib.sha256(i.read_bytes()).hexdigest()}
                                       for i in indices]})
    # SAM M5 is the MD5 of the uppercase sequence without whitespace.
    name = None
    digest = None
    def finish():
        if name in required_m5 and digest.hexdigest() != required_m5[name]:
            raise ValueError(f"FASTA/BAM reference sequence M5 mismatch: {name}")
    with fasta_handle as handle:
        for line in handle:
            if line.startswith(">"):
                finish()
                name = line[1:].split()[0]
                digest = hashlib.md5()
            elif digest is not None:
                digest.update("".join(line.split()).upper().encode("ascii"))
        finish()
    return {"bam_inputs": identities, "reference_sequence_check":
            "partial_or_absent_M5; matching_reference_provenance_required" if missing_m5 else "M5_verified"}


def calculate_depth(paths, destination, threads, environment=None, conda="conda"):
    destination = Path(destination)
    receipt = destination.with_suffix(".sha256")
    if receipt.exists():
        if not destination.is_file() or hashlib.sha256(destination.read_bytes()).hexdigest() != receipt.read_text().strip():
            raise ValueError(f"cached BAM depth was modified: {destination}")
        return
    temporary = destination.with_suffix(".tmp.tsv")
    temporary.unlink(missing_ok=True)
    command = tool_command(["jgi_summarize_bam_contig_depths", "--outputDepth", str(temporary),
                            "--percentIdentity", "97", "--minMapQual", "0", "--maxEdgeBases", "75",
                            *[str(p) for p in paths]], environment, conda)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    with destination.with_suffix(".log").open("w") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=env)
    if result.returncode or not temporary.is_file() or not temporary.stat().st_size:
        raise ValueError(f"JGI coverage failed; see {destination.with_suffix('.log')}")
    temporary.replace(destination)
    receipt.write_text(hashlib.sha256(destination.read_bytes()).hexdigest()+"\n")
