#!/usr/bin/env python3
"""Install one user-facing environment with private, conflicting tool runtimes."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

VERSION = "0.40"
SOURCE = Path(__file__).resolve().parents[1]


def run(command: list[str]) -> None:
    print("[install] " + " ".join(command), flush=True)
    subprocess.run(command, check=True)


def install(prefix: Path, conda: str, gpu: str) -> None:
    prefix = prefix.expanduser().resolve()
    conda = str(Path(shutil.which(conda) or conda).resolve())
    marker = prefix / ".ebm-installing.json"
    settings = {"version": VERSION, "gpu": gpu, "conda": conda}
    if prefix.exists():
        if not marker.is_file() or json.loads(marker.read_text()) != settings:
            raise ValueError(f"use a new installation directory: {prefix}")
    else:
        prefix.mkdir(parents=True)
        marker.write_text(json.dumps(settings))

    package = prefix / "share/eukbinminer"
    logs = prefix / "share/eukbinminer/install"
    logs.mkdir(parents=True, exist_ok=True)
    cache = package / "install-cache"
    for name in ("conda", "pip", "xdg", "tmp"):
        (cache / name).mkdir(parents=True, exist_ok=True)
    # Keep downloads and build files off quota-limited home directories.
    os.environ.update(CONDA_PKGS_DIRS=str(cache / "conda"),
                      PIP_CACHE_DIR=str(cache / "pip"),
                      XDG_CACHE_HOME=str(cache / "xdg"), TMPDIR=str(cache / "tmp"),
                      CONDA_REGISTER_ENVS="false")
    completed = logs / "steps.json"
    steps = set(json.loads(completed.read_text())) if completed.exists() else set()
    if not (prefix / "bin/python").is_file():
        steps.discard("core")

    def step(name: str, commands: list[list[str]]) -> None:
        if name in steps:
            return
        for command in commands:
            run(command)
        steps.add(name)
        completed.write_text(json.dumps(sorted(steps), indent=2) + "\n")

    def create(target: Path, specs: list[str]) -> list[str]:
        # A failed transaction can leave package records without their files.
        records = [json.loads(p.read_text()) for p in (target / "conda-meta").glob("*.json")]
        action = "install" if records else "create"
        repair = []
        if records:
            repair = [f"{r['name']}={r['version']}={r['build']}" for r in records]
        return [conda, action, "--yes", "--prefix", str(target), "--override-channels",
                "-c", "conda-forge", "-c", "bioconda",
                *(["--force-reinstall"] if records else []), *specs, *repair]

    step("core", [create(prefix, ["python=3.12", "pip", "seqkit", "miniprot=0.18"])])
    python = str(prefix / "bin/python")
    torch_index = "cu124" if gpu == "T" else "cpu"
    step("pytorch", [[python, "-m", "pip", "install", "torch==2.6.0", "--index-url",
                      f"https://download.pytorch.org/whl/{torch_index}"]])
    step("neural_binners", [[python, "-m", "pip", "install", "vamb==5.0.4", "remag==0.4.4"],
                            [python, "-m", "pip", "check"]])

    specs = {
        "concoct": ["concoct=1.1.0", "biopython", "samtools>=1.17"],
        "metabat2": ["metabat2=2.18_23_gc869c52", "samtools>=1.17"],
        "eukcc": ["eukcc=2.1.3", "taxonkit", "csvtk"],
        "busco": ["python=3.10", "busco=6"],
    }
    environments = {"CONDA_EXE": conda}
    for tool, packages in specs.items():
        target = package / "environments" / tool
        step(tool, [create(target, packages)])
        environments[tool.upper() + "_ENV"] = str(target)

    # Launch every executable in its own environment, including BAM helpers.
    probes = {
        "concoct": [["concoct", "--version"], ["samtools", "--version"],
                    ["concoct_coverage_table.py", "--help"]],
        "metabat2": [["metabat2", "--help"], ["samtools", "--version"],
                     ["jgi_summarize_bam_contig_depths", "--help"]],
        "eukcc": [["eukcc", "--version"]],
        "busco": [["busco", "--version"]],
    }
    for tool, commands in probes.items():
        for command in commands:
            result = subprocess.run([conda, "run", "--no-capture-output", "-p",
                environments[tool.upper() + "_ENV"], *command], capture_output=True, text=True)
            (logs / (command[0] + "-" + tool + ".txt")).write_text(result.stdout + result.stderr)
            # MetaBAT help exits 1 in some builds; missing libraries must still fail.
            valid_help = command[0] in {"metabat2", "jgi_summarize_bam_contig_depths"}
            if result.returncode and not (valid_help and result.returncode == 1
                    and "Usage" in result.stdout + result.stderr):
                raise RuntimeError(f"{tool} validation failed: {result.stderr[-1000:]}")
    for command in [[python, "-m", "pip", "check"],
                    [str(prefix / "bin/remag"), "--version"],
                    [str(prefix / "bin/vamb"), "--help"],
                    [str(prefix / "bin/seqkit"), "version"]]:
        run(command)
    run([python, "-c", "import torch; assert torch.version.cuda is not None" if gpu == "T"
         else "import torch; assert torch.version.cuda is None"])

    for name in ("scripts", "envs"):
        shutil.copytree(SOURCE / name, package / name, dirs_exist_ok=True)
    for name in ("README.md", "LICENSE", "CITATION.cff"):
        shutil.copy2(SOURCE / name, package / name)
    (package / "installation.json").write_text(json.dumps({
        **settings, "environments": environments,
    }, indent=2) + "\n")
    shutil.copy2(SOURCE / "bin/eukbinminer", prefix / "bin/eukbinminer")
    (prefix / "bin/eukbinminer").chmod(0o755)
    # Activation selects the exact interpreter and native executables together.
    run([conda, "run", "--no-capture-output", "-p", str(prefix), "eukbinminer", "--version"])
    marker.unlink()
    shutil.rmtree(cache)
    print(f"\nInstalled EukBinMiner {VERSION}. Activate with:\nconda activate {prefix}")


def main() -> int:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--help", action="help")
    parser.add_argument("--prefix", type=Path, required=True)
    parser.add_argument("--gpu", choices=("T", "F"), default="F")
    parser.add_argument("--conda", default=os.environ.get("CONDA_EXE", "conda"), help=argparse.SUPPRESS)
    args = parser.parse_args()
    if platform.system() != "Linux" or platform.machine() not in {"x86_64", "AMD64"}:
        parser.error("the bundled installer requires Linux x86-64")
    try:
        install(args.prefix, args.conda, args.gpu)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"Installation failed: {error}\nFix the reported error and repeat the same command to resume.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
