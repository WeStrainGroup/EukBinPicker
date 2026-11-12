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

VERSION = "0.3"
SOURCE = Path(__file__).resolve().parents[1]


def run(command: list[str]) -> None:
    print("[install] " + " ".join(command), flush=True)
    subprocess.run(command, check=True)


def install(prefix: Path, conda: str) -> None:
    prefix = prefix.expanduser().resolve()
    conda = str(Path(shutil.which(conda) or conda).resolve())
    state = prefix.with_name("." + prefix.name + ".ebp-install")
    marker = state / "settings.json"
    settings = {"version": VERSION, "runtime": "cu124", "conda": conda}
    package = prefix / "share/eukbinpicker"
    installed = package / "installation.json"
    record = marker if marker.is_file() else installed
    if prefix.exists() or state.exists():
        if not record.is_file() or any(json.loads(record.read_text()).get(k) != v
                                       for k, v in settings.items()):
            raise ValueError(f"use a new installation directory: {prefix}")
    state.mkdir(parents=True, exist_ok=True)
    marker.write_text(json.dumps(settings))
    logs = state / "logs"
    logs.mkdir(exist_ok=True)
    cache = state / "cache"
    for name in ("conda", "pip", "xdg", "tmp"):
        (cache / name).mkdir(parents=True, exist_ok=True)
    # Temporary files stay outside the environment being created.
    os.environ.update(CONDA_PKGS_DIRS=str(cache / "conda"),
                      PIP_CACHE_DIR=str(cache / "pip"),
                      XDG_CACHE_HOME=str(cache / "xdg"), TMPDIR=str(cache / "tmp"),
                      CONDA_REGISTER_ENVS="false")
    completed = logs / "steps.json"
    previous = completed if completed.exists() else package / "install/steps.json"
    steps = set(json.loads(previous.read_text())) if previous.exists() else set()
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
    step("pytorch", [[python, "-m", "pip", "install", "torch==2.6.0+cu124", "--index-url",
                      "https://download.pytorch.org/whl/cu124"]])
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
                    ["cut_up_fasta.py", "--help"],
                    ["concoct_coverage_table.py", "--help"]],
        "metabat2": [["metabat2", "--help"], ["samtools", "--version"],
                     ["jgi_summarize_bam_contig_depths", "--help"]],
        "eukcc": [["eukcc", "--version"]],
        "busco": [["busco", "--version"]],
    }
    for tool, commands in probes.items():
        for command in commands:
            result = subprocess.run([conda, "run", "--no-capture-output", "-p",
                environments[tool.upper() + "_ENV"], *command], capture_output=True, text=True,
                cwd=logs)
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
    run([python, "-c", "import torch; assert torch.version.cuda == '12.4'; "
         "assert torch.ones(2, device='cpu').sum().item() == 2"])

    for name in ("scripts", "envs", "docs"):
        shutil.copytree(SOURCE / name, package / name, dirs_exist_ok=True)
    for name in ("README.md", "LICENSE", "CITATION.cff"):
        shutil.copy2(SOURCE / name, package / name)
    (package / "installation.json").write_text(json.dumps({
        **settings, "environments": environments,
    }, indent=2) + "\n")
    shutil.copy2(SOURCE / "bin/eukbinpicker", prefix / "bin/eukbinpicker")
    (prefix / "bin/eukbinpicker").chmod(0o755)
    # Activation selects the exact interpreter and native executables together.
    run([conda, "run", "--no-capture-output", "-p", str(prefix), "eukbinpicker", "--version"])
    shutil.copytree(logs, package / "install", dirs_exist_ok=True)
    shutil.rmtree(state)
    print(f"\nInstalled EukBinPicker {VERSION}.")


def main() -> int:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--help", action="help")
    parser.add_argument("--prefix", type=Path, required=True)
    parser.add_argument("--conda", default=os.environ.get("CONDA_EXE", "conda"), help=argparse.SUPPRESS)
    args = parser.parse_args()
    if platform.system() != "Linux" or platform.machine() not in {"x86_64", "AMD64"}:
        parser.error("the bundled installer requires Linux x86-64")
    try:
        install(args.prefix, args.conda)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"Installation failed: {error}\nFix the reported error and repeat the same command to resume.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
