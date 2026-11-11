#!/usr/bin/env python3
"""Launch native ReMAG with an explicit CPU/CUDA policy, before its imports."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import runpy
import shutil
import sys

from run_seeded_vamb import console_interpreter


def main() -> None:
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument("--gpu", choices=("T", "F"), required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("a ReMAG console command is required")
    executable = shutil.which(command[0])
    if executable is None:
        parser.error(f"ReMAG executable not found: {command[0]}")
    interpreter = console_interpreter(Path(executable))
    if interpreter and interpreter.resolve() != Path(sys.executable).resolve():
        os.execv(str(interpreter), [str(interpreter), str(Path(__file__).resolve()), *sys.argv[1:]])

    if args.gpu == "F":
        os.environ["CUDA_VISIBLE_DEVICES"] = ""
    import torch

    if args.gpu == "F":
        # Native ReMAG automatically tries CUDA then Apple MPS. Keep CPU explicit.
        torch.cuda.is_available = lambda: False
        mps = getattr(getattr(torch, "backends", None), "mps", None)
        if mps is not None:
            mps.is_available = lambda: False
    elif not torch.cuda.is_available():
        parser.error("--gpu T requires an available CUDA device in the ReMAG environment")
    print(f"[EukBinMiner] ReMAG device: {'cuda' if args.gpu == 'T' else 'cpu'}", flush=True)
    sys.argv = [executable, *command[1:]]
    runpy.run_path(executable, run_name="__main__")


if __name__ == "__main__":
    main()
