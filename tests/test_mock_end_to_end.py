from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
CLI = PROJECT / "bin" / "eukbinpicker"

MOCK_TOOL = r"""#!/usr/bin/env python3
import csv
import gzip
import json
import os
import sys
from pathlib import Path


name = Path(sys.argv[0]).name
if os.environ.get("MOCK_CHECK_COMPACT"):
    # Observe the real output tree while an external tool is running.
    for parent in Path.cwd().parents:
        if (parent / '0_log/run.json').exists():
            assert {p.name for p in (parent / '0_log').iterdir()} == {'run.json', 'run.log'}
            assert not list(parent.rglob('.complete.json'))
            assert not (parent / '0_log/work').exists()
            active = [p.name for p in (parent / '.work').iterdir() if p.is_dir() and p.name not in {'inputs', 'cache'}]
            assert len(active) <= 1, active
            break
args = sys.argv[1:]


def value(flag):
    return args[args.index(flag) + 1]


def records(path):
    observed = []
    identifier = None
    sequence = []
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        lines = handle.read().splitlines()
    for raw in lines:
        if raw.startswith(">"):
            if identifier is not None:
                observed.append((identifier, "".join(sequence)))
            identifier = raw[1:].split()[0]
            sequence = []
        elif raw:
            sequence.append(raw.strip())
    if identifier is not None:
        observed.append((identifier, "".join(sequence)))
    return observed


def write_fasta(path, observed):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "".join(f">{identifier}\n{sequence}\n" for identifier, sequence in observed),
        encoding="utf-8",
    )


def binner_records(source, binner):
    observed = records(source)
    if os.environ.get("MOCK_BIN_MODE") != "half":
        return observed
    if binner == "remag":
        return [row for row in observed if row[0] in {"c1", "c2"}]
    if binner == "metabat2":
        return [row for row in observed if row[0] in {"c1", "c3"}]
    return observed


if name == "conda":
    if (
        len(args) < 5
        or args[:2] != ["run", "--no-capture-output"]
        or args[2] not in {"-n", "-p"}
    ):
        raise SystemExit(f"unexpected mock conda invocation: {args!r}")
    os.execvp(args[4], args[4:])
elif name == "seqkit":
    input_list = Path(value("-X"))
    print("file\tnum_seqs\tsum_len\tN50\tmax_len")
    for raw_path in input_list.read_text(encoding="utf-8").splitlines():
        observed = records(raw_path)
        lengths = sorted((len(sequence) for _, sequence in observed), reverse=True)
        total = sum(lengths)
        running = 0
        n50 = 0
        for length in lengths:
            running += length
            if running >= total / 2:
                n50 = length
                break
        print(f"{Path(raw_path).resolve()}\t{len(lengths)}\t{total}\t{n50}\t{max(lengths)}")
elif name == "remag":
    if "--version" in args:
        print("remag 0.4.4")
    else:
        expected_seed = os.environ.get("MOCK_EXPECTED_SEED")
        if expected_seed is not None and value("--random-seed") != expected_seed:
            raise SystemExit("unexpected ReMAG random seed")
        if int(value("--min-bin-size")) < 1:
            raise SystemExit("invalid REMAG export threshold")
        source = args[0]
        out = Path(value("-o"))
        coverage_files = [
            args[index + 1]
            for index, argument in enumerate(args)
            if argument == "-c"
        ]
        expected_coverage_count = os.environ.get(
            "MOCK_EXPECTED_REMAG_COVERAGE_COUNT"
        )
        if (
            expected_coverage_count is not None
            and len(coverage_files) != int(expected_coverage_count)
        ):
            raise SystemExit("unexpected ReMAG coverage file count")
        out.mkdir(parents=True, exist_ok=True)
        (out / "mock_coverage_inputs.json").write_text(
            json.dumps(
                [
                    Path(path).read_text(encoding="utf-8").splitlines()
                    for path in coverage_files
                ]
            ),
            encoding="utf-8",
        )
        write_fasta(out / "bins" / "remag-bin.fa", binner_records(source, "remag"))
elif name == "metabat2":
    if "-s" in args:
        raise SystemExit("MetaBAT2 native bin minimum overridden")
    if "--help" in args:
        print("MetaBAT 2.18")
    else:
        expected_seed = os.environ.get("MOCK_EXPECTED_SEED")
        if expected_seed is not None and value("--seed") != expected_seed:
            raise SystemExit("unexpected MetaBAT2 random seed")
        table = Path(value("-a")).read_text().splitlines()
        if "--cvExt" in args:
            if "-var" in table[0] or "contigLen" in table[0] or "NA" in table[1:]:
                raise SystemExit("invalid native mean-only table")
        elif "-var" not in table[0] or any("NA" in row.split("\t") for row in table[1:]):
            raise SystemExit("invalid observed variance table")
        source = value("-i")
        prefix = Path(value("-o"))
        if "--onlyLabel" not in args:
            raise SystemExit("expected membership-only MetaBAT2 output")
        Path(str(prefix) + ".1").write_text("".join(
            identifier + " total_depth=10\n" for identifier, _ in binner_records(source, "metabat2")))
elif name in {
    "concoct",
    "cut_up_fasta.py",
    "merge_cutup_clustering.py",
    "extract_fasta_bins.py",
}:
    if name == "concoct" and "--version" in args:
        print("CONCOCT 1.1.0")
    elif os.environ.get("MOCK_CONCOCT_SUCCESS") == "1":
        if name == "cut_up_fasta.py":
            source = args[0]
            bed = Path(value("-b"))
            bed.parent.mkdir(parents=True, exist_ok=True)
            bed_rows = []
            for identifier, sequence in records(source):
                bed_rows.append(
                    f"{identifier}\t0\t{len(sequence)}\t{identifier}\n"
                )
                print(f">{identifier}\n{sequence}")
            bed.write_text("".join(bed_rows), encoding="utf-8")
        elif name == "concoct":
            coverage = Path(value("--coverage_file"))
            header = coverage.read_text(encoding="utf-8").splitlines()[0]
            expected_columns = os.environ.get("MOCK_EXPECTED_COVERAGE_COLUMNS")
            if expected_columns is not None:
                observed_columns = header.split("\t")[1:]
                if observed_columns != expected_columns.split(","):
                    raise SystemExit("unexpected CONCOCT coverage columns")
            if any(flag in args for flag in ("--no_cov_normalization", "--no_total_coverage", "--clusters")):
                raise SystemExit("CONCOCT native defaults were overridden")
            out = Path(value("-b"))
            out.mkdir(parents=True, exist_ok=True)
            threshold = value("--length_threshold")
            identifiers = [
                row.split("\t")[0]
                for row in coverage.read_text(encoding="utf-8").splitlines()[1:]
            ]
            (out / f"clustering_gt{threshold}.csv").write_text(
                "contig_id,cluster_id\n"
                + "".join(f"{identifier},0\n" for identifier in identifiers),
                encoding="utf-8",
            )
        elif name == "merge_cutup_clustering.py":
            print(Path(args[0]).read_text(encoding="utf-8"), end="")
        else:
            source = args[0]
            output_path = Path(value("--output_path"))
            write_fasta(output_path / "0.fa", records(source))
    else:
        print("intentional mock CONCOCT failure", file=sys.stderr)
        raise SystemExit(17)
elif name == "eukcc":
    if "--help" in args:
        print("EukCC 2.1.3")
    else:
        out = Path(value("--out"))
        out.mkdir(parents=True, exist_ok=True)
        bin_path = args[-1]
        if os.environ.get("MOCK_FAIL_QC_SOURCE") and os.environ["MOCK_FAIL_QC_SOURCE"] in bin_path:
            print("synthetic QC failure", file=sys.stderr)
            raise SystemExit(23)
        if os.environ.get("MOCK_UNASSESSABLE"):
            print("No placement marker genes found", file=sys.stderr)
            raise SystemExit(1)
        with gzip.open(out / "savestate.json.gz", "wt") as handle:
            json.dump(json.dumps({"marker_set": {"profiles": ["marker1"]}, "scmg_data": [{"query": "marker1", "target": "protein1"}]}), handle)
        if "metabat2" in bin_path:
            completeness, contamination = "91", "1"
        else:
            completeness, contamination = "95", "2"
        (out / "quality.csv").write_text(
            "completeness,contamination,ncbi_lng\n"
            f"{completeness},{contamination},Eukaryota;4751\n",
            encoding="utf-8",
        )
elif name == "busco":
    if "--version" in args:
        print("BUSCO 6")
    else:
        Path("busco_cwd_marker.log").write_text("native BUSCO starts a log in cwd\n")
        out = Path(value("--out_path"))
        bin_id = value("--out")
        (out / 'full_table.tsv').write_text('Busco id\tStatus\tSequence\nmarker1\tComplete\tprotein1\n')
        out.mkdir(parents=True, exist_ok=True)
        bin_path = value("--in")
        if "metabat2" in bin_path:
            complete, duplicated = "98.0", "1.0"
        else:
            complete, duplicated = "80.0", "2.0"
        single = str(float(complete) - float(duplicated))
        (out / f"short_summary.specific.{bin_id}.txt").write_text(
            f"C:{complete}%[S:{single}%,D:{duplicated}%],"
            f"F:1.0%,M:{99.0 - float(complete)}%,n:100\n",
            encoding="utf-8",
        )
elif name == "vamb":
    if "--version" in args:
        print("Vamb 5.0.4")
    else:
        expected_seed = os.environ.get("MOCK_EXPECTED_SEED")
        if expected_seed is not None and value("--seed") != expected_seed:
            raise SystemExit("unexpected VAMB random seed")
        mode = args[1]
        out = Path(value("--outdir"))
        source = value("--fasta")
        write_fasta(out / "bins" / f"{mode}-bin.fna", records(source))
        (out / "clusters.tsv").write_text("cluster\tcontig\n", encoding="utf-8")
        (out / "log.txt").write_text("mock VAMB complete\n", encoding="utf-8")
else:
    print(f"unsupported mock tool: {name}", file=sys.stderr)
    raise SystemExit(99)
"""


class MockEndToEndTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.tools = self.root / "tools"
        self.tools.mkdir()
        for name in (
            "conda",
            "seqkit",
            "remag",
            "metabat2",
            "concoct",
            "cut_up_fasta.py",
            "merge_cutup_clustering.py",
            "extract_fasta_bins.py",
            "eukcc",
            "busco",
            "vamb",
        ):
            tool = self.tools / name
            tool.write_text(MOCK_TOOL, encoding="utf-8")
            tool.chmod(0o755)

        self.fasta = self.root / "contigs.fa"
        sequences = {
            "c1": "ACGT" * 25,
            "c2": "A" * 99 + "C",
            "c3": "A" * 98 + "CC",
        }
        self.fasta.write_text(
            "".join(
                f">{identifier}\n{sequence}\n"
                for identifier, sequence in sequences.items()
            ),
            encoding="utf-8",
        )
        self.depth = self.root / "depth.tsv"
        self.depth.write_text(
            "contig_id\tdepth\nc1\t3\nc2\t2\nc3\t1\n",
            encoding="utf-8",
        )
        self.eukcc_db = self.root / "eukcc_db"
        self.eukcc_db.mkdir()
        (self.eukcc_db / "db.info").write_text("mock\n", encoding="utf-8")
        self.busco_db = self.root / "fungi_odb12"
        (self.busco_db / "hmms").mkdir(parents=True)
        (self.busco_db / "dataset.cfg").write_text("name=fungi\n", encoding="utf-8")
        self.environment = dict(os.environ)
        self.environment["PATH"] = (
            str(self.tools) + os.pathsep + self.environment.get("PATH", "")
        )
        mock_modules = self.root / "mock_python_modules"
        mock_modules.mkdir()
        (mock_modules / "numpy.py").write_text(
            "class _Random:\n"
            "    def seed(self, value):\n"
            "        self.value = value\n"
            "random = _Random()\n",
            encoding="utf-8",
        )
        (mock_modules / "torch.py").write_text(
            "def manual_seed(value):\n"
            "    return value\n"
            "class _Cuda:\n"
            "    @staticmethod\n"
            "    def is_available():\n"
            "        return False\n"
            "cuda = _Cuda()\n",
            encoding="utf-8",
        )
        self.environment["PYTHONPATH"] = os.pathsep.join(
            filter(
                None,
                (
                    str(mock_modules),
                    self.environment.get("PYTHONPATH", ""),
                ),
            )
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def run_cli(
        self,
        *arguments: str,
        mock_bin_mode: str | None = None,
        expected_seed: int | None = None,
        expected_coverage_columns: tuple[str, ...] | None = None,
        concoct_success: bool = False,
    ) -> subprocess.CompletedProcess[str]:
        environment = dict(self.environment)
        if mock_bin_mode is not None:
            environment["MOCK_BIN_MODE"] = mock_bin_mode
        if expected_seed is not None:
            environment["MOCK_EXPECTED_SEED"] = str(expected_seed)
        if expected_coverage_columns is not None:
            environment["MOCK_EXPECTED_COVERAGE_COLUMNS"] = ",".join(
                expected_coverage_columns
            )
            environment["MOCK_EXPECTED_REMAG_COVERAGE_COUNT"] = str(
                len(expected_coverage_columns)
            )
        if concoct_success:
            environment["MOCK_CONCOCT_SUCCESS"] = "1"
        result = subprocess.run(
            [sys.executable, str(CLI), *arguments],
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=environment,
        )

        if result.returncode and "--output" in arguments:
            output = Path(arguments[arguments.index("--output") + 1])
            for log in sorted(output.rglob("*.log")):
                result.stderr += f"\n--- {log.relative_to(output)} ---\n"
                result.stderr += log.read_text(errors="replace")[-4000:]
        return result

    def base_run(self, output: Path) -> list[str]:
        return [
            "--sample",
            "sample1",
            "--contigs",
            str(self.fasta),
            "--coverage",
            str(self.depth),
            "--output",
            str(output),
            "--min-retain-size",
            "1",
            "--min-selection-size",
            "1",
            "--min-input-length",
            "1",
            "--threads",
            "1",
        ]


if __name__ == "__main__":
    unittest.main()
