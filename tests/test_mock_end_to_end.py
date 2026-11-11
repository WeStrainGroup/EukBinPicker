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
import json
import os
import sys
from pathlib import Path


name = Path(sys.argv[0]).name
args = sys.argv[1:]


def value(flag):
    return args[args.index(flag) + 1]


def records(path):
    observed = []
    identifier = None
    sequence = []
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
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
    if "--min-bin-size" in args or "-m" in args:
        raise SystemExit("ReMAG native defaults overridden")
    if "--version" in args:
        print("remag 0.4.4")
    else:
        expected_seed = os.environ.get("MOCK_EXPECTED_SEED")
        if expected_seed is not None and value("--random-seed") != expected_seed:
            raise SystemExit("unexpected ReMAG random seed")
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
        write_fasta(
            Path(str(prefix) + ".1.fa"),
            binner_records(source, "metabat2"),
        )
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
        out = Path(value("--out_path"))
        bin_id = value("--out")
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
            "--min-retain-length",
            "1",
            "--min-selection-length",
            "1",
            "--min-input-length",
            "1",
            "--threads",
            "1",
        ]

    def test_configured_split_conda_environments_reach_backends(self) -> None:
        output = self.root / "split-conda-environments"
        configured = {
            "CONDA_EXE": str(self.tools / "conda"),
            "REMAG_ENV": str(self.root / "envs" / "mock-remag"),
            "CONCOCT_ENV": str(self.root / "envs" / "mock-concoct"),
            "METABAT2_ENV": str(self.root / "envs" / "mock-metabat2"),
            "EUKCC_ENV": str(self.root / "envs" / "mock-eukcc"),
            "BUSCO_ENV": str(self.root / "envs" / "mock-busco"),
        }
        for suffix, value in configured.items():
            self.environment[f"EUKBINPICKER_{suffix}"] = value
        completed = self.run_cli(
            *self.base_run(output),
            "--binners",
            "remag",
            "--qc-tools",
            "all",
            "--qc-decision",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
            "--busco-db",
            str(self.busco_db),
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        backend = (output / "0_log/config" / "backend.env").read_text(encoding="utf-8")
        for suffix, value in configured.items():
            self.assertIn(f"{suffix}={value}\n", backend)

    def test_partial_binner_failure_busco_decision_and_thresholds(self) -> None:
        output = self.root / "busco-threshold-050"
        completed = self.run_cli(
            *self.base_run(output),
            "--binners",
            "remag,concoct,metabat2",
            "--qc-tools",
            "all",
            "--qc-decision",
            "busco",
            "--eukcc-db",
            str(self.eukcc_db),
            "--busco-db",
            str(self.busco_db),
            "--min-overlap",
            "0.50",
            mock_bin_mode="half",
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        summary = json.loads((output / "0_log/summary.json").read_text(encoding="utf-8"))
        self.assertEqual(summary["tool"], "EukBinPicker")
        self.assertEqual(summary["completed_binners"], ["remag", "metabat2"])
        self.assertEqual(summary["representatives"], 1)
        self.assertEqual(summary["alternatives"], 1)
        with (output / "0_log/selection" / "representatives.tsv").open(
            encoding="utf-8", newline=""
        ) as handle:
            representative = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(representative["decision_qc"], "busco")
        self.assertEqual(
            representative["contamination_source"],
            "busco_duplicated",
        )
        self.assertEqual(representative["qc_pass"], "yes")
        self.assertEqual(representative["contamination"], "1")
        self.assertEqual(representative["quality_score"], "93")
        self.assertNotIn("redundancy", representative)
        self.assertNotIn("decision_quality", representative)
        self.assertNotIn("eukcc_quality", representative)
        self.assertEqual(representative["selected_binner"], "metabat2")
        normalized_depth = (
            output / "0_log/work" / "input" / "depth.jgi.tsv"
        ).read_text(encoding="utf-8").splitlines()[0]
        self.assertEqual(
            normalized_depth,
            "contigName\tcontigLen\ttotalAvgDepth\t"
            "sample1.bam\tsample1.bam-var",
        )
        mags = sorted((output / "3_selection/bins").glob("*.fa"))
        self.assertEqual(len(mags), 1)

        strict_output = self.root / "busco-threshold-051"
        strict = self.run_cli(
            *self.base_run(strict_output),
            "--binners",
            "remag,metabat2",
            "--qc-tools",
            "busco",
            "--qc-decision",
            "busco",
            "--busco-db",
            str(self.busco_db),
            "--min-overlap",
            "0.51",
            mock_bin_mode="half",
        )
        self.assertEqual(strict.returncode, 0, strict.stderr)
        strict_summary = json.loads(
            (strict_output / "0_log/summary.json").read_text(encoding="utf-8")
        )
        self.assertEqual(strict_summary["representatives"], 2)
        self.assertEqual(strict_summary["alternatives"], 0)

    def test_repeated_mapping_depths_reach_all_four_binners(self) -> None:
        alpha = self.root / "alpha.tsv"
        beta = self.root / "beta.tsv"
        alpha.write_text(
            "contig_id\tdepth\nc1\t3\nc2\t2\nc3\t1\n",
            encoding="utf-8",
        )
        beta.write_text(
            "contig_id\tdepth\nc1\t30\nc2\t20\nc3\t10\n",
            encoding="utf-8",
        )
        output = self.root / "multi-coverage-all-binners"
        completed = self.run_cli(
            "--sample",
            "sample1",
            "--contigs",
            str(self.fasta),
            "--coverage",
            str(alpha),
            "--coverage",
            str(beta),
            "--output",
            str(output),
            "--binners",
            "remag,concoct,metabat2,vamb",
            "--qc-tools",
            "eukcc",
            "--qc-decision",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
            "--min-retain-length",
            "1",
            "--min-selection-length",
            "1",
            "--min-input-length",
            "1",
            "--threads",
            "1",
            expected_coverage_columns=("alpha", "beta"),
            concoct_success=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)

        summary = json.loads((output / "0_log/summary.json").read_text(encoding="utf-8"))
        failed_diagnostics = "\n".join(
            path.read_text(encoding="utf-8", errors="replace")
            for path in (output / "0_log/work/results" / "standalone" / "sample1" / "failed").glob(
                "**/*"
            )
            if path.is_file()
        )
        self.assertEqual(
            summary["completed_binners"],
            ["remag", "concoct", "metabat2", "vamb"],
            completed.stderr
            + (output / "0_log/logs" / "binning.classic.log").read_text(
                encoding="utf-8"
            )
            + failed_diagnostics,
        )
        run_metadata = json.loads(
            (output / "0_log" / "run.json").read_text(encoding="utf-8")
        )
        self.assertEqual(run_metadata["depth_source_count"], 2)
        self.assertEqual(run_metadata["coverage_columns"], ["alpha", "beta"])
        self.assertEqual(run_metadata["coverage_column_count"], 2)

        normalized = (
            output / "0_log/work" / "input" / "depth.jgi.tsv"
        ).read_text(encoding="utf-8").splitlines()
        self.assertEqual(
            normalized[0],
            "contigName\tcontigLen\ttotalAvgDepth\talpha\talpha-var\t"
            "beta\tbeta-var",
        )
        self.assertEqual(normalized[1], "c1\t100\t33\t3\tNA\t30\tNA")

        preflight = output / "0_log/work/results" / "standalone" / "sample1" / "preflight"
        self.assertEqual(
            (preflight / "union.concoct.coverage.tsv")
            .read_text(encoding="utf-8")
            .splitlines()[0],
            "contig_id\talpha\tbeta",
        )
        preflight_metadata = json.loads(
            (preflight / "preflight.json").read_text(encoding="utf-8")
        )
        self.assertEqual(preflight_metadata["depth"]["coverage_column_count"], 2)

        branch_root = (
            output
            / "0_log/work/results"
            / "standalone"
            / "sample1"
            / "results"
            / "union"
        )
        remag_inputs = json.loads(
            (branch_root / "remag" / "mock_coverage_inputs.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            remag_inputs,
            [
                ["c1\t3", "c2\t2", "c3\t1"],
                ["c1\t30", "c2\t20", "c3\t10"],
            ],
        )
        self.assertEqual(
            (branch_root / "concoct" / "fragment_coverage.tsv")
            .read_text(encoding="utf-8")
            .splitlines()[0],
            "contig_id\talpha\tbeta",
        )
        self.assertEqual(
            (output / "0_log/work/inputs" / "r000001" / "abundance.tsv")
            .read_text(encoding="utf-8")
            .splitlines()[0],
            "contigname\talpha\tbeta",
        )
        vamb_metadata = json.loads(
            (output / "0_log/work/inputs" / "r000001" / "input.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(vamb_metadata["coverage_columns"], ["alpha", "beta"])
        self.assertEqual(vamb_metadata["coverage_column_count"], 2)

    def test_single_depth_concoct_keeps_native_defaults(self) -> None:
        output = self.root / "single-depth-concoct"
        completed = self.run_cli(
            *self.base_run(output),
            "--binners",
            "concoct",
            "--qc-tools",
            "eukcc",
            "--qc-decision",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
            expected_coverage_columns=("sample1.bam",),
            concoct_success=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        summary = json.loads((output / "0_log/summary.json").read_text(encoding="utf-8"))
        self.assertEqual(summary["completed_binners"], ["concoct"])

    def test_multi_column_matrix_reaches_all_four_binners(self) -> None:
        matrix = self.root / "mapped-samples.coverage.tsv"
        matrix.write_text(
            "contig_id\tdonor_A\tdonor_B\tdonor_C\n"
            "c1\t3\t30\t300\n"
            "c2\t2\t20\t200\n"
            "c3\t1\t10\t100\n",
            encoding="utf-8",
        )
        output = self.root / "multi-column-all-binners"
        completed = self.run_cli(
            "--sample",
            "sample1",
            "--contigs",
            str(self.fasta),
            "--coverage",
            str(matrix),
            "--output",
            str(output),
            "--binners",
            "remag,concoct,metabat2,vamb",
            "--qc-tools",
            "eukcc",
            "--qc-decision",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
            "--min-retain-length",
            "1",
            "--min-selection-length",
            "1",
            "--min-input-length",
            "1",
            "--threads",
            "1",
            expected_coverage_columns=("donor_A", "donor_B", "donor_C"),
            concoct_success=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        summary = json.loads((output / "0_log/summary.json").read_text(encoding="utf-8"))
        self.assertEqual(
            summary["completed_binners"],
            ["remag", "concoct", "metabat2", "vamb"],
        )
        normalized = (
            output / "0_log/work" / "input" / "depth.jgi.tsv"
        ).read_text(encoding="utf-8").splitlines()
        self.assertEqual(
            normalized[0],
            "contigName\tcontigLen\ttotalAvgDepth\tdonor_A\tdonor_A-var\t"
            "donor_B\tdonor_B-var\tdonor_C\tdonor_C-var",
        )
        self.assertEqual(normalized[1], "c1\t100\t333\t3\tNA\t30\tNA\t300\tNA")
        preflight = output / "0_log/work/results" / "standalone" / "sample1" / "preflight"
        self.assertEqual(
            (preflight / "union.concoct.coverage.tsv")
            .read_text(encoding="utf-8")
            .splitlines()[0],
            "contig_id\tdonor_A\tdonor_B\tdonor_C",
        )
        self.assertEqual(
            (output / "0_log/work/inputs" / "r000001" / "abundance.tsv")
            .read_text(encoding="utf-8")
            .splitlines()[0],
            "contigname\tdonor_A\tdonor_B\tdonor_C",
        )




    def test_zero_qc_candidates_is_a_successful_empty_result(self) -> None:
        output = self.root / "zero-candidates"
        arguments = self.base_run(output)
        qc_index = arguments.index("--min-selection-length") + 1
        arguments[qc_index] = "1000000"
        completed = self.run_cli(
            *arguments,
            "--binners",
            "metabat2",
            "--qc-tools",
            "all",
            "--qc-decision",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
            "--busco-db",
            str(self.busco_db),
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        summary = json.loads((output / "0_log/summary.json").read_text(encoding="utf-8"))
        self.assertEqual(summary["qc_candidates"], 0)
        self.assertEqual(summary["representatives"], 0)
        self.assertEqual(list((output / "3_selection/bins").glob("*.fa")), [])
        self.assertTrue((output / "3_selection/bins" / "manifest.tsv").is_file())

    def test_empty_contig_input_is_a_successful_tool_free_result(self) -> None:
        empty_fasta = self.root / "empty.fa"
        empty_fasta.write_bytes(b"")
        unrelated_depth = self.root / "empty.depth.tsv"
        unrelated_depth.write_text(
            "contigName\tcontigLen\ttotalAvgDepth\n"
            "unrelated\t1000\t1\n",
            encoding="utf-8",
        )
        output = self.root / "empty-input"
        completed = self.run_cli(
            "--sample",
            "empty1",
            "--contigs",
            str(empty_fasta),
            "--coverage",
            str(unrelated_depth),
            "--output",
            str(output),
            "--binners",
            "all",
            "--qc-tools",
            "all",
            "--qc-decision",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
            "--busco-db",
            str(self.busco_db),
            "--random-seed",
            "6311",
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        summary = json.loads(
            (output / "0_log/summary.json").read_text(encoding="utf-8")
        )
        run_metadata = json.loads(
            (output / "0_log" / "run.json").read_text(encoding="utf-8")
        )
        self.assertTrue(summary["empty_input"])
        self.assertEqual(summary["completed_binners"], [])
        self.assertNotIn("taxonomy", summary)
        self.assertEqual(summary["discovered_bins"], 0)
        self.assertEqual(summary["representatives"], 0)
        self.assertEqual(run_metadata["contigs"], 0)
        self.assertFalse((output / "0_log/logs").exists())
        branch_rows = (
            output
            / "0_log/work/results"
            / "standalone"
            / "empty1"
            / "branch_status.tsv"
        ).read_text(encoding="utf-8").splitlines()
        self.assertEqual(len(branch_rows), 5)
        self.assertTrue(
            all(
                "not_run_empty_input" in row
                for row in branch_rows[1:]
            )
        )
        self.assertTrue((output / "0_log/selection" / "representatives.tsv").is_file())
        self.assertTrue((output / "3_selection/bins" / "manifest.tsv").is_file())

    def test_contigs_below_minimum_are_a_tool_free_empty_result(self) -> None:
        output = self.root / "below-minimum-input"
        arguments = self.base_run(output)
        minimum_index = arguments.index("--min-input-length") + 1
        arguments[minimum_index] = "101"
        completed = self.run_cli(
            *arguments,
            "--binners",
            "all",
            "--qc-tools",
            "all",
            "--qc-decision",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
            "--busco-db",
            str(self.busco_db),
            "--random-seed",
            "6311",
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        summary = json.loads(
            (output / "0_log/summary.json").read_text(encoding="utf-8")
        )
        run_metadata = json.loads(
            (output / "0_log" / "run.json").read_text(encoding="utf-8")
        )
        self.assertTrue(summary["empty_input"])
        self.assertEqual(
            summary["empty_input_reason"],
            "no_input_contigs_at_minimum_bp",
        )
        self.assertEqual(summary["completed_binners"], [])
        self.assertEqual(summary["qc_tools"], ["eukcc", "busco"])
        self.assertEqual(summary["random_seed"], 6311)
        self.assertEqual(summary["discovered_bins"], 0)
        self.assertEqual(summary["representatives"], 0)
        self.assertEqual(run_metadata["contigs"], 3)
        self.assertEqual(run_metadata["eligible_contigs"], 0)
        self.assertEqual(run_metadata["eligible_bp"], 0)
        self.assertFalse((output / "0_log/logs").exists())
        branch_rows = (
            output
            / "0_log/work/results"
            / "standalone"
            / "sample1"
            / "branch_status.tsv"
        ).read_text(encoding="utf-8").splitlines()
        self.assertEqual(len(branch_rows), 5)
        self.assertTrue(
            all(
                "not_run_empty_input" in row
                for row in branch_rows[1:]
            )
        )
        self.assertTrue((output / "0_log/.complete.json").is_file())
        self.assertTrue((output / "3_selection/bins" / "manifest.tsv").is_file())

    def test_direct_cli_help_and_removed_options(self) -> None:
        help_result = self.run_cli("--help")
        self.assertEqual(help_result.returncode, 0, help_result.stderr)
        self.assertIn("--min-input-length", help_result.stdout)
        self.assertIn("--max-selection-contamination", help_result.stdout)
        self.assertNotIn("--kunpeng-db", help_result.stdout)
        self.assertIn("--random-seed INT", help_result.stdout)
        self.assertIn("usage: eukbinpicker", help_result.stdout)
        self.assertIn("default: 6311", help_result.stdout)
        for removed in (
            "--input-min-contig-exclusive",
            "--high-completeness",
            "--medium-completeness",
            "--seqkit",
            "--vamb",
            "--remag-env",
            "--kraken-db",
            "--kraken-memory-mapping",
        ):
            self.assertNotIn(removed, help_result.stdout.split())

        old_run = self.run_cli("run")
        self.assertNotEqual(old_run.returncode, 0)
        removed_taxonomy = self.run_cli("taxonomy")
        self.assertNotEqual(removed_taxonomy.returncode, 0)

    def test_random_seed_is_unified_and_recorded(self) -> None:
        output = self.root / "random-seed"
        seed = 7311
        completed = self.run_cli(
            *self.base_run(output),
            "--binners",
            "remag,metabat2,vamb",
            "--qc-tools",
            "eukcc",
            "--qc-decision",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
            "--random-seed",
            str(seed),
            expected_seed=seed,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)

        run_metadata = json.loads(
            (output / "0_log" / "run.json").read_text(encoding="utf-8")
        )
        summary = json.loads(
            (output / "0_log/summary.json").read_text(encoding="utf-8")
        )
        self.assertEqual(run_metadata["random_seed"], seed)
        self.assertEqual(summary["random_seed"], seed)
        backend_environment = (
            output / "0_log/config" / "backend.env"
        ).read_text(encoding="utf-8")
        self.assertIn(f"RANDOM_SEED={seed}\n", backend_environment)

        branch_root = (
            output
            / "0_log/work/results"
            / "standalone"
            / "sample1"
            / "results"
            / "union"
        )
        for binner in ("remag", "metabat2", "vamb"):
            fingerprint = json.loads(
                (branch_root / binner / "fingerprint.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(
                fingerprint["fields"]["random_seed"],
                str(seed),
            )

        conflicting_reuse = self.run_cli(
            *self.base_run(output),
            "--binners",
            "remag,metabat2,vamb",
            "--qc-tools",
            "eukcc",
            "--qc-decision",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
            "--random-seed",
            "8311",
        )
        self.assertEqual(conflicting_reuse.returncode, 2)
        self.assertIn(
            "completed output random seed differs",
            conflicting_reuse.stderr,
        )

        invalid = self.run_cli(
            "--output", str(self.root / "invalid-seed"),
            "--sample",
            "sample1",
            "--contigs",
            str(self.fasta),
            "--coverage",
            str(self.depth),
            "--binners",
            "metabat2",
            "--qc-tools",
            "eukcc",
            "--qc-decision",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
            "--random-seed",
            "0",
        )
        self.assertEqual(invalid.returncode, 2)
        self.assertIn("--random-seed must be positive", invalid.stderr)

    def test_quality_thresholds_are_inclusive(self) -> None:
        accepted_output = self.root / "quality-boundary-accepted"
        accepted = self.run_cli(
            *self.base_run(accepted_output),
            "--binners",
            "metabat2",
            "--qc-tools",
            "busco",
            "--qc-decision",
            "busco",
            "--busco-db",
            str(self.busco_db),
            "--min-selection-completeness",
            "98",
            "--max-selection-contamination",
            "1",
        )
        self.assertEqual(accepted.returncode, 0, accepted.stderr)
        accepted_summary = json.loads(
            (accepted_output / "0_log/summary.json").read_text(encoding="utf-8")
        )
        self.assertEqual(accepted_summary["representatives"], 1)

        rejected_output = self.root / "quality-boundary-rejected"
        rejected = self.run_cli(
            *self.base_run(rejected_output),
            "--binners",
            "metabat2",
            "--qc-tools",
            "busco",
            "--qc-decision",
            "busco",
            "--busco-db",
            str(self.busco_db),
            "--min-selection-completeness",
            "98.1",
            "--max-selection-contamination",
            "1",
        )
        self.assertEqual(rejected.returncode, 0, rejected.stderr)
        rejected_summary = json.loads(
            (rejected_output / "0_log/summary.json").read_text(encoding="utf-8")
        )
        self.assertEqual(rejected_summary["representatives"], 0)

        eukcc_output = self.root / "eukcc-contamination-boundary"
        eukcc = self.run_cli(
            *self.base_run(eukcc_output),
            "--binners",
            "metabat2",
            "--qc-tools",
            "eukcc",
            "--qc-decision",
            "eukcc",
            "--eukcc-db",
            str(self.eukcc_db),
            "--min-selection-completeness",
            "91",
            "--max-selection-contamination",
            "1",
        )
        self.assertEqual(eukcc.returncode, 0, eukcc.stderr)
        with (eukcc_output / "0_log/selection" / "representatives.tsv").open(
            encoding="utf-8",
            newline="",
        ) as handle:
            representative = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(
            representative["contamination_source"],
            "eukcc_contamination",
        )
        self.assertEqual(representative["contamination"], "1")


if __name__ == "__main__":
    unittest.main()
