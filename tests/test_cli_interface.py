"""Public input routing and inclusive selection defaults, using synthetic data."""
import gzip
import io
import json
import unittest
from unittest.mock import patch

from test_multi_coverage import MODULE as cli
import test_mock_end_to_end as mock


class CliInterfaceTests(unittest.TestCase):
    setUp = mock.MockEndToEndTests.setUp
    tearDown = mock.MockEndToEndTests.tearDown
    run_cli = mock.MockEndToEndTests.run_cli
    base_run = mock.MockEndToEndTests.base_run

    def test_only_long_help_is_supported_and_subcommands_are_rejected(self):
        help_result = self.run_cli("--help")
        self.assertEqual(help_result.returncode, 0, help_result.stderr)
        self.assertNotIn("[-h]", help_result.stdout)
        self.assertNotIn("doctor", help_result.stdout)
        self.assertNotIn("validate", help_result.stdout)
        for removed in ("-h", "doctor", "validate"):
            with self.subTest(removed=removed):
                output = self.root / "rejected-command"
                result = self.run_cli(*self.base_run(output), removed)
                self.assertEqual(result.returncode, 2)
                self.assertIn(f"unrecognized arguments: {removed}", result.stderr)
                self.assertFalse(output.exists())

    def test_normal_run_reports_invalid_inputs_and_databases_before_binning(self):
        for case in ("coverage", "database"):
            with self.subTest(case=case):
                output = self.root / f"bad-{case}"
                if case == "coverage":
                    self.depth.write_text("contig_id\tdepth\nc1\t-1\nc2\t2\nc3\t1\n")
                    database = self.eukcc_db
                    message = "finite and non-negative"
                else:
                    self.depth.write_text("contig_id\tdepth\nc1\t1\nc2\t2\nc3\t1\n")
                    database = self.root / "missing-db"
                    message = "EukCC needs a valid --eukcc-db"
                result = self.run_cli(*self.base_run(output), "--binners", "metabat2",
                                      "--eukcc-db", str(database))
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertIn(message, result.stderr)
                self.assertNotIn("Traceback", result.stderr)
                self.assertFalse((output / "0_log/logs/binning.classic.log").exists())

    def mock_cuda(self):
        modules = self.root / "mock_python_modules"
        (modules / "torch.py").write_text('''import os
from types import SimpleNamespace
__version__ = "mock"
version = SimpleNamespace(cuda="mock")
def manual_seed(seed): pass
cuda = SimpleNamespace(
    is_available=lambda: os.environ.get("CUDA_VISIBLE_DEVICES") != "",
    get_device_name=lambda index: "mock GPU",
    manual_seed_all=lambda seed: None)
backends = SimpleNamespace(mps=SimpleNamespace(is_available=lambda: True))
''')
        package = modules / "remag"
        package.mkdir()
        (package / "__init__.py").write_text("")
        (package / "utils.py").write_text('''import torch
from types import SimpleNamespace
def get_torch_device():
    return SimpleNamespace(type="cuda" if torch.cuda.is_available() else "mps" if torch.backends.mps.is_available() else "cpu")
''')
        self.environment["CUDA_VISIBLE_DEVICES"] = "3"
        self.environment["DEVICE_TRACE"] = str(self.root / "devices.jsonl")
        for name in ("remag", "vamb"):
            tool = self.tools / name
            script = tool.read_text()
            script = script.replace('args = sys.argv[1:]', '''args = sys.argv[1:]
if "--version" not in args:
    from remag.utils import get_torch_device
    device = get_torch_device().type if name == "remag" else ("cuda" if "--cuda" in args else "cpu")
    with open(os.environ["DEVICE_TRACE"], "a") as handle:
        handle.write(json.dumps(dict(tool=name, device=device, visible=os.environ.get("CUDA_VISIBLE_DEVICES"))) + "\\n")
''')
            tool.write_text(script)

    def test_gpu_switch_controls_both_tools_and_preserves_allocated_device(self):
        self.mock_cuda()
        for gpu, device, visible in (("F", "cpu", ""), ("T", "cuda", "3")):
            with self.subTest(gpu=gpu):
                output = self.root / f"gpu-{gpu}"
                result = self.run_cli(*self.base_run(output), "--binners", "remag,vamb",
                    "--eukcc-db", str(self.eukcc_db), "--gpu", gpu)
                self.assertEqual(result.returncode, 0, result.stderr)
                traces = [json.loads(line) for line in (self.root / "devices.jsonl").read_text().splitlines()][-2:]
                self.assertEqual({row["tool"] for row in traces}, {"remag", "vamb"})
                self.assertTrue(all(row["device"] == device and row["visible"] == visible for row in traces), traces)
                config = json.loads((output / "0_log/run.json").read_text())
                self.assertEqual(config["inputs"]["gpu"], gpu)
                self.assertEqual(set(config["inputs"]["cuda_devices"]), {"remag", "vamb"} if gpu == "T" else set())
        changed = self.run_cli(*self.base_run(self.root / "gpu-F"), "--binners", "remag,vamb",
            "--eukcc-db", str(self.eukcc_db), "--gpu", "T")
        self.assertNotEqual(changed.returncode, 0)
        self.assertIn("parameters or software differ", changed.stderr)

    def test_vamb_only_gpu_requires_available_cuda_before_output(self):
        self.mock_cuda()
        self.environment["CUDA_VISIBLE_DEVICES"] = ""
        output = self.root / "no-gpu"
        result = self.run_cli(*self.base_run(output), "--binners", "vamb",
            "--eukcc-db", str(self.eukcc_db), "--gpu", "T")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("VAMB", result.stderr)
        self.assertFalse(output.exists())

    def test_defaults_keep_low_scoring_but_valid_qc(self):
        # Below the previous 50% completeness and above 10% contamination.
        tool = self.tools / "eukcc"
        tool.write_text(tool.read_text().replace('"91", "1"', '"0", "100"'))
        output = self.root / "default-selection"
        args = self.base_run(output)
        i = args.index("--min-retain-length")
        del args[i:i + 2]
        result = self.run_cli(*args, "--binners", "metabat2",
                              "--eukcc-db", str(self.eukcc_db))
        self.assertEqual(result.returncode, 0, result.stderr)
        record = json.loads((output / "0_log/run.json").read_text())
        summary = record["summary"]
        self.assertEqual(record["inputs"]["qc_tools"], ["eukcc"])
        self.assertEqual(summary["representatives"], 1)
        self.assertFalse((output / "2_qc/busco").exists())

    def test_multiple_files_and_custom_columns_keep_sample_pairing(self):
        self.depth.write_text("id\tmean\tvar\nc1\t3\t7\nc2\t2\t5\nc3\t1\t0\n")
        second = self.root / "second.tsv.gz"
        with gzip.open(second, "wt") as handle:
            handle.write("id\tmean\tvar\nc1\t4\t8\nc2\t0\t0\nc3\t2\t1\n")
        output = self.root / "custom-columns"
        result = self.run_cli("--output", str(output), "--sample", "s", "--contigs", str(self.fasta),
            "--coverage", str(self.depth), str(second),
            "--cov-contig-column", "id", "--cov-mean-column", "mean",
            "--cov-variance-column", "var", "--eukcc-db", str(self.eukcc_db))
        self.assertEqual(result.returncode, 0, result.stderr)
        plan = json.loads((output / "0_log/run.json").read_text())
        self.assertEqual(plan["inputs"]["coverage_input"], "table")
        self.assertEqual(plan["inputs"]["coverage_columns"], ["depth", "second"])
        self.assertEqual(plan["inputs"]["variance_available"], [True, True])

    def test_bam_signature_without_extension_routes_to_bam_validation(self):
        bam = self.root / "alignments.data"
        with gzip.open(bam, "wb") as handle:
            handle.write(b"BAM\x01synthetic")
        args = cli.build_parser().parse_args(["--sample", "s", "--contigs", str(self.fasta),
            "--coverage", str(bam), "--output", str(self.root / "out")])
        cli.validate_args(args)
        self.assertEqual(args.bam, [bam.resolve()])
        self.assertIsNone(args.depth)

    def test_old_options_are_rejected_and_new_limits_are_checked(self):
        base = self.base_run(self.root / "out")
        for option in ("--study", "--depth", "--bam", "--decision-qc", "--min-contig-bp",
                       "--retain-min-bin-bp", "--qc-min-bin-bp", "--min-completeness",
                       "--max-contamination", "--overlap-threshold", "--coverage-columns", "--remag-gpu"):
            with self.subTest(option=option), patch("sys.stderr", io.StringIO()):
                with self.assertRaises(SystemExit):
                    cli.build_parser().parse_args([*base, option, "1"])
        for option, value in (("--min-retain-length", "-1"),
                              ("--min-selection-length", "-1"),
                              ("--min-selection-completeness", "nan"),
                              ("--max-selection-contamination", "101"),
                              ("--min-overlap", "0")):
            with self.subTest(option=option):
                args = cli.build_parser().parse_args([*base, option, value])
                with self.assertRaises(ValueError):
                    cli.validate_args(args)
