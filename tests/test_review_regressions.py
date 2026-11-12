from __future__ import annotations
import csv
import fcntl
import json
import shutil
import unittest
import os
import subprocess
from unittest.mock import patch
from pathlib import Path
from test_multi_coverage import MODULE as cli
import test_mock_end_to_end as mock
from selection import fasta_contigs, Candidate, shared_overlap
from qc_parsers import parse_summary
from depth_schema import identify_depth_schema


class ReviewRegressionTests(unittest.TestCase):
    setUp = mock.MockEndToEndTests.setUp
    tearDown = mock.MockEndToEndTests.tearDown
    run_cli = mock.MockEndToEndTests.run_cli
    base_run = mock.MockEndToEndTests.base_run

    def arguments(self, output):
        return [*self.base_run(output), "--binners", "metabat2",
                "--qc-tools", "eukcc", "--qc-decision", "eukcc",
                "--eukcc-db", str(self.eukcc_db)]

    def test_remag_gpu_probes_the_backend_environment(self):
        info = {"device": "cuda", "torch_version": "test", "cuda_version": "12.4", "gpu_name": "mock GPU"}
        for environment, selector in (("remag-gpu", "-n"), ("/env/remag-gpu", "-p")):
            with patch.dict(os.environ, {"EUKBINPICKER_REMAG_ENV": environment,
                                        "EUKBINPICKER_CONDA_EXE": "/mock/conda"}, clear=True):
                with patch.object(cli.subprocess, "run", return_value=subprocess.CompletedProcess(
                        [], 0, json.dumps(info), "")) as run:
                    self.assertEqual(cli.require_remag_cuda(), info)
                    self.assertEqual(run.call_args.args[0][:6],
                                     ["/mock/conda", "run", "--no-capture-output", selector, environment, "python"])
        with patch.dict(os.environ, {}, clear=True), patch.object(cli.shutil, "which", return_value=str(self.tools / "remag")):
            with patch.object(cli.subprocess, "run", return_value=subprocess.CompletedProcess(
                    [], 0, json.dumps(info), "")) as run:
                cli.require_remag_cuda()
                self.assertEqual(run.call_args.args[0][:2], [cli.sys.executable, "-c"])

    def test_remag_gpu_rejects_cpu_and_probe_errors(self):
        for result in (subprocess.CompletedProcess([], 2, '{"device":"cpu"}', ""),
                       subprocess.CompletedProcess([], 0, '{"device":"cpu"}', ""),
                       subprocess.CompletedProcess([], 1, "", "No module named torch")):
            with patch.object(cli.subprocess, "run", return_value=result), patch.object(cli.shutil, "which", return_value=str(self.tools / "remag")):
                with self.assertRaises(cli.UserError):
                    cli.require_remag_cuda()
        with patch.object(cli.subprocess, "run", side_effect=subprocess.TimeoutExpired("probe", 60)), patch.object(cli.shutil, "which", return_value=str(self.tools / "remag")):
            with self.assertRaises(cli.UserError):
                cli.require_remag_cuda()


    def test_gpu_requires_a_supported_binner(self):
        args = cli.build_parser().parse_args([
            *self.arguments(self.root / "gpu-wrong-binner"), "--gpu", "T"])
        with self.assertRaisesRegex(cli.UserError, "requires remag or vamb"):
            cli.validate_args(args)


    def test_removed_taxonomy_options_are_rejected(self):
        for extra in (["--binners", "taxvamb"], ["--taxonomy", "unused.tsv"],
                      ["--kunpeng-db", "unused"]):
            result = self.run_cli(*self.arguments(self.root / "removed"), *extra)
            self.assertEqual(result.returncode, 2, result.stderr)
        result = self.run_cli("--help")
        self.assertNotIn("--taxonomy", result.stdout)
        self.assertNotIn("--kunpeng-db", result.stdout)


    def test_distinct_contig_ids_preserve_identical_sequences(self):
        left = self.root / "left.fa"
        right = self.root / "right.fa"
        left.write_text(">a\nACGT\n>b\nACGT\n")
        right.write_text(">c\nACGT\n")
        lc, rc = fasta_contigs(left), fasta_contigs(right)
        self.assertEqual(len(lc), 2)
        a = Candidate({}, 90, 0, 2, 4, 90, lc, 8)
        b = Candidate({}, 80, 0, 1, 4, 80, rc, 4)
        self.assertEqual(shared_overlap(a, b), (0, 0.0))
        right.write_text(">b\nACGT\n")
        b.contig_bp = fasta_contigs(right)
        self.assertEqual(shared_overlap(a, b), (4, 1.0))

    def test_invalid_quality_cannot_pass(self):
        for completeness, contamination in (("101", "0"), ("90", "-1"), ("nan", "0"), ("90", "inf")):
            passed, score = cli.qc_metrics(completeness, contamination,
                min_completeness=50, max_contamination=10, contamination_penalty=5)
            self.assertEqual((passed, score), ("no", ""))

    def test_jgi_alias_is_applied(self):
        schema = identify_depth_schema(
            ["contigName", "contigLen", "totalAvgDepth", "sample.bam", "sample.bam-var"],
            "assembly", single_coverage_name="mapping_A")
        self.assertEqual(schema.coverage_names, ("mapping_A",))
        self.assertEqual(schema.source_coverage_names, ("sample.bam",))

    def test_invalid_busco_summary_is_rejected(self):
        summary = self.root / "summary.txt"
        for text in ("C:101%[S:100%,D:1%],F:0%,M:0%,n:100",
                     "C:90%[S:80%,D:1%],F:5%,M:5%,n:100",
                     "C:90%[S:89%,D:1%],F:0%,M:0%,n:100"):
            summary.write_text(text)
            self.assertIsNone(parse_summary(summary))
