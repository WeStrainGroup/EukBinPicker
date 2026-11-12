"""Generated and imported candidates share one QC and selection pass."""
import csv
import json
import unittest
import test_import_bins as imports
import test_mock_end_to_end as mocks


class CombinedInputsTests(unittest.TestCase):
    setUp = mocks.MockEndToEndTests.setUp
    tearDown = mocks.MockEndToEndTests.tearDown
    run_cli = mocks.MockEndToEndTests.run_cli
    bins = imports.ImportBinsTests.bins

    def arguments(self, *sources, binners="metabat2", minimum="1"):
        return ["--sample", "sample1", "--contigs", str(self.fasta),
                "--coverage", str(self.depth), "--bins", *map(str, sources),
                "--binners", binners, "--min-input-length", minimum,
                "--min-selection-length", "1", "--eukcc-db", str(self.eukcc_db),
                "--output", str(self.root / "out")]

    def summary(self):
        return json.loads((self.root / "out/0_log/summary.json").read_text())

    def test_joint_qc_selection_and_safe_reuse(self):
        source = self.bins("previous")
        args = self.arguments(source)
        result = self.run_cli(*args, expected_seed=6311)
        self.assertEqual(result.returncode, 0, result.stderr)
        summary = self.summary()
        self.assertEqual(summary["input_mode"], "combined")
        self.assertEqual(summary["native_binner_status"], {"metabat2": "completed"})
        self.assertEqual(summary["completed_native_binners"], ["metabat2"])
        self.assertEqual(summary["imported_sources"], ["previous"])
        self.assertEqual(summary["qc_candidates"], 2)
        self.assertEqual(summary["representatives"], 1)
        self.assertEqual(summary["alternatives"], 1)
        out = self.root / "out"
        with (out / "3_selection/representatives.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(rows[0]["binner"], "metabat2")
        self.assertEqual((out / "1_raw/previous/one.fa").read_bytes(), (source / "one.fa").read_bytes())
        self.assertTrue(list((out / "1_raw/metabat2").glob("*.fa")))
        again = self.run_cli(*args)
        self.assertEqual(again.returncode, 0, again.stderr)
        self.assertIn("already complete", again.stderr)
        (source / "one.fa").write_text(">c1\nACGT\n")
        changed = self.run_cli(*args)
        self.assertNotEqual(changed.returncode, 0)
        self.assertIn("differ", changed.stderr)

    def test_colliding_source_labels_and_repeated_arguments(self):
        a, b = self.bins("a/metabat2"), self.bins("b/metabat2")
        result = self.run_cli(*self.arguments(a), "--bins", str(b))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.summary()["completed_binners"], ["metabat2", "metabat2_2", "metabat2_3"])
        self.assertEqual(self.summary()["qc_candidates"], 3)
        for label, source in [("metabat2_2", a), ("metabat2_3", b)]:
            self.assertEqual((self.root / "out/1_raw" / label / "one.fa").read_bytes(), (source / "one.fa").read_bytes())

    def test_zero_native_eligible_contigs_still_assesses_whole_imports(self):
        source = self.bins("previous", ">short_not_in_binner_input\nACGT\n")
        result = self.run_cli(*self.arguments(source, minimum="1500"))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.summary()["completed_native_binners"], [])
        self.assertEqual(self.summary()["native_binner_status"], {"metabat2": "skipped_empty_input"})
        self.assertEqual(self.summary()["qc_candidates"], 1)
        self.assertEqual(self.summary()["representatives"], 1)
        self.assertEqual((self.root / "out/1_raw/previous/one.fa").read_bytes(), (source / "one.fa").read_bytes())

    def test_native_failure_does_not_discard_imports(self):
        source = self.bins("previous")
        result = self.run_cli(*self.arguments(source, binners="concoct"))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("all requested binners failed", result.stderr)
        self.assertEqual(self.summary()["native_binner_status"], {"concoct": "failed"})
        self.assertEqual(self.summary()["completed_native_binners"], [])
        self.assertEqual(self.summary()["completed_binners"], ["previous"])
        self.assertEqual(self.summary()["representatives"], 1)

    def test_busco_and_eukcc_assess_both_sources(self):
        source = self.bins("previous")
        result = self.run_cli(*self.arguments(source), "--qc-tools", "all",
                              "--qc-decision", "busco", "--busco-db", str(self.busco_db))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.summary()["qc_candidates"], 2)
        self.assertEqual(self.summary()["decision_qc"], "busco")
        self.assertEqual(len(list((self.root / "out/2_qc/busco").glob("*/busco_cwd_marker.log"))), 2)
        self.assertEqual(self.summary()["representatives"], 1)

    def test_invalid_import_rejected_before_native_run(self):
        source = self.bins("bad", ">c1\nACGT\n>c1\nACGT\n")
        result = self.run_cli(*self.arguments(source))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("duplicate", result.stderr)
        self.assertFalse((self.root / "out/0_log/run.json").exists())

    def test_at_least_one_input_required(self):
        result = self.run_cli("--sample", "sample1", "--output", str(self.root / "out"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("provide --coverage, --bins, or both", result.stderr)
        self.assertNotIn("Traceback", result.stderr)
