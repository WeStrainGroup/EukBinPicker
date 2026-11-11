"""Public output layout and full QC provenance, exercised through the CLI."""
import csv
from pathlib import Path
import unittest

import test_mock_end_to_end as mock_workflow


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class OutputLayoutTests(unittest.TestCase):
    def setUp(self):
        self.case = mock_workflow.MockEndToEndTests()
        self.case.setUp()
        self.addCleanup(self.case.tearDown)

    def test_public_tables_include_both_qc_tools_and_final_paths(self):
        case = self.case
        output = case.root / "layout"
        args = [*case.base_run(output), "--binners", "remag,metabat2",
                "--qc-tools", "all", "--eukcc-db", str(case.eukcc_db),
                "--busco-db", str(case.busco_db)]
        result = case.run_cli(*args)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual({p.name for p in output.iterdir() if not p.name.startswith(".")},
                         {"0_log", "1_raw", "2_qc", "3_selection"})
        qc = rows(output / "2_qc/summary.tsv")
        self.assertEqual({r["binner"] for r in qc}, {"remag", "metabat2"})
        for row in qc:
            self.assertEqual(row["eukcc_status"], "ok")
            self.assertEqual(row["busco_status"], "ok")
            self.assertTrue(Path(row["bin_path"]).is_file())
            self.assertIn("/1_raw/", row["bin_path"])
            self.assertGreater(int(row["longest_contig"]), 0)
        reps = rows(output / "3_selection/representatives.tsv")
        self.assertTrue(reps)
        represented = []
        for row in reps:
            self.assertTrue(Path(row["bin_path"]).is_file())
            self.assertIn("/3_selection/bins/", row["bin_path"])
            self.assertTrue(row["selection_reason"])
            self.assertIn("eukcc_completeness", row)
            self.assertIn("busco_complete", row)
            represented.extend(row["represented_bins"].split(";"))
        self.assertCountEqual(represented, [r["bin_id"] for r in qc])
        resumed = case.run_cli(*args)
        self.assertEqual(resumed.returncode, 0, resumed.stderr)
        raw = next((output / "1_raw").glob("*/*"))
        raw.write_text(raw.read_text() + "\n")
        changed = case.run_cli(*args)
        self.assertNotEqual(changed.returncode, 0)
        self.assertIn("retained bins or public summaries changed", changed.stderr)

    def test_bins_below_qc_length_are_retained(self):
        case = self.case
        output = case.root / "retained"
        args = case.base_run(output)
        args[args.index("--min-selection-length") + 1] = "100000"
        result = case.run_cli(*args, "--binners", "metabat2", "--eukcc-db", str(case.eukcc_db))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(list((output / "1_raw/metabat2").glob("*")))
        self.assertEqual(rows(output / "2_qc/summary.tsv"), [])
        self.assertEqual(rows(output / "3_selection/representatives.tsv"), [])
