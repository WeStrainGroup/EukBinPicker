"""Exercise native argument delivery, conflict errors and output reuse."""
import json
import shlex
import unittest
from test_multi_coverage import MODULE as cli

import test_mock_end_to_end as mock


class BinnerParameterTests(unittest.TestCase):
    setUp = mock.MockEndToEndTests.setUp
    tearDown = mock.MockEndToEndTests.tearDown
    base_run = mock.MockEndToEndTests.base_run
    run_cli = mock.MockEndToEndTests.run_cli

    def test_parameters_reach_all_binners_as_literal_arguments(self):
        trace = self.root / "trace"
        trace.mkdir()
        self.environment["MOCK_NATIVE_TRACE"] = str(trace)
        for name in ("remag", "concoct", "metabat2", "vamb"):
            tool = self.tools / name
            source = tool.read_text().replace('args = sys.argv[1:]', '''args = sys.argv[1:]
if "--version" not in args and "--help" not in args:
    (Path(os.environ["MOCK_NATIVE_TRACE"]) / (name + ".json")).write_text(json.dumps(args))
''')
            if name == "concoct":
                source = source.replace('if any(flag in args for flag in ("--no_cov_normalization", "--no_total_coverage", "--clusters")):', 'if False:')
            tool.write_text(source)
        marker = self.root / "must-not-exist"
        literal = f"sample boundary; $(touch {marker}) `touch {marker}`"
        parameters = {
            "remag": ["--epochs", "200", "--batch-size", "1024"],
            "concoct": ["--clusters", "50", "--no_cov_normalization"],
            "metabat2": ["--maxP", "90", "--minS", "70"],
            "vamb": ["-e", "200", "-n", "128", "128", "-o", literal],
        }
        output = self.root / "all-binners"
        arguments = [*self.base_run(output), "--eukcc-db", str(self.eukcc_db)]
        for name, tokens in parameters.items():
            arguments += [f"--{name}-parameters", shlex.join(tokens)]
        result = self.run_cli(*arguments, concoct_success=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        metadata = json.loads((output / "0_log/run.json").read_text())
        self.assertEqual(metadata["inputs"]["binner_parameters"], parameters)
        self.assertEqual(set(metadata['summary']['native_binner_status'].values()), {'completed'})
        for name, tokens in parameters.items():
            received = json.loads((trace / f"{name}.json").read_text())
            self.assertEqual(received[-len(tokens):], tokens)
            commands = [c['argv'] for c in metadata['commands'] if c['stage'] == name]
            self.assertTrue(any(command[-len(tokens):] == tokens for command in commands))
        self.assertFalse(marker.exists())
        reused = self.run_cli(*arguments, concoct_success=True)
        self.assertEqual(reused.returncode, 0, reused.stderr)
        self.assertIn("already complete", reused.stderr)
        changed = self.run_cli(*arguments, "--remag-parameters=--epochs 201", concoct_success=True)
        self.assertNotEqual(changed.returncode, 0)
        self.assertIn("parameters or software differ", changed.stderr)

    def test_quoted_single_flags_do_not_require_equals(self):
        base = ["--sample", "x", "--contigs", "x.fa", "--coverage", "x.tsv", "--output", "out"]
        for name, value in (("concoct", "--no_cov_normalization"), ("vamb", "-o")):
            parsed = cli.parse_arguments([*base, f"--{name}-parameters", value])
            self.assertEqual(getattr(parsed, f"{name}_parameters"), value)

    def test_missing_native_value_does_not_consume_next_public_option(self):
        result = self.run_cli(*self.base_run(self.root / "missing"),
            "--concoct-parameters", "--threads", "8")
        self.assertEqual(result.returncode, 2)
        self.assertIn("expected one argument", result.stderr)

    def test_managed_options_and_invalid_parameter_strings_fail_before_output(self):
        cases = (
            ("remag", "--threads=8"), ("remag", "-o/tmp/elsewhere"),
            ("remag", "--coverage other.tsv"), ("remag", "--filter-only"),
            ("concoct", "--length_threshold 500"), ("concoct", "--base=/tmp/elsewhere"),
            ("metabat2", "-t8"), ("metabat2", "--cvExt"), ("metabat2", "-vt8"),
            ("vamb", "--cuda"), ("vamb", "--minfasta 1"), ("vamb", "--out /tmp/elsewhere"),
            ("remag", '--epochs "200'), ("remag", "-- --output elsewhere"),
            ("concoct", "--help"), ("concoct", "-v"),
        )
        for name, text in cases:
            with self.subTest(name=name, text=text):
                output = self.root / "rejected"
                result = self.run_cli(*self.base_run(output), "--eukcc-db", str(self.eukcc_db),
                                      f"--{name}-parameters={text}")
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertIn(f"--{name}-parameters", result.stderr)
                self.assertFalse(output.exists())

    def test_parameters_for_unselected_binner_are_not_silently_ignored(self):
        result = self.run_cli(*self.base_run(self.root / "unselected"), "--binners", "vamb",
            "--eukcc-db", str(self.eukcc_db), "--remag-parameters=--epochs 200")
        self.assertEqual(result.returncode, 2)
        self.assertIn("requires remag in --binners", result.stderr)
