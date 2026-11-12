"""External bins use the same QC and selector, without invoking any binner."""
import csv
import gzip
import json
from pathlib import Path
import unittest

import test_mock_end_to_end as mocks
from test_multi_coverage import MODULE as cli
from import_bins import discover_bin_sets, inspect_bin_sets


class ImportBinsTests(unittest.TestCase):
    setUp = mocks.MockEndToEndTests.setUp
    tearDown = mocks.MockEndToEndTests.tearDown
    run_cli = mocks.MockEndToEndTests.run_cli

    def bins(self, name, sequence=None):
        directory = self.root / name
        directory.mkdir(parents=True)
        (directory / 'one.fa').write_text(sequence or self.fasta.read_text())
        return directory

    def arguments(self, *sources):
        return ['--sample', 'sample1',
                '--bins', *map(str, sources), '--eukcc-db', str(self.eukcc_db),
                '--min-selection-length', '1', '--output', str(self.root / 'out')]

    def test_full_workflow_preserves_bins_and_skips_binners(self):
        a, b = self.bins('custom'), self.bins('metabat2')
        # Deliberately omit binning executables and all coverage. Short contigs
        # must survive the default 1500 bp threshold unchanged in import mode.
        for name in ('remag', 'metabat2', 'concoct', 'vamb'):
            (self.tools / name).unlink()
        original = self.fasta.read_text()
        self.fasta.unlink()
        args = self.arguments(a, b)
        result = self.run_cli(*args)
        self.assertEqual(result.returncode, 0, result.stderr)
        out = self.root / 'out'
        summary = json.loads((out / '0_log/summary.json').read_text())
        self.assertEqual(summary['input_mode'], 'bins')
        self.assertEqual(summary['qc_candidates'], 2)
        self.assertEqual(summary['representatives'], 1)
        self.assertEqual(summary['alternatives'], 1)
        self.assertEqual((out / '1_raw/custom/one.fa').read_bytes(), (a / 'one.fa').read_bytes())
        self.assertFalse((out / '0_log/work/input/depth.jgi.tsv').exists())
        metadata = json.loads((out / '0_log/run.json').read_text())
        self.assertIsNone(metadata['contigs_path'])
        self.assertIsNone(metadata['contigs_sha256'])
        self.assertEqual(metadata['contigs'], 3)
        with (out / '3_selection/representatives.tsv').open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(rows[0]['binner'], 'metabat2')
        again = self.run_cli(*args)
        self.assertEqual(again.returncode, 0, again.stderr)
        self.assertIn('already complete', again.stderr)
        (a / 'one.fa').write_text('>c1\n'+'ACGT'*25+'\n')
        self.assertEqual(len((out / '1_raw/custom/one.fa').read_text()), len(original))
        changed = self.run_cli(*args)
        self.assertNotEqual(changed.returncode, 0)
        self.assertIn('differ', changed.stderr)

    def test_retention_and_qc_thresholds_are_separate(self):
        a = self.bins('custom')
        (a / 'small.fa').write_text('>c1\n'+'ACGT'*25+'\n')
        result = self.run_cli(*self.arguments(a), '--min-retain-length', '200', '--min-selection-length', '301')
        self.assertEqual(result.returncode, 0, result.stderr)
        out = self.root / 'out'
        self.assertTrue((out / '1_raw/custom/one.fa').exists())
        self.assertFalse((out / '1_raw/custom/small.fa').exists())
        summary = json.loads((out / '0_log/summary.json').read_text())
        self.assertEqual(summary['discovered_bins'], 2)
        self.assertEqual(summary['qc_candidates'], 0)
        self.assertEqual(summary['representatives'], 0)

    def test_coverage_mode_still_requires_contigs(self):
        result = self.run_cli('--sample', 'sample1', '--coverage', str(self.depth),
                              '--output', str(self.root / 'out'))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('--contigs is required with --coverage', result.stderr)
        self.assertNotIn('Traceback', result.stderr)

    def test_rejects_malformed_fasta(self):
        for text, error in [('>c1\nACGT\n>c1\nACGT\n', 'duplicate'), ('', 'no FASTA')]:
            with self.subTest(error=error):
                a = self.root / 'bad'
                a.mkdir(exist_ok=True)
                (a / 'one.fa').write_text(text)
                result = self.run_cli(*self.arguments(a))
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(error, result.stderr)
                self.assertFalse((self.root / 'out/0_log/run.json').exists())

    def test_source_names_gzip_and_filename_collisions(self):
        a, b = self.bins('tool_a/bins'), self.bins('tool_b/bins')
        sources = discover_bin_sets([a, b])
        self.assertEqual(list(sources), ['bins', 'bins_2'])
        (b / 'one.fa').unlink()
        with gzip.open(b / 'one.fna.gz', 'wt') as handle:
            handle.write(self.fasta.read_text())
        self.assertEqual(len(inspect_bin_sets(discover_bin_sets([a, b]))[0]), 2)
        with self.assertRaisesRegex(ValueError, 'duplicate bin directory'):
            discover_bin_sets([a, a])
        (a / 'one.fasta').write_text(self.fasta.read_text())
        with self.assertRaisesRegex(ValueError, 'duplicate identifiers'):
            discover_bin_sets([a])

    def test_bins_only_rejects_binning_options_or_coverage_without_contigs(self):
        a = self.bins('custom')
        for options in [('--contigs', str(self.fasta)), ('--coverage', str(self.depth)), ('--gpu', 'T'),
                        ('--binners', 'vamb'), ('--cov-contig-column', 'id')]:
            with self.subTest(options=options):
                result = self.run_cli(*self.arguments(a), *options)
                self.assertNotEqual(result.returncode, 0)

    def test_gzip_import_busco_and_completed_output_verification(self):
        a = self.bins('custom')
        (a / 'one.fa').unlink()
        with gzip.open(a / 'one.fa.gz', 'wt') as handle:
            handle.write(self.fasta.read_text())
        args = [*self.arguments(a), '--qc-tools', 'busco', '--qc-decision', 'busco',
                '--busco-db', str(self.busco_db)]
        result = self.run_cli(*args)
        self.assertEqual(result.returncode, 0, result.stderr)
        retained = self.root / 'out/1_raw/custom/one.fa'
        self.assertEqual(retained.read_text(), self.fasta.read_text())
        self.assertEqual(len(list((self.root / 'out/2_qc/busco').glob('*/busco_cwd_marker.log'))), 1)
        self.assertEqual(self.run_cli(*args).returncode, 0)
        retained.write_text('>c1\nACGT\n')
        changed = self.run_cli(*args)
        self.assertNotEqual(changed.returncode, 0)
        self.assertIn('modified', changed.stderr)

    def test_repeated_bins_argument_and_same_source_not_deduplicated(self):
        a = self.bins('custom')
        (a / 'two.fa').write_text(self.fasta.read_text())
        b = self.bins('other')
        args = self.arguments(a)
        result = self.run_cli(*args, '--bins', str(b))
        self.assertEqual(result.returncode, 0, result.stderr)
        summary = json.loads((self.root / 'out/0_log/summary.json').read_text())
        self.assertEqual(summary['qc_candidates'], 3)
        self.assertEqual(summary['representatives'], 2)
