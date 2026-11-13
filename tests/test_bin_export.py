"""Size-based export preserves whole-bin membership without tiny FASTA files."""
from pathlib import Path
import tempfile
import unittest
from test_multi_coverage import MODULE as api
from binners import export_memberships


class BinExportTests(unittest.TestCase):
    def test_inclusive_whole_bin_size_and_original_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / 'input.fa'
            text = '>a original header\n' + 'A' * 50000 + '\n>b\n' + 'C' * 50000 + '\n'
            fasta.write_text(text + '>small\nACGT\n')
            bins = root / 'bins'; bins.mkdir()
            export_memberships(api, fasta, [('a', '1'), ('b', '1'), ('small', '2')], bins, 100000)
            self.assertEqual([p.name for p in bins.iterdir()], ['1.fa'])
            self.assertEqual((bins / '1.fa').read_text(), text)

    def test_invalid_memberships_fail_before_output(self):
        for rows in [[('missing','1')], [('a','1'), ('a','2')], [('a','../escape')]]:
            with self.subTest(rows=rows), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp); fasta = root / 'input.fa'; fasta.write_text('>a\nACGT\n')
                bins = root / 'bins'; bins.mkdir()
                with self.assertRaises(ValueError):
                    export_memberships(api, fasta, rows, bins, 0)
                self.assertEqual(list(bins.iterdir()), [])
