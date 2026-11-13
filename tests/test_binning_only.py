"""Binning-only output, safe continuation and two-stage equivalence."""
import csv
import json
import unittest

from test_mock_end_to_end import MockEndToEndTests
import test_cli_interface as cli_tests


class BinningOnlyTests(unittest.TestCase):
    setUp = MockEndToEndTests.setUp
    tearDown = MockEndToEndTests.tearDown
    run_cli = MockEndToEndTests.run_cli
    base_run = MockEndToEndTests.base_run
    mock_cuda = cli_tests.CliInterfaceTests.mock_cuda

    def arguments(self, output, binners='remag,vamb', retain=1):
        args = self.base_run(output)
        at = args.index('--min-selection-size')
        del args[at:at+2]
        args[args.index('--min-retain-size')+1] = str(retain)
        return [*args, '--binners', binners, '--binning-only', 'T']

    def state(self, output):
        return json.loads((output/'0_log/run.json').read_text())

    def test_gpu_binning_needs_no_qc_and_creates_only_candidates_and_logs(self):
        self.mock_cuda()
        self.environment['MOCK_CHECK_COMPACT'] = '1'
        # Observe the output while native tools run, not just after cleanup.
        for name in ('remag', 'vamb', 'seqkit'):
            tool = self.tools/name
            tool.write_text(tool.read_text().replace(
                "if (parent / '0_log/run.json').exists():",
                "if (parent / '0_log/run.json').exists():\n"
                "            assert not (parent / '2_qc').exists()\n"
                "            assert not (parent / '3_selection').exists()"))
        for tool in ('eukcc', 'busco'):
            (self.tools/tool).unlink()
        self.environment['EUKBINPICKER_EUKCC_DB'] = '/missing/database'
        self.environment['EUKBINPICKER_BUSCO_DB'] = '/missing/lineage'
        output = self.root/'gpu'
        result = self.run_cli(*self.arguments(output), '--gpu', 'T')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual({p.name for p in output.iterdir()}, {'0_log', '1_raw'})
        self.assertEqual({p.name for p in (output/'0_log').iterdir()}, {'run.json', 'run.log'})
        state = self.state(output)
        self.assertEqual(state['status'], 'binning_completed')
        self.assertTrue(state['inputs']['binning_only'])
        self.assertEqual(state['inputs']['qc_tools'], [])
        self.assertIsNone(state['inputs']['decision_qc'])
        self.assertIsNone(state['inputs']['qc_min_bin_bp'])
        self.assertEqual(state['summary']['qc_status'], 'skipped_binning_only')
        self.assertEqual(state['summary']['selection_status'], 'skipped_binning_only')
        self.assertNotIn('representatives', state['summary'])
        self.assertEqual(set(state['stages']), {'remag', 'vamb'})
        self.assertFalse({'eukcc', 'busco'} & state['inputs']['software'].keys())
        self.assertFalse(state['partial'])
        self.assertEqual(state['summary']['retained_bins'], 2)
        traces = [json.loads(line) for line in (self.root/'devices.jsonl').read_text().splitlines()]
        self.assertEqual({row['tool'] for row in traces}, {'remag', 'vamb'})
        self.assertTrue(all(row['device'] == 'cuda' for row in traces))

    def test_retention_filters_native_and_imported_whole_bins(self):
        source = self.root/'previous'; source.mkdir()
        (source/'large.fa').write_text(self.fasta.read_text())
        (source/'small.fa').write_text('>small\nACGT\n')
        original = {p.name:p.read_bytes() for p in source.iterdir()}
        for cutoff, expected in ((300, 2), (301, 0), (600000, 0)):
            output = self.root/f'retain_{cutoff}'
            result = self.run_cli(*self.arguments(output, 'remag', cutoff), '--bins', str(source))
            self.assertEqual(result.returncode, 0, result.stderr)
            state = self.state(output)
            self.assertEqual(len(state['bins']), expected)
            self.assertTrue(all(row['bp'] >= cutoff for row in state['bins']))
            self.assertEqual(len(list((output/'1_raw').rglob('*.fa'))), expected)
            self.assertEqual(state['summary']['imported_sources'], ['previous'])
        self.assertEqual({p.name:p.read_bytes() for p in source.iterdir()}, original)

    def test_same_mode_resume_verifies_outputs_and_rejects_mode_switch(self):
        output = self.root/'out'; args = self.arguments(output)
        first = self.run_cli(*args)
        self.assertEqual(first.returncode, 0, first.stderr)
        original = self.state(output)
        again = self.run_cli(*args)
        self.assertEqual(again.returncode, 0, again.stderr)
        self.assertIn('already complete', again.stderr)
        self.assertEqual(self.state(output)['commands'], original['commands'])
        switched = self.run_cli(*self.base_run(output), '--binners', 'remag,vamb', '--eukcc-db', str(self.eukcc_db))
        self.assertNotEqual(switched.returncode, 0)
        self.assertIn('use a new output directory', switched.stderr)
        candidate = next((output/'1_raw').rglob('*.fa'))
        candidate.write_text(candidate.read_text()+'N\n')
        corrupt = self.run_cli(*args)
        self.assertNotEqual(corrupt.returncode, 0)
        self.assertIn('missing or modified', corrupt.stderr)

    def test_inapplicable_options_and_bins_only_are_rejected_before_output(self):
        output = self.root/'bad'
        for option,value in (('--qc-tools','all'), ('--eukcc-db',str(self.eukcc_db)),
                             ('--min-selection-size','1'), ('--min-overlap','0.5')):
            result = self.run_cli(*self.arguments(output), f'{option}={value}')
            self.assertEqual(result.returncode, 2)
            self.assertIn('skips QC and selection; omit', result.stderr)
            self.assertFalse(output.exists())
        source = self.root/'existing'; source.mkdir(); (source/'one.fa').write_text(self.fasta.read_text())
        result = self.run_cli('--sample','s','--bins',str(source),'--output',str(output),'--binning-only','T')
        self.assertEqual(result.returncode, 2)
        self.assertIn('requires --contigs and --coverage', result.stderr)
        self.assertFalse(output.exists())

    def test_partial_failure_is_explicit_and_all_failed_is_an_error(self):
        output = self.root/'partial'
        result = self.run_cli(*self.arguments(output, 'remag,concoct'))
        self.assertEqual(result.returncode, 0, result.stderr)
        state = self.state(output)
        self.assertEqual(state['status'], 'binning_completed')
        self.assertTrue(state['partial'])
        self.assertEqual(state['summary']['native_binner_status'], {'remag':'completed', 'concoct':'failed'})
        self.assertTrue((output/'.work').is_dir())
        self.assertFalse((output/'2_qc').exists())
        failed = self.root/'failed'
        result = self.run_cli(*self.arguments(failed, 'concoct'))
        self.assertEqual(result.returncode, 2, result.stderr)
        self.assertEqual(self.state(failed)['status'], 'failed')
        self.assertNotIn('summary', self.state(failed))

    def test_empty_filtered_input_completes_without_qc_or_native_commands(self):
        output = self.root/'empty'
        args = self.arguments(output)
        args[args.index('--min-input-length')+1] = '1000'
        result = self.run_cli(*args)
        self.assertEqual(result.returncode, 0, result.stderr)
        state = self.state(output)
        self.assertEqual(state['commands'], [])
        self.assertEqual(state['summary']['retained_bins'], 0)
        self.assertEqual(set(state['summary']['native_binner_status'].values()), {'skipped_empty_input'})
        self.assertEqual({p.name for p in output.iterdir()}, {'0_log','1_raw'})

    def test_gpu_then_cpu_import_matches_default_full_workflow(self):
        self.mock_cuda()
        first, second, full = (self.root/name for name in ('gpu','cpu','full'))
        result = self.run_cli(*self.arguments(first), '--gpu', 'T', expected_seed=6311)
        self.assertEqual(result.returncode, 0, result.stderr)
        sources = [first/'1_raw'/name for name in ('remag','vamb')]
        result = self.run_cli(*self.base_run(second), '--binners', 'concoct,metabat2', '--bins', *map(str,sources),
                              '--binning-only','F','--eukcc-db',str(self.eukcc_db), concoct_success=True, expected_seed=6311)
        self.assertEqual(result.returncode, 0, result.stderr)
        result = self.run_cli(*self.base_run(full), '--gpu', 'T', '--eukcc-db',str(self.eukcc_db),
                              concoct_success=True, expected_seed=6311)
        self.assertEqual(result.returncode, 0, result.stderr)
        a, b = self.state(second), self.state(full)
        for state in (a,b):
            self.assertEqual(state['status'], 'completed')
            self.assertFalse(state['inputs']['binning_only'])
        self.assertEqual(a['summary']['imported_sources'], ['remag','vamb'])
        self.assertEqual(set(a['summary']['native_binner_status']), {'concoct','metabat2'})
        key = lambda rows: sorted((r['source'],r['original_bin_id'],r['sha256'],r['bp']) for r in rows)
        self.assertEqual(key(a['bins']), key(b['bins']))
        for row in a['bins']:
            if row['source'] in ('remag','vamb'):
                self.assertTrue(row['source_path'].startswith(str((first/'1_raw'/row['source']).resolve())))
        # Compare stable scientific fields, not output paths or reuse provenance.
        def selected(output):
            with (output/'3_selection/representatives.tsv').open() as handle:
                return [{k:r[k] for k in ('bin_id','source','original_bin_id','completeness','contamination','quality_score')}
                        for r in csv.DictReader(handle,delimiter='\t')]
        self.assertEqual(selected(second), selected(full))
        for p in (second/'3_selection/bins').glob('*.fa'):
            self.assertEqual(p.read_bytes(), (full/'3_selection/bins'/p.name).read_bytes())
        self.assertEqual(self.state(first)['status'], 'binning_completed')


if __name__ == '__main__':
    unittest.main()
