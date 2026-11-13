"""Regression checks for scoped resources and in-memory reuse, using synthetic inputs."""
import json
import os
import sys
import types
import unittest
from pathlib import Path
from unittest.mock import patch

import test_workflow as fixtures
from test_multi_coverage import MODULE as cli
from binners import prepare
from workflow import Run


class EfficiencyTests(unittest.TestCase):
    setUp = fixtures.WorkflowTests.setUp
    tearDown = fixtures.WorkflowTests.tearDown
    run_cli = fixtures.WorkflowTests.run_cli
    base_run = fixtures.WorkflowTests.base_run
    arguments = fixtures.WorkflowTests.arguments
    state = fixtures.WorkflowTests.state
    table = fixtures.WorkflowTests.table
    source = fixtures.WorkflowTests.source

    def test_jgi_half_budget_is_local_and_other_libraries_keep_full_budget(self):
        keys = ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS')
        original = {key: os.environ.get(key) for key in keys}
        for threads, expected in ((1, 1), (3, 1), (8, 4)):
            with self.subTest(threads=threads):
                out = self.root / f'threads_{threads}'
                out.mkdir(); (out / '0_log').mkdir()
                args = types.SimpleNamespace(output=out, threads=threads, random_seed=6311, gpu='F',
                    contigs=self.fasta, bam=[self.root / 'synthetic.bam'], sample='s', min_input_length=1,
                    _contig_lengths=cli.fasta_lengths(self.fasta), _input_states={}, _normalized_depth=None)
                run = Run(cli, args, {})
                run.state = {'commands': []}
                run.lock_fd = os.open(out, os.O_RDONLY)
                execute = run.execute
                observed = []
                def jgi(label, command, **kwargs):
                    observed.append(command)
                    target = Path(command[command.index('--outputDepth') + 1])
                    text = 'contigName\tcontigLen\ttotalAvgDepth\tA\tA-var\n'
                    text += ''.join(f'{name}\t{bp}\t2.5\t2.5\t0.75\n' for name, bp in args._contig_lengths.items())
                    program = 'from pathlib import Path; Path(' + repr(str(target)) + ').write_text(' + repr(text) + ')'
                    return execute(label, [sys.executable, '-c', program], cwd=kwargs['cwd'], env=kwargs['env'])
                try:
                    run.execute = jgi
                    prepared = prepare(run)
                    self.assertEqual([r[3:] for r in prepared[3]], [['2.5', '0.75']] * 3)
                    command = observed[0]
                    for option, value in (('--percentIdentity', '97'), ('--minMapQual', '0'), ('--maxEdgeBases', '75')):
                        self.assertEqual(command[command.index(option)+1], value)
                    self.assertEqual(run.state['commands'][0]['thread_environment'],
                                     {key: str(expected if key == 'OMP_NUM_THREADS' else threads) for key in keys})
                    result = execute('after', [sys.executable, '-c',
                        'import os,json; print(json.dumps({key:os.environ[key] for key in '+repr(keys)+'}))'], capture=True)
                    self.assertEqual(json.loads(result.stdout), dict.fromkeys(keys, str(threads)))
                    self.assertFalse((out/'.work/inputs/jgi.tsv').exists())
                finally:
                    os.close(run.lock_fd)
        self.assertEqual({key: os.environ.get(key) for key in keys}, original)

    def test_tables_are_validated_once_and_reused_without_text_roundtrip(self):
        args = cli.build_parser().parse_args(self.arguments('metabat2'))
        cli.validate_args(args)
        with patch.object(cli, 'normalized_depth_rows', wraps=cli.normalized_depth_rows) as normalize:
            inputs = cli.validate_inputs(args)
            header, rows = args._normalized_depth
            expected_text = '\t'.join(header)+'\n'+''.join('\t'.join(row)+'\n' for row in rows)
            self.assertEqual(inputs['normalized_depth_bytes'], len(expected_text.encode()))
            run = Run(cli, args, inputs)
            with patch.object(cli, 'fasta_lengths', side_effect=AssertionError('reparsed FASTA')), patch.object(cli, 'normalized_depth_inputs', side_effect=AssertionError('text roundtrip')):
                fasta, actual_header, schema, actual_rows = prepare(run)
            self.assertEqual(normalize.call_count, 1)
            self.assertEqual((actual_header, actual_rows), (header, rows))
            self.assertEqual(fasta.read_text(), self.fasta.read_text())
            self.assertIsNone(args._normalized_depth)

    def test_changed_table_is_not_silently_served_from_memory(self):
        args = cli.build_parser().parse_args(self.arguments('metabat2'))
        cli.validate_args(args); inputs = cli.validate_inputs(args)
        self.depth.write_text(self.depth.read_text() + 'new\t4\n')
        with self.assertRaisesRegex(ValueError, 'input changed after validation'):
            prepare(Run(cli, args, inputs))

    def test_native_reference_fingerprints_are_read_once_before_binning(self):
        import workflow
        args = cli.build_parser().parse_args(self.arguments())
        cli.validate_args(args); inputs = cli.validate_inputs(args)
        run = Run(cli, args, inputs)
        run.out.mkdir(); (run.out/'0_log').mkdir(); (run.work/'inputs').mkdir(parents=True)
        run.state = {'stages': {}, 'bins': [], 'commands': []}
        references = []
        def bins(task, source, stage, prepared):
            result = stage/'bins'; result.mkdir(); return result
        def inventory(source, paths, stage, reference):
            references.append(reference); return []
        run.inventory = inventory
        with patch.object(workflow, 'prepare', return_value=(self.fasta, None, None, None)), patch.object(workflow, 'bin_contigs', side_effect=bins), patch.object(workflow, 'fasta_contigs', wraps=workflow.fasta_contigs) as fingerprints:
            run.sources()
        self.assertEqual(fingerprints.call_count, 1)
        self.assertEqual(len(references), 4)
        self.assertTrue(all(reference is references[0] for reference in references))
        self.assertEqual(len(references[0]), 3)

    def test_exact_qc_reuse_keeps_source_rows_and_evidence_separate(self):
        a = self.source('a'); b = self.source('b')
        partial = self.source('partial', text='>c1\n'+'ACGT'*25+'\n')
        reformatted = self.source('reformatted', text=self.fasta.read_text().replace('>c1', '>c1 description'))
        args = ['--sample','s','--bins',str(a),str(b),str(partial),str(reformatted),
                '--output',str(self.root/'out'),'--min-retain-size','1','--min-selection-size','1',
                '--qc-tools','all','--eukcc-db',str(self.eukcc_db),'--busco-db',str(self.busco_db)]
        result = self.run_cli(*args)
        self.assertEqual(result.returncode, 0, result.stderr)
        state = self.state(); rows = {r['source']:r for r in self.table('2_qc/summary.tsv')}
        self.assertEqual(len(rows), 4)
        for tool in ('eukcc','busco'):
            self.assertEqual(sum(c['stage'].startswith(tool+':') for c in state['commands']), 3)
            self.assertEqual(rows['b'][tool+'_reused_from'], 'a/one')
            self.assertEqual(state['stages'][tool+':b/one']['reused_from'], 'a/one')
            self.assertTrue(all(not rows[source][tool+'_reused_from'] for source in ('a','partial','reformatted')))
            self.assertNotEqual(rows['a'][tool+'_evidence'], rows['b'][tool+'_evidence'])
            self.assertTrue((self.root/'out'/rows['b'][tool+'_evidence']).is_dir())
        again = self.run_cli(*args)
        self.assertEqual(again.returncode, 0, again.stderr)
        evidence = self.root/'out'/rows['b']['eukcc_evidence']/'eukcc.csv'
        evidence.write_text('damaged')
        self.assertNotEqual(self.run_cli(*args).returncode, 0)

    def test_unassessable_qc_is_not_reused(self):
        self.environment['MOCK_UNASSESSABLE'] = '1'
        a, b = self.source('a'), self.source('b')
        result = self.run_cli('--sample','s','--bins',str(a),str(b),'--output',str(self.root/'out'),
                              '--min-retain-size','1','--min-selection-size','1','--eukcc-db',str(self.eukcc_db))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(sum(c['stage'].startswith('eukcc:') for c in self.state()['commands']), 2)
        self.assertTrue(all(row['eukcc_status'] == 'no_placement_markers' for row in self.table('2_qc/summary.tsv')))

    def test_failed_qc_is_retried_before_any_duplicate_can_reuse_it(self):
        a, b = self.source('fail_first_source'), self.source('duplicate')
        args = ['--sample','s','--bins',str(a),str(b),'--output',str(self.root/'out'),
                '--min-retain-size','1','--min-selection-size','1','--eukcc-db',str(self.eukcc_db)]
        self.environment['MOCK_FAIL_QC_SOURCE'] = 'fail_first_source'
        self.assertNotEqual(self.run_cli(*args).returncode, 0)
        self.environment.pop('MOCK_FAIL_QC_SOURCE')
        result = self.run_cli(*args)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(sum(c['stage']=='eukcc:fail_first_source/one' for c in self.state()['commands']), 2)
        self.assertEqual(self.state()['stages']['eukcc:duplicate/one']['reused_from'], 'fail_first_source/one')

    def test_concoct_skips_only_unused_original_matrix_output(self):
        result = self.run_cli(*self.arguments('concoct'), concoct_success=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        command = next(c['argv'] for c in self.state()['commands'] if '--composition_file' in c['argv'])
        self.assertIn('--no_original_data', command)
        self.assertNotIn('--no_cov_normalization', command)
        self.assertNotIn('--no_total_coverage', command)
