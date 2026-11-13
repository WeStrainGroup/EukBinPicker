"""Public contracts: compact live output, provenance, safe recovery and selection."""
import csv
import fcntl
import json
import os
import types
from pathlib import Path
import unittest

from test_mock_end_to_end import MockEndToEndTests
from test_multi_coverage import MODULE as cli
from selection import select_rows


class WorkflowTests(unittest.TestCase):
    setUp = MockEndToEndTests.setUp
    tearDown = MockEndToEndTests.tearDown
    run_cli = MockEndToEndTests.run_cli
    base_run = MockEndToEndTests.base_run

    def arguments(self, binners='all'):
        return [*self.base_run(self.root / 'out'), '--eukcc-db', str(self.eukcc_db), '--binners', binners]

    def state(self):
        return json.loads((self.root / 'out/0_log/run.json').read_text())

    def table(self, relative):
        with (self.root / 'out' / relative).open() as handle:
            return list(csv.DictReader(handle, delimiter='\t'))

    def source(self, folder='previous', name='one.fa', text=None):
        path = self.root / folder
        path.mkdir(parents=True, exist_ok=True)
        (path / name).write_text(text or self.fasta.read_text())
        return path

    def test_all_four_native_and_compact_live_and_final_outputs(self):
        self.environment['MOCK_CHECK_COMPACT'] = '1'
        result = self.run_cli(*self.arguments(), concoct_success=True, expected_seed=6311)
        self.assertEqual(result.returncode, 0, result.stderr)
        out = self.root / 'out'
        self.assertEqual({p.name for p in out.iterdir()}, {'0_log','1_raw','2_qc','3_selection'})
        self.assertEqual({p.name for p in (out/'0_log').iterdir()}, {'run.json','run.log'})
        state = self.state()
        self.assertEqual(state['status'], 'completed')
        self.assertEqual(set(state['summary']['native_binner_status'].values()), {'completed'})
        self.assertEqual(state['summary']['qc_candidates'], 4)
        self.assertEqual(state['summary']['representatives'], 1)
        reps = self.table('3_selection/representatives.tsv')
        self.assertEqual(reps[0]['bin_id'], 'bin_1')
        self.assertEqual(reps[0]['source'], 'metabat2')
        self.assertEqual((out/'3_selection/bins/bin_1.fa').read_bytes(), self.fasta.read_bytes())
        self.assertEqual(len(self.table('3_selection/selection.tsv')), 4)
        self.assertEqual({r['representative'] for r in self.table('3_selection/selection.tsv')}, {'bin_1'})
        for row in self.table('2_qc/summary.tsv'):
            evidence = out/'2_qc/evidence'/row['source']/row['bin_id']/'eukcc'
            self.assertTrue((evidence/'eukcc.csv').is_file())
            self.assertTrue((evidence/'markers.json').is_file())
        again = self.run_cli(*self.arguments(), concoct_success=True)
        self.assertEqual(again.returncode, 0, again.stderr)
        self.assertIn('already complete', again.stderr)

    def test_bins_only_colliding_names_and_distinct_evidence(self):
        a, b = self.source('a/custom'), self.source('b/custom')
        result = self.run_cli('--sample','sample1','--bins',str(a),str(b),'--output',str(self.root/'out'),
                              '--min-retain-size','1','--min-selection-size','1','--qc-tools','all','--qc-decision','busco',
                              '--eukcc-db',str(self.eukcc_db),'--busco-db',str(self.busco_db))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.state()['summary']['imported_sources'], ['custom','custom_2'])
        for source in ('custom','custom_2'):
            p=self.root/'out/2_qc/evidence'/source/'one'
            self.assertTrue((p/'eukcc/eukcc.csv').is_file())
            self.assertTrue((p/'busco/short_summary.txt').is_file())
            self.assertTrue((p/'busco/full_table.tsv.gz').is_file())
        self.assertEqual(self.state()['summary']['native_binner_status'], {})
        # Preserve the established source/bin tie-break despite the new display IDs.
        self.assertEqual(self.table('3_selection/representatives.tsv')[0]['source'], 'custom_2')

    def test_concoct_integer_coverage_is_written_as_float(self):
        from binners import bin_contigs
        from depth_schema import identify_depth_schema
        stage = self.root / 'adapter'; stage.mkdir()
        header = ['contigName', 'contigLen', 'totalAvgDepth', 'A', 'A-var']
        schema = identify_depth_schema(header, 'sample1', require_jgi=True)
        args = types.SimpleNamespace(binner_parameters={'concoct': []}, min_input_length=1,
                                     threads=4, random_seed=6311, bam=[])
        class StopAfterCoverage(Exception):
            pass
        def execute(label, command, **kwargs):
            if command[0] == 'cut_up_fasta.py':
                kwargs['stdout'].write_text('>x\nACGT\n>y\nACGT\n')
                (stage / 'contigs_10K.bed').write_text('x\t0\t4\tx\ny\t0\t4\ty\n')
            elif command[0] == 'concoct':
                raise StopAfterCoverage
        run = types.SimpleNamespace(args=args, api=cli, execute=execute)
        rows = [['x', '4', '20', '20', '1'], ['y', '4', '0', '0', '0']]
        with self.assertRaises(StopAfterCoverage):
            bin_contigs(run, 'concoct', stage, (self.fasta, header, schema, rows))
        lines = (stage / 'fragment_coverage.tsv').read_text().splitlines()
        self.assertEqual(lines, ['contig\tA', 'x\t20.0', 'y\t0.0'])

    def test_retention_and_qc_thresholds_remain_separate(self):
        a=self.source()
        self.source('previous','small.fa','>small\nACGT\n')
        result=self.run_cli('--sample','sample1','--bins',str(a),'--output',str(self.root/'out'),
                            '--min-retain-size','5','--min-selection-size','301','--eukcc-db',str(self.eukcc_db))
        self.assertEqual(result.returncode,0,result.stderr)
        rows=self.table('2_qc/summary.tsv')
        self.assertEqual(len(rows),1)
        self.assertEqual(rows[0]['qc_status'],'below_min_selection_size')
        self.assertEqual(self.table('3_selection/selection.tsv')[0]['reason'],'below_min_selection_size')
        self.assertEqual(self.state()['summary']['representatives'],0)

    def test_default_retention_is_inclusive_and_zero_keeps_small_imports(self):
        source = self.root / 'sizes'
        source.mkdir()
        for size in (99999, 100000, 100001):
            (source / f'{size}.fa').write_text(f'>c{size}\n' + 'A' * size + '\n')
        originals = {p.name: p.read_bytes() for p in source.iterdir()}
        for name, extra, expected in (('out', [], {'100000.fa', '100001.fa'}),
                                      ('all', ['--min-retain-size', '0'], set(originals))):
            result = self.run_cli('--sample', 's', '--bins', str(source),
                                  '--output', str(self.root / name),
                                  '--eukcc-db', str(self.eukcc_db), *extra)
            self.assertEqual(result.returncode, 0, result.stderr)
            retained = self.root / name / '1_raw/sizes'
            self.assertEqual({p.name for p in retained.iterdir()}, expected)
        self.assertEqual({p.name: p.read_bytes() for p in source.iterdir()}, originals)
        args = cli.build_parser().parse_args(['--sample', 's', '--bins', str(source),
                                             '--output', str(self.root / 'unused')])
        self.assertEqual(args.min_retain_size, 100000)
        self.assertEqual(args.min_selection_size, 500000)

    def test_vamb_export_uses_retention_minimum_without_changing_clustering(self):
        from binners import bin_contigs
        from depth_schema import identify_depth_schema
        schema = identify_depth_schema(['contigName', 'contigLen', 'totalAvgDepth', 'A', 'A-var'], 's', require_jgi=True)
        for minimum, expected in ((0, '1'), (100000, '100000')):
            stage = self.root / f'export-{minimum}'
            stage.mkdir()
            commands = []
            def execute(label, command, **kwargs):
                commands.append(command)
                (Path(command[command.index('--outdir') + 1]) / 'bins').mkdir(parents=True)
            args = types.SimpleNamespace(binner_parameters={'vamb': []}, min_input_length=1500,
                min_retain_size=minimum, threads=4, random_seed=6311, gpu='F')
            api = types.SimpleNamespace(SCRIPTS=self.root, executable=lambda *args: 'vamb')
            run = types.SimpleNamespace(args=args, api=api, execute=execute, state={'stages': {'vamb': {}}})
            rows = [['c1', '2000', '1', '1', '0'], ['c2', '2000', '2', '2', '0']]
            bin_contigs(run, 'vamb', stage, (self.fasta, [], schema, rows))
            command = commands[0]
            self.assertEqual(command[command.index('--minfasta') + 1], expected)
            self.assertEqual(command[command.index('-m') + 1], '1500')
            self.assertEqual(command[command.index('--seed') + 1], '6311')

    def test_partial_binner_failure_preserves_successful_sources(self):
        result=self.run_cli(*self.arguments('concoct,metabat2'))
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertTrue(self.state()['partial'])
        self.assertEqual(self.state()['summary']['native_binner_status'],{'concoct':'failed','metabat2':'completed'})
        self.assertFalse((self.root/'out/.work/binning/metabat2').exists())
        self.assertTrue((self.root/'out/.work/binning/concoct').exists())
        self.assertEqual(self.state()['summary']['representatives'],1)

    def test_all_binners_failed_is_not_a_success(self):
        result=self.run_cli(*self.arguments('concoct'))
        self.assertNotEqual(result.returncode,0)
        self.assertEqual(self.state()['status'],'failed')
        self.assertIn('all requested binners failed',result.stderr)
        again=self.run_cli(*self.arguments('concoct'),concoct_success=True)
        self.assertEqual(again.returncode,0,again.stderr)
        self.assertEqual(self.state()['status'],'completed')

    def test_combined_empty_native_input_keeps_imports_intact(self):
        a=self.source(text='>outside_input\nACGT\n')
        result=self.run_cli(*self.arguments('metabat2'),'--min-input-length','1500','--bins',str(a))
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(self.state()['summary']['native_binner_status'],{'metabat2':'skipped_empty_input'})
        self.assertEqual(self.state()['summary']['representatives'],1)
        self.assertEqual((self.root/'out/1_raw/previous/one.fa').read_text(),'>outside_input\nACGT\n')

    def test_failed_qc_resumes_without_rerunning_completed_work(self):
        self.environment['MOCK_FAIL_QC_SOURCE']='previous'
        a=self.source()
        args=[*self.arguments('metabat2'),'--bins',str(a)]
        result=self.run_cli(*args)
        self.assertNotEqual(result.returncode,0)
        before=self.state()
        self.assertEqual(before['status'],'failed')
        self.assertEqual(before['stages']['eukcc:previous/one']['status'],'failed')
        self.environment.pop('MOCK_FAIL_QC_SOURCE')
        result=self.run_cli(*args)
        self.assertEqual(result.returncode,0,result.stderr)
        commands=self.state()['commands']
        self.assertEqual(sum(c['stage']=='metabat2' for c in commands),1)
        self.assertEqual(sum(c['stage'].startswith('eukcc:metabat2/') for c in commands),1)
        self.assertEqual(sum(c['stage']=='eukcc:previous/one' for c in commands),2)
        self.assertFalse((self.root/'out/.work').exists())

    def test_unassessable_qc_has_explicit_status_and_no_representative(self):
        self.environment['MOCK_UNASSESSABLE']='1'
        result=self.run_cli(*self.arguments('metabat2'))
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(self.table('2_qc/summary.tsv')[0]['eukcc_status'],'no_placement_markers')
        self.assertEqual(self.table('3_selection/selection.tsv')[0]['reason'],'qc_unassessable')
        self.assertEqual(self.state()['summary']['representatives'],0)

    def test_reuse_checks_parameters_files_and_evidence(self):
        args=self.arguments('metabat2')
        result=self.run_cli(*args)
        self.assertEqual(result.returncode,0,result.stderr)
        changed=self.run_cli(*args,'--random-seed','6312')
        self.assertNotEqual(changed.returncode,0)
        self.assertIn('differ',changed.stderr)
        p=next((self.root/'out/2_qc/evidence').rglob('eukcc.csv'))
        original=p.read_bytes();p.write_text('modified')
        corrupt=self.run_cli(*args)
        self.assertNotEqual(corrupt.returncode,0)
        self.assertIn('modified',corrupt.stderr)
        p.write_bytes(original)
        (self.root/'out/3_selection/bins/bin_1.fa').unlink()
        self.assertNotEqual(self.run_cli(*args).returncode,0)

    def test_directory_lock_blocks_concurrent_writers(self):
        out=self.root/'out';out.mkdir()
        fd=os.open(out,os.O_RDONLY)
        try:
            fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB)
            result=self.run_cli(*self.arguments('metabat2'))
            self.assertNotEqual(result.returncode,0)
            self.assertIn('another process',result.stderr)
        finally:os.close(fd)

    def test_input_boundary_and_multiple_coverage_dimensions(self):
        self.fasta.write_text('>short\nAAA\n>boundary\nACGT\n>long\nAAAAA\n')
        self.depth.write_text('contigName\tcontigLen\ttotalAvgDepth\tA\tA-var\tB\tB-var\nshort\t3\t3\t1\t2\t2\t3\nboundary\t4\t5\t2\t3\t3\t4\nlong\t5\t7\t3\t4\t4\t5\n')
        result=self.run_cli(*self.arguments(),'--min-input-length','4',concoct_success=True,
                            expected_coverage_columns=('A','B'))
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(self.state()['summary']['native_binner_status'],dict.fromkeys(cli.BINNERS,'completed'))
        for p in (self.root/'out/1_raw').rglob('*.fa'):
            self.assertNotIn('>short',p.read_text())
            self.assertIn('>boundary',p.read_text())

    def test_selection_chain_same_source_and_shared_bp(self):
        rows=[]
        for name,source,ids,c in [('a','remag',('c1','c2'),95),('b','metabat2',('c2','c3'),90),('c','concoct',('c3','c4'),85)]:
            path=self.root/f'{name}.fa';path.write_text(''.join(f'>{i}\nACGT\n' for i in ids))
            rows.append(dict(bin_id=name,sample_rank=1,study='s',sample='s',binner=source,bin_path=str(path),
                        decision_qc='eukcc',contamination_source='eukcc_contamination',qc_pass='yes',
                        completeness=c,contamination=0,contigs=2,n50=4,bp=8))
        reps,alts=select_rows(rows,0.5,5,50,10,1,'s','s')
        self.assertEqual([r['representative_bin_id'] for r in reps],['a','c'])
        self.assertEqual(alts[0]['alternative_bin_id'],'b')
        reps,_=select_rows([rows[0],{**rows[1],'binner':'remag'}],0.5,5,50,10,1,'s','s')
        self.assertEqual([r['shared_contig_fraction'] for r in reps],['0.500000','0.500000'])

    def test_source_named_inputs_does_not_collide_with_work_directory(self):
        a=self.source('inputs')
        result=self.run_cli(*self.arguments('metabat2'),'--bins',str(a))
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertTrue((self.root/'out/1_raw/inputs/one.fa').is_file())
        self.assertFalse((self.root/'out/.work').exists())

    def test_invalid_native_bin_does_not_leak_into_public_candidates(self):
        a=self.source()
        p=self.tools/'metabat2'
        p.write_text(p.read_text().replace('binner_records(source, "metabat2")', '[("foreign", "AAAA")]'))
        result=self.run_cli(*self.arguments('metabat2'),'--bins',str(a))
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(self.state()['summary']['native_binner_status'],{'metabat2':'failed'})
        self.assertFalse((self.root/'out/1_raw/metabat2').exists())
        self.assertEqual({r['source'] for r in self.table('2_qc/summary.tsv')},{'previous'})

    def test_gzip_import_and_native_source_name_collision(self):
        import gzip
        a=self.root/'metabat2';a.mkdir()
        with gzip.open(a/'same.fa.gz','wt') as h:h.write(self.fasta.read_text())
        result=self.run_cli(*self.arguments('metabat2'),'--bins',str(a))
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(self.state()['summary']['imported_sources'],['metabat2_2'])
        self.assertEqual((self.root/'out/1_raw/metabat2_2/same.fa').read_text(),self.fasta.read_text())
        rows=self.table('2_qc/summary.tsv')
        self.assertIn('same.fa.gz',[r['original_filename'] for r in rows])

    def test_termination_records_interruption_and_resumes(self):
        import signal,subprocess,time,sys
        tool=self.tools/'eukcc'
        tool.write_text(tool.read_text().replace('        bin_path = args[-1]',
            '        if os.environ.get("MOCK_SLEEP"):\n            import time; time.sleep(30)\n        bin_path = args[-1]'))
        env=dict(self.environment,MOCK_SLEEP='1')
        process=subprocess.Popen([sys.executable,str(cli.__file__),*self.arguments('metabat2')],env=env,
                                 stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        try:
            deadline=time.monotonic()+15
            while time.monotonic()<deadline:
                if (self.root/'out/0_log/run.json').exists():
                    state=self.state()
                    if any(c['stage'].startswith('eukcc:') for c in state['commands']):break
                time.sleep(0.05)
            else:self.fail('QC did not start')
            process.send_signal(signal.SIGTERM)
            _,err=process.communicate(timeout=15)
            self.assertEqual(process.returncode,130,err.decode())
            self.assertEqual(self.state()['status'],'interrupted')
            result=self.run_cli(*self.arguments('metabat2'))
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertEqual(sum(c['stage']=='metabat2' for c in self.state()['commands']),1)
        finally:
            if process.poll() is None:process.kill();process.communicate()
