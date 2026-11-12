"""One run record, one log, and short-lived work for each active stage."""
import csv
import fcntl
import gzip
import io
import math
import importlib.metadata
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import signal
from datetime import datetime, timezone

from bam_inputs import tool_command
from binners import prepare, bin_contigs
from quality import assess
from bin_files import is_fasta, original_bin_id, safe_component
from selection import fasta_contigs, select_rows


def now():
    return datetime.now(timezone.utc).isoformat()


class Run:
    def __init__(self, api, args, inputs):
        self.api, self.args, self.inputs = api, args, inputs
        self.out = args.output.resolve()
        self.work = self.out / '.work'
        self.record = self.out / '0_log/run.json'
        self.log = self.out / '0_log/run.log'
        self.last_output = ''
        self.state = {}
        self.env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONHASHSEED=str(args.random_seed),
                        OMP_NUM_THREADS=str(args.threads), MKL_NUM_THREADS=str(args.threads))
        if args.gpu == 'F':
            self.env['CUDA_VISIBLE_DEVICES'] = ''

    def save(self):
        self.state['updated_at'] = now()
        self.api.atomic_text(self.record, json.dumps(self.state, indent=2, sort_keys=True) + '\n')

    def message(self, text):
        self.api.note(text)
        with self.log.open('a') as handle:
            handle.write(f'[{now()}] {text}\n')

    def execute(self, label, command, *, tool=None, cwd=None, stdout=None, capture=False, check=True, env=None):
        command = list(map(str, command))
        if tool:
            command = tool_command(command, self.api.environment_value(tool.upper() + '_ENV'),
                                   self.api.environment_value('CONDA_EXE') or 'conda')
        work = cwd or self.work
        work.mkdir(parents=True, exist_ok=True)
        cache = self.work / 'cache'
        cache.mkdir(exist_ok=True)
        environment = dict(self.env, TMPDIR=str(work), TMP=str(work), TEMP=str(work),
                           XDG_CACHE_HOME=str(cache), MPLCONFIGDIR=str(cache / 'matplotlib'),
                           NUMBA_CACHE_DIR=str(cache / 'numba'), TORCH_HOME=str(cache / 'torch'))
        environment.update(env or {})
        entry = {'stage': label, 'argv': command, 'cwd': str(work), 'started_at': now()}
        self.state['commands'].append(entry)
        self.save()
        self.message(f'{label}: {shlex.join(command)}')
        with self.log.open('ab') as log:
            offset = log.tell()
            target = stdout.open('wb') if stdout else None
            try:
                process = subprocess.Popen(command, cwd=work, env=environment, stderr=log,
                    stdout=target if target else (subprocess.PIPE if capture else log), start_new_session=True, pass_fds=(self.lock_fd,))
                try:
                    output, _ = process.communicate()
                except BaseException:
                    # Stop the whole native process group before releasing the output lock.
                    try:
                        os.killpg(process.pid, signal.SIGTERM)
                    except ProcessLookupError:
                        pass
                    try:
                        process.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
                    entry.update(exit_code=process.returncode, finished_at=now(), interrupted=True)
                    self.save()
                    raise
                result = subprocess.CompletedProcess(command, process.returncode, output)
            finally:
                if target:
                    target.close()
        with self.log.open('rb') as handle:
            handle.seek(max(offset, self.log.stat().st_size - 128000))
            self.last_output = handle.read().decode(errors='replace')
        entry.update(exit_code=result.returncode, finished_at=now())
        self.save()
        if check and result.returncode:
            raise ValueError(f'{label} exited with {result.returncode}; see 0_log/run.log')
        return result

    def checksums(self, paths):
        return {str(p.relative_to(self.out)): self.api.sha256_file(p) for p in sorted(paths) if p.is_file()}

    def verify(self, checksums):
        for name, expected in checksums.items():
            path = self.out / name
            if not path.is_file() or self.api.sha256_file(path) != expected:
                raise ValueError(f'recorded result is missing or modified: {name}; restore it or use a new output directory')

    def new_stage(self, name):
        stage = self.work / name
        if stage.exists():
            shutil.rmtree(stage)
        stage.mkdir(parents=True)
        return stage

    def inventory(self, source, paths, stage, native_reference=None):
        paths = sorted(paths)
        names = [safe_component(original_bin_id(path)) for path in paths]
        if len(set(names)) != len(names):
            raise ValueError(f'duplicate bin filenames within source {source}')
        rows, seen = [], set()
        # SeqKit accepts a file list; write one per active source, not per bin.
        if not paths:
            return rows
        listing = stage / 'seqkit_inputs.txt'
        listing.write_text(''.join(str(p.resolve()) + '\n' for p in paths))
        result = self.execute(f'stats:{source}', ['seqkit', 'stats', '-a', '-T', '-j', str(self.args.threads),
                                                '-X', str(listing)], cwd=stage, capture=True)
        stats = {Path(r['file']).resolve(): r for r in csv.DictReader(io.StringIO(result.stdout.decode()), delimiter='\t')}
        if set(stats) != {p.resolve() for p in paths}:
            raise ValueError(f'SeqKit result does not match source {source}')
        listing.unlink()
        destination = self.out / '1_raw' / source
        if destination.exists():
            shutil.rmtree(destination)
        destination.mkdir(parents=True)
        for path, identifier in zip(paths, names):
            contigs = fasta_contigs(path)
            if native_reference is not None:
                if any(key not in native_reference for key in contigs):
                    raise ValueError(f'{source}/{identifier}: native bin contains a foreign or changed contig')
                if seen.intersection(contigs):
                    raise ValueError(f'{source}: a contig occurs in more than one native bin')
                seen.update(contigs)
            lengths = sorted(contigs.values(), reverse=True)
            bp, total, n50 = sum(lengths), 0, 0
            for size in lengths:
                total += size
                if total * 2 >= bp:
                    n50 = size
                    break
            stat = stats[path.resolve()]
            values = tuple(int(stat[k].replace(',', '')) for k in ('num_seqs', 'sum_len', 'N50', 'max_len'))
            if values != (len(lengths), bp, n50, max(lengths)):
                raise ValueError(f'SeqKit/FASTA statistics disagree: {path}')
            if bp < self.args.min_retain_length:
                continue
            target = destination / (path.name[:-3] if path.name.lower().endswith('.gz') else path.name)
            if path.name.lower().endswith('.gz'):
                with gzip.open(path, 'rb') as src, target.open('wb') as dst:
                    shutil.copyfileobj(src, dst)
            elif native_reference is not None:
                shutil.move(str(path), target)
            else:
                shutil.copy2(path, target)
            rows.append(dict(source=source, original_bin_id=identifier, original_filename=path.name, sample=self.args.sample,
                             source_path=str(path) if native_reference is None else '',
                             bin_path=str(target.relative_to(self.out)), sha256=self.api.sha256_file(target),
                             contigs=len(lengths), bp=bp, n50=n50, longest_contig=max(lengths)))
        return rows

    def inputs_identity(self):
        api, args = self.api, self.args
        return {**self.inputs, 'version': api.VERSION, 'threads': args.threads,
                'cuda_devices': getattr(args, 'cuda_devices', {}),
                'contigs_path': str(args.contigs.resolve()) if args.contigs else None,
                'contigs_sha256': api.sha256_file(args.contigs) if args.contigs else None,
                'depth_inputs': [{'path': str(x.path), 'sha256': api.sha256_file(x.path)}
                                 for x in (api.parse_depth_inputs(args.depth) if args.depth else ())],
                'environments': {name: api.environment_value(name) for name in
                    ('CONDA_EXE', 'REMAG_ENV', 'CONCOCT_ENV', 'METABAT2_ENV', 'EUKCC_ENV', 'BUSCO_ENV')},
                'code_sha256': {p.name: api.sha256_file(p) for p in sorted(api.SCRIPTS.glob('*.py'))},
                'software': self.software_identity(),
                'cli_sha256': api.sha256_file(Path(api.__file__))}

    def software_identity(self):
        records = {}
        tools = (*self.args.binner_list, *self.args.qc_list, 'seqkit')
        for name in tools:
            environment = self.api.environment_value(name.upper() + '_ENV')
            prefix = Path(environment) if environment and '/' in environment else None
            executable = prefix / 'bin' / name if prefix else Path(shutil.which(name) or name)
            record = {'environment': environment, 'executable': str(executable)}
            if executable.is_file():
                record['sha256'] = self.api.sha256_file(executable)
            if prefix and (prefix / 'conda-meta/history').is_file():
                record['conda_history_sha256'] = self.api.sha256_file(prefix / 'conda-meta/history')
                for metadata in (prefix / 'conda-meta').glob(name + '-*.json'):
                    package = json.loads(metadata.read_text())
                    if package.get('name') == name:
                        record['version'] = package.get('version')
            if 'version' not in record:
                try:
                    record['version'] = importlib.metadata.version(name)
                except importlib.metadata.PackageNotFoundError:
                    record['version'] = 'not_available'
            records[name] = record
        return records

    def start(self):
        identity = self.inputs_identity()
        if self.record.is_file():
            self.state = json.loads(self.record.read_text())
            if self.state.get('inputs') != identity:
                raise ValueError('output inputs, parameters or software differ; use a new output directory')
            if self.state.get('status') == 'completed':
                self.verify(self.state['checksums'])
                expected = set(self.state['checksums'])
                actual = {str(p.relative_to(self.out)) for folder in ('1_raw', '2_qc', '3_selection')
                          for p in (self.out / folder).rglob('*') if p.is_file()}
                if actual != expected:
                    raise ValueError('completed output contains missing or unexpected files; use a new output directory')
                self.message(f'already complete: {self.out}')
                return False
            for state in self.state['stages'].values():
                if state.get('status') == 'completed':
                    self.verify(state.get('checksums', {}))
        elif any(self.out.iterdir()):
            raise ValueError('output lacks matching run metadata; use a new output directory')
        else:
            self.state = dict(schema_version=4, tool='EukBinPicker', version=self.api.VERSION,
                              inputs=identity, stages={}, commands=[], bins=[], started_at=now())
        (self.out / '0_log').mkdir(exist_ok=True)
        for name in ('1_raw', '2_qc', '3_selection'):
            (self.out / name).mkdir(exist_ok=True)
        self.work.mkdir(exist_ok=True)
        self.state.update(status='running', error=None)
        self.save()
        return True

    def sources(self):
        args = self.args
        prepared = None
        for source in (*args.binner_list, *args.bin_sources):
            previous = self.state['stages'].get(source, {})
            if previous.get('status') in {'completed', 'skipped_empty_input'}:
                continue
            native = source in args.binner_list
            if native and not self.inputs['eligible_contigs']:
                self.state['stages'][source] = {'status': 'skipped_empty_input'}
                self.save()
                continue
            stage = self.new_stage('binning/' + source)
            self.state['stages'][source] = {'status': 'running'}
            self.save()
            try:
                if native:
                    if prepared is None:
                        prepared = prepare(self)
                    bins = bin_contigs(self, source, stage, prepared)
                    paths = [p for p in bins.iterdir() if is_fasta(p)]
                    reference = fasta_contigs(prepared[0])
                else:
                    paths, reference = args.bin_sources[source], None
                rows = self.inventory(source, paths, stage, reference)
                self.state['bins'] = [r for r in self.state['bins'] if r['source'] != source] + rows
                self.state['stages'][source].update(status='completed', retained_bins=len(rows),
                    checksums=self.checksums((self.out / '1_raw' / source).glob('*')))
                self.save()
                shutil.rmtree(stage)
            except (OSError, ValueError) as error:
                incomplete = self.out / '1_raw' / source
                if incomplete.exists():
                    shutil.move(str(incomplete), stage / 'incomplete_bins')
                self.state['stages'][source].update(status='failed', error=str(error))
                self.save()
                if not native:
                    raise
                self.message(f'{source} failed: {error}; continuing with other sources')
        if prepared is not None:
            shutil.rmtree(self.work / 'inputs')
        binning = self.work / 'binning'
        if binning.exists() and not any(binning.iterdir()):
            binning.rmdir()
        completed = [s for s in args.binner_list if self.state['stages'][s]['status'] == 'completed']
        if args.binner_list and self.inputs['eligible_contigs'] and not completed and not args.bin_sources:
            raise ValueError('all requested binners failed; see 0_log/run.log')
        if args.binner_list and self.inputs['eligible_contigs'] and not completed and args.bin_sources:
            self.message('all requested binners failed; continuing with imported bin sets')

    def qc(self):
        for row in self.state['bins']:
            if row['bp'] < self.args.min_selection_length:
                row['qc_status'] = 'below_min_selection_length'
                continue
            for tool in self.args.qc_list:
                key = f'{tool}:{row["source"]}/{row["original_bin_id"]}'
                previous = self.state['stages'].get(key, {})
                if previous.get('status') == 'completed':
                    continue
                stage = self.new_stage('qc')
                self.state['stages'][key] = {'status': 'running'}
                self.save()
                try:
                    row.update(assess(self, row, tool, stage))
                    evidence = self.out / '2_qc/evidence' / row['source'] / row['original_bin_id'] / tool
                    self.state['stages'][key] = dict(status='completed',
                        checksums=self.checksums(evidence.rglob('*')) if evidence.exists() else {})
                    self.save()
                    shutil.rmtree(stage)
                except (OSError, ValueError) as error:
                    row[f'{tool}_status'] = 'failed'
                    self.state['stages'][key] = dict(status='failed', error=str(error))
                    self.save()
                    raise
            decision = self.args.qc_decision
            row['qc_status'] = row.get(f'{decision}_status', 'not_assessed')
        self.save()

    def export(self):
        api, args = self.api, self.args
        rows = self.state['bins']
        eligible = []
        for row in rows:
            tool = args.qc_decision
            c = row.get(f'{tool}_completeness' if tool == 'eukcc' else 'busco_complete', '')
            d = row.get(f'{tool}_contamination' if tool == 'eukcc' else 'busco_duplicated', '')
            passed, score = api.qc_metrics(c, d, min_completeness=args.min_selection_completeness,
                           max_contamination=args.max_selection_contamination, contamination_penalty=args.contamination_penalty)
            if row.get('qc_status') == 'ok':
                try:
                    valid = (math.isfinite(float(c)) and 0 <= float(c) <= 100
                             and math.isfinite(float(d)) and float(d) >= 0)
                except (TypeError, ValueError):
                    valid = False
                if not valid:
                    row['qc_status'] = 'invalid_metrics'
            if row.get('qc_status') != 'ok':
                passed, score = 'no', ''
            row.update(decision_qc=tool, completeness=c, contamination=d, quality_score=score, qc_pass=passed)
            eligible.append({**row, 'bin_id': row['source'] + '/' + row['original_bin_id'],
                'sort_id': row['source'] + '__' + row['original_bin_id'],
                'binner': row['source'], 'sample_rank': 1, 'study': 'standalone',
                'bin_path': str(self.out / row['bin_path']),
                'contamination_source': 'eukcc_contamination' if tool == 'eukcc' else 'busco_duplicated'})
        reps, alternatives = select_rows(eligible, args.min_overlap, args.contamination_penalty,
                                        args.min_selection_completeness, args.max_selection_contamination,
                                        1, 'standalone', args.sample)
        by_key = {r['source'] + '/' + r['original_bin_id']: r for r in rows}
        renamed = {r['representative_bin_id']: f'bin_{i}' for i, r in enumerate(reps, 1)}
        rejected = {r['alternative_bin_id']: r for r in alternatives}
        destination = self.out / '3_selection/bins'
        if destination.exists():
            shutil.rmtree(destination)
        destination.mkdir()
        representatives = []
        for rep in reps:
            key = rep['representative_bin_id']
            original = by_key[key]
            new_id = renamed[key]
            target = destination / f'{new_id}.fa'
            shutil.copy2(self.out / original['bin_path'], target)
            represented = [key, *(a['alternative_bin_id'] for a in alternatives if a['representative_bin_id'] == key)]
            representatives.append({**original, 'bin_id': new_id, 'bin_path': str(target.relative_to(self.out)),
                'shared_contig_bp': rep['shared_contig_bp'], 'shared_contig_fraction': rep['shared_contig_fraction'],
                'represented_bins': ';'.join(represented)})
        decisions = []
        for key, row in by_key.items():
            alt = rejected.get(key)
            reason = ('selected' if key in renamed else alt['reason'] if alt else
                      'below_min_selection_length' if row['bp'] < args.min_selection_length else
                      'qc_unassessable' if row.get('qc_status') != 'ok' else 'failed_quality_thresholds')
            matches = alt['matching_representatives'].split(';') if alt else ([key] if key in renamed else [])
            decisions.append(dict(source=row['source'], bin_id=row['original_bin_id'], reason=reason,
                representative=renamed.get(alt['representative_bin_id'] if alt else key, ''),
                overlap=alt['overlap'] if alt else '', shared_bp=alt['shared_bp'] if alt else '',
                matching_representatives=';'.join(renamed[k] for k in matches)))
        columns = ['source', 'bin_id', 'original_filename', 'sample', 'contigs', 'bp', 'n50', 'longest_contig', 'sha256', 'bin_path',
                   'qc_status', 'decision_qc', 'completeness', 'contamination', 'quality_score', 'qc_pass']
        metrics = sorted({k for r in rows for k in r if k.startswith(('eukcc_', 'busco_'))})
        api.write_tsv(self.out / '2_qc/summary.tsv', [*columns, *metrics],
                      ({**r, 'bin_id': r['original_bin_id']} for r in rows))
        api.write_tsv(self.out / '3_selection/representatives.tsv',
            ['bin_id', 'source', 'original_bin_id', *columns[2:], *metrics,
             'shared_contig_bp', 'shared_contig_fraction', 'represented_bins'], representatives)
        api.write_tsv(self.out / '3_selection/selection.tsv',
            ['source', 'bin_id', 'reason', 'representative', 'overlap', 'shared_bp', 'matching_representatives'], decisions)
        self.state['summary'] = dict(input_mode=args.input_mode, retained_bins=len(rows),
            qc_candidates=sum(r['bp'] >= args.min_selection_length for r in rows),
            representatives=len(reps), alternatives=len(alternatives),
            native_binner_status={s: self.state['stages'][s]['status'] for s in args.binner_list},
            imported_sources=list(args.bin_sources))
        self.state['checksums'] = self.checksums(p for name in ('1_raw', '2_qc', '3_selection')
                                               for p in (self.out / name).rglob('*'))
        self.verify(self.state['checksums'])
        failed = any(s['status'] == 'failed' for s in self.state['stages'].values())
        if not failed:
            shutil.rmtree(self.work)
        # Failed native stages remain available for diagnosis, even if other sources completed.
        self.state.update(status='completed', completed_at=now(), partial=failed)
        self.save()
        self.message(f'complete: {len(reps)} representatives, {len(alternatives)} alternatives')


def interrupted(signum, frame):
    raise KeyboardInterrupt(f'received signal {signum}')


def run(api, args):
    api.validate_args(args)
    inputs = api.validate_inputs(args)
    args.cuda_devices = {}
    if args.gpu == 'T' and inputs['eligible_contigs']:
        for name in args.binner_list:
            if name in {'remag', 'vamb'}:
                args.cuda_devices[name] = getattr(api, f'require_{name}_cuda')()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    # Lock the directory inode, not a JSON file that is atomically replaced.
    descriptor = os.open(args.output, os.O_RDONLY)
    old_handler = signal.signal(signal.SIGTERM, interrupted)
    try:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError('another process owns output directory')
        task = Run(api, args, inputs)
        task.lock_fd = descriptor
        if not task.start():
            return
        try:
            task.sources()
            task.qc()
            task.export()
        except BaseException as error:
            task.state.update(status='interrupted' if isinstance(error, (KeyboardInterrupt, SystemExit)) else 'failed',
                              error=str(error) or type(error).__name__)
            task.save()
            task.message(f'{task.state["status"]}: {task.state["error"]}; rerun the same command to resume')
            raise
    finally:
        signal.signal(signal.SIGTERM, old_handler)
        os.close(descriptor)
