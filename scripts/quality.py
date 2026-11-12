"""Run one QC tool and retain its small, interpretable evidence files."""
import gzip
import json
import shutil

from qc_parsers import parse_summary, METRIC_COLUMNS
from qc_parsers import quality_row


def keep_file(source, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def eukcc_result(stage, evidence, rc, log_text):
    reports = []
    for path in sorted(stage.rglob('*.csv')):
        row, _ = quality_row(path)
        if row is not None:
            reports.append((path, row))
    if rc == 0:
        if len(reports) != 1:
            raise ValueError(f'EukCC produced {len(reports)} quality reports; expected exactly one')
        path, metrics = reports[0]
        keep_file(path, evidence / 'eukcc.csv')
        result = {'eukcc_status': 'ok', **{f'eukcc_{key}': value for key, value in metrics.items()}}
    else:
        # Only recognize native unassessable states with their specific evidence.
        diagnostic = log_text
        for path in stage.rglob('*'):
            if path.is_file() and path.suffix in {'.log', '.stdout', '.stout'}:
                diagnostic += path.read_text(errors='replace')[-100000:]
        if 'No placement marker genes found' in diagnostic:
            status = 'no_placement_markers'
        elif any(message in diagnostic for message in ('No marker gene set could be found',
                 'Could not identify a single suitable marker gene set')):
            status = 'no_suitable_marker_set'
        elif ('Metaeuk produced empty fasta.' in diagnostic
              and any(p.stat().st_size == 0 for p in stage.rglob('*_metaeuk_cleaned.faa'))
              and any(p.stat().st_size == 0 for p in stage.rglob('*_metaeuk.codon.fas')) and not reports):
            status = 'metaeuk_empty_proteins'
        elif ('Found 1 marker genes, placing them in the tree using epa-ng' in diagnostic
              and any(p.stat().st_size == 0 for p in stage.rglob('hmmalign.stdout'))
              and not list(stage.rglob('query.fasta')) and not list(stage.rglob('epa-ng.stout'))):
            status = 'placement_alignment_failed'
        else:
            raise ValueError(f'EukCC failed (exit {rc}); see 0_log/run.log')
        result = {'eukcc_status': status}
    # EukCC's state contains the selected marker set and the hits used for QC.
    for path in sorted(stage.glob('*.json.gz')):
        with gzip.open(path, 'rt') as handle:
            state = json.load(handle)
        if isinstance(state, str):
            state = json.loads(state)
        if isinstance(state, dict) and 'marker_set' in state:
            detail = {key: state[key] for key in ('marker_set', 'scmg_data', 'quality') if key in state}
            evidence.mkdir(parents=True, exist_ok=True)
            (evidence / 'markers.json').write_text(json.dumps(detail, indent=2) + '\n')
            break
    for pattern in ('scmg_marker_table.csv.gz', 'found_markers_*_GA.csv', 'found_markers_*.tbl'):
        for path in sorted(stage.rglob(pattern)):
            keep_file(path, evidence / path.name)
    return result


def busco_result(stage, evidence, rc):
    if rc:
        raise ValueError(f'BUSCO failed (exit {rc}); see 0_log/run.log')
    reports = [(path, parse_summary(path)) for path in sorted(stage.rglob('short_summary*.txt'))]
    reports = [(p, row) for p, row in reports if row is not None]
    if not reports or len({tuple(row[k] for k in METRIC_COLUMNS) for _, row in reports}) != 1:
        raise ValueError('BUSCO summary is missing or ambiguous')
    path, metrics = min(reports, key=lambda item: (not item[0].name.startswith('short_summary.specific.'),
                                                 len(item[0].parts), str(item[0])))
    keep_file(path, evidence / 'short_summary.txt')
    tables = sorted(stage.rglob('full_table.tsv'))
    if len(tables) > 1:
        raise ValueError('BUSCO produced multiple full tables for the specified lineage')
    if tables:
        with tables[0].open('rb') as src, gzip.open(evidence / 'full_table.tsv.gz', 'wb') as dst:
            shutil.copyfileobj(src, dst)
    return {'busco_status': 'ok', **{f'busco_{k}': v for k, v in metrics.items()}}


def assess(run, row, tool, stage):
    args = run.args
    output = stage / 'native'
    output.mkdir()
    evidence = run.out / '2_qc/evidence' / row['source'] / row['original_bin_id'] / tool
    if evidence.exists():
        shutil.rmtree(evidence)
    source = run.out / row['bin_path']
    reference = run.state['inputs']['qc_references'][f'{tool}_db']
    if tool == 'eukcc':
        command = ['eukcc', 'single', '--db', reference, '--out', str(output),
                   '--threads', str(args.threads), str(source)]
    else:
        # No --tar: archive creation adds I/O; retain only summary and marker table.
        command = ['busco', '--in', str(source), '--out', 'result', '--out_path', str(output),
                   '--mode', 'genome', '--lineage_dataset', reference, '--cpu', str(args.threads), '--offline']
    result = run.execute(f'{tool}:{row["source"]}/{row["original_bin_id"]}', command,
                         tool=tool, cwd=stage, check=False)
    if tool == 'eukcc':
        metrics = eukcc_result(output, evidence, result.returncode, run.last_output)
    else:
        metrics = busco_result(output, evidence, result.returncode)
    if evidence.is_dir():
        metrics[f'{tool}_evidence'] = str(evidence.relative_to(run.out))
    return metrics
