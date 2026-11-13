"""Native binner commands; only tool-required input adapters are written."""
import csv
from pathlib import Path
import sys
import tempfile

from depth_schema import identify_depth_schema


def prepare(run):
    api, args = run.api, run.args
    folder = run.work / 'inputs'
    folder.mkdir(parents=True, exist_ok=True)
    lengths = api.fasta_lengths(args.contigs, allow_empty=True)
    fasta = folder / 'contigs.fa'
    # One shared, length-filtered FASTA for all four binners.
    with api.open_text(args.contigs) as src, fasta.open('w') as dst:
        keep = False
        for line in src:
            if line.startswith('>'):
                keep = lengths[line[1:].split()[0]] >= args.min_input_length
            if keep:
                dst.write(line)
    sources = api.parse_depth_inputs(args.depth) if args.depth else ()
    if args.bam:
        jgi = folder / 'jgi.tsv'
        run.execute('coverage', ['jgi_summarize_bam_contig_depths', '--outputDepth', str(jgi),
                    '--percentIdentity', '97', '--minMapQual', '0', '--maxEdgeBases', '75',
                    *map(str, args.bam)], tool='metabat2', cwd=folder)
        sources = (api.DepthInput(jgi),)
    text, _ = api.normalized_depth_inputs(sources, lengths, args.sample,
                  **({} if args.bam else api.depth_column_options(args)))
    lines = list(csv.reader(text.splitlines(), delimiter='\t'))
    schema = identify_depth_schema(lines[0], args.sample, require_jgi=True)
    rows = [row for row in lines[1:] if lengths[row[schema.id_index]] >= args.min_input_length]
    (folder / 'jgi.tsv').unlink(missing_ok=True)
    return fasta, lines[0], schema, rows


def write_table(path, header, rows):
    with path.open('w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        if header is not None:
            writer.writerow(header)
        writer.writerows(rows)


def export_memberships(api, fasta, memberships, bins, minimum):
    """Write only whole bins meeting the size cutoff; keep original contig records."""
    lengths = api.fasta_lengths(fasta)
    assignments, sizes = {}, {}
    for contig, group in memberships:
        if contig not in lengths:
            raise ValueError(f'unknown contig in bin membership: {contig}')
        if contig in assignments:
            raise ValueError(f'duplicate contig in bin membership: {contig}')
        if not group or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-' for c in group) or group in {'.', '..'}:
            raise ValueError(f'invalid bin name: {group!r}')
        assignments[contig] = group
        sizes[group] = sizes.get(group, 0) + lengths[contig]
    selected = {group for group, size in sizes.items() if size >= minimum}
    # One open output at a time; no sequence cache or files for discarded bins.
    output = None
    try:
        with api.open_text(fasta) as source:
            for line in source:
                if line.startswith('>'):
                    if output is not None:
                        output.close()
                    group = assignments.get(line[1:].split()[0])
                    output = (bins / f'{group}.fa').open('a') if group in selected else None
                if output is not None:
                    output.write(line)
    finally:
        if output is not None:
            output.close()


def bin_contigs(run, name, stage, prepared):
    args, api = run.args, run.api
    fasta, header, schema, rows = prepared
    bins = stage / 'bins'
    bins.mkdir()
    means = {row[schema.id_index]: [row[i] for i in schema.coverage_indices] for row in rows}
    extra = args.binner_parameters[name]
    length, threads, seed = map(str, (args.min_input_length, args.threads, args.random_seed))
    if name in {'remag', 'vamb'} and len(rows) < 2:
        return bins
    if name == 'metabat2':
        coverage = stage / 'coverage.tsv'
        variances = all(run.state['inputs']['variance_available'])
        write_table(coverage, header if variances else ['contigName', *schema.coverage_names],
                    rows if variances else ([key, *values] for key, values in means.items()))
        labels = stage / 'labels'
        labels.mkdir()
        prefix = f'{args.sample}.metabat2'
        command = ['metabat2', '-i', str(fasta), '-a', str(coverage), '-o', str(labels / prefix),
                   '-m', length, '--seed', seed, '-t', threads, '--onlyLabel']
        run.execute(name, [*command, *([] if variances else ['--cvExt']), *extra], tool=name, cwd=stage)
        def memberships():
            for path in sorted(labels.glob(prefix + '.*')):
                if not path.name[len(prefix) + 1:].isdigit():
                    continue
                with path.open() as handle:
                    for line in handle:
                        if line.strip():
                            yield line.split()[0], path.name
        export_memberships(api, fasta, memberships(), bins, args.min_retain_size)
    elif name == 'remag':
        coverage_args = []
        for i, _ in enumerate(schema.coverage_names):
            coverage = stage / f'coverage_{i + 1}.tsv'
            write_table(coverage, None, ([key, values[i]] for key, values in means.items()))
            coverage_args += ['-c', str(coverage)]
        command = ['python' if api.environment_value('REMAG_ENV') else sys.executable,
                   str(api.SCRIPTS / 'run_remag.py'), '--gpu', args.gpu, '--', 'remag', str(fasta),
                   *coverage_args, '-o', str(stage), '--skip-bacterial-filter',
                   '--min-bin-size', str(max(1, args.min_retain_size)),
                   '--min-contig-length', length, '--random-seed', seed, '-t', threads, *extra]
        run.execute(name, command, tool=name, cwd=stage)
    elif name == 'concoct':
        fragments, bed = stage / 'contigs_10K.fa', stage / 'contigs_10K.bed'
        run.execute(name, ['cut_up_fasta.py', str(fasta), '-c', '10000', '-o', '0', '--merge_last',
                          '-b', str(bed)], tool=name, cwd=stage, stdout=fragments)
        if sum(bp >= args.min_input_length for bp in api.fasta_lengths(fragments).values()) < 2:
            return bins
        coverage = stage / 'fragment_coverage.tsv'
        if args.bam:
            run.execute(name, ['concoct_coverage_table.py', str(bed), *map(str, args.bam)],
                        tool=name, cwd=stage, stdout=coverage)
        else:
            with bed.open() as handle:
                fragments_rows = [line.rstrip().split('\t') for line in handle if line.strip()]
            # CONCOCT adds fractional pseudocounts in place; keep numeric columns floating point.
            write_table(coverage, ['contig', *schema.coverage_names],
                        ([row[3], *(float(value) for value in means[row[0]])] for row in fragments_rows))
        native = stage / 'concoct'
        native.mkdir()
        threshold = str(args.min_input_length - 1)
        run.execute(name, ['concoct', '--composition_file', str(fragments), '--coverage_file', str(coverage),
                          '--length_threshold', threshold, '--threads', threads, '--seed', seed,
                          '-b', str(native) + '/', *extra], tool=name, cwd=stage)
        merged = stage / 'clustering_merged.csv'
        run.execute(name, ['merge_cutup_clustering.py', str(native / f'clustering_gt{threshold}.csv')],
                    tool=name, cwd=stage, stdout=merged)
        with merged.open() as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames != ['contig_id', 'cluster_id']:
                raise ValueError('invalid CONCOCT membership header')
            export_memberships(api, fasta,
                               ((row['contig_id'], str(int(row['cluster_id']))) for row in reader),
                               bins, args.min_retain_size)
    elif name == 'vamb':
        indices = [i for i in range(len(schema.coverage_names)) if any(float(v[i]) > 0 for v in means.values())]
        if not indices:
            raise ValueError('all coverage columns are zero; cannot run VAMB')
        coverage = stage / 'abundance.tsv'
        write_table(coverage, ['contigname', *(schema.coverage_names[i] for i in indices)],
                    ([key, *(values[i] for i in indices)] for key, values in means.items()))
        run.state['stages'][name]['coverage_columns'] = [schema.coverage_names[i] for i in indices]
        # VAMB insists on a new output directory and a short Unix socket path.
        native = stage / 'vamb'
        command = [sys.executable, str(api.SCRIPTS / 'run_seeded_vamb.py'), '--seed', seed, '--',
                   api.executable('vamb', 'VAMB'), 'bin', 'default', '--outdir', str(native),
                   '--fasta', str(fasta), '--abundance_tsv', str(coverage), '-m', length,
                   '-p', threads, '--seed', seed, '--minfasta', str(max(1, args.min_retain_size)), *(['--cuda'] if args.gpu == 'T' else []), *extra]
        socket_root = Path('/dev/shm') if Path('/dev/shm').is_dir() else stage
        with tempfile.TemporaryDirectory(prefix='eukbin-', dir=socket_root) as socket_dir:
            run.execute(name, command, cwd=stage, env={'TMPDIR': socket_dir, 'TMP': socket_dir, 'TEMP': socket_dir})
        bins = native / 'bins'
    else:
        raise ValueError(f'unknown binner: {name}')
    if not bins.is_dir():
        raise ValueError(f'{name} produced no bins directory')
    return bins
