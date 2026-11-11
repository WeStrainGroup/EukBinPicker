# EukBinMiner

Current version: **0.3.0**.

EukBinMiner is a standalone command-line workflow for mining,
quality-controlling, and selecting non-redundant eukaryotic bins from one
sample's contigs plus one or more mapped-sample depth/coverage dimensions. It is distributed as a standalone Python/Bash CLI.

## Workflow

```text
contigs + one or more depth/coverage dimensions
  -> one or more binners (default: all)
  -> seqkit bin statistics
  -> EukCC and/or BUSCO
  -> within-sample, cross-binner overlap comparison
  -> best representative
  -> non-redundant representative-bin FASTA files
```

The default and `--binners all` run four binners: `remag`, `concoct`,
`metabat2`, and `vamb`. They do not require taxonomy or Kun-peng.
A failed binner does not stop the workflow when another requested binner
produced a usable branch. Bins emitted by the same binner are never compared
to one another.

Two bins from different binners are considered overlapping when

```text
shared_contig_bp / min(bin_a_bp, bin_b_bp) >= overlap_threshold
```

The default threshold is `0.50`. EukCC is the only QC tool run by default.
Bins at least 500,000 bp with successful QC and valid completeness/contamination
values enter selection. Default quality limits are completeness at least 0%
and contamination at most 100%; these do not imply high-quality MAGs.
Eligible representatives are ordered by quality score
(`completeness - contamination_penalty * contamination`), completeness,
contamination, contig count, N50, and bin ID.

EukCC supplies its native contamination percentage. When BUSCO controls the
decision, EukBinMiner uses BUSCO duplicated percentage as a contamination
proxy. These measurements are useful for the same filtering role but are not
biologically identical.

Reproducible binning uses one public seed, `--random-seed` (default `6311`).
The parsed value is supplied to ReMAG, CONCOCT, MetaBAT2, and VAMB
and is recorded in run metadata and every binner branch fingerprint.

## Calculation settings

By default, the four binners retain their native algorithm settings, matching
standalone use. The agreed exceptions are the inclusive 1,500-bp input minimum,
seed 6311, and ReMAG skipping its own classifier so that ECM defines the input.
CONCOCT uses the standard 10-kb non-overlapping fragments with the last fragment
merged; BAM mode computes native fragment coverage. VAMB exports nonempty FASTA
files while retaining native training and bin splitting. The 500-kb selection
minimum does not replace each tool's native output threshold. All four native
parameter strings are empty by default.

Version 0.3.0 changes the public interface and EBM selection defaults described
below. It does not rewrite frozen benchmark inputs, code, device choices or
results. CPU/GPU choice must be matched explicitly between independent and EBM
runs when assessing candidate equivalence.

| Stage | Applied behavior |
| --- | --- |
| Input | Retain contigs at least 1,500 bp by default. Every retained contig must have coverage in each supplied mapping dimension. |
| ReMAG | Native metagenomics defaults, bacterial filtering skipped, seed 6311; native 500-kb output minimum. The supplied contigs define the candidate input set. |
| CONCOCT | Cut into 10-kb fragments with the last fragment merged. BAM input uses native fragment coverage; table input approximates fragments by the parent mean. Use the native 400-cluster default without adaptive overrides. Fewer than two eligible fragments produce zero bins. |
| CONCOCT coverage | Native normalization and total-coverage feature for every coverage dimensionality. |
| MetaBAT2 | Native 200-kb output minimum. Preserves observed JGI mean/variance pairs. Mean-only input uses native `--cvExt`; missing variance is never filled with zero. |
| VAMB | Default training schedule and native bin splitting, seed applied before model import. All-zero coverage columns are removed and recorded. |
| Raw bins | Additional retention minimum defaults to 0 (no extra size filter) and cannot lower native output thresholds. VAMB exports nonempty FASTA at 1 bp because native default emits assignments only. Bins at least 500,000 bp enter QC by default. |
| EukCC | Native completeness and contamination estimates. Verified no-marker or empty-protein outcomes retain an unassessed status. |
| BUSCO | Genome mode with the supplied lineage, offline operation and archived intermediate files. Duplicated BUSCO percentage supplies the decision proxy when selected. |
| Selection | Valid QC required; completeness ≥0%, contamination/proxy ≤100%, score = completeness − 5 × contamination/proxy. Compare only different binners within the current assembly. Unassessed bins stay recorded but cannot be ranked. |

Selection proceeds in descending quality order and compares each candidate
against the representatives already retained. A candidate is suppressed only
when it directly overlaps a retained representative. Each suppressed candidate
records its directly overlapping representative and any ambiguity.
The contig-ID and normalized-sequence pair defines a shared source contig.

## Installation

Clone the source and put `bin/` on `PATH`:

```bash
git clone https://github.com/WeStrainGroup/EukBinMiner.git
export PATH="$PWD/EukBinMiner/bin:${PATH}"
eukbinminer --version
```

EukBinMiner discovers tools from the current `PATH`. The combined environment
specification is `envs/eukbinminer.yml`. The earlier EukCC 2.1.3 / BUSCO 6
dependency solve found incompatible MetaEuk requirements, so this combined
configuration still requires validation. On systems that keep tools in separate named Conda
environments, set `EUKBINMINER_REMAG_ENV`, `EUKBINMINER_CONCOCT_ENV`,
`EUKBINMINER_METABAT2_ENV`, `EUKBINMINER_EUKCC_ENV`, and
`EUKBINMINER_BUSCO_ENV`. The corresponding backend is then run with
`conda run --no-capture-output -n ENV`; an absolute environment path uses
`conda run -p PREFIX`. VAMB and seqkit remain resolved from
`PATH`. `EUKBINMINER_CONDA_EXE` may be an absolute Conda executable. Legacy
`EUKMAGMINER_*_ENV` names are also accepted. The normal run checks its inputs and reference databases before binning,
and reports concrete input, dependency or device errors when encountered.

The repository contains a local-development EukBinMiner recipe.
The source uses CC BY-NC-SA 4.0. [GitHub](https://github.com/WeStrainGroup/EukBinMiner). The Conda recipe and combined environment are development
artifacts. A combined EukCC/BUSCO environment still needs a validated packaging solution.
The per-tool specifications are in `envs/`.

The old `eukmagminer` executable remains as a deprecated compatibility alias
for existing scripts. New commands, documentation, summaries, and package
metadata use `EukBinMiner`/`eukbinminer`.

## CPU and GPU execution

`--gpu F` is the default and runs every tool on CPU, even inside a GPU
allocation. EBM hides CUDA devices from child processes and disables ReMAG's
native CUDA/MPS auto-selection in its launcher before importing ReMAG.

`--gpu T` uses native CUDA for every selected ReMAG and VAMB branch. VAMB
receives its native `--cuda` flag for training and clustering. CONCOCT,
MetaBAT2 and the current QC workflow remain on CPU. ReMAG and VAMB execute
sequentially, so one allocated GPU can serve both. At least one of these two
binners must be selected when requesting GPU mode.

The flag does not allocate GPUs or install CUDA. A CUDA-enabled PyTorch build
and a visible NVIDIA GPU are required in each selected tool's environment.
EBM checks both applicable environments before creating the output directory
and fails if CUDA is unavailable.
The scheduler's `CUDA_VISIBLE_DEVICES` restriction is preserved in GPU mode.
Requested mode and CUDA device/build details are recorded in `config/run.json`;
CPU/GPU outputs cannot be silently reused across different modes.

For a separate ReMAG environment, set `EUKBINMINER_REMAG_ENV`; VAMB uses its
console-script Python environment. Neither bundled environment specification
guarantees a CUDA-enabled PyTorch build.

```bash
export EUKBINMINER_REMAG_ENV=/path/to/validated-remag-gpu
eukbinminer --sample sample1 --contigs euk_contigs.fna \
  --coverage depth.tsv --output sample1_bins \
  --eukcc-db /path/to/eukcc_db --gpu T
```

GPU support does not imply a fixed end-to-end speedup. Input processing and QC
still take CPU time, and CPU/GPU bins need not be identical even with the same
seed. Software routing is covered by mock regression tests; benchmark timing
requires representative runs on real hardware.

## Reference databases

Pass the unpacked database directories directly:

```text
--eukcc-db /path/to/eukcc_db
--busco-db /path/to/fungi_odb12
```

The equivalent environment variables are `EUKBINMINER_EUKCC_DB` and
`EUKBINMINER_BUSCO_DB`. `--busco-db` must point to the unpacked BUSCO lineage
directory containing `dataset.cfg` and `hmms/`.

## Usage

Use `eukbinminer --help` to list parameters or `eukbinminer --version` to
print the version. There are no subcommands and no `-h` alias. Run EBM directly;
invalid arguments, inputs, databases and unavailable requested CUDA produce
specific error messages and a nonzero exit code.

One invocation handles one assembly. Organize runs with `--output`; there is
no `--study` argument. The shared backend uses a fixed internal `standalone`
group in intermediate manifests and paths.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--coverage FILE [FILE ...]` | Required | Coverage tables or indexed BAMs; repeatable, with no mixing of BAMs and tables. |
| `--cov-contig-column` | Auto | Contig ID column in a custom table. |
| `--cov-mean-column` | Auto | Comma-separated sample mean columns in order. |
| `--cov-variance-column` | Optional | Observed variance columns paired with the mean columns in the same order. |
| `--gpu` | `F` | `T` uses CUDA for selected ReMAG/VAMB; `F` runs all tools on CPU. |
| `--remag-parameters` | Empty | Native ReMAG algorithm options, supplied as one quoted string. |
| `--concoct-parameters` | Empty | Native CONCOCT algorithm options, supplied as one quoted string. |
| `--metabat2-parameters` | Empty | Native MetaBAT2 algorithm options, supplied as one quoted string. |
| `--vamb-parameters` | Empty | Native VAMB `bin default` options, supplied as one quoted string. |
| `--qc-tools` | `eukcc` | Run EukCC, BUSCO or both (`all`). |
| `--qc-decision` | `eukcc` | Executed QC tool used for eligibility and ranking. |
| `--min-input-length` | `1500` | Minimum contig length entering binning, in bp. |
| `--min-retain-length` | `0` | Additional raw-bin size filter, in bp; 0 keeps all nonempty native outputs. |
| `--min-selection-length` | `500000` | Minimum retained-bin length entering QC and selection, in bp. |
| `--min-selection-completeness` | `0` | Minimum completeness for selection, in percent. |
| `--max-selection-contamination` | `100` | Maximum contamination/proxy for selection, in percent. |
| `--contamination-penalty` | `5` | Penalty in completeness minus penalty times contamination. |
| `--min-overlap` | `0.50` | Minimum shared-contig fraction that suppresses a lower-ranked cross-binner candidate. |

All thresholds are inclusive. `selection` denotes the representative-selection
stage. Selection length must be at least retention length. Quality limits do
not admit missing, failed or non-finite QC results. No aliases for the old CLI
parameter names are provided.

```bash
eukbinminer \
  --sample SAMPLE_001 \
  --contigs contigs.fa \
  --coverage depth.tsv \
  --output eukbinminer-output \
  --binners all \
  --qc-tools eukcc \
  --qc-decision eukcc \
  --min-input-length 1500 \
  --min-retain-length 0 \
  --min-selection-length 500000 \
  --min-selection-completeness 0 \
  --max-selection-contamination 100 \
  --contamination-penalty 5 \
  --min-overlap 0.50 \
  --random-seed 6311 \
  --eukcc-db /path/to/eukcc_db
```

For differential coverage from reads of several samples mapped back to the
same contig assembly, repeat `--coverage`. `NAME=TSV` is recommended for generic
two-column files because it records the intended sample dimension explicitly:

```bash
eukbinminer \
  --sample ASSEMBLY_SAMPLE \
  --contigs contigs.fa \
  --coverage sample_A=sample_A.depth.tsv \
  --coverage sample_B=sample_B.depth.tsv \
  --coverage sample_C=sample_C.depth.tsv \
  --output eukbinminer-output \
  --eukcc-db /path/to/eukcc_db
```

The equivalent single-matrix form is:

```text
contig_id  sample_A  sample_B  sample_C
contig_1   2.5       0.0       8.1
contig_2   4.0       1.2       3.3
```

and is supplied once as `--coverage coverage.tsv`.

`--qc-tools` accepts `eukcc`, `busco`, `eukcc,busco`, or `all`; it defaults to
`eukcc`. Selecting `--qc-decision busco` also requires `--qc-tools busco` or `all`.
`--qc-decision` selects which executed QC tool controls final eligibility and
representative ordering. Contigs enter binning when their length is greater
than or equal to `--min-input-length` (default 1500 bp).
`--random-seed` accepts a positive integer and defaults to `6311`.
If that inclusive filter leaves no eligible contigs, the run writes a
checksum-bearing empty result without running binners or QC; it records all requested binners as
`not_run_empty_input` rather than claiming that they executed.

Each depth input may be a JGI table, a generic two-column table, or a generic
matrix containing an ID column, optional `contigLen`, and one or more numeric
coverage columns. `--coverage` may be repeated for independent mapping outputs;
every supplied file must contain every FASTA contig. A `NAME=TSV` label is
valid for a single-coverage file. Without an explicit label, repeated generic
files use a meaningful coverage header when present and otherwise infer the
dimension name from the filename. Input order and column order are preserved.

Observed variance is retained exactly. When no variance is supplied, the internal
canonical table records `NA`; MetaBAT2 receives only IDs and means through its
native `--cvExt` interface. This delegates the mean-only model to MetaBAT2 and
does not create observed variances. Mixed observed/missing variance across
samples is rejected when MetaBAT2 is selected: supply all variance columns or
omit MetaBAT2. Generic `totalAvgDepth` is the sum of coverage dimensions; a single
JGI input retains its supplied total. VAMB, ReMAG and CONCOCT use the sample
means. VAMB records and drops all-zero dimensions. ReMAG receives one headerless
TSV per dimension through repeated `-c` options.

For custom tables, select the ID, mean columns and optional paired variance
columns explicitly; unrelated columns are ignored:

```bash
eukbinminer --sample assembly --contigs contigs.fa --coverage custom.tsv \
  --cov-contig-column sequence_id --cov-mean-column A_mean,B_mean \
  --cov-variance-column A_var,B_var --output output \
  --qc-tools eukcc --eukcc-db /path/to/eukcc_db
```

Omit `--cov-variance-column` if unavailable. The lists pair by position. FASTA
provides contig lengths; any supplied `contigLen` is checked. Standard JGI tables
need no column selectors. The same selectors apply to every repeated table.

Alternatively, supply **one coordinate-sorted, indexed BAM per sample**:

```bash
eukbinminer --sample assembly --contigs eukaryotic_contigs.fa \
  --coverage mapping/sample_A.bam mapping/sample_B.bam \
  --output output --qc-tools eukcc --eukcc-db /path/to/eukcc_db
```

`--coverage` detects BAM by its file signature or `.bam` extension; other inputs
are parsed as tables (plain or gzip-compressed TSV). All inputs in one run must
be BAMs or all must be tables. Column selectors and `NAME=` aliases apply only
to tables. BAM mode requires `samtools` and
`jgi_summarize_bam_contig_depths` on PATH (or in `EUKBINMINER_METABAT2_ENV`), even
when MetaBAT2 is not selected. Every BAM needs a nonempty `.bai` or `.csi` index,
a unique basename, coordinate-sort declaration and the same ordered reference
ID/length dictionary. BAMs may map against the complete assembly before ECM:
all input FASTA IDs and lengths must match that reference. EBM checks sequence
MD5 against BAM `M5` tags where available. Without those tags, IDs and lengths
cannot establish sequence identity; the user must ensure the FASTA came from
the mapped assembly. BAMs and indexes must remain immutable during a run.

EBM reuses these BAMs without remapping. Native JGI computes contig mean and
variance with 97% identity, minimum mapping quality 0 and native edge exclusion
(maximum 75 bases). CONCOCT cuts the input internally and obtains fragment means
using its native `concoct_coverage_table.py` (`samtools bedcov` divided by interval
length). These tools have different native read/position filtering semantics;
fragment means are not expected to reproduce JGI means exactly.

With table-only input, CONCOCT fragments inherit their parent's mean as an
explicit approximation. Means and variances cannot recover positional coverage.
CONCOCT normalization, total coverage and cluster count use native defaults
for both single and multiple coverage dimensions.

The normal run checks BAM inputs and indexes before calculating coverage.
Resume records resolved BAM paths, sizes, modification times and
header hashes plus index hashes, rather than hashing every BAM. Generated JGI
coverage has its own checksum. This detects routine input changes but is not a
cryptographic checksum of the BAM alignment payload.

FASTA IDs, depth IDs, and coverage names must be unique. Rows must have exactly
the declared number of fields. Every supplied coverage, observed variance, and
`totalAvgDepth` value must be finite and non-negative, and any supplied
`contigLen` must match the FASTA. Duplicate depth IDs are fatal even when rows
are identical. Production adapters may canonicalize read-only legacy inputs,
but the public CLI always validates the canonical copy strictly. The resolved
coverage names/count, source files, total-depth policy, and checksums are
recorded in validation and run metadata.

## Native binner parameters

Use each tool's original option names inside its `--TOOL-parameters` string.
An empty string retains the current defaults. For example:

```bash
eukbinminer --sample sample1 --contigs euk_contigs.fa \
  --coverage coverage.tsv --output sample1_bins \
  --eukcc-db /path/to/eukcc_db \
  --remag-parameters="--epochs 200 --batch-size 1024" \
  --concoct-parameters="--clusters 200" \
  --metabat2-parameters="--maxP 90" \
  --vamb-parameters="-e 200 -n 128 128"
```

Use `=` and quote the whole string, including a single flag. Use separate short
options and values (`-e 200`); combined short flags are not supported. Quoting
inside the string preserves multiword values. Parameters are passed as literal
arguments, without shell expansion or command execution.

Input paths, output paths, threads, seeds, CUDA, common input-length settings and
coverage-format switches are managed by EBM. Conflicting options (including
their aliases) are rejected; set the corresponding EBM option instead. Native
help/version and workflow-changing options such as ReMAG `--filter-only` are
also rejected. Parameters for an unselected binner produce an error. Other
option names and values are validated by the native tool, and a native error
fails that branch with its original error recorded in the logs. Independently
successful branches can still complete under EBM's existing failure policy.

Resolved argument lists are stored in `config/run.json` under `binner_parameters`.
Each executed binner branch also saves `command.txt` with its final command;
these record the actual execution paths, which may be staging paths. Changing
parameters requires a new output directory, preventing stale results from
being silently reused. Native overrides also supersede the corresponding
algorithm defaults described above.

## Important output

- `mags/`: copied non-redundant representative-bin FASTA files and
  `manifest.tsv` (the path is retained for output-schema compatibility)
- `selection/representatives.tsv`: selected representatives and their metrics,
  including `shared_contig_bp` and `shared_contig_fraction` (0 to 1). These count
  original contigs also present in any other final retained bin from the same
  assembly, including the same binner. Each shared contig is counted once by
  length; discarded alternatives do not contribute. Identity uses original ID
  plus normalized sequence, not sequence similarity. This diagnostic does not
  change selection, QC scores or bin contents; it is not a contamination estimate.
- `selection/alternatives.tsv`: suppressed cross-binner alternatives
- `qc/`: seqkit, EukCC, BUSCO, and decision tables
- `results/`: retained successful raw binner branches
- `summary.json`: run-level machine-readable summary
- `config/run.json`: resolved inputs, parameters, checksums, and random seed
- `.complete.json`: checksum-bearing completion marker

EukBinMiner never deletes raw bins. Reusing an output directory is allowed
only when its validated completion state is reusable; stale or partial branch
outputs are rejected rather than overwritten.

Public QC result fields are `completeness`, `contamination`,
`contamination_source`, `qc_pass`, and `quality_score`.
`contamination_source` is `eukcc_contamination` or `busco_duplicated`.

## Development

Run the local tests with:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -v
```

The suite exercises mock external tools, real file conversion, QC parsing,
overlap selection, failure isolation, rerun validation and installed layouts.
Historical real-data runs use frozen versions. The changes in this source
release are validated by software regression tests.

## Review changes

- Version 0.3.0: unified `--coverage`, consistent input/retention/selection
  names, EukCC-only default QC, quality limits of 0%/100%, and `--gpu T/F`
  controlling native ReMAG/VAMB CUDA. Four `--TOOL-parameters` options forward
  native algorithm settings. Previous CLI parameter names, auxiliary subcommands
  and the `-h` alias are removed.

- Version 0.2.0: add BAM or table coverage input, explicit custom column selection,
  native CONCOCT fragment coverage from BAM and native MetaBAT2 mean-only input.
  Missing variance is no longer fabricated. Current CAMI benchmark uses frozen
  pre-0.2.0 code and is not retroactively updated: its standalone CONCOCT uses
  BAM fragment coverage, whereas EBM inherits contig means. That benchmark
  compares full implementations, not only the final bin-selection strategy.

- Version 0.1.4: add --remag-gpu to require native ReMAG CUDA, with environment
  checks and recorded device/build information.
- Version 0.1.3: report final-bin shared contig length and fraction without
  changing selection or recovering additional bins.

- Version 0.1.2: default and `all` select only ReMAG, CONCOCT, MetaBAT2 and VAMB.
  TaxVAMB, taxonomy options and classifier dependencies have been removed.
  QC and overlap settings are unchanged.

- Shared overlap uses the original contig ID plus normalized sequence identity.
  Distinct source IDs with identical sequences retain their individual lengths.
- Reuse checks the input checksums, configuration, CLI/backend code and final
  output checksums. An output-directory lock prevents concurrent writers.
  Reference database paths must point to immutable datasets.
- Invalid QC metrics and inconsistent BUSCO summaries are rejected.
  EukCC bins with verified empty MetaEuk proteins remain explicitly unassessed.
- VAMB drops coverage columns that are zero across every selected contig and
  records their names. A branch with no informative coverage fails.
- A failed binner allows independently successful binners to finish.
  Shell-stage failures propagate before a branch can be published.
- Installed distributions include backend environment specifications.

## License

CC BY-NC-SA 4.0. See [LICENSE](LICENSE).
External tools and reference datasets retain their respective licenses.
