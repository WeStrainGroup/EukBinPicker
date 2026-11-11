# EukBinMiner

**Recover eukaryotic genome bins with complementary binning methods.**

Version **0.40** · [GitHub](https://github.com/WeStrainGroup/EukBinMiner)

## Introduction

Genome-resolved metagenomics provides access to uncultivated eukaryotes, as illustrated by [Tara Oceans](https://pmc.ncbi.nlm.nih.gov/articles/PMC9903769/). However, fragmented assemblies, uneven coverage and diverse genome architectures make recovery challenging. Different binning methods can recover different genomes from the same sample. Following the complementary-binner rationale of [DAS Tool](https://github.com/cmks/DAS_Tool), EukBinMiner (EBM) evaluates multiple candidate sets to select representative eukaryotic bins.

EBM combines [REMAG](https://github.com/danielzmbp/remag), [CONCOCT](https://github.com/BinPro/CONCOCT), [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/) and [VAMB](https://github.com/RasmussenLab/vamb). These methods use contig sequences and coverage without requiring prokaryotic marker genes or taxonomy for binning. They offer complementary approaches: contrastive learning, Gaussian mixture models, graph clustering and variational autoencoders. REMAG specifically targets eukaryotes; the other methods are general metagenomic binners applied here to eukaryotic contigs.

A typical workflow is to assemble quality-controlled reads, map reads to the complete assembly, and identify eukaryotic contigs with [EukContigMiner](https://github.com/WeStrainGroup/EukContigMiner). Supply the resulting eukaryotic FASTA and the full assembly's coverage table or BAM files to EBM. Coverage may come from the source sample alone or from multiple samples mapped to that assembly. **Each run processes one assembly.**

```text
Eukaryotic contigs + full-assembly coverage tables or BAMs
                          ↓
             REMAG · CONCOCT · MetaBAT2 · VAMB
                          ↓
                 EukCC or BUSCO assessment
                          ↓
          Rank overlapping candidates within the assembly
                          ↓
   Select the highest-ranked cross-binner representative
               for each putative genome
```

## Installation

On **Linux x86-64**, with Conda and curl installed, one command installs and configures EBM, all four binners, EukCC, BUSCO and the BAM-processing tools:

```bash
curl -fsSL https://raw.githubusercontent.com/WeStrainGroup/EukBinMiner/v0.40/install.sh | bash -s -- --prefix "$HOME/conda-envs/eukbinminer"
```

Then activate the environment:

```bash
conda activate "$HOME/conda-envs/eukbinminer"
eukbinminer --version
```

The installer manages conflicting dependencies in private environments automatically; no manual tool-path configuration is needed. Use a new installation directory. An interrupted installation can be resumed with the same command.

Tool versions: REMAG 0.4.4, CONCOCT 1.1.0, MetaBAT2 2.18, VAMB 5.0.4, EukCC 2.1.3 and BUSCO 6.

For GPU support, add `--gpu T` to the installation command. This installs CUDA 12.4-enabled PyTorch for REMAG and VAMB; a compatible NVIDIA driver and an allocated GPU are required when running with `--gpu T`. CPU installation and CPU execution are the defaults.

Reference databases are downloaded separately: provide the unpacked [EukCC database](https://github.com/Finn-Lab/EukCC) with `--eukcc-db`, or a [BUSCO lineage](https://busco.ezlab.org/) with `--busco-db`. These paths can also be set through `EUKBINMINER_EUKCC_DB` and `EUKBINMINER_BUSCO_DB`.

## Quick start

```bash
eukbinminer \
  --sample sample1 \
  --contigs euk_contigs.fa \
  --coverage assembly_depth.tsv \
  --eukcc-db /path/to/eukcc_db \
  --output sample1_ebm
```

By default, EBM runs all four binners on CPU with 4 threads, retains contigs ≥1,500 bp, and assesses bins ≥500,000 bp with EukCC. Selected FASTA files are written to `sample1_ebm/bins/`.

## Coverage input

`--coverage` accepts tables or BAM files, including multiple files or repeated `--coverage` arguments. All dimensions must refer to the same assembly; do not mix BAMs and tables in one run.

**Full-assembly coverage is supported.** EBM selects the rows matching the input eukaryotic FASTA, preserving its contig order and coverage values. Other contigs are excluded and counted in run metadata. Every supplied FASTA contig must have coverage; missing values are not replaced with zero.

### Coverage tables

[Standard JGI depth tables](https://bitbucket.org/berkeleylab/metabat/) are recognized automatically. Generic mean-coverage matrices and custom schemas are also supported. Tables must be tab-separated, with a header; plain and gzip-compressed files are accepted.

For example, `custom.tsv` may contain:

```text
sequence_id	A_mean	B_mean	A_var	B_var
contig_1	12.4	8.1	2.0	1.4
contig_2	0.0	3.2	0.0	0.6
contig_3	5.7	2.1	0.8	0.3
```

```bash
eukbinminer \
  --sample sample1 \
  --contigs euk_contigs.fa \
  --coverage custom.tsv \
  --cov-contig-column sequence_id \
  --cov-mean-column A_mean,B_mean \
  --cov-variance-column A_var,B_var \
  --eukcc-db /path/to/eukcc_db \
  --output sample1_ebm
```

If the FASTA contains only `contig_1` and `contig_3`, only those rows enter binning. Contig IDs must be unique, values finite and non-negative, and any supplied contig lengths must match the FASTA. Input tables are validated before subsetting.

Variance columns are optional and pair with mean columns in the same order. MetaBAT2 uses observed variance when available and its native `--cvExt` mode for mean-only input; variance must be present for all dimensions or none. REMAG, CONCOCT and VAMB use mean coverage.

### BAM files

```bash
eukbinminer \
  --sample sample1 \
  --contigs euk_contigs.fa \
  --coverage sample_A.bam sample_B.bam \
  --eukcc-db /path/to/eukcc_db \
  --output sample1_ebm
```

Provide one coordinate-sorted BAM per mapping sample, each with a `.bai` or `.csi` index, unique basename, and the same ordered reference IDs and lengths. The eukaryotic FASTA must be a matching subset of the mapped assembly. EBM reuses these alignments without remapping reads.

The installer includes SAMtools and `jgi_summarize_bam_contig_depths`. JGI produces contig means and variances; CONCOCT computes native fragment coverage from the BAMs. With table input, CONCOCT approximates fragment coverage using each parent contig's mean.

## Parameters

### Input and output

| Parameter | Default | Description |
| --- | --- | --- |
| `--sample` | Required | Assembly/sample label: letters, digits, `.`, `_` or `-`. |
| `--contigs` | Required | Eukaryotic contig FASTA with unique IDs. |
| `--coverage` | Required | One or more coverage tables or indexed BAMs; repeatable. |
| `--output` | Required | Run output directory. |
| `--cov-contig-column` | Auto | Contig ID column in custom tables. |
| `--cov-mean-column` | Auto | Comma-separated mean columns, in sample order. |
| `--cov-variance-column` | None | Observed variance columns paired with mean columns. |

Column selectors apply only to tables, and the same selectors apply to every supplied table. Single-coverage tables can be named explicitly with `--coverage sample_A=sample_A.tsv`.

### Binning and resources

| Parameter | Default | Description |
| --- | --- | --- |
| `--binners` | `all` | All four tools or a comma-separated subset of `remag,concoct,metabat2,vamb`. |
| `--threads` | `4` | Threads supplied to each tool. |
| `--random-seed` | `6311` | Shared positive seed for binning. |
| `--gpu` | `F` | `T`: CUDA for selected REMAG/VAMB; `F`: CPU. |
| `--remag-parameters` | Empty | Quoted native REMAG options. |
| `--concoct-parameters` | Empty | Quoted native CONCOCT options. |
| `--metabat2-parameters` | Empty | Quoted native MetaBAT2 options. |
| `--vamb-parameters` | Empty | Quoted native VAMB options. |

REMAG and VAMB run sequentially, so one allocated GPU can serve both. CONCOCT, MetaBAT2 and QC remain on CPU. `--gpu T` checks CUDA availability in the selected GPU-capable tools; it does not request hardware from a scheduler.

### Quality assessment

| Parameter | Default | Description |
| --- | --- | --- |
| `--qc-tools` | `eukcc` | `eukcc`, `busco`, `eukcc,busco` or `all`. |
| `--qc-decision` | `eukcc` | Executed QC tool used for filtering and ranking. |
| `--eukcc-db` | Environment | Unpacked EukCC database; required when EukCC is selected. |
| `--busco-db` | Environment | Unpacked BUSCO lineage containing `dataset.cfg` and `hmms/`. |

To rank with BUSCO, include it in `--qc-tools` and set `--qc-decision busco`. BUSCO duplicated percentage is used as a contamination proxy; it is not equivalent to EukCC's contamination estimate.

### Filtering and selection

| Parameter | Default | Description |
| --- | --- | --- |
| `--min-input-length` | `1500` | Minimum input contig length, in bp. |
| `--min-retain-length` | `0` | Additional bin-retention minimum, in bp; 0 disables it. |
| `--min-selection-length` | `500000` | Minimum bin length entering QC and selection, in bp. |
| `--min-selection-completeness` | `0` | Minimum completeness, in %. |
| `--max-selection-contamination` | `100` | Maximum contamination or BUSCO proxy, in %. |
| `--contamination-penalty` | `5` | Contamination penalty in the ranking score. |
| `--min-overlap` | `0.50` | Cross-binner overlap threshold. |

All thresholds are inclusive. Selection length must be at least retention length. Valid QC estimates are required even with the permissive 0%/100% defaults; selected bins are not necessarily high-quality MAGs.

### Information

| Parameter | Description |
| --- | --- |
| `--help` | Show all command-line options. |
| `--version` | Print the installed version. |

### Customizing a binner

Use the tool's native option names inside one quoted string:

```bash
eukbinminer \
  --sample sample1 \
  --contigs euk_contigs.fa \
  --coverage assembly_depth.tsv \
  --eukcc-db /path/to/eukcc_db \
  --output sample1_custom \
  --concoct-parameters "--clusters 200" \
  --metabat2-parameters "--maxP 90"
```

No `=` is required. Input/output paths, threads, seeds, devices and shared input-length settings are managed by EBM; use their public EBM options instead. Native options and actual commands are recorded with the run.

## How selection works

Candidates passing the length and QC thresholds are ranked by:

```text
score = completeness − contamination_penalty × contamination
```

EBM visits candidates in descending quality order. A candidate is excluded if it overlaps an already selected bin from a different binner:

```text
overlap(A, B) = shared_contig_bp / min(length(A), length(B))
```

The default cutoff is 0.50. This selects the highest-ranked representative among competing cross-binner candidates within the assembly. Shared contigs must match in original ID and normalized sequence; overlap is a proxy for competing genome reconstructions, not a taxonomic identity test. Bins from the same binner are not compared. EBM preserves complete candidate bins without merging, splitting or refinement.

Native algorithm defaults are retained, with the shared input threshold and seed. REMAG's classifier is skipped. CONCOCT uses standard 10-kb fragments, native normalization and 400 clusters. MetaBAT2 and REMAG retain their native 200-kb and 500-kb output minima. VAMB retains native training and bin splitting, with nonempty FASTA export enabled. The retention setting does not lower these native minima.

## Outputs

| Path | Contents |
| --- | --- |
| `bins/` | Selected bin FASTA files and `manifest.tsv`. |
| `selection/representatives.tsv` | Selected bins, QC metrics and shared-contig statistics. |
| `selection/alternatives.tsv` | Excluded alternatives and corresponding representatives. |
| `qc/` | Bin statistics and QC results. |
| `results/` | Retained binner results. |
| `config/run.json` | Resolved inputs, parameters and provenance. |
| `summary.json` | Run summary and binner outcomes. |

`shared_contig_fraction` is the fraction of a selected bin's length also present in any other final selected bin, including those from the same binner. Each shared contig is counted once; excluded alternatives do not contribute. This diagnostic does not change selection and is not a contamination estimate.

A failed binner can be recorded while other successful branches complete; inspect `summary.json` for per-binner outcomes. Reuse requires matching inputs, parameters, software and verified outputs. Use a new output directory when changing a run.

## License

EukBinMiner is distributed under [CC BY-NC-SA 4.0](LICENSE). External tools and databases retain their respective licenses.
