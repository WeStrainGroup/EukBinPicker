# EukBinPicker

**Recover eukaryotic genome bins with complementary binning methods.**

Version **0.1** · [GitHub](https://github.com/WeStrainGroup/EukBinPicker)

## Introduction

While challenging, genome-resolved metagenomics helps reveal the taxonomic and functional diversity of environmental eukaryotes and their ecological roles. Different binning methods can recover different genomes from the same sample. Following the complementary-binner rationale of [DAS Tool](https://github.com/cmks/DAS_Tool), EukBinPicker (EBP) evaluates multiple candidate sets to select representative eukaryotic bins.

EBP combines [REMAG](https://github.com/danielzmbp/remag), [CONCOCT](https://github.com/BinPro/CONCOCT), [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/) and [VAMB](https://github.com/RasmussenLab/vamb). These methods use contig sequences and coverage without requiring additional info like prokaryotic marker genes or taxonomy for binning and offer complementary approaches: contrastive learning, Gaussian mixture models, graph clustering and variational autoencoders.

A typical workflow is to assemble quality-controlled reads, map reads to the complete assembly, and identify eukaryotic contigs with [EukContigMiner](https://github.com/WeStrainGroup/EukContigMiner) or other similar tools. Supply the resulting eukaryotic FASTA and related coverage table or BAM files to EBP. Coverage may come from the source sample alone or from multiple samples mapped to that assembly.

```txt
Eukaryotic contigs + Coverage tables or BAMs
                          ↓
             CONCOCT · MetaBAT2 · VAMB · REMAG
                          ↓
                 EukCC or BUSCO assessment
                          ↓
          Rank overlapping candidates within the assembly
                          ↓
   Select the highest-ranked cross-binner representative
               for each putative genome
```

## Installation

Install on **Linux x86-64** from the [WeStrainGroup Conda channel](https://anaconda.org/westraingroup/eukbinpicker):

```bash
conda create -n eukbinpicker --override-channels -c westraingroup -c conda-forge eukbinpicker=0.1
conda activate eukbinpicker
```

Alternatively, install directly from GitHub with Conda and curl available:

```bash
curl -fsSL https://raw.githubusercontent.com/WeStrainGroup/EukBinPicker/v0.1/install.sh | bash -s -- --prefix "$HOME/conda-envs/eukbinpicker"
conda activate "$HOME/conda-envs/eukbinpicker"
```

Tool versions: REMAG 0.4.4, CONCOCT 1.1.0, MetaBAT2 2.18, VAMB 5.0.4, EukCC 2.1.3 and BUSCO 6. PyTorch 2.6.0 with CUDA 12.4 supports both execution modes.

Reference databases are downloaded separately: provide the unpacked [EukCC database](https://github.com/Finn-Lab/EukCC) with `--eukcc-db`, or a [BUSCO lineage](https://busco.ezlab.org/) with `--busco-db`. These paths can also be set through `EUKBINPICKER_EUKCC_DB` and `EUKBINPICKER_BUSCO_DB`.

## Quick start

```bash
eukbinpicker \
  --sample sample_1 \
  --contigs euk_contigs.fa \
  --coverage assembly_depth.tsv \
  --eukcc-db /path/to/eukcc_db \
  --output sample_1_ebp
```

By default, EBP runs all four binners on CPU with 4 threads, retains contigs ≥1,500 bp, and assesses bins ≥500,000 bp with EukCC. Selected FASTA files are written to `sample_1_ebp/3_selection/bins/`.

## Coverage input

`--coverage` accepts tables or BAM files, including multiple files or repeated `--coverage` arguments. All dimensions must refer to the same assembly; do not mix BAMs and tables in one run.

**Full-assembly coverage is supported.** EBP selects the rows matching the input eukaryotic FASTA, preserving its contig order and coverage values. Other contigs are excluded and counted in run metadata. Every supplied FASTA contig must have coverage; missing values are not replaced with zero.

### Coverage tables

[Standard JGI depth tables](https://bitbucket.org/berkeleylab/metabat/) are recognized automatically. Generic mean-coverage matrices and custom schemas are also supported. Tables must be tab-separated, with a header; plain and gzip-compressed files are accepted.

For example, `custom.tsv` may contain:

```txt
sequence_id	A_mean	B_mean	A_var	B_var
contig_1	12.4	8.1	2.0	1.4
contig_2	0.0	3.2	0.0	0.6
contig_3	5.7	2.1	0.8	0.3
```

```bash
eukbinpicker \
  --sample sample_1 \
  --contigs euk_contigs.fa \
  --coverage custom.tsv \
  --cov-contig-column sequence_id \
  --cov-mean-column A_mean,B_mean \
  --cov-variance-column A_var,B_var \
  --eukcc-db /path/to/eukcc_db \
  --output sample_1_ebp
```

If the FASTA contains only `contig_1` and `contig_3`, only those rows enter binning. Contig IDs must be unique, values finite and non-negative, and any supplied contig lengths must match the FASTA. Input tables are validated before subsetting.

Variance columns are optional and pair with mean columns in the same order. MetaBAT2 uses observed variance when available and its native `--cvExt` mode for mean-only input; variance must be present for all dimensions or none. REMAG, CONCOCT and VAMB use mean coverage.

### BAM files

```bash
eukbinpicker \
  --sample sample_1 \
  --contigs euk_contigs.fa \
  --coverage sample_A.bam sample_B.bam \
  --eukcc-db /path/to/eukcc_db \
  --output sample_1_ebp
```

Provide one coordinate-sorted BAM per mapping sample, each with a `.bai` or `.csi` index, unique basename, and the same ordered reference IDs and lengths. The eukaryotic FASTA must be a matching subset of the mapped assembly. EBP reuses these alignments without remapping reads.

The installer includes SAMtools and `jgi_summarize_bam_contig_depths`. JGI produces contig means and variances; CONCOCT computes native fragment coverage from the BAMs. With table input, CONCOCT approximates fragment coverage using each parent contig's mean.

## Parameters

### Input and output

| Parameter | Default | Description |
| --- | --- | --- |
| `--sample` | Required | Assembly/sample label: letters, digits, `.`, `_` or `-`. |
| `--contigs` | Required | Eukaryotic contig FASTA with unique IDs. |
| `--coverage` | Required | One or more coverage tables or indexed BAMs; repeatable. |
| `--cov-contig-column` | Auto | Contig ID column name in custom tables. |
| `--cov-mean-column` | Auto | Mean-coverage column names in custom tables, comma-separated. |
| `--cov-variance-column` | Auto | Coverage-variance column names in custom tables, comma-separated and paired with means. |
| `--output` | Required | Run output directory. |

Standard formats are detected automatically. For custom schemas, specify the contig and mean columns together; variance columns are optional. Selectors apply only to tables, and the same selectors apply to every supplied table. Coverage dimensions use table headers, or filenames for multiple generic or custom single-coverage tables.

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

Native algorithm defaults are retained, with the shared input threshold and seed. REMAG's classifier is skipped. CONCOCT uses standard 10-kb fragments, native normalization and 400 clusters. MetaBAT2 and REMAG retain their native 200-kb and 500-kb output minima. VAMB retains native training and bin splitting, with nonempty FASTA export enabled. The retention setting does not lower these native minima.

**Customizing a binner.** Use the tool's native option names inside one quoted string:

```bash
eukbinpicker \
  --sample sample_1 \
  --contigs euk_contigs.fa \
  --coverage assembly_depth.tsv \
  --eukcc-db /path/to/eukcc_db \
  --output sample_1_custom \
  --concoct-parameters "--clusters 200" \
  --metabat2-parameters "--maxP 90"
```

Input/output paths, threads, seeds, devices and shared input-length settings are managed by EBP; use their public EBP options instead. Native options and actual commands are recorded with the run.

### Quality assessment

| Parameter | Default | Description |
| --- | --- | --- |
| `--qc-tools` | `eukcc` | `eukcc`, `busco`, `eukcc,busco` or `all`. |
| `--qc-decision` | `eukcc` | Executed QC tool used for filtering and ranking. |
| `--eukcc-db` | Environment | Unpacked EukCC database; required when EukCC is selected. |
| `--busco-db` | Environment | Unpacked BUSCO lineage containing `dataset.cfg` and `hmms/`. |

To rank with BUSCO, include it in `--qc-tools` and set `--qc-decision busco`. BUSCO duplicated percentage is used as a contamination proxy.

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

## How selection works

Candidates passing the length and QC thresholds are ranked by:

```txt
score = completeness − contamination_penalty × contamination
```

EBP visits candidates in descending quality order. A candidate is excluded if it overlaps an already selected bin from a different binner:

```txt
overlap(A, B) = shared_contig_bp / min(length(A), length(B))
```

The default cutoff is 0.50. This selects the highest-ranked representative among competing cross-binner candidates within the assembly. Shared contigs must match in original ID and normalized sequence; overlap is a proxy for competing genome reconstructions, not a taxonomic identity test. Bins from the same binner are not compared. EBP preserves complete candidate bins without merging, splitting or refinement.

## Outputs

| Path | Contents |
| --- | --- |
| `0_log/` | `run.json`, `summary.json`, readable logs and internal work records. |
| `1_raw/<binner>/` | Original bin FASTA files passing the retention threshold, grouped by binner. |
| `2_qc/eukcc/`, `2_qc/busco/` | Results from each requested QC tool. |
| `2_qc/summary.tsv` | One row per assessed bin: originating binner, SeqKit statistics and all requested QC metrics. |
| `3_selection/representatives.tsv` | One row per selected representative: selection reason, SeqKit and QC statistics, shared-contig statistics and represented bins. |
| `3_selection/bins/` | Selected FASTA files and their checksum manifest. |

`represented_bins` includes the representative itself and excluded candidates assigned to it. When a candidate matches multiple representatives, its primary assignment follows the selection algorithm; all matches remain available in `0_log/selection/alternatives.tsv`.

`shared_contig_fraction` is the fraction of a selected bin's length also present in any other final selected bin, including those from the same binner. Each shared contig is counted once; excluded alternatives do not contribute. This diagnostic does not change selection and is not a contamination estimate.

A failed binner can be recorded while other successful branches complete; inspect `0_log/summary.json` for per-binner outcomes. Reuse requires matching inputs, parameters, software and verified outputs. Use a new output directory when changing a run.

## License

EukBinPicker is distributed under [CC BY-NC-SA 4.0](https://github.com/WeStrainGroup/EukBinPicker/blob/main/LICENSE). External tools and databases retain their respective licenses.
