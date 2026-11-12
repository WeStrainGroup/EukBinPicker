# EukBinPicker

**Recover eukaryotic genome bins with complementary binning methods.**

Version **0.3** · [GitHub](https://github.com/WeStrainGroup/EukBinPicker)

## Introduction

Genome-resolved metagenomics helps reveal eukaryotic diversity and function, but genome recovery remains challenging. A typical workflow assembles quality-controlled reads, maps reads to the complete assembly, and identifies eukaryotic contigs with [EukContigMiner](https://github.com/WeStrainGroup/EukContigMiner) or similar tools. These contigs and their coverage profiles are then used to recover metagenome-assembled genomes (MAGs).

Binning performance varies across samples and environments. Following the complementary-binner rationale of [DAS Tool](https://github.com/cmks/DAS_Tool), EukBinPicker (EBP) evaluates candidate bins from multiple sources and selects representatives based on quality estimates and shared contigs.

EBP integrates [CONCOCT](https://github.com/BinPro/CONCOCT), [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/), [VAMB](https://github.com/RasmussenLab/vamb) and [REMAG](https://github.com/danielzmbp/remag), which use contig sequences and coverage without requiring prokaryotic marker genes or taxonomic assignments for binning. It also accepts existing bin sets, alone or alongside new binning. All candidates undergo the same quality assessment and selection workflow.

![EukBinPicker workflow](docs/ebp_workflow.png)

## Installation

Install with Conda on Linux x86-64:

```bash
conda create -n eukbinpicker --override-channels -c westraingroup -c conda-forge eukbinpicker=0.3
conda activate eukbinpicker
```

Or install from GitHub with Conda and curl available:

```bash
curl -fsSL https://raw.githubusercontent.com/WeStrainGroup/EukBinPicker/v0.3/install.sh | bash -s -- --prefix "$HOME/conda-envs/eukbinpicker"
conda activate "$HOME/conda-envs/eukbinpicker"
```

Download reference databases separately. Provide an unpacked [EukCC database](https://github.com/Finn-Lab/EukCC) with `--eukcc-db`, or a [BUSCO lineage](https://busco.ezlab.org/) with `--busco-db`. Alternatively, set `EUKBINPICKER_EUKCC_DB` or `EUKBINPICKER_BUSCO_DB`.

## Quick start

### Bin contigs

```bash
eukbinpicker \
  --sample sample_1 \
  --contigs euk_contigs.fa \
  --coverage assembly_depth.tsv \
  --eukcc-db /path/to/eukcc_db \
  --output sample_1_ebp
```

### Select from existing bins

```bash
eukbinpicker \
  --sample sample_1 \
  --bins /path/to/bins_A /path/to/bins_B \
  --eukcc-db /path/to/eukcc_db \
  --output sample_1_selection
```

### Combine new and existing bins

```bash
eukbinpicker \
  --sample sample_1 \
  --contigs euk_contigs.fa \
  --coverage assembly_depth.tsv \
  --bins /path/to/bins_A /path/to/bins_B \
  --eukcc-db /path/to/eukcc_db \
  --output sample_1_combined
```

Binning defaults to all four tools, CPU, 4 threads and contigs ≥1,500 bp. Bins ≥500,000 bp enter EukCC assessment. Selected FASTA files are written to `3_selection/bins/` within the output directory.

## Inputs

### Contigs and existing bins

Supply eukaryotic contigs for binning and eukaryotic bins for import. EBP does not classify inputs; REMAG's classifier is skipped. All inputs must come from the same assembly and retain original contig IDs.

`--contigs` is required with `--coverage`. Coverage may represent the source sample alone or multiple samples mapped to that assembly. Full-assembly coverage is supported: EBP selects rows matching the supplied FASTA and preserves their values. Every supplied FASTA contig must have coverage, including those below the length threshold; missing values are not filled with zero.

`--bins` accepts one directory per source, including tools or parameter settings outside the bundled binners. Each directory must directly contain `.fa`, `.fna` or `.fasta` files, optionally gzip-compressed. Multiple directories and repeated `--bins` arguments are supported; subdirectories are not searched.

Directory names label imported sources. Conflicting names receive numeric suffixes, such as `metabat2_2`. EBP records paths and checksums, validates FASTA structure, and trusts assembly provenance. Imported bins need no coverage and remain intact, including contigs below `--min-input-length`. Retained bins are copied to `1_raw/`; gzip files are decompressed and source files are unchanged.

With `--bins` alone, omit contigs, coverage and binner options. With combined input, generated and imported bins retain separate source labels and enter one shared QC and selection workflow.

### Coverage tables

`--coverage` accepts one or more tab-separated tables with headers, either plain or gzip-compressed. Standard JGI depth tables and mean-coverage matrices are detected automatically. Multiple files and repeated `--coverage` arguments are supported; all must refer to the same assembly.

For custom columns, provide the contig and mean column names together. Variance columns are optional and pair with means in the same order. The same selectors apply to every input table. Coverage dimensions use table headers, or filenames for multiple single-coverage generic/custom tables.

For example, `custom.tsv`:

```txt
sequence_id	A_mean	B_mean	A_var	B_var
contig_1	12.4	8.1	2.0	1.4
contig_2	0.0	3.2	0.0	0.6
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

Contig IDs must be unique, values finite and non-negative, and any supplied lengths must match the FASTA. Tables are validated before subsetting.

MetaBAT2 uses observed variance when available and native `--cvExt` mode for mean-only input; variance must cover all dimensions or none when MetaBAT2 is selected. REMAG, CONCOCT and VAMB use means. With table input, CONCOCT approximates fragment coverage using each parent contig's mean.

### BAM files

Use BAM files in place of coverage tables:

```bash
eukbinpicker \
  --sample sample_1 \
  --contigs euk_contigs.fa \
  --coverage sample_A.bam sample_B.bam \
  --eukcc-db /path/to/eukcc_db \
  --output sample_1_ebp
```

Provide one coordinate-sorted BAM per mapping sample, each with a `.bai` or `.csi` index, a unique basename, and the same ordered reference IDs and lengths. The input FASTA must be a matching subset of the mapped assembly. Do not mix BAMs and tables in one run.

EBP reuses alignments without remapping. The installer includes SAMtools and `jgi_summarize_bam_contig_depths`; JGI calculates contig means and variances, while CONCOCT calculates fragment coverage directly from BAMs.

## Parameters

Default values are shown below. **—** means no preset value; input requirements are stated in the description.

### Input and output

| Parameter | Default | Description |
| --- | --- | --- |
| `--sample` | — | Required sample label: letters, digits, `.`, `_` or `-`. |
| `--contigs` | — | Eukaryotic FASTA with unique IDs; required with `--coverage`. |
| `--coverage` | — | Coverage tables or indexed BAMs; multiple files and repeated arguments accepted. |
| `--bins` | — | Existing bin directories; multiple directories and repeated arguments accepted. |
| `--cov-contig-column` | — | Contig ID column in custom tables. |
| `--cov-mean-column` | — | Comma-separated mean columns in custom tables. |
| `--cov-variance-column` | — | Comma-separated variance columns, paired with means. |
| `--output` | — | Required output directory. |

Provide `--coverage`, `--bins`, or both. Coverage-column selectors apply only to tables; omit them for automatic format detection.

### Binning and resources

| Parameter | Default | Description |
| --- | --- | --- |
| `--binners` | `all` | All four tools or a comma-separated subset of `remag,concoct,metabat2,vamb`. |
| `--threads` | `4` | Threads supplied to each tool. |
| `--random-seed` | `6311` | Shared positive seed for binning. |
| `--gpu` | `F` | `T`: CUDA for selected REMAG/VAMB; `F`: CPU. |
| `--remag-parameters` | — | Quoted native REMAG options. |
| `--concoct-parameters` | — | Quoted native CONCOCT options. |
| `--metabat2-parameters` | — | Quoted native MetaBAT2 options. |
| `--vamb-parameters` | — | Quoted native VAMB options. |

EBP retains native algorithm defaults with a shared input threshold and seed. CONCOCT uses 10-kb fragments, native normalization and 400 clusters. MetaBAT2 and REMAG retain their native 200-kb and 500-kb output minima. VAMB retains native training and bin splitting, with nonempty FASTA export enabled. `--min-retain-length` does not lower native output minima.

To customize a binner, pass native options in one quoted string:

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

Use EBP options for input/output paths, threads, seeds, devices and input length. Native options and executed commands are recorded with the run.

REMAG and VAMB run sequentially and can share one allocated GPU. CONCOCT, MetaBAT2 and QC use CPU. `--gpu T` checks CUDA availability; it does not request hardware from a scheduler.

### Quality assessment

| Parameter | Default | Description |
| --- | --- | --- |
| `--qc-tools` | `eukcc` | `eukcc`, `busco`, `eukcc,busco` or `all`. |
| `--qc-decision` | `eukcc` | Executed QC tool used for filtering and ranking. |
| `--eukcc-db` | — | Unpacked EukCC database; required when EukCC is selected. |
| `--busco-db` | — | Unpacked BUSCO lineage with `dataset.cfg` and `hmms/`; required when BUSCO is selected. |

Database paths may also come from the environment variables listed under Installation. To rank with BUSCO, include it in `--qc-tools` and set `--qc-decision busco`. Its duplicated percentage is used as a contamination proxy.

### Filtering and selection

| Parameter | Default | Description |
| --- | --- | --- |
| `--min-input-length` | `1500` | Minimum contig length for binning, in bp; imported bins remain intact. |
| `--min-retain-length` | `0` | Additional bin-retention minimum, in bp; `0` disables it. |
| `--min-selection-length` | `500000` | Minimum bin length entering QC and selection, in bp. |
| `--min-selection-completeness` | `0` | Minimum completeness, in %. |
| `--max-selection-contamination` | `100` | Maximum contamination or BUSCO proxy, in %. |
| `--contamination-penalty` | `5` | Contamination penalty in the ranking score. |
| `--min-overlap` | `0.50` | Overlap threshold between candidates from different sources. |

Thresholds are inclusive. Selection length must be at least retention length. Valid QC estimates are required even with the permissive completeness and contamination defaults; selected bins are not necessarily high-quality MAGs.

### Information

| Parameter | Description |
| --- | --- |
| `--help` | Show command-line options. |
| `--version` | Print the installed version. |

## How selection works

Candidates passing the length and QC thresholds are ranked by:

```txt
score = completeness − contamination_penalty × contamination
```

EBP visits candidates in descending quality order. A candidate is excluded if its overlap with an already selected bin from a different source reaches `--min-overlap` (default: 0.50):

```txt
overlap(A, B) = shared_contig_bp / min(length(A), length(B))
```

A source is one native binner or one imported directory. Bins within a source are not compared. Shared contigs must match in original ID and normalized sequence; overlap identifies competing reconstructions, not taxonomic identity. EBP preserves complete candidate bins without merging, splitting or refinement.

## Outputs

| Path | Contents |
| --- | --- |
| `0_log/` | Input metadata, run summary, commands and internal records. |
| `1_raw/` | Retained bin FASTA files, grouped by source. |
| `2_qc/eukcc/`, `2_qc/busco/` | Results from each requested QC tool. |
| `2_qc/summary.tsv` | Source, SeqKit statistics and QC metrics for each assessed bin. |
| `3_selection/representatives.tsv` | Selected representatives, quality metrics, shared-contig statistics and represented bins. |
| `3_selection/bins/` | Selected FASTA files and their checksum manifest. |

`represented_bins` lists each representative and its excluded alternatives. Candidates matching multiple representatives receive one primary assignment; all matches are recorded in `0_log/selection/alternatives.tsv`.

Check `0_log/summary.json` for completed native binners, imported sources and selection counts. Imported candidates can still be processed if native branches fail or no contigs pass the input-length threshold. Reusing an output directory requires matching inputs, parameters, software and intact outputs; otherwise, use a new directory.

## License

EukBinPicker is distributed under [CC BY-NC-SA 4.0](https://github.com/WeStrainGroup/EukBinPicker/blob/main/LICENSE). External tools and databases retain their respective licenses.
