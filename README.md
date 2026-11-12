# EukBinPicker

**Recover eukaryotic genome bins with complementary binning methods.**

## Introduction

Genome-resolved metagenomics helps reveal eukaryotic diversity and function, but genome recovery remains challenging. A typical workflow assembles quality-controlled reads, maps reads to the complete assembly, and identifies eukaryotic contigs with [EukContigMiner](https://github.com/WeStrainGroup/EukContigMiner) or similar tools. These contigs and their coverage profiles are then used to recover metagenome-assembled genomes (MAGs).

Binning performance varies across samples and environments. Following the complementary-binner rationale of [DAS Tool](https://github.com/cmks/DAS_Tool), EukBinPicker (EBP) evaluates candidate bins from multiple sources and selects representatives based on quality estimates and shared contigs.

EBP integrates [CONCOCT](https://github.com/BinPro/CONCOCT), [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/), [VAMB](https://github.com/RasmussenLab/vamb) and [REMAG](https://github.com/danielzmbp/remag), which use contig sequences and coverage without requiring prokaryotic marker genes or taxonomic assignments for binning. It also accepts existing bin sets, alone or alongside new binning. All candidates undergo the same quality assessment and selection workflow.

![EukBinPicker workflow](docs/ebp_workflow.png)

## Installation

Install with Conda on Linux x86-64:

```bash
conda create -n eukbinpicker --override-channels -c westraingroup -c conda-forge eukbinpicker=0.3.1
conda activate eukbinpicker
```

Or install from GitHub with Conda and curl available:

```bash
curl -fsSL https://raw.githubusercontent.com/WeStrainGroup/EukBinPicker/v0.3.1/install.sh | bash -s -- --prefix "$HOME/conda-envs/eukbinpicker"
conda activate "$HOME/conda-envs/eukbinpicker"
```

Download reference databases separately. Provide an unpacked [EukCC database](https://github.com/Finn-Lab/EukCC) with `--eukcc-db`, or a [BUSCO lineage](https://busco.ezlab.org/) with `--busco-db`. Alternatively, set `EUKBINPICKER_EUKCC_DB` or `EUKBINPICKER_BUSCO_DB`.

## Quick start

### From contigs

Bin eukaryotic contigs using their coverage profiles, then select representative bins.

```bash
eukbinpicker \
  --sample sample_1 \
  --contigs euk_contigs.fa \
  --coverage assembly_depth.tsv \
  --eukcc-db /path/to/eukcc_db \
  --output sample_1_ebp
```

### From existing bin sets

Select representative bins from existing bin sets.

```bash
eukbinpicker \
  --sample sample_1 \
  --bins /path/to/bins_A /path/to/bins_B \
  --eukcc-db /path/to/eukcc_db \
  --output sample_1_ebp
```

### From both

Run new binning and select representatives jointly with existing bins from the same assembly.

```bash
eukbinpicker \
  --sample sample_1 \
  --contigs euk_contigs.fa \
  --coverage assembly_depth.tsv \
  --bins /path/to/bins_A /path/to/bins_B \
  --eukcc-db /path/to/eukcc_db \
  --output sample_1_ebp
```

Selected FASTA files are written to `3_selection/bins/` within the output directory.

## Coverage inputs

`--coverage`: coverage tables **or** indexed BAMs. Do not mix the two formats. Map reads from the source sample or multiple samples to the same assembly. Full-assembly coverage is accepted; EBP extracts the supplied contigs. Every input must include all FASTA contig IDs, including those below the length threshold. Zero coverage is valid; missing entries are not.

### Tables

Use tab-separated tables with headers, optionally gzip-compressed. JGI depth tables and mean-coverage matrices are detected automatically. IDs must be unique, values finite and non-negative, and supplied contig lengths must match the FASTA.

For custom headers, specify the contig and mean columns. For this table:

```txt
sequence_id	A_mean	B_mean	A_var	B_var
contig_1	12.4	8.1	2.0	1.4
contig_2	0.0	3.2	0.0	0.6
```

Add these options to the binning command:

```txt
--cov-contig-column sequence_id
--cov-mean-column A_mean,B_mean
--cov-variance-column A_var,B_var
```

Variance columns are optional and follow the mean-column order. If using MetaBAT2, provide variance for all coverage dimensions or none. Column settings apply to every table.

### BAMs

Use one coordinate-sorted BAM per mapping sample, each with a `.bai` or `.csi` index. BAM basenames must be unique; reference IDs, lengths and order must match across BAMs. The contig FASTA must match a subset of these references.

Example: `--coverage sample_A.bam sample_B.bam`. EBP calculates coverage from the alignments without remapping reads.

## Parameters

Default values are shown below. **—** means no preset value; input requirements are stated in the description.

### Input and output

| Parameter | Default | Description |
| --- | --- | --- |
| `--sample` | — | Required sample label: letters, digits, `.`, `_` or `-`. |
| `--contigs` | — | Eukaryotic FASTA with unique IDs; required with `--coverage`. |
| `--coverage` | — | Coverage tables or indexed BAMs; multiple files and repeated arguments accepted. |
| `--bins` | — | One directory per source, one FASTA per bin (`.fa`, `.fna`, `.fasta`, optionally `.gz`); no nested directories. Multiple inputs accepted. |
| `--cov-contig-column` | — | Contig ID column in custom tables. |
| `--cov-mean-column` | — | Comma-separated mean columns in custom tables. |
| `--cov-variance-column` | — | Comma-separated variance columns, paired with means. |
| `--output` | — | Required output directory. |

Provide `--coverage`, `--bins`, or both. Use eukaryotic inputs from the same assembly with original contig IDs. Imported bins remain intact; assembly provenance is assumed. Coverage-column selectors apply only to tables; omit them for automatic detection.

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

EBP retains native algorithm defaults with a shared input threshold and seed. CONCOCT uses 10-kb fragments, native normalization and 400 clusters. MetaBAT2 and REMAG retain their native 200-kb and 500-kb output minima. VAMB retains native training and bin splitting, with nonempty FASTA export enabled. `--min-retain-length` does not lower native output minima. REMAG's classifier is skipped: supply contigs already identified as eukaryotic.

MetaBAT2 uses observed variance or native mean-only coverage; REMAG, CONCOCT and VAMB use means. CONCOCT calculates fragment coverage from BAMs, or approximates it from parent-contig means when given tables.

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

EBP considers candidates in descending quality order. If a candidate overlaps an already selected, higher-ranked bin from a different source within the same assembly, and the overlap reaches `--min-overlap` (default: 0.50), that selected bin represents the candidate. Otherwise, the candidate is retained as another representative.

```txt
overlap(A, B) = shared_contig_bp / min(length(A), length(B))
```

A source is one native binner or one imported directory. Bins within a source are not compared. Shared contigs must match in original ID and normalized sequence; overlap is a proxy for competing reconstructions, not a taxonomic identity test. EBP preserves complete candidate bins without merging, splitting or refinement.

## Outputs

| Path | Contents |
| --- | --- |
| `0_log/` | Input metadata, run summary, commands and internal records. |
| `1_raw/` | Retained bin FASTA files, grouped by source. |
| `2_qc/` | Per-tool QC results and `summary.tsv` with source, SeqKit statistics and QC metrics for each assessed bin. |
| `3_selection/` | `bins/` with selected FASTAs and checksums; `representatives.tsv` with quality metrics, shared-contig statistics and represented bins. |

Check `0_log/summary.json` for source and selection counts. The `native_binner_status` field distinguishes `completed`, `failed` and `skipped_empty_input`; completion does not imply that a binner produced bins. Logs are in `0_log/logs/`. Imported candidates can still be processed if native branches fail or no contigs pass the input-length threshold. Reusing an output directory requires matching inputs, parameters, software and intact outputs; otherwise, use a new directory.

## License

EukBinPicker is distributed under [CC BY-NC-SA 4.0](https://github.com/WeStrainGroup/EukBinPicker/blob/main/LICENSE). External tools and databases retain their respective licenses.
