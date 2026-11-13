# EukBinPicker

Recover eukaryotic genome bins with complementary binning methods.

## Introduction

Genome-resolved metagenomics helps reveal eukaryotic diversity and function, but genome recovery remains challenging. A typical workflow assembles quality-controlled reads, maps reads to the complete assembly, and identifies eukaryotic contigs with [EukContigMiner](https://github.com/WeStrainGroup/EukContigMiner) or similar tools. These contigs and their coverage profiles are then used to recover eukaryotic MAGs.

The binning performance of existing tools varies across samples and environments. Following the complementary-binner rationale of [DAS Tool](https://github.com/cmks/DAS_Tool), EukBinPicker (EBP) evaluates candidate bins from multiple sources and selects representatives based on quality estimates and shared contigs.

EBP integrates [CONCOCT](https://github.com/BinPro/CONCOCT), [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/), [VAMB](https://github.com/RasmussenLab/vamb) and [REMAG](https://github.com/danielzmbp/remag), which use contig sequences and coverage without requiring additional information, such as prokaryotic marker genes or taxonomic assignments, for binning. It also accepts existing bin sets, alone or alongside new binning. All candidates undergo the same quality assessment and selection workflow.

![EukBinPicker workflow](docs/ebp_workflow.png)

## Installation

Install with Conda:

```bash
conda create -n eukbinpicker --override-channels -c westraingroup -c conda-forge eukbinpicker=0.42
conda activate eukbinpicker
```

Or install a tagged release from GitHub with Git and Conda available:

```bash
git clone --depth 1 --branch v0.42 https://github.com/WeStrainGroup/EukBinPicker.git
cd EukBinPicker
bash install.sh --prefix "$HOME/conda-envs/eukbinpicker"
conda activate "$HOME/conda-envs/eukbinpicker"
```

Download reference databases separately. Provide an unpacked [EukCC database](https://github.com/Finn-Lab/EukCC) with `--eukcc-db`, or a [BUSCO lineage](https://busco.ezlab.org/) with `--busco-db`. Alternatively, set `EUKBINPICKER_EUKCC_DB` or `EUKBINPICKER_BUSCO_DB`.

## Quick start

Bin eukaryotic contigs using their coverage profiles, then select representative bins.

```bash
eukbinpicker \
  --sample sample_1 \
  --contigs euk_contigs.fa \
  --coverage assembly_depth.tsv \
  --eukcc-db /path/to/eukcc_db \
  --output sample_1_ebp
```

Select representative bins from existing bin sets.

```bash
eukbinpicker \
  --sample sample_1 \
  --bins /path/to/bins_A /path/to/bins_B \
  --eukcc-db /path/to/eukcc_db \
  --output sample_1_ebp
```

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

Use coverage from one or more samples mapped to the same complete assembly. EBP extracts the rows matching the supplied eukaryotic FASTA; all its contigs must have coverage entries. Zero coverage is valid, but missing entries are not.

**BAM input.** Pass `--coverage sample_A.bam sample_B.bam`. Use coordinate-sorted, indexed BAMs with unique basenames and the same reference IDs, lengths and order. For example, map one sample's clean paired-end reads with Bowtie2 and prepare the BAM with SAMtools:

```bash
bowtie2-build assembly.fa assembly_index
bowtie2 -x assembly_index \
  -1 sample_A.clean_R1.fq.gz -2 sample_A.clean_R2.fq.gz \
  --very-sensitive --no-unal --seed 6311 -p 4 \
  | samtools sort -@ 4 -o sample_A.bam -
samtools index sample_A.bam
```

For additional samples, repeat the mapping, sorting and indexing with the same assembly index. EBP uses these alignments without remapping reads.

**JGI tables.** Alternatively, generate a coverage table from the BAMs, with one mean/variance pair per sample:

```bash
jgi_summarize_bam_contig_depths \
--outputDepth assembly_depth.tsv \
sample_A.bam sample_B.bam
```

Pass `--coverage assembly_depth.tsv`; EBP recognizes the [JGI format](https://bitbucket.org/berkeleylab/metabat/src/master/) automatically. For single-sample coverage, use one BAM.

**Custom tables.** You can also supply a tab-separated table with a header, unique contig IDs and finite, non-negative coverage values. For columns named `sequence_id`, `A_mean`, `B_mean`, `A_var` and `B_var`, add:

```bash
--cov-contig-column sequence_id \
--cov-mean-column A_mean,B_mean \
--cov-variance-column A_var,B_var
```

Variance is optional; if provided, pair it with the mean columns in the same order. Only MetaBAT2 uses variance, requiring it for all samples or none; mean-only input uses its native `--cvExt` mode. The other three binners use mean coverage.

### Binning

| Parameter | Default | Description |
| --- | --- | --- |
| `--min-input-length` | `1500` | Minimum contig length for binning, in bp; imported bins remain intact. |
| `--binners` | `all` | All four tools or a comma-separated subset of `remag,concoct,metabat2,vamb`. |
| `--threads` | `4` | Threads supplied to each tool. |
| `--random-seed` | `6311` | Shared positive seed for binning. |
| `--gpu` | `F` | `T`: CUDA for selected REMAG/VAMB; other tools remain on CPU. `F`: CPU for all tools. |
| `--remag-parameters` | — | Quoted native REMAG options. |
| `--concoct-parameters` | — | Quoted native CONCOCT options. |
| `--metabat2-parameters` | — | Quoted native MetaBAT2 options. |
| `--vamb-parameters` | — | Quoted native VAMB options. |

EBP retains native algorithm defaults, with a shared `--min-input-length` and random seed. MetaBAT2 requires a contig-length threshold of at least 1,500 bp; VAMB requires at least 250 bp. REMAG and CONCOCT use the supplied threshold. Keep `--min-input-length` at 1,500 bp or above when running all four binners. This threshold applies to contigs, whereas `--min-retain-size` applies to whole bins.

EBP applies `--min-retain-size` before writing bin FASTAs: REMAG and VAMB use native export thresholds; CONCOCT uses bin membership lists. MetaBAT2 receives the same value through `--minClsSize`, which also affects cluster processing. The default matches MetaBAT2's native setting; `0` is passed through to disable its size cutoff.

REMAG's classifier is skipped because the supplied contigs are expected to have already been identified as eukaryotic. CONCOCT calculates fragment coverage from BAMs, or approximates it from parent-contig means when given tables.

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

### Quality assessment and selection

| Parameter | Default | Description |
| --- | --- | --- |
| `--min-retain-size` | `200000` | Minimum bin size to retain, in bp; `0` disables this filter. |
| `--min-selection-size` | `500000` | Minimum bin size entering QC and selection, in bp. |
| `--qc-tools` | `eukcc` | `eukcc`, `busco`, `eukcc,busco` or `all`. |
| `--qc-decision` | `eukcc` | Executed QC tool used for filtering and ranking. |
| `--eukcc-db` | — | Unpacked EukCC database; required when EukCC is selected. |
| `--busco-db` | — | Unpacked BUSCO lineage with `dataset.cfg` and `hmms/`; required when BUSCO is selected. |
| `--min-selection-completeness` | `0` | Minimum completeness, in %. |
| `--max-selection-contamination` | `100` | Maximum contamination or BUSCO proxy, in %. |
| `--contamination-penalty` | `5` | Contamination penalty in the ranking score. |
| `--min-overlap` | `0.50` | Overlap threshold between candidates from different sources. |

The retention threshold applies to both newly generated and imported bins. Only whole bins meeting it enter `1_raw/`; imported source files remain unchanged.

All retained bins meeting `--min-selection-size` undergo QC. By default, bins with valid QC estimates enter representative selection under the permissive 0% completeness and 100% contamination thresholds. Use `--min-selection-completeness` and `--max-selection-contamination` to set stricter quality requirements.

To rank with BUSCO, include it in `--qc-tools` and set `--qc-decision busco`. BUSCO's duplicated percentage is used as a contamination proxy.

### Information

| Parameter | Description |
| --- | --- |
| `--help` | Show command-line options. |
| `--version` | Print the installed version. |

## How selection works

Candidates passing the size and QC thresholds are ranked by:

```txt
score = completeness − contamination_penalty × contamination
```

Ties are resolved by higher completeness, lower contamination, fewer contigs, higher N50, then lexicographic order of `source__bin_id`.

EBP considers candidates in descending quality order. If a candidate overlaps an already selected, higher-ranked bin from a different source within the same assembly, and the overlap reaches `--min-overlap` (default: 0.50), that selected bin represents the candidate. Otherwise, the candidate is retained as another representative.

```txt
overlap(A, B) = shared_contig_bp / min(size(A), size(B))
```

A source is one native binner or one imported directory. Bins within a source are not compared. Shared contigs must match in original ID and normalized sequence; overlap is a proxy for competing reconstructions, not a taxonomic identity test. EBP preserves complete candidate bins without merging, splitting or refinement.

## Outputs

```text
sample_1_ebp/
├── 0_log/
│   ├── run.json                 # Inputs, commands, status and result checksums
│   └── run.log                  # Progress, warnings and tool output
├── 1_raw/
│   └── <source>/
│       └── <original_bin>.fa    # Retained candidate bins
├── 2_qc/
│   ├── summary.tsv              # Statistics, QC metrics and status for every retained bin
│   └── evidence/
│       └── <source>/
│           └── <bin_id>/
│               ├── eukcc/       # Quality report and marker evidence, when available
│               └── busco/       # Short summary and compressed full marker table
└── 3_selection/
    ├── representatives.tsv      # Final bin IDs, original sources and quality statistics
    ├── selection.tsv            # Candidate decisions and matching representatives
    └── bins/
        ├── bin_1.fa
        ├── bin_2.fa
        └── ...
```

Final bins are numbered in selection order. `representatives.tsv` maps each new ID to its source, original bin ID and filename, and includes quality metrics, shared-contig statistics and represented candidates. `selection.tsv` records why each candidate was selected or excluded, its assigned representative and any additional matching representatives. Candidate references use `source/bin_id`.

QC evidence is grouped by source and bin; only the requested QC tools are included. `summary.tsv` distinguishes bins below the QC size threshold, unassessable bins and valid quality estimates.

`run.json` records workflow status, per-stage outcomes, actual commands and checksums. A completed run may contain failed native binners; check `partial` and `summary.native_binner_status`. Valid imported bins can still be selected if native binning fails or no contigs pass the input-length threshold.

Temporary files are limited to `.work/`. Each completed stage releases its working files after its results have been checked; failed-stage files are retained for diagnosis. After an interrupted or failed run, repeat the same command to reuse verified completed stages and rerun unfinished stages. Reuse requires matching inputs, parameters, software and intact results. Use a new output directory when these change, including when upgrading from an earlier version.

## Contact

Xinyu Wang: [wangxinyu30@westlake.edu.cn](mailto:wangxinyu30@westlake.edu.cn)

## License

EukBinPicker is distributed under [CC BY-NC-SA 4.0](https://github.com/WeStrainGroup/EukBinPicker/blob/main/LICENSE). External tools and databases retain their respective licenses.
