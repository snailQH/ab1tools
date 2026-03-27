# ab1tools

**The comprehensive toolkit for AB1 chromatogram files** — generate, visualize, analyze, convert, compare, trim, and extract. Like [samtools](https://github.com/samtools/samtools) for BAM and [bedtools](https://github.com/arq5x/bedtools2) for BED, ab1tools fills the ecosystem gap for the AB1/ABIF chromatogram format.

ab1tools bridges next-generation sequencing (NGS) data with the Sanger chromatogram ecosystem. It converts BAM alignments from **any sequencing platform** (Illumina, Oxford Nanopore, PacBio, Ion Torrent, etc.) into realistic 4-channel AB1 files compatible with **SnapGene**, **Chromas**, **FinchTV**, and **ApE**.

## 1. Key Capabilities

| # | Capability | Commands | Description |
|---|-----------|----------|-------------|
| 1 | **Generate** | `from-seq`, `single`, `batch` | Create AB1 from BAM pileup or FASTA/FASTQ — unprecedented NGS-to-AB1 bridge |
| 2 | **Visualize** | `plot-ab1`, `smart` | Export publication-quality chromatogram PNGs with region selection and hetero highlighting |
| 3 | **Analyze** | `call`, `smart`, `stats` | Dual-mode SDVC variant detection (LOD <1% with control, ~3% without), quality metrics |
| 4 | **Convert** | `convert` | AB1 to FASTA, FASTQ, JSON, or CSV — bidirectional format bridge |
| 5 | **Inspect** | `view` | Display AB1 metadata, trace statistics, GC content, channel stats |
| 6 | **Compare** | `compare` | Signal correlation, base concordance, mismatch positions between two AB1 files |
| 7 | **Process** | `trim`, `extract` | Quality-based trimming (Mott's algorithm), sub-region extraction into new AB1 |

## 2. Installation

```bash
pip install -e .
```

Dependencies: pysam, numpy, matplotlib, biopython (Python >= 3.9)

## 3. The 12 Subcommands

### 3.1 `from-seq` — Generate AB1 from FASTA/FASTQ (no BAM required)

```bash
ab1tools from-seq sequence.fasta -o output/
ab1tools from-seq consensus.fastq -o output/ --name my_sample
```

### 3.2 `single` — Convert a single BAM + consensus sample

```bash
ab1tools single --bam alignment.bam --consensus consensus.fastq -o output/
```

### 3.3 `batch` — Batch convert all barcode samples

```bash
ab1tools batch /path/to/Analysis/ -o output/ --noise 5 --phasing 0.1 --decay 0.5
```

### 3.4 `smart` — Hetero detection + targeted PNGs + CSV report

```bash
ab1tools smart /path/to/Analysis/ -o output/ --noise 5 --phasing 0.1 --decay 0.5
```

Detects heterozygous sites (minor allele >= 5%), generates targeted chromatogram PNGs with orange-highlighted variant regions, and outputs `heterozygous_sites.csv`.

### 3.5 `plot-ab1` — Export custom region PNG from existing AB1

```bash
ab1tools plot-ab1 sample.ab1 --start 100 --end 200 --bases-per-row 50 -o region.png
```

### 3.6 `call` — SDVC variant detection from AB1 (dual-mode)

```bash
# Mode 1: Control-free (any AB1 file, LOD ~2-3%)
ab1tools call sample.ab1 --threshold 0.05

# Mode 2: Control-enhanced (with matched control, LOD <1%)
ab1tools call sample.ab1 --control wild_type.ab1 --threshold 0.03

# VCF output for pipeline integration
ab1tools call sample.ab1 --format vcf -o variants.vcf

# With annotated chromatogram plot
ab1tools call sample.ab1 --plot
```

The SDVC (Signal Deconvolution-based Variant Calling) algorithm features dynamic thresholding, Bayesian prior support, and multi-position consistency enhancement.

### 3.7 `view` — Display AB1 metadata and trace statistics

```bash
ab1tools view sample.ab1
ab1tools view sample.ab1 --json
```

Output: sample name, sequence length, GC content, base composition, per-channel signal statistics (mean, max, min, std).

### 3.8 `convert` — Convert AB1 to other formats

```bash
ab1tools convert sample.ab1 -f fasta -o sample.fasta
ab1tools convert sample.ab1 -f fastq -o sample.fastq
ab1tools convert sample.ab1 -f json -o sample.json
ab1tools convert sample.ab1 -f csv -o sample.csv
```

### 3.9 `stats` — Trace-level quality statistics

```bash
ab1tools stats sample.ab1
ab1tools stats sample.ab1 --format csv -o quality.csv
ab1tools stats sample.ab1 --format json
```

Output: mean/median quality, Q20/Q30 percentages, SNR, peak height distribution.

### 3.10 `compare` — Compare two AB1 files

```bash
ab1tools compare sample.ab1 control.ab1
ab1tools compare sample.ab1 replicate.ab1 --format json -o comparison.json
ab1tools compare sample.ab1 control.ab1 --plot
```

Output: base concordance %, per-channel Pearson correlation, mismatch positions.

### 3.11 `trim` — Quality-based trace trimming

```bash
ab1tools trim sample.ab1 -o trimmed.ab1 --quality-threshold 20
ab1tools trim sample.ab1 --quality-threshold 30 --min-length 100
```

Uses Mott's algorithm (Ewing et al., 1998) to find the optimal high-quality region.

### 3.12 `extract` — Extract sub-region into new AB1

```bash
ab1tools extract sample.ab1 --start 100 --end 200 -o region.ab1
ab1tools extract sample.ab1 --start 3000 --end 3050
```

Extracts any base range into a self-contained AB1 file with rebased peak positions. Like `samtools view -b file.bam chr:start-end` for AB1.

## 4. Signal Parameters

Apply to `from-seq`, `single`, `batch`, and `smart` modes:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--spacing` | 10 | Points between peaks |
| `--sigma` | 2.0 | Gaussian peak width |
| `--scale` | 1024 | Peak amplitude |
| `--noise` | 0 | Random noise level (5 recommended for realism) |
| `--phasing` | 0.0 | Tailing/drag effect (0.1 recommended) |
| `--decay` | 0.0 | Signal decay (0.5 = mild Sanger-like) |
| `--min-mapq` | 20 | Minimum mapping quality filter |
| `--no-plot` | off | Skip PNG generation (AB1 only) |

**Recommended realism:** `--noise 5 --phasing 0.1 --decay 0.5`

## 5. Platform Support

ab1tools is **platform-agnostic** — it works with any sequencing platform that produces BAM files:

| Platform | Tested | Notes |
|----------|--------|-------|
| Oxford Nanopore (MinION, PromethION) | Yes | Primary development platform |
| PacBio (Revio, Sequel II) | Compatible | Standard BAM format |
| Illumina (MiSeq, NovaSeq, etc.) | Compatible | Standard BAM format |
| Ion Torrent (S5, Genexus) | Compatible | Standard BAM format |
| MGI/DNBSEQ | Compatible | Standard BAM format |

The `from-seq` and `convert` modes accept any FASTA/FASTQ regardless of source.

## 6. Output Formats

| Format | Command | Description |
|--------|---------|-------------|
| **AB1** | from-seq, single, batch, smart, trim, extract | ABIF binary, compatible with SnapGene/Chromas/FinchTV/ApE |
| **PNG** | plot-ab1, smart, call --plot | Multi-panel chromatograms, configurable resolution |
| **CSV** | smart, call, stats, convert | Variant reports, quality metrics, per-position data |
| **VCF** | call --format vcf | VCF v4.2 with AF, CONF, SNR fields |
| **FASTA** | convert -f fasta | Standard sequence output |
| **FASTQ** | convert -f fastq | Sequence + quality scores |
| **JSON** | view --json, convert -f json, stats --format json | Structured data for programmatic use |

## 7. Architecture

```
ab1tools/
├── signal.py          # BAM pileup → base frequencies → 4-channel Gaussian traces
├── variant_caller.py  # SDVC dual-mode variant detection algorithm
├── abif_writer.py     # ABIF binary format writer
├── abif_reader.py     # ABIF binary format reader
├── visualize.py       # Multi-panel chromatogram PNG rendering
├── hetero.py          # Heterozygous site detection and reporting
└── cli.py             # CLI entry point with 12 subcommands
```

## 8. Demo Data

```
demo_data/
├── input/                  # 8 barcode samples (epi2me wf-amplicon format)
└── output/                 # Test outputs from all modes
    ├── batch/              # batch mode (8 AB1 + PNGs)
    ├── smart/              # smart mode (AB1 + hetero PNGs + CSV)
    ├── call-csv/           # SDVC variant calling (CSV + LOD analysis)
    ├── call-vcf/           # VCF output
    ├── call-plot/          # Annotated variant chromatograms
    ├── view/               # Metadata inspection
    ├── convert/            # Format conversion (FASTA, FASTQ, CSV)
    ├── stats/              # Quality statistics
    ├── compare/            # File comparison reports
    ├── trim/               # Quality-trimmed AB1
    └── extract/            # Extracted sub-regions
```

## 9. Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| [pysam](https://github.com/pysam-developers/pysam) | >= 0.21 | BAM reading and pileup |
| [numpy](https://numpy.org/) | >= 1.24 | Signal array operations |
| [matplotlib](https://matplotlib.org/) | >= 3.7 | Chromatogram rendering |
| [biopython](https://biopython.org/) | >= 1.81 | FASTA/FASTQ parsing |

## 10. Citation

If you use ab1tools in your research, please cite:

> ab1tools: a comprehensive toolkit for Sanger chromatogram generation, analysis, and NGS bridging. (Manuscript in preparation)

## 11. License

MIT License
