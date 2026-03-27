# ab1tools Wiki

**ab1tools** is the first comprehensive command-line toolkit for AB1 Sanger chromatogram files. It fills an ecosystem gap in bioinformatics: while BAM has [samtools](https://github.com/samtools/samtools), VCF has [bcftools](https://github.com/samtools/bcftools), and BED has [bedtools](https://github.com/arq5x/bedtools2), AB1/ABIF chromatogram files have had no unified CLI toolkit — until now.

## What ab1tools Does

ab1tools provides **12 subcommands** organized into 6 functional categories:

| Category | Commands | Purpose |
|----------|----------|---------|
| [**Generate**](Generate-AB1-Files) | `from-seq`, `single`, `batch` | Create AB1 chromatograms from NGS data or sequences |
| [**Visualize**](Visualize-Chromatograms) | `plot-ab1`, `smart` | Export publication-quality chromatogram PNGs |
| [**Analyze**](Analyze-Variants) | `call`, `smart`, `stats` | Variant detection (SDVC algorithm), quality metrics |
| [**Convert**](Convert-Formats) | `convert` | AB1 to FASTA, FASTQ, JSON, or CSV |
| [**Inspect**](Inspect-AB1-Files) | `view` | Display metadata, trace statistics, channel info |
| [**Process**](Process-AB1-Files) | `compare`, `trim`, `extract` | Compare, quality-trim, or extract sub-regions |

## Quick Start

```bash
# Install
pip install -e .

# Generate AB1 from a FASTA sequence
ab1tools from-seq my_sequence.fasta -o output/

# Generate AB1 from BAM alignment (any platform)
ab1tools single --bam aligned.bam --consensus consensus.fasta -o output/

# View AB1 metadata
ab1tools view output/my_sample.ab1

# Detect variants
ab1tools call output/my_sample.ab1

# Extract a region
ab1tools extract output/my_sample.ab1 --start 100 --end 500 -o region.ab1
```

## Platform Support

ab1tools works with **any sequencing platform** that produces BAM files:

- **Oxford Nanopore** (MinION, PromethION, Flongle)
- **PacBio** (Revio, Sequel II/IIe — HiFi reads)
- **Illumina** (MiSeq, NextSeq, NovaSeq, iSeq)
- **Ion Torrent** (S5, Genexus)
- **MGI/DNBSEQ**, **Element**, **Ultima**, and others

The `from-seq` and `convert` modes accept any FASTA/FASTQ regardless of source platform.

## Wiki Pages

- [Installation and Setup](Installation)
- [Generate AB1 Files](Generate-AB1-Files) — `from-seq`, `single`, `batch`
- [Visualize Chromatograms](Visualize-Chromatograms) — `plot-ab1`, `smart`
- [Analyze Variants (SDVC)](Analyze-Variants) — `call`, `stats`
- [Convert Formats](Convert-Formats) — `convert`
- [Inspect AB1 Files](Inspect-AB1-Files) — `view`
- [Process AB1 Files](Process-AB1-Files) — `compare`, `trim`, `extract`
- [Signal Parameters Guide](Signal-Parameters)
- [AB1 Format Reference](AB1-Format-Reference)
- [Use Cases and Workflows](Use-Cases)
- [FAQ](FAQ)
