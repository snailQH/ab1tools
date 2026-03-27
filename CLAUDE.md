# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**ab1tools** is a toolkit for generating, analyzing, and visualizing AB1 Sanger chromatogram files from various sequencing inputs. It bridges Nanopore amplicon sequencing data with the Sanger chromatogram ecosystem, enabling researchers to use familiar tools like SnapGene, Chromas, FinchTV, and ApE to inspect Nanopore consensus results.

Key capabilities:
1. **BAM → AB1 conversion** — realistic 4-channel traces from Nanopore BAM pileup frequencies, with mixed peaks at variant sites
2. **FASTA/FASTQ → AB1 conversion** — generate AB1 from standalone sequences (no BAM required), for plasmid assemblies or synthetic constructs
3. **Heterozygous site detection** — auto-detect variants, generate targeted PNGs with highlighted mixed-peak regions, plus CSV report
4. **Custom region PNG export** — read existing AB1 files and export PNGs for any region with user-defined base range and row width
5. **Batch processing** — auto-detect and process all barcode samples from epi2me wf-amplicon output in one command

## Commands

```bash
# Install
pip install -e .

# 1. From sequence file (FASTA or FASTQ, no BAM required)
ab1tools from-seq sequence.fasta -o output/
ab1tools from-seq consensus.fastq -o output/ --name my_sample

# 2. Single sample (BAM + consensus)
ab1tools single --bam sample.bam --consensus consensus.fastq -o output/

# 3. Batch (auto-detects barcodeXX/ subdirectories)
ab1tools batch /path/to/Analysis/ -o output/

# 4. Smart mode (hetero detection + targeted PNGs)
ab1tools smart /path/to/Analysis/ -o output/ --noise 5 --phasing 0.1 --decay 0.5

# 5. Plot AB1 (export custom region PNG from existing AB1)
ab1tools plot-ab1 sample.ab1 --start 100 --end 200 --bases-per-row 50 -o region.png
```

## Architecture

```
ab1tools/
├── signal.py        # BAM pileup → base frequencies → 4-channel Gaussian traces
├── abif_writer.py   # ABIF binary format writer (DATA9-12, PBAS1, PLOC1, FWO_, SMPL1)
├── abif_reader.py   # ABIF binary format reader (for plot-ab1 mode)
├── visualize.py     # Multi-panel chromatogram PNG with region support
├── hetero.py        # Heterozygous site detection and CSV reporting
└── cli.py           # CLI entry point with 5 subcommands
```

**Signal pipeline (BAM mode):** `extract_base_frequencies()` → `generate_traces()` → `write_ab1()` + `plot_chromatogram()`

**Signal pipeline (from-seq mode):** `read_sequence()` → `sequence_to_base_frequencies()` → `generate_traces()` → `write_ab1()` + `plot_chromatogram()`

**Plot pipeline (plot-ab1 mode):** `read_ab1()` → `plot_chromatogram(base_range=...)`

Key parameters: SPACING=10, SIGMA=2.0, SCALE=1024, NOISE=0, PHASING=0.0.

## AB1 Format Notes

ABIF binary channel order is **GATC** (DATA9=G, DATA10=A, DATA11=T, DATA12=C). Traces are big-endian int16, peak positions are big-endian uint16. The 128-byte header contains the root directory entry pointing to tag directory entries at end of file. Data blocks are 4-byte aligned.

## Demo Data

`demo_data/input/` contains 8 barcode samples from epi2me wf-amplicon with structure: `Analysis/barcodeXX/alignments/*.bam` + `Analysis/barcodeXX/consensus/consensus.fastq`, plus standalone FASTA/FASTQ in `sequences/`.

`demo_data/output/` contains test output from all 5 modes.

## Development Workflow

Follow these 6 steps for every new feature or release:

### Step 1 — Implement
- Add new functions/modules to the `ab1tools/` package (one file per concern)
- Add CLI subcommands in `cli.py` with clear `--help` descriptions
- Keep imports and function signatures clean

### Step 2 — Test All Modes
- Run comprehensive tests covering **ALL** modes, not just the new feature
- Test multiple scenarios per feature (different inputs, edge cases, parameter combos)
- Validate AB1 files (check ABIF header magic bytes + version)
- Count all outputs (AB1 files, PNGs, CSV rows) for the report

### Step 3 — Organize Files
- Test input data: `demo_data/input/` (Analysis/, raw/, sequences/)
- Test output: `demo_data/output/` with subdirectories per mode
- Test output is committed to the repo so reviewers can see results
- Test reports: `demo_data/test_report_v{X.Y.Z}.md`

### Step 4 — Write Test Report
- File name: `demo_data/test_report_v{X.Y.Z}.md`
- Include: summary table, test environment tree, detailed per-test sections with exact commands, AB1 validation results, output size summary, version changelog, conclusion
- Use numbered sections throughout
- Keep historical reports (don't delete old ones)

### Step 5 — Update Documentation
- `README.md`: numbered sections (## 1., ## 2., ### 3.1, ### 3.2, etc.)
- `CLAUDE.md`: update architecture section when adding new modules
- Bump version in both `ab1tools/__init__.py` and `pyproject.toml`

### Step 6 — Release
- Commit: `"Release vX.Y.Z: {summary}"` with feature list
- Push to **both** remotes: `origin` (GitHub) and `gitlab`
- Create GitHub release: `gh release create vX.Y.Z`
- Create GitLab release: via REST API (curl)
- Release notes must cite the test report file (e.g., "See `demo_data/test_report_v1.2.0.md`")
- Include all subcommands summary in release notes

## Technical Reference

`docs/ideas.md` contains the original design discussion with code examples (Chinese text, English code).
