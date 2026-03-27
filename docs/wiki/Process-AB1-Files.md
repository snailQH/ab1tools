# Process AB1 Files

ab1tools provides three processing operations for AB1 files: comparison, quality trimming, and region extraction.

---

## `compare` — Compare two AB1 files

Compute signal-level similarity and base call concordance between two AB1 files.

```bash
# Compare sample vs control
ab1tools compare sample.ab1 control.ab1

# Save report
ab1tools compare sample.ab1 replicate.ab1 -o comparison.txt

# JSON output for programmatic use
ab1tools compare sample.ab1 control.ab1 --format json -o report.json

# Overlay plot (4-channel dual-trace PNG)
ab1tools compare sample.ab1 control.ab1 --plot -o comparison.txt
```

### Output metrics

| Metric | Description |
|--------|-------------|
| Base concordance | Percentage of matching base calls at each position |
| Mismatch count | Number of positions with different base calls |
| Channel correlation | Pearson r for each of the 4 trace channels (A, C, G, T) |
| Mean correlation | Average of the 4 channel correlations |
| Mismatch positions | List of discordant positions with the bases from each file |

### Example output

```
Comparing: barcode01 vs bc01_fasta
  Compared bases: 3788

Base concordance: 99.95% (3786/3788)
Mismatches:       2

Signal correlation (Pearson r):
  A: 0.9994    C: 0.9995    G: 0.9997    T: 0.9997
  Mean: 0.9996

Mismatch positions:
  pos 1916: A → C
  pos 3016: G → T
```

### Use cases

- **BAM vs pure sequence:** compare `single` output (with variants) vs `from-seq` output (pure) to find variant positions
- **Replicate consistency:** compare technical replicates to assess reproducibility
- **Cross-platform validation:** compare AB1 generated from Nanopore vs Illumina BAM for the same sample
- **Pre/post editing:** compare CRISPR-edited sample AB1 with wild-type AB1
- **Real vs simulated:** compare a real Sanger AB1 with ab1tools-generated AB1 for signal fidelity assessment
- **Analogous to:** `diff`, `bcftools isec`

---

## `trim` — Quality-based trace trimming

Trim low-quality regions from the start and end of an AB1 file using Mott's algorithm (Ewing et al., 1998), the same approach used in the phred base-caller.

```bash
ab1tools trim sample.ab1 -o trimmed.ab1
ab1tools trim sample.ab1 --quality-threshold 30
ab1tools trim sample.ab1 --quality-threshold 20 --min-length 100
```

### How Mott's algorithm works

1. Compute per-position quality scores from trace signal ratios
2. Calculate `adjusted[i] = Q[i] - threshold` for each position
3. Find the contiguous region that **maximizes** the sum of adjusted scores (Kadane's algorithm)
4. This region is the optimal high-quality segment

Positions with Q above the threshold contribute positively; positions below contribute negatively. The algorithm automatically finds the best balance.

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--quality-threshold` | 20 | Minimum Q score (positions below this are "bad") |
| `--min-length` | 50 | Minimum length of trimmed region (keeps full sequence if too short) |

### Output

A new AB1 file containing only the trimmed region, with rebased peak positions.

```
  Input:    sample.ab1 (3788 bp)
  Trimmed:  bases 50–3650 (3601 bp)
  Removed:  49 bp (left) + 138 bp (right)
  Mean Q:   28.3 (trimmed region)
  Output:   sample_trimmed.ab1
```

### Use cases

- **Pre-processing** before variant calling to remove unreliable end regions
- **Sanger data cleanup** for real AB1 files from sequencing facilities
- **Quality filtering** in batch processing pipelines
- **Analogous to:** `seqtk trimfq`, `cutadapt`, phred/phrap trimming

---

## `extract` — Extract sub-region into new AB1

Cut any base range from an AB1 file into a new, self-contained AB1 file.

```bash
# Extract bases 100-200
ab1tools extract sample.ab1 --start 100 --end 200 -o region.ab1

# Extract first 500 bases
ab1tools extract sample.ab1 --end 500

# Extract from position 3000 to end
ab1tools extract sample.ab1 --start 3000

# Auto-named output: sample_100-200.ab1
ab1tools extract sample.ab1 --start 100 --end 200
```

### What gets extracted

The new AB1 file is fully self-contained:
- **Trace data:** all 4 channels for the specified region (with padding)
- **Base calls:** subsequence of the original
- **Peak positions:** rebased to start from 0
- **Quality scores:** re-estimated from trace signals
- **Sample name:** appended with `_start-end`

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--start` / `-s` | 1 | Start position (1-based, inclusive) |
| `--end` / `-e` | end | End position (1-based, inclusive) |
| `--output` / `-o` | auto | Output path (default: `<input>_<start>-<end>.ab1`) |

### Use cases

- **Focus on variant region:** extract the area around a SNP or editing site for detailed analysis
- **Split long amplicons:** divide a 4kb amplicon into 500bp segments for easier visualization
- **Share specific regions:** send only the relevant portion of a chromatogram to collaborators
- **Feed into `ab1tools call`:** run variant detection on a specific region only
- **Analogous to:** `samtools view -b file.bam chr:start-end`, `bedtools getfasta`

### Example workflow: extract + call

```bash
# Extract region around known CRISPR target site
ab1tools extract edited_sample.ab1 --start 450 --end 550 -o target_region.ab1

# Run SDVC on extracted region
ab1tools call target_region.ab1 --threshold 0.03

# Visualize the region
ab1tools plot-ab1 target_region.ab1 --bases-per-row 50 --dpi 300
```
