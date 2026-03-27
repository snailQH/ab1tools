# Convert Formats

ab1tools provides bidirectional format conversion between AB1 and common bioinformatics formats.

---

## `convert` — AB1 to FASTA, FASTQ, JSON, or CSV

```bash
ab1tools convert sample.ab1 -f fasta -o sample.fasta
ab1tools convert sample.ab1 -f fastq -o sample.fastq
ab1tools convert sample.ab1 -f json -o sample.json
ab1tools convert sample.ab1 -f csv -o sample.csv
```

### Output formats

#### FASTA
Standard sequence output, 80-character line wrapping.
```
>sample_name
ACTCAGCGTTGGAACTGGGTGGAGGGGGACTGGGATGGCTCCGTCCTTACAGGCCCAGTCC...
```

**Use case:** Extract the base-called sequence from AB1 for BLAST, alignment, or assembly input.

#### FASTQ
Sequence with quality scores estimated from trace peak heights.
```
@sample_name
ACTCAGCGTTGGAACTGGGTGGAG...
+
<<<<<<<<<<<<<<<<<<<<<<<<<...
```

**Use case:** Feed AB1-derived sequences into quality-aware pipelines (variant callers, trimmers) that expect FASTQ input.

#### JSON
Full structured data including traces, peaks, and metadata — machine-readable for custom analysis scripts.
```json
{
  "sample_name": "barcode01",
  "sequence": "ACTCAGCGTTGG...",
  "length": 3788,
  "peak_positions": [0, 10, 20, ...],
  "traces": {
    "A": [0, 0, 2, 5, 20, 100, 512, 1024, 512, 100, ...],
    "C": [...],
    "G": [...],
    "T": [...]
  }
}
```

**Use case:** Import AB1 data into Python/R/JavaScript for custom visualization or analysis. Build web-based chromatogram viewers. Feed into machine learning pipelines.

#### CSV
Per-position table with base calls and 4-channel peak intensities.
```csv
position,base,peak_pos,A,C,G,T
1,A,0,1024,2,0,0
2,C,10,2,1024,0,2
3,T,20,0,2,0,1024
```

**Use case:** Open in Excel for manual review. Plot signal profiles in R/ggplot2. Compute custom statistics. Feed into statistical analysis workflows.

### Direction of conversion

```
FASTA/FASTQ → AB1:  use ab1tools from-seq
BAM + ref   → AB1:  use ab1tools single/batch
AB1 → FASTA:        use ab1tools convert -f fasta
AB1 → FASTQ:        use ab1tools convert -f fastq
AB1 → JSON:         use ab1tools convert -f json
AB1 → CSV:          use ab1tools convert -f csv
AB1 → PNG:          use ab1tools plot-ab1
```

This provides the **complete bidirectional bridge** between AB1 and all common bioinformatics formats.
