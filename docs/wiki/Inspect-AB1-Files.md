# Inspect AB1 Files

## `view` — Display AB1 metadata and trace statistics

Quickly inspect the contents and quality of any AB1 file without opening a GUI viewer.

```bash
ab1tools view sample.ab1
```

### Example output

```
File:            sample.ab1
Sample:          barcode01
Sequence length: 3788 bp
Trace length:    37880 points
Peak spacing:    10
GC content:      56.6%
Base composition: A=765 C=1125 G=1018 T=880

Channel statistics:
  Chan     Mean     Max     Min      Std
     A    103.1    1024       0    251.7
     C    152.3    1024       0    294.1
     G    137.8    1024       0    283.1
     T    118.8    1024       0    267.1
```

### JSON output

```bash
ab1tools view sample.ab1 --json
```

Returns a JSON object suitable for programmatic parsing:

```json
{
  "file": "sample.ab1",
  "sample_name": "barcode01",
  "sequence_length": 3788,
  "trace_length": 37880,
  "peak_spacing": 10,
  "base_composition": {"A": 765, "C": 1125, "G": 1018, "T": 880, "N": 0},
  "gc_content_%": 56.6,
  "channel_stats": {
    "A": {"mean": 103.1, "max": 1024, "min": 0, "std": 251.7},
    ...
  }
}
```

### Use cases

- **Quick QC:** check file integrity, sequence length, GC content before detailed analysis
- **Batch inspection:** script `ab1tools view --json` across hundreds of files and parse results
- **Debug:** verify that ab1tools-generated AB1 files have correct trace data
- **File inventory:** catalog a directory of AB1 files with their properties
- **Analogous to:** `samtools view -H` for BAM header inspection
