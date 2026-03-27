# Generate AB1 Files

ab1tools can generate realistic Sanger-like AB1 chromatogram files from three types of input. This is the toolkit's most distinctive capability — **no other published tool converts NGS data to AB1 format**.

---

## `from-seq` — Generate AB1 from FASTA or FASTQ

The simplest mode: convert any sequence into an AB1 file with clean, single-peak chromatogram traces.

```bash
ab1tools from-seq sequence.fasta -o output/
ab1tools from-seq consensus.fastq -o output/ --name my_sample
ab1tools from-seq assembly.fasta -o output/ --noise 5 --phasing 0.1 --decay 0.5
```

### When to use `from-seq`

**1. De novo assembly results (no reference)**

After assembling reads with tools like Flye, Canu, Hifiasm, or SPAdes, you have a consensus FASTA but no BAM alignment. Use `from-seq` to convert the assembly into AB1 for visual inspection in SnapGene or Chromas.

```bash
# Nanopore de novo assembly with Flye
flye --nano-hq reads.fastq --out-dir assembly/
ab1tools from-seq assembly/assembly.fasta -o ab1_output/ --name my_assembly

# PacBio HiFi assembly with Hifiasm
hifiasm -o asm reads.hifi.fq.gz
ab1tools from-seq asm.bp.p_ctg.fa -o ab1_output/

# Illumina SPAdes assembly
spades.py -1 R1.fq.gz -2 R2.fq.gz -o spades_out/
ab1tools from-seq spades_out/contigs.fasta -o ab1_output/
```

**2. Plasmid verification**

After whole-plasmid sequencing (e.g., via Plasmidsaurus or in-house Nanopore), convert the consensus to AB1 for inspection in SnapGene alongside the plasmid map.

```bash
ab1tools from-seq plasmid_consensus.fasta -o ab1_output/ --name pUC19_verified
```

**3. Synthetic gene sequences**

For ordered synthetic genes or gene fragments, create AB1 files for archival or integration into existing LIMS workflows that expect Sanger data.

```bash
ab1tools from-seq synthetic_gene.fasta -o ab1_output/ --name GFP_optimized
```

**4. Reference sequences for comparison**

Generate a "clean" AB1 from a reference sequence to use as a control in `ab1tools compare` or `ab1tools call --control`.

```bash
ab1tools from-seq reference.fasta -o control/ --name wild_type_control
```

**5. Teaching and simulation**

Create example chromatograms for courses, demonstrations, or algorithm benchmarking.

```bash
# Clean trace (ideal)
ab1tools from-seq example.fasta -o teaching/ --name clean_example

# Realistic trace with noise and decay
ab1tools from-seq example.fasta -o teaching/ --name realistic_example \
    --noise 5 --phasing 0.1 --decay 0.5
```

### Output

- `<name>.ab1` — AB1 chromatogram file
- `<name>.png` — Full-sequence chromatogram PNG (unless `--no-plot`)

---

## `single` — Convert a single BAM + consensus sample

Generates AB1 from aligned sequencing reads. This mode extracts **per-position base frequencies** from BAM pileup data, creating chromatogram traces where **mixed positions show proportional multi-channel peaks** — just like real Sanger data at heterozygous or variant sites.

```bash
ab1tools single --bam aligned.bam --consensus consensus.fasta -o output/
```

### When to use `single`

**1. Reference-based alignment (reads mapped to a reference genome)**

After aligning NGS reads to a reference genome with minimap2, bwa, or bowtie2, use `single` to generate AB1. Variant sites will show mixed peaks proportional to their allele frequency.

```bash
# Nanopore reads aligned to reference
minimap2 -a reference.fa reads.fastq | samtools sort -o aligned.bam
samtools index aligned.bam
# Extract consensus
samtools consensus aligned.bam -o consensus.fasta
# Generate AB1
ab1tools single --bam aligned.bam --consensus consensus.fasta -o output/
```

```bash
# Illumina reads aligned to reference
bwa mem reference.fa R1.fq.gz R2.fq.gz | samtools sort -o aligned.bam
samtools index aligned.bam
samtools consensus aligned.bam -o consensus.fasta
ab1tools single --bam aligned.bam --consensus consensus.fasta -o output/
```

**2. De novo assembly + read mapping back to consensus**

After assembling, map the original reads back to the assembly to capture variant frequencies. This is the recommended workflow for amplicon sequencing.

```bash
# Step 1: Assemble
flye --nano-hq reads.fastq --out-dir assembly/

# Step 2: Map reads back to assembly
minimap2 -a assembly/assembly.fasta reads.fastq | samtools sort -o aligned.bam
samtools index aligned.bam

# Step 3: Generate AB1 with variant information
ab1tools single --bam aligned.bam --consensus assembly/assembly.fasta -o output/ \
    --noise 5 --phasing 0.1 --decay 0.5
```

This workflow preserves allele frequency information in the chromatogram, so heterozygous positions show proportional mixed peaks.

**3. Amplicon sequencing (epi2me wf-amplicon output)**

The epi2me wf-amplicon workflow produces BAM + consensus FASTQ per barcode:

```bash
ab1tools single \
    --bam Analysis/barcode01/alignments/barcode01.aligned.sorted.bam \
    --consensus Analysis/barcode01/consensus/consensus.fastq \
    -o output/ --name barcode01
```

**4. CRISPR editing analysis**

After amplicon sequencing of CRISPR-edited samples, generate AB1 to visualize editing efficiency as mixed chromatogram peaks.

```bash
# Align edited amplicon reads to the target reference
minimap2 -a target_gene.fa edited_reads.fastq | samtools sort -o edited.bam
samtools index edited.bam

# Generate AB1 — editing sites will show multi-peak patterns
ab1tools single --bam edited.bam --consensus target_gene.fa -o crispr_output/ \
    --name KRAS_edited
```

**5. Tumor/somatic variant visualization**

From targeted panel or amplicon sequencing of tumor samples:

```bash
# Illumina amplicon panel aligned to reference
ab1tools single --bam tumor_panel.bam --consensus target_regions.fa -o output/ \
    --name tumor_sample_01
```

Variant allele frequencies (e.g., 15% KRAS G12D) will appear as proportional secondary peaks in the chromatogram.

### Signal parameters for realism

```bash
# Clean (default) — sharp, uniform peaks
ab1tools single --bam aligned.bam --consensus consensus.fa -o output/

# Realistic Sanger-like — recommended for publication figures
ab1tools single --bam aligned.bam --consensus consensus.fa -o output/ \
    --noise 5 --phasing 0.1 --decay 0.5

# Strong decay (simulating old sequencer or long read degradation)
ab1tools single --bam aligned.bam --consensus consensus.fa -o output/ \
    --decay 2.0 --noise 10
```

### Output

- `<name>.ab1` — AB1 with pileup-derived mixed peaks
- `<name>.png` — Full chromatogram PNG (unless `--no-plot`)

---

## `batch` — Batch convert all barcode samples

Auto-detects and processes all barcode directories from an epi2me wf-amplicon analysis in one command.

```bash
ab1tools batch /path/to/Analysis/ -o output/
ab1tools batch /path/to/Analysis/ -o output/ --noise 5 --phasing 0.1 --decay 0.5
ab1tools batch /path/to/Analysis/ -o output/ --no-plot  # AB1 only, faster
```

### Expected input structure

```
Analysis/
├── barcode01/
│   ├── alignments/barcode01.aligned.sorted.bam  (+ .bai index)
│   └── consensus/consensus.fastq
├── barcode02/
│   ├── alignments/barcode02.aligned.sorted.bam
│   └── consensus/consensus.fastq
└── ...
```

### When to use `batch`

- Processing an entire Nanopore amplicon sequencing run
- Generating AB1 files for all samples in a multiplexed experiment
- Creating a library of chromatograms for archival or downstream analysis

### Output

One AB1 file (+ optional PNG) per barcode sample in the output directory.

---

## Key Concept: How Base Frequencies Become Chromatogram Traces

The signal generation pipeline works the same for all modes:

```
Input (BAM pileup or sequence)
    ↓
Per-position base frequencies: {A: 0.7, C: 0.0, G: 0.3, T: 0.0}
    ↓
4-channel Gaussian peak modeling:
    S_b(x) = Σ f_i(b) · G(x; μ_i, σ) + ε(x)
    ↓
Optional effects: phasing, noise, signal decay
    ↓
AB1 binary file (ABIF format)
```

- **Pure positions** (e.g., A=100%): single sharp peak in A channel
- **Heterozygous positions** (e.g., A=70%, G=30%): two overlapping peaks with proportional heights
- **Multi-allelic positions** (e.g., A=50%, G=30%, T=20%): three peaks visible

This is how real Sanger chromatograms work — ab1tools faithfully reproduces this behavior from NGS data.
