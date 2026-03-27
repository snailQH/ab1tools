# Use Cases and Workflows

Real-world workflows showing how ab1tools fits into bioinformatics pipelines.

---

## 1. Nanopore Amplicon Sequencing (Primary Use Case)

### Full workflow: sequencing → analysis → chromatogram

```bash
# Step 1: Sequencing (Oxford Nanopore MinION/PromethION)
# Basecalling with Dorado or Guppy → FASTQ

# Step 2: Analysis with epi2me wf-amplicon
nextflow run epi2me-labs/wf-amplicon \
    --fastq reads/ --reference reference.fa --out_dir Analysis/

# Step 3: Generate AB1 chromatograms for all barcodes
ab1tools batch Analysis/ -o ab1_output/ --noise 5 --phasing 0.1 --decay 0.5

# Step 4: Smart mode — auto-detect variants
ab1tools smart Analysis/ -o smart_output/ --noise 5 --phasing 0.1 --decay 0.5

# Step 5: Inspect specific variant
ab1tools extract smart_output/barcode01.ab1 --start 3000 --end 3050 -o variant_region.ab1
ab1tools call variant_region.ab1 --threshold 0.03 --plot
```

---

## 2. Illumina Amplicon Panel (Clinical)

### Targeted gene panel → AB1 for variant review

```bash
# Step 1: Align reads to reference
bwa mem -t 8 reference.fa R1.fq.gz R2.fq.gz | samtools sort -o aligned.bam
samtools index aligned.bam

# Step 2: Generate consensus
samtools consensus aligned.bam -o consensus.fasta

# Step 3: Generate AB1 — variant sites show proportional mixed peaks
ab1tools single --bam aligned.bam --consensus consensus.fasta \
    -o output/ --name patient_001 --noise 5 --phasing 0.1

# Step 4: Detect variants with SDVC
ab1tools call output/patient_001.ab1 --threshold 0.03 --format vcf -o variants.vcf

# Step 5: Visual review in SnapGene
# Open output/patient_001.ab1 in SnapGene for clinical review
```

### Clinical validation: virtual Sanger confirmation

Many clinical labs require Sanger confirmation of NGS variants (AMP 2023 guidelines). ab1tools can generate a "virtual Sanger" trace from the same NGS data:

```bash
# Generate AB1 from the panel BAM
ab1tools single --bam clinical_panel.bam --consensus reference.fa \
    -o validation/ --name BRCA1_patient_001

# Extract the region around the variant
ab1tools extract validation/BRCA1_patient_001.ab1 --start 1450 --end 1550 \
    -o validation/BRCA1_variant_region.ab1

# Generate high-res chromatogram for report
ab1tools plot-ab1 validation/BRCA1_variant_region.ab1 \
    --bases-per-row 50 --dpi 300 -o validation/BRCA1_variant_figure.png
```

---

## 3. CRISPR Editing Efficiency Analysis

### Amplicon sequencing → editing efficiency visualization

```bash
# Step 1: Sequence the edited amplicon (any platform)
minimap2 -a target_gene.fa edited_reads.fastq | samtools sort -o edited.bam
samtools index edited.bam

# Step 2: Generate AB1 — editing sites show mixed peaks
ab1tools single --bam edited.bam --consensus target_gene.fa \
    -o crispr/ --name sgRNA_01_edited

# Step 3: Generate control AB1 from wild-type reference
ab1tools from-seq target_gene.fa -o crispr/ --name wild_type --no-plot

# Step 4: SDVC with control for maximum sensitivity
ab1tools call crispr/sgRNA_01_edited.ab1 \
    --control crispr/wild_type.ab1 --threshold 0.02 --plot

# Step 5: Compare edited vs wild-type
ab1tools compare crispr/sgRNA_01_edited.ab1 crispr/wild_type.ab1 \
    --plot -o crispr/editing_comparison.txt
```

---

## 4. PacBio HiFi Amplicon Analysis

```bash
# Step 1: Align HiFi reads to reference
pbmm2 align reference.fa hifi_reads.bam aligned.bam --sort

# Step 2: Generate consensus
samtools consensus aligned.bam -o consensus.fasta

# Step 3: Generate AB1
ab1tools single --bam aligned.bam --consensus consensus.fasta \
    -o output/ --name pacbio_sample
```

---

## 5. De Novo Assembly QC

### Verify assembly by mapping reads back and generating chromatograms

```bash
# Step 1: Assemble
flye --nano-hq reads.fastq --out-dir assembly/

# Step 2: Map reads back to assembly
minimap2 -a assembly/assembly.fasta reads.fastq | samtools sort -o realigned.bam
samtools index realigned.bam

# Step 3: Generate AB1 with variant info
ab1tools single --bam realigned.bam --consensus assembly/assembly.fasta \
    -o assembly_qc/ --name assembly_check

# Step 4: Check for assembly errors (unexpected variants)
ab1tools call assembly_qc/assembly_check.ab1 --threshold 0.05 -o assembly_qc/errors.csv

# Step 5: Quality statistics
ab1tools stats assembly_qc/assembly_check.ab1 --format csv -o assembly_qc/quality.csv
```

---

## 6. Plasmid Verification (Replacing Sanger)

### Nanopore whole-plasmid sequencing → AB1 for SnapGene

As services like Plasmidsaurus replace traditional Sanger sequencing for plasmid verification, researchers still need chromatogram visualization:

```bash
# Step 1: Receive plasmid consensus from Plasmidsaurus or in-house Nanopore
# (FASTA file of the assembled plasmid sequence)

# Step 2: Generate AB1 for SnapGene
ab1tools from-seq plasmid_consensus.fasta -o output/ --name pUC19_clone_A3

# Step 3: Open in SnapGene alongside the plasmid map
# SnapGene aligns the AB1 chromatogram to the map automatically
```

---

## 7. Teaching and Education

### Create controlled example chromatograms for courses

```bash
# Clean homozygous trace
echo ">example" > /tmp/homo.fa
echo "ACGTACGTACGTACGTACGT" >> /tmp/homo.fa
ab1tools from-seq /tmp/homo.fa -o teaching/ --name homozygous_example

# Simulate a heterozygous position (manually create mixed-freq BAM)
# Or use existing demo data:
ab1tools smart demo_data/input/Analysis/ -o teaching/hetero_examples/

# Generate traces at different noise levels for comparison
for noise in 0 5 10 20; do
    ab1tools from-seq /tmp/homo.fa -o teaching/ \
        --name "noise_${noise}" --noise $noise --phasing 0.1 --decay 0.5
done
```

---

## 8. Benchmarking Variant Callers

### Generate ground-truth AB1 files with known variant frequencies

```bash
# Step 1: Create controlled data with known variants
# (Use the simulation scripts in docs/research/phase4_benchmarks/)

# Step 2: Run ab1tools SDVC
ab1tools call simulated_20pct_variant.ab1 --threshold 0.03 -o sdvc_results.csv

# Step 3: Submit to TIDE for comparison
# Upload simulated AB1 to https://tide.nki.nl

# Step 4: Compare results
# ab1tools-generated AB1 files have known ground-truth frequencies,
# making them ideal for benchmarking any chromatogram analysis tool
```

---

## 9. Batch Processing Pipeline (Shell Script)

```bash
#!/bin/bash
# Process all AB1 files in a directory

INPUT_DIR="sanger_data/"
OUTPUT_DIR="results/"
mkdir -p "$OUTPUT_DIR"

for ab1 in "$INPUT_DIR"/*.ab1; do
    name=$(basename "$ab1" .ab1)
    echo "Processing: $name"

    # View metadata
    ab1tools view "$ab1" --json > "$OUTPUT_DIR/${name}_meta.json"

    # Quality stats
    ab1tools stats "$ab1" --format csv -o "$OUTPUT_DIR/${name}_stats.csv"

    # Variant calling
    ab1tools call "$ab1" --threshold 0.05 -o "$OUTPUT_DIR/${name}_variants.csv"

    # Convert to FASTA
    ab1tools convert "$ab1" -f fasta -o "$OUTPUT_DIR/${name}.fasta"
done

echo "Done. Results in $OUTPUT_DIR/"
```

---

## 10. Forensic / Mixed Sample Analysis

### Visualize contributor mixtures as chromatogram peaks

```bash
# Align STR/SNP amplicon reads to reference
minimap2 -a str_reference.fa mixture_reads.fastq | samtools sort -o mixture.bam
samtools index mixture.bam

# Generate AB1 — mixture ratios visible as peak height ratios
ab1tools single --bam mixture.bam --consensus str_reference.fa \
    -o forensic/ --name evidence_01

# Detect mixed positions
ab1tools call forensic/evidence_01.ab1 --threshold 0.10 --format csv
```
