# SARS-CoV-2 Variant Discovery & Lineage Assignment

> Adapted from the [Galaxy Training Network SARS-CoV-2 variant discovery tutorial](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/sars-cov-2-variant-discovery/tutorial.html), converted from the graphical Galaxy workflow to an equivalent command-line pipeline.

Analysis of ARTIC-protocol paired-end amplicon sequencing data from a
SARS-CoV-2 positive patient sample, to call intra-host variants, build a
consensus genome, and assign a PANGO lineage.

Sample used: **SRR17054503** (one sample from a 16-sample Omicron-era batch;
this run processes it individually as a first pass before scaling to the
full batch).

---

## Pipeline

```
FASTQ (paired-end)
  → fastp (QC, default settings)
  → BWA-MEM (mapping to NC_045512.2 reference)
  → ivar trim (remove ARTIC v4 primer sequences)
  → lofreq call (sensitive variant calling)
  → bcftools reheader (fix missing ##contig header — see notes)
  → SnpEff (functional annotation, NC_045512.2 database)
  → ivar consensus (build full consensus genome)
  → Pangolin (PANGO lineage assignment)
```

---

## 1. Reference Preparation

```bash
wget "https://www.ebi.ac.uk/ena/browser/api/fasta/MN908947.3?download=true" -O ref/sars_cov2_ref.fasta

# Rename header to the NCBI RefSeq accession SnpEff's database expects
sed -i '1s/.*/>NC_045512.2/' ref/sars_cov2_ref.fasta

samtools faidx ref/sars_cov2_ref.fasta
bwa index ref/sars_cov2_ref.fasta

wget https://zenodo.org/record/5888324/files/ARTIC_nCoV-2019_v4.bed -O ref/artic_v4_primers.bed
```

## 2. QC

```bash
fastqc fastq/SRR17054503_R1.fastq.gz fastq/SRR17054503_R2.fastq.gz -o fastqc_before/
```

No adapter contamination detected. The flagged "Per base sequence content"
and "Sequence Duplication Levels" warnings are expected and normal for
amplicon data (primer sequences create a fixed start, and amplicon
duplication is inherent to the protocol) — not a quality problem.

## 3. Trimming (fastp, default settings)

```bash
fastp -i fastq/SRR17054503_R1.fastq.gz -I fastq/SRR17054503_R2.fastq.gz \
      -o fastq/SRR17054503_R1.trimmed.fastq.gz -O fastq/SRR17054503_R2.trimmed.fastq.gz \
      --thread 4 --html fastp_report.html --json fastp_report.json
```

## 4. Mapping

```bash
bwa mem -t 4 ref/sars_cov2_ref.fasta \
    fastq/SRR17054503_R1.trimmed.fastq.gz fastq/SRR17054503_R2.trimmed.fastq.gz \
    > BAMs/SRR17054503.sam

samtools view -b BAMs/SRR17054503.sam > BAMs/SRR17054503.bam
samtools sort -@ 4 -o BAMs/SRR17054503.sorted.bam BAMs/SRR17054503.bam
samtools index BAMs/SRR17054503.sorted.bam
```

Result: **99.34% mapped**, 99.03% properly paired — expected to be very
high given the small, concentrated viral genome.

## 5. Primer Trimming (ivar)

```bash
ivar trim -i BAMs/SRR17054503.sorted.bam \
    -b ref/artic_v4_primers.bed \
    -p BAMs/SRR17054503.primertrimmed \
    -q 15 -m 30 -s 4

samtools sort -@ 4 -o BAMs/SRR17054503.primertrimmed.sorted.bam BAMs/SRR17054503.primertrimmed.bam
samtools index BAMs/SRR17054503.primertrimmed.sorted.bam
```

## 6. Variant Calling (lofreq)

```bash
lofreq call -f ref/sars_cov2_ref.fasta \
    -o vcf/SRR17054503.vcf \
    BAMs/SRR17054503.primertrimmed.sorted.bam
```

### Note: missing `##contig` header

`lofreq call` does not always write a `##contig` line to the VCF header,
which then breaks downstream tools (SnpEff, bcftools) that need to confirm
the chromosome name/length. Fixed with:

```bash
bcftools reheader -f ref/sars_cov2_ref.fasta.fai -o vcf/SRR17054503.fixed.vcf vcf/SRR17054503.vcf
```

## 7. Annotation (SnpEff)

```bash
snpEff -v NC_045512.2 vcf/SRR17054503.fixed.vcf > vcf/SRR17054503.annotated.vcf
```

## 8. Summary Table

```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/ANN\n' vcf/SRR17054503.annotated.vcf | \
awk -F'\t' '{split($6,a,"|"); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"a[4]"\t"a[2]"\t"a[3]}' > vcf/SRR17054503_summary.txt
```

Numerous missense variants were found across the **S (Spike)** gene — the
hallmark of the Omicron lineage's heavy mutation load in the region
responsible for immune evasion and transmissibility. Allele frequencies
ranged from near-fixed (AF ≈ 0.99, present in nearly all viral copies) to
very low (AF ≈ 0.01–0.1, minor variant populations), demonstrating
lofreq's sensitivity to low-frequency intra-host variation.

## 9. Consensus Genome

```bash
samtools mpileup -aa -A -d 0 -Q 0 BAMs/SRR17054503.primertrimmed.sorted.bam | \
ivar consensus -p vcf/SRR17054503_consensus -q 20 -t 0.5 -m 10
```

Genome quality:
```
Total length: 29,873 bp
N count:      801 bp (2.68%)
Resolved:     97.32%
```

A low N-fraction (well under Pangolin's 30% QC threshold) indicates high
sequencing depth and coverage across nearly the entire genome.

## 10. Lineage Assignment (Pangolin)

Pangolin required a dedicated conda environment due to an `iqtree` version
incompatibility inside its internal snakemake pipeline when installed
alongside other tools:

```bash
mamba create -n pangolin_env -c bioconda -c conda-forge pangolin
conda activate pangolin_env
pangolin vcf/SRR17054503_consensus.fa --outfile vcf/SRR17054503_lineage.csv
```

---

## Result

```
lineage:        BA.1
scorpio_call:   Omicron (BA.1-like)  (support: 0.88)
qc_status:      pass
ambiguity_score: 0.04
Alt alleles:    51
Usher placements: BA.1 (1/1)
```

The sample was confidently classified as **Omicron BA.1**, confirmed
independently by two internal classification systems (pangoLEARN/UShER
placement and the scorpio constellation caller), consistent with the batch
being drawn from an early Omicron surveillance dataset.

---

## Files

- `vcf/SRR17054503_summary.txt` — annotated variant table
- `vcf/SRR17054503_consensus.fa` — full consensus genome
- `vcf/SRR17054503_lineage.csv` — Pangolin lineage assignment report
- `fastp_report.html` — QC/trimming report

---

## Next step

This run processes a single sample as a first pass. The tutorial's full
value comes from batch-processing all 16 samples together and comparing
lineages/variant patterns across the cohort — a natural extension of this
pipeline using the same steps in a loop or script per sample.
