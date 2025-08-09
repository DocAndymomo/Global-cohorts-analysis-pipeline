# Global Cohorts Analysis Pipeline

## Introduction
The Global Cohorts Analysis Pipeline enables the systematic recovery and high-precision profiling of microbial communities from human genomic and transcriptomic sequencing data, with minimized contamination. This pipeline integrates two-step host sequence removal, robust species annotation, and multi-center contamination filtering, making it highly suitable for large-scale, cross-cohort microbiome studies. By combining deep host depletion, high-resolution taxonomic classification, and multi-cohort contamination filtering, the pipeline achieves lower contamination rates and greater analytical precision than conventional metagenomic workflows.

---

## Overview

### Main Features
- Two-step host removal using both GRCh38.p14 and CHM13 v2.0 reference genomes to maximize removal of host sequences.
- High-accuracy species annotation with Kraken2 + Bracken using the PlusPF database.
- Cross-cohort contamination removal using a multi-center contamination database and decontam filtering.
- Final microbiome abundance tables generated from high-confidence reads for downstream statistical analyses.

---

## Step-by-Step Pipeline

### Step 1 – Raw Data Conversion (`01-Raw2fastq`)
Convert raw sequencing output (e.g., BCL format) into FASTQ files, ensuring proper sample demultiplexing and metadata integrity.

### Step 2 – Quality Control (`02-fastp`)
Perform read trimming and quality filtering using fastp (v0.23.4) to remove low-quality bases, adapters, and sequencing artifacts.

### Step 3 – First Host Removal (GRCh38) (`03-bowtie2_grch`)
Align reads to the human reference genome GRCh38.p14 using bowtie2 (v2.5.2) to remove human-origin sequences.

### Step 4 – Second Host Removal (CHM13) (`04-bowtie2_chm`)
Align the remaining reads to the CHM13 v2.0 telomere-to-telomere (T2T) human reference using bowtie2 (v2.5.2) for deeper host sequence depletion.

### Step 5 – Taxonomic Classification & Quantitation (`05-kraken_beacken`)
Classify microbial reads using Kraken2 (v2.1.3) with the PlusPF (2024.1.12) database.  
Refine species-level abundance estimates using Bracken for more accurate microbial quantitation.

### Step 6 – Contamination Removal (`06-Decontam`)
Use a multi-center contamination database and decontam filtering to remove background contaminants, ensuring cross-cohort comparability.

### Step 7 – Microbiome Data Extraction (`07-kraken_tools`)
Extract high-confidence microbial reads and generate final microbiome abundance tables using KrakenTools, ready for downstream statistical and ecological analyses.

---

## Advantages Over Conventional Methods
- Reduced Contamination — Multi-step host removal and cohort-specific contaminant filtering.
- Improved Precision — Dual-host genome alignment combined with Bracken refinement for species-level accuracy.
- Scalable — Optimized for integration of datasets from multiple global cohorts.
- Reproducible — Modular structure with version-controlled tools and databases.
