# Earwig_genome_project
This repository contains all the scripts used to assemble and annotate the Earwig genome.
The pipeline is presented in three parts:

1. Genome assembly
2. Denovo repeat library
3. Genome annotation
---
## 1. Genome assembly

Genome is assembled using linked reads from 10x chromium and long reads from Oxford nanopore.
Long and linked reads were individually assembled and then merged together. After multiple iterations of scaffolding, gapclosing, and haplotigs and contaminants removal, assembly was polished with mRNA-seq reads to obtain final assembly. Schematic representation in figure below (Created with BioRender.com).

![Alt text](Earwig_assembly_pipeline.png?raw=true "Title")

Workflow:
1. [Linked read assembly](Linked_reads_only_assembly.md)
2. [Long read assembly](Long_read_assembly.md)
3. [Merging two assemblies](merging_assemblies.md)
4. Further processing with long reads
5. Further processing with linked reads
6. Processing with RNA-seq reads
---
## 2. Denovo repeat library
A comprehensive denovo repeat library is prepared for the assembled genome. It was used for repeat content analysis, repeat masking and as input for annotation pipeline.

Workflow:
1. De novo repeat identification
      1. Repeatmoduler
      2. LTRharvest & LTRdigest
      3. TransposonPSI
      4. Mite-tracker
2. External repeat database
      1. Sine database
      2. Dfam database
      3. Repbase database
3. Concatenating, filtering, and classifying repeats
      1. RepeatClassifier
4. Repeat masking the genome
      1. RepeatMasker

---

## 3. Genome annotation
Maker2 pipeline is used for genome annotation. mRNA-seq data is denovo assembled using Trinity. Other relavant publicly available datasets were downloaded and used as input.

1. Denovo mRNA-seq assembly
        1. Trimmomatic
        2. Trinity
2. Downloading publicly available datasets
3. Configuring and running Maker2 pipeline
