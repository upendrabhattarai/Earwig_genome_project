# Earwig_genome_project
This repository contains all the scripts used to assemble and annotate the Earwig genome.
The pipeline is presented in three parts:

1. Genome assembly
2. Denovo repeat library
3. Genome annotation
---
1. Genome assembly
---
Genome is assembled using linked reads from 10x chromium and long reads from Oxford nanopore.
Long and linked reads were individually assembled and then merged together. After multiple iterations of scaffolding, gapclosing, and haplotigs and contaminants removal, assembly was polished with mRNA-seq reads to obtain final assembly. Schematic representation in figure below (Created with BioRender.com).

