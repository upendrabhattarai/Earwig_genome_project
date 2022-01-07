# 1. Repeat library preparation
## i. RepeatModuler
Repeatmodeler uses RECON, RepeatScout and Tandem Repeats Finder (TRF) to build a repeat library for an assembly.
RECON is good in finding less conserved elements in the genome while RepeatScout finds highly conserved repetitive elements [See more](https://www.pnas.org/content/117/17/9451).

Repeatmoduler gets confused with unusual characters in the sequence header so we renamed them using following oneliner.
```
awk '/^>/{print ">C" ++i; next}{print}' EW_assembly.fasta > EW_assembly_renamed.fa #This renames sequence header as C1, C2, C3 and so on.
```

`Script for RepeatModuler`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --job-name rep.ew
#SBATCH --mem=30G
#SBATCH --time=7-00:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module purge
module load RepeatModeler/2.0.2-Miniconda3

# First, we need to build the database
BuildDatabase -name earwig -engine ncbi ref.fa

# Second, run repeatmoduler
RepeatModeler -engine ncbi -pa 10 -database earwig
```
This will produce consensi.fa.classified, a classified repeat library

