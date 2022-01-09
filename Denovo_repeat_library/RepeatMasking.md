# Repeatmasking
we used the prepared repeat library to repeat mask the genome.


`script used for repeatmasker`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 30
#SBATCH --job-name rep.M.ew
#SBATCH --mem=30G
#SBATCH --time=70:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread



RepeatMasker -pa $SLURM_NTASKS -lib EW_Rep_CombinedLibrary.lib.minlen50.nr.classified.filtered.fa \
EW_assembly.fasta --gff -dir ./
```
