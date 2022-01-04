#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --job-name quast
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load QUAST/5.0.2-gimkl-2018b

quast.py -t 10 --eukaryote --large --conserved-genes-finding \
genome.assembly.fasta \
-o quast
