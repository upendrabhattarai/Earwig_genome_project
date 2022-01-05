## 1. RagTag
Before processing with raw linked-read. We used the Supernova produced assembly to scaffold the assembly. We used [RagTag](https://github.com/malonge/RagTag) for that purpurse.

`Script to run RagTag`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --job-name Ragtag.EW
#SBATCH --mem=10G
#SBATCH --time=05:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH="/path/to/the/installation/directory/of/RagTag:$PATH"  # Installation path, conda installation available

ragtag.py scaffold \
                path/to/the/supernova/assembly/EW_10xSN.fasta \
                path/to/the/gapclosed/assembly/after/LR_Gapcloser/assembly.fasta
 ```
 ## 2.ARBitR
 We used raw linked read data to scaffold the assembly with [ARBitR](https://github.com/markhilt/ARBitR)
 Before that, we had to process the assembly file and map linked-reads to the processed assembly using 'Longranger' [Scripts here](longranger.sh). 
 The resulting `Bam` file was used as input for ARBitR.
 
 `Script to run ARBitR`
 ```
 #!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --job-name Arbitr.ew
#SBATCH --mem=5G
#SBATCH --time=05:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH="/path/to/installation/directory/of/ARBitR/src:$PATH"
export PATH="/path/to/conda/environment/where/dependencies/are/installed:$PATH" # dependencies were installed in conda environment

arbitr.py -i ragtag.scaffold.fasta \                                    # This is the output from the RagTag above
           -o output.arbitr.scaffolds \
           path/to/bam/file/from/longranger/outs/possorted_bam.bam      # Bamfile produced by longranger-align
```

