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
 Before that, we had to process the assembly file and map linked-reads to the processed assembly using 'Longranger' [Scripts here](longranger-align.md). 
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
## 3. Arks & Links
We further used Arks pipeline to scaffold the assembly. We had to process the linked read with longranger-basics to prepare an input [Scripts here](longranger-basic.md)

`Script to run Arks pipeline`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 16
#SBATCH --job-name arks.EW
#SBATCH --mem=25G
#SBATCH --time=18:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load BWA/0.7.17-GCC-9.2.0
module load SAMtools/1.13-GCC-9.2.0
module load BEDTools/2.29.2-GCC-9.2.0
module load LINKS/1.8.7-GCC-9.2.0

export PATH="/Path/to/installation/directory/arks-1.0.4/Arks:$PATH"
export PATH="/Path/to/installation/directory/arks-1.0.4/Examples:$PATH"

arks-make arks draft=assembly \    # assembly file to scaffold should be in working directory (or symlink) should have the name supplied to draft with .fasta extension
                reads=barcoded \    # barcoded.fastq.gz should be in the working directory (or symlink)
                threads=16
```

