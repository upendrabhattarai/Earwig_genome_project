## Longranger pipeline
We used longranger to format the assembly and align the linked reads to it to get the aligned bam file. which was used as input for ARBitR

`Script to run longranger`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --job-name longR_EW_align
#SBATCH --mem=50G
#SBATCH --time=72:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH=/path/to/installation/folder/of/ongranger-2.2.2:$PATH

# First, we have to format the assembly before alignment, this will create a folder called [refdata-assembly] with assembly information
longranger mkref path/to/assembly/file/to/align/sequences/assembly.fasta                  

# Second, aligning the linked reads to the formatted assembly
longranger align --id=Earwig \
                  --fastqs=path/to/the/folder/containing/all/the/linked/reads \
                  --reference=path/to/the/refdata-assembly/folder/created/above
```
It will create a folder `Earwig` as the name of `--id` we provided in the script. Final output is the bam file in: `./Earwig/outs/possorted_bam.bam`
