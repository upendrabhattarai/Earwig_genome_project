## Processing with RNA-seq Reads

## 1. Rascaf
We have two paired-end RNA-seq data sets, mRNA-seq and total RNA-seq. Both data were already quality filtered using `Trimmomatic`.
So we will leverage both of them to scaffold the assembly using Rascaf.
To do so, we mapped each RNA-seq data to the assembly (output from the Arks pipeline) using `Hisat2` 

`Script to run hisat2`
```
#!/bin/bash -e

#SBATCH --job-name=hisat.mRNA-seq
#SBATCH --account=uoo02752
#SBATCH --nodes 1 
#SBATCH --cpus-per-task 1 
#SBATCH --ntasks 10
#SBATCH --mem=40G
#SBATCH --time=20:00:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-type=All
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load HISAT2/2.2.1-gimkl-2020a
module load Python/3.9.5-gimkl-2020a
module load SAMtools/1.13-GCC-9.2.0

# First, index the assembly
hisat2-build -p 10 path/to/assembly/output/from/Arks/step/assembly.fasta \
                    Hisat2.RNA.EW

# Second, run the alignment
hisat2 -x Hisat2.RNA.EW -1 path/to/First/pair/of/mRNA-seq/data/mRNA_merged_EW_R1_001_trim.fastq \
                        -2 path/to/Sirst/pair/of/mRNA-seq/data//mRNA_merged_EW_R2_001_trim.fastq \
                        -S path/to/output/EW_mRNA_alignment.sam
                        
samtools view -bS path/to/output/EW_mRNA_alignment.sam > path/to/output/EW_mRNA_alignment.bam
samtools sort path/to/output/EW_mRNA_alignment.bam -o path/to/output/EW_mRNA_alignment_sorted.bam
```

`Script to run Rascaf`
```
#!/bin/bash -e

#SBATCH --job-name=Rascaf.EW.mRNA
#SBATCH --account=uoo02752
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --mem=5G
#SBATCH --time=04:00:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-type=All
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH="/path/to/installation/directory/of/rascaf:$PATH"

rascaf -b path/to/sorted/bam/file/from/above/step/EW_mRNA_alignment_sorted.bam \
        -f path/to/assembly/output/from/Arks/step/assembly.fasta \                # This is the assembly to scaffold, output from Arks pipeline.
        -o EW_mRNA_scaffold
```

Hisat2 and Rascaf step for total RNA-seq data were run as above. That gave us two main output files
One from mRNA-scaffolding is: `EW_mRNA_scaffold.out`
and the other from tRNA-scaffolding is: `EW_tRNA_scaffold.out`

Finally, Rascaf combined these two outputs to get a final scaffolded assembly.
```
#!/bin/bash -e

#SBATCH --job-name=Rascaf.merge
#SBATCH --account=uoo02752
#SBATCH --nodes 1 
#SBATCH --cpus-per-task 1 
#SBATCH --ntasks 10
#SBATCH --mem=5G
#SBATCH --time=02:00:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-type=All
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH="/path/to/installation/directory/of/rascaf:$PATH"

rascaf-join -r EW_tRNA_scaffold.out -r EW_mRNA_scaffold.out -o EW_RNA_scaffold
```
This gave us the scaffolded assembly in fasta format `EW_RNA_scaffold.fa`
