# Braker pipeline

We used Braker to train Augustus with the RNA-seq reads and predictions from Genemark-ES.
RNA-seq reads were aligned with `STAR` and sorted and indexed with samtools

`Script for alignment with STAR`
```
#!/bin/bash -e

#SBATCH --nodes 1 
#SBATCH --cpus-per-task 1 
#SBATCH --ntasks=10
#SBATCH --job-name star.ew
#SBATCH --mem=50G 
#SBATCH --time=04:00:00 
#SBATCH --account=uoo02752 
#SBATCH --output=%x_%j.out 
#SBATCH --error=%x_%j.err 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=bhaup057@student.otago.ac.nz 
#SBATCH --hint=nomultithread

module load STAR/2.7.9a-GCC-9.2.0
module load SAMtools/1.13-GCC-9.2.0
All_R1="Comma separated list of R1 files"
All_R2="Comma searated list of R2 files"

# Step-1
STAR --runMode genomeGenerate --runThreadN 8 --genomeSAindexNbases 13 --genomeDir index --genomeFastaFiles path/to/assembly/EW_assembly.fasta

# Step-2
STAR --twopassMode Basic --genomeDir index --runThreadN 10 \
--readFilesIn $All_R1 $All_R2 \
--readFilesCommand zcat \
--outFileNamePrefix RNA_ --outSAMtype BAM Unsorted

# Step-3
samtools sort -T ./ -m 4G --threads 10 -o RNA.sorted.bam RNA_Aligned.out.bam
samtools index RNA.sorted.bam
```

We then ran `Braker` with input from the genemark and star aligner. 
Braker can also run genemark internally and use the result to train augustus. However, as we ran `Genemark-ES` separately, we will use that output in braker to speed up the process.
To train Augustus, we have augustus config file downloaded locally at `/nesi/nobackup/uoo02752/nematode/busco_downloads/MyAugustusConfig_Ew` we will provide this path for training.
We also need to provide path to the bin directory and script directory from augustus installation folder. We have augustus available as module in `Nesi` and that is where we got these paths.
Loading Augustus as module didn't take care of the path for bin and script directory for some reason. Braker also utilizes CBD_tools, which we installed using `conda`.
We will be using soft masked genome using the repeat library we prepared as an input for braker.


`Script to run braker`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 25
#SBATCH --job-name braker.rnaseq
#SBATCH --mem=150G
#SBATCH --time=1-00:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module purge
module load BRAKER/2.1.6-gimkl-2020a-Perl-5.30.1-Python-3.8.2
export AUGUSTUS_CONFIG_PATH=/nesi/nobackup/uoo02752/nematode/busco_downloads/MyAugustusConfig_Ew
export AUGUSTUS_BIN_PATH=/opt/nesi/CS400_centos7_bdw/AUGUSTUS/3.3.3-gimkl-2020a/bin
export AUGUSTUS_SCRIPTS_PATH=/nesi/nobackup/uoo02752/earwig/earwig_annotation/scripts
export CDBTOOLS_PATH=/nesi/nobackup/uoo02752/nematode/bin/miniconda3/bin

braker.pl --species=ForAuricularia --genome=path/to/soft-masked/assembly/EW_assembly.fasta.masked \
--geneMarkGtf=/path/to/output/from/Genemrk-ES_output/genemark.gtf \
--bam=path/to/STAR/alignment/of/rna-seq/reads/rna_seq_alignment/RNA.sorted.bam --softmasking --cores=25
```
This will train Augustus and the trained data will be in species folder inside the Augustus config path we provided in the script above.
Which we can use as an input to Maker. Braker will also produce annotation `braker.gtf` file. That can also be used in Maker.
We can convert `braker.gtf` to `gff3` file using [AGAT toolkit](https://github.com/NBISweden/AGAT)

