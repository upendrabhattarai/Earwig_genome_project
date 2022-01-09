# Processing mRNA-seq data
We processed our mRNA-seq data through quality control and denovo assembly before using them in the maker pipeline for genome annotation.
We used Trimmomatic for filtering out adapters and low quality reads and Trinity for denovo assembly. Below are the scripts we used
## Quality control
`Script: Trimmomatic`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --job-name trimmRNA.EW
#SBATCH --mem=50G
#SBATCH --time=5:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread


mkdir trim2
mkdir unpaired2
module load Trimmomatic/0.38-Java-1.8.0_144
for f in $(<names)                             # We prepared a text files: "names" with the list of all the fastq files to use in this loop
do
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar PE -phred33 -threads 10 "${f}_R1_001.fastq.gz" "${f}_R2_001.fastq.gz" \
"trim2/${f}_R1_001_trim.fastq.gz" "unpaired2/${f}_R1_001_trim_unpaired.fastq.gz" \
"trim2/${f}_R2_001_trim.fastq.gz" "unpaired2/${f}_R2_001_trim_unparied.fastq.gz" \
ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 HEADCROP:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:35
done
```
We used `FastQC` before and after Trimmomatic to visualize the data quality and make sure there are no adapters and low quality reads left after Trimmomatic.


## Denovo assembly

`Script: Trinity`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --job-name=tri.mRNAEW
#SBATCH --account=uoo02752   # your NeSI project code
#SBATCH --time=30:00:00       # maximum run time
#SBATCH --ntasks=1            # always 1
#SBATCH --cpus-per-task=18    # number of threads to use for Trinity
#SBATCH --mem=300G            # maximum memory available to Trinity
#SBATCH --hint=nomultithread  # disable hyper-threading
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz

module load Trinity/2.11.0-gimkl-2020a

Trinity \
--CPU ${SLURM_CPUS_PER_TASK} --max_memory 300G \
--seqType fq \
--left "comma separated list of all *_R1_001_trim.fastq.gz from trim folder after Trimmomatic above" \
--right "comma separated list of all *_R2_001_trim.fastq.gz from trim folder after Trimmomatic above" \
--min_contig_length 200 \
--normalize_max_read_cov 200 \
--SS_lib_type RF \
--output trinity_out
```
