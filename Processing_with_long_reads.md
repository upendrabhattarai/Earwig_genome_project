## 1. Purgehaplotigs

We ran [Purgehaplotigs](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/) on the merged assembly to remove haplotigs.
Because of the low coverage of long reads, peaks overlapped so we setup the parameters for the coverage steps to only pick contigs to produce a haploid representation, but not merge or collapse sequences. as discussed in this [issue](https://bitbucket.org/mroachawri/purge_haplotigs/issues/69/low-coverage-dataset-unsure-of-parameters)

`Scripts for Purgehaplotigs`

```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --partition=large,bigmem
#SBATCH --job-name purgehap.ew
#SBATCH --mem=25G
#SBATCH --time=10:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load SAMtools/1.12-GCC-9.2.0
module load minimap2/2.20-GCC-9.2.0
module load BEDTools/2.29.2-GCC-9.2.0

# Commands below are run in sequences in the same working directory

# First, mapping longreads back to assembly
minimap2 -t 10 -ax map-ont 
          path/to/the/merged/assembly/EW_Flye_SN_merge.fasta \
          path/to/the/filtered/long-read/data/EW_nanopore_merged_filtered_porechop.fastq.gz \
          --secondary=no | samtools sort -m 5G -o aligned.bam -T tmp.ali
          
# Second, create histogram with purgehaplotigs
export PATH="/nesi/nobackup/uoo02752/.conda/envs/purge_haplotigs_env/bin:$PATH" # path to installation of purge_haplotigs software, we installed it in conda env.
purge_haplotigs hist -t 10 -b aligned.bam                                       # This bam file is produced in the step above
          -g path/to/the/merged/assembly/EW_Flye_SN_merge.fasta 

# Third, coverage step with custom parameters, still needs path to purgehaplotigs installation as in second step
purge_haplotigs cov -i aligned.bam.gencov -l 1 -m 11 -h 190 -o coverage_stats.csv

# Fourth, purging step
purge_haplotigs purge \
                -g path/to/the/merged/assembly/EW_Flye_SN_merge.fasta \
                -c coverage_stats.csv -b aligned.bam -d
```
From the Purgehaplotigs, we will get `curated.fasta` as an output.
---
## 2. Rails & Cobbler

We ran [Rails & Cobbler](https://github.com/bcgsc/RAILS) with a wrapper script `runRAILSminimapSTREAM.sh` using minimap as the aligner.
It first aligns the reads with the assembly and runs `cobbler` for gap-filling. Then the scaffolding is done with `RAILS`. 
For some reason, it requires all the executable scripts during installation in the working directory. 
So, I softlinked all those scripts in the working directory as:

```
ln -s path/to/installation/folder/* /path/to/Rails&Cobbler/working/directory
```

`Script used for Rails & Cobbler pipeline`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --job-name Rails&Cobbler.EW
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load Perl/5.30.1-GCC-9.2.0
module load minimap2
module load SAMtools/1.13-GCC-9.2.0

sh runRAILSminimapSTREAM.sh curated.fasta \                                                           # curated.fasta is the output from purgehaplotigs
                            path/to/the/filtered/long-read/data/EW_nanopore_merged_filtered_porechop.fastq.gz\
                            -d 250 -i 0.85 -g 500 -l 2 ont \                                      
                            /scale_wlg_persistent/filesets/opt_nesi/CS400_centos7_bdw/SAMtools/1.13-GCC-9.2.0/bin/samtools -t 10
```
---
## 3. LRScaf

We ran multiple iterations of [LRScaf](https://github.com/shingocat/lrscaf) in series to scaffold the genome. We checked the assembly quality after each round using [Quast](quast.sh). After five rounds there was no significant improvement in scaffolding, hence the fifth round scaffolded assembly was further processed.

`Script to run LRscaf`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --job-name lrscaf.EW
#SBATCH --mem=30G
#SBATCH --time=05:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load minimap2/2.20-GCC-9.2.0

minimap2 -t 10 \
          path/to/output/from/rails&cobbler/assembly.fasta \
          path/to/the/filtered/long-read/data/EW_nanopore_merged_filtered_porechop.fastq.gz > ./aln.mm

export PATH="/nesi/nobackup/uoo02752/bin/lrscaf/target/:$PATH"                  # Installation directory

java -Xms40g -Xmx40g -jar /nesi/nobackup/uoo02752/bin/lrscaf/target/LRScaf-1.1.11.jar \
          --contig path/to/output/from/rails&cobbler/assembly.fasta \
          --alignedFile aln.mm -t mm -p 10 --output path/to/output/directory/scaffolds1
```
---
## 4. LR_Gapcloser

We ran multiple iterations of [LR_Gapcloser](https://github.com/CAFS-bioinformatics/LR_Gapcloser) to close gaps in the assembly. 
LR_Gapcloser takes fasta file as input so had to convert longread fastq file in to fasta file. [Quast](quast.sh) was run on output of each iteration. Eighth was selected as there was not much improvement after that.

`Script to run LR_Gapcloser`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --job-name lr-gapEW
#SBATCH --mem=25G
#SBATCH --time=72:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load BWA/0.7.17-gimkl-2017a
export PATH=/nesi/nobackup/uoo02752/bin/LR_Gapcloser/src/:$PATH       # Installation path

sh LR_Gapcloser.sh -i path/to/scaffolded/assembly/from/lrscaf/scaffolds.fasta \
                    -l path/to/the/filtered/long-read/data/EW_nanopore_merged_filtered_porechop.fasta \
                    -s n -t 10 -r 10                                  #-r 15 makes it run for 10 iterations
```
---
