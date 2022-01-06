To Finalize our assembly, we ran the folloing pipeline:

## 1. Purgehaplotigs

`Script used to run Purgehaplotig`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --job-name Purgehap.EW
#SBATCH --mem=40G
#SBATCH --time=06:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load SAMtools/1.12-GCC-9.2.0
module load minimap2/2.20-GCC-9.2.0
module load BEDTools/2.29.2-GCC-9.2.0

# Steps below are run in sequentially
# First, mapping
minimap2 -t 10 -ax map-ont path/to/EW_RNA_scaffold.fa \                            # this is the final output from Rascaf step
                   path/to/EW_nanopore_merged_filtered_porechop.fasta \
                    --secondary=no | samtools sort -m 5G -o aligned.bam -T tmp.ali

# Second, create histogram with purgehaplotigs
export PATH="/nesi/nobackup/uoo02752/.conda/envs/purge_haplotigs_env/bin:$PATH"       # path to installation of purge_haplotigs software, we installed it in conda env.
purge_haplotigs hist -b aligned.bam -g EW_RNA_scaffold.fa -t 10

# Third, coverage step with custom parameters, still needs path to purgehaplotigs installation as in second step
purge_haplotigs cov -i aligned.bam.gencov -l 2 -m 15 -h 190 -o coverage_stats.csv

# Fourth, purging step
purge_haplotigs purge -g EW_RNA_scaffold.fa -c coverage_stats.csv -b aligned.bam
```

## 2. RagTag
We than scaffolded the assembly with the discarded haplotigs using RagTag. This step increased the N50 of the assembly significantly.

`Scritp for RagTag`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --job-name ragtag.nem
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH="/path/to/the/installation/of/RagTag:$PATH"

ragtag.py scaffold curated.haplotigs.fasta \        # This is the haplotigs removed by Purgehaplotigs
                    curated.fasta                    # This is the final output from Purgehaplotigs
```
We ran RagTag again for the output from this round using the same `curated.haplotigs.fasta` to scaffold again.

## 3. Blobtools
We used Blobtools to remove contaminants from our assembly.
Carriedout installation following [this page](https://blobtoolkit.genomehubs.org/install/).

### 3.1 Coverage
We got the coverage data `bam file` by mapping long reads to the assembly.

`Script for coverage`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --job-name coverage.EW
#SBATCH --mem=40G
#SBATCH --time=07:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load minimap2/2.18-GCC-9.2.0
module load SAMtools/1.12-GCC-9.2.0

minimap2 -ax map-ont -t 10 path/to/assembly/from/ragtag/above/EW_Assembly.fasta \
                      path/to/EW_nanopore_merged_filtered_porechop.fasta | samtools sort -@10 -O BAM -o EW_Assembly.bam -
```

### 3.2 Blastn hits
We blasted the assembly against `nt` database following the Blobtools manual [How to get the dtabase](https://blobtoolkit.genomehubs.org/install/#databases)

`Script for blastn`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 10
#SBATCH --ntasks 1
#SBATCH --job-name Blastn.EW
#SBATCH --mem=2G
#SBATCH --time=1-00:00:00
#SBATCH --account=uoo02772
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load BLAST/2.10.0-GCC-9.2.0
export BLASTDB='/path/to/the/downloaded/nt/database/nt'

blastn -db nt \
        -query path/to/assembly/from/ragtag/above/EW_Assembly.fasta \
        -outfmt "6 qseqid staxids bitscore std" \
        -max_target_seqs 10 \
        -max_hsps 1 \
        -evalue 1e-25 \
        -num_threads 10 \
        -out results/EW_Assembly.ncbi.blastn.out
```
### 3.3 Diamond.blastx
We blasted the assembly against the uniprot database. We can download and format the database following the instructions [here](https://blobtoolkit.genomehubs.org/install/#databases)
`Script for Diamond.blastx`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 20
#SBATCH --job-name blastx.EW
#SBATCH --mem=20G
#SBATCH --time=2-00:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load DIAMOND/2.0.6-GCC-9.2.0


diamond blastx \
              --query path/to/assembly/from/ragtag/above/EW_Assembly.fasta \
              --db path/to/the/downloaded/and/formatted/uniprot/database/reference_proteomes.dmnd \
              --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
              --sensitive \
              --max-target-seqs 1 \
              --evalue 1e-25 \
              --threads 20 \
              --memory-limit 100 \
              --out EW_Assembly.diamond.blastx.out
```
### 3.4 Database
We then created blobdatabase, added blastn, diamond.blastx, and coverage info, and then filtered.
taxdump dataset was downloaded as described [here](https://blobtoolkit.genomehubs.org/install/#databases)
We created a text file called `Earwig_Assembly.yaml`, used in the script below.
`Earwig_Assembly.yaml` included
```
assembly:
  accession: 000_000000000.0
  alias: M_negrescens
  bioproject: 00000000
  biosample: 00000000000
  record_type: scaffolds
taxon:
  name: Forficula auricularia
  ```

`Script to create, add, and filter database`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --job-name blob.create
#SBATCH --mem=10G
#SBATCH --time=01:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH="/nesi/nobackup/uoo02752/.conda/envs/btk_env/:$PATH"                      #for the dependencies installed as suggested in blobmanual
export PATH="/nesi/nobackup/uoo02752/nematode/bin/blobtoolkit.new/blobtools2:$PATH"
export PATH="/nesi/nobackup/uoo02752/.conda/envs/btk_env/lib/python3.6/site-packages:$PATH"

# First create the database

blobtools create \
       --fasta path/to/assembly/from/ragtag/above/EW_Assembly.fasta \
       --meta Earwig_Assembly.yaml \
       --taxid 13068 \
       --taxdump path/to/downloaded/taxdump/database/taxdump \  
       Earwig_Assembly

# Second, add other info init.

blobtools add \
       --hits path/to/blastn/hits/EW_Assembly.ncbi.blastn.out \
       --hits path/to/diamond.blastx/hits/EW_Assembly.diamond.blastx.out \
       --cov path/to//coverage/bam/file/EW_Assembly.bam \
       --taxrule bestsumorder \
       --taxdump path/to/downloaded/taxdump/database/taxdump \
       Earwig_Assembly

# Third step is to filter 
blobtools filter \
        --query-string "length--Min=1000&EW_Assembly_cov--Min=5.00" \
        --fasta path/to/assembly/from/ragtag/above/EW_Assembly.fasta \
        --output path/to/output/folder/filter \
        Earwig_Assembly

# Fourth step to inverse filter, to get all those scaffolds, those were discarded in above process
blobtools filter \
        --query-string "length--Min=1000&EW_Assembly_cov--Min=5.00" \
        --fasta path/to/assembly/from/ragtag/above/EW_Assembly.fasta \
        --output path/to/output/folder/filter.inverse \
        --invert \
        Earwig_Assembly
```
From above script we got two output
1. The filtered assembly `EW_Assembly.filtered.fasta` which is in `path/to/output/folder/filter`
2. The other one is the inverse filtered file `EW_Assembly.filtered.inverse.fasta` which is in `path/to/output/folder/filter.inverse`

## 4. RagTag
Now we will use `EW_Assembly.filtered.inverse.fasta` to scaffold `EW_Assembly.filtered.fasta` with RagTag.

`Script for RagTag`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --job-name ragtag.EW
#SBATCH --mem=6G
#SBATCH --time=01:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH="/nesi/nobackup/uoo02752/nematode/bin/miniconda3/bin:$PATH"

ragtag.py scaffold path/to/output/folder/filter.inverse/EW_Assembly.filtered.inverse.fasta \
                  path/to/output/folder/filter/EW_Assembly.filtered.fasta
```
We ran RagTag four times in series, each time having the output from the last ragtag and `EW_Assembly.filtered.inverse.fasta` as inputs.
We also noticed that we need to rename the sequence header on the ragtag output fasta file. We used the following oneliner to do so.
```
awk '/^>/{print ">Scaff_" ++i; next}{print}' Ragtag_assembly.fasta > Ragtag_assembly_header_renamed.fasta
```
## 5. Pilon
This is the final step. We used the mRNA-seq data to finally polish the assembly with Pilon. we ran two iterations of this.

`Script for pilong`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --job-name pilon_EW
#SBATCH --mem=180G
#SBATCH --time=30:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load Bowtie2/2.4.1-GCC-9.2.0
module load SAMtools/1.12-GCC-9.2.0
module load Python/3.9.5-gimkl-2020a
module load Pilon/1.24-Java-15.0.2

#First step is to map sequences to assembly with bowtie
bowtie2-build path/to/assembly/from/ragtag/above/EW_assembly.fasta 
              Earwig
              
bowtie2 -p 10 --local -x Earwig \
              -1 path/to/R1/of/mRNA-seq/data/EW_all_merged_R1_trim.fq \
               -2 path/to/R2/of/mRNA-seq/data/EW_all_merged_R2_trim.fq | samtools sort > EW_assembly.fasta.bam
samtools index EW_assembly.fasta.bam EW_assembly.fasta.bai

# Second is to run Pilon
java -Xmx180G -jar $EBROOTPILON/pilon.jar --genome path/to/assembly/from/ragtag/above/EW_assembly.fasta \
              --frags EW_assembly.fasta.bam --fix snps,indels \
              --output path/to/output/EW_assembly.pilon \
              --gapmargin 1 --mingap 10000000 --threads 10 --verbose --changes \
              2>Pilon.stderr.txt 1>Pilon.stdout.txt
```
