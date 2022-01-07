# 1. Repeat library preparation
## ii. LTRharvest & LTRdigest
We need to install GenomeTools first [See](http://genometools.org/)
Finds LTR retrotransposons in the genome.

`Script for LTRharvest`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --job-name ltrpipeline
#SBATCH --mem=30G
#SBATCH --time=20:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH=$PATH:/path/to/genometools/gt/bin
export PATH=$PATH:/path/to/genometools/gt_libgenometool_lib_header/bin
module load HMMER/3.2.1-gimkl-2018b

# Run following commands sequentially
# First, create index
gt suffixerator -dna -db EW_assembly_renamed.fa -lcp -ssp -suf -tis -des -lossless

# Second, find candidate elements
gt ltrharvest -index EW_assembly_renamed.fa -tabout no -seqids -md5 > EW_assembly_renamed.ltrh.gff3

# Third, quick check of results
gt stat EW_assembly_renamed.ltrh.gff3 > stat.EW_assembly_renamed.gff3.out

# Fourth, Sort the gff3 file produced
gt gff3 -sort EW_assembly_renamed.ltrh.gff3 > EW_assembly_renamed.ltrh.sorted.gff3
```
We got the datasets from GyDB_collection to use for LTRdigest
```
wget http://gydb.org/gydbModules/collection/collection/db/GyDB_collection.zip
unzip GyDB_collection.zip
```
Additionally, the LTRdigest paper by `Steinbiss et al (2009)` gives a list of LTR retrotransposon Pfam domains in tables B1 and B2.
For example, PF03732 = Retrotrans_gag; we can download the HMM from Pfam using:
```
wget http://pfam.xfam.org/family/PF03732/hmm
```
Convert it to HMMER2 format (which ltrdigest needs) and rename it
```
hmmconvert -2 hmm > PF03732.hmm
```
Move `PF03732.hmm` to `GyDB_collection/profile` to keep them together with other `*.hmm` files
Now annotate candidates with PHMM hits using LTRDigest

`Script for LTRDigest`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --job-name ltrdigest
#SBATCH --mem=30G
#SBATCH --time=20:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH=$PATH:/path/to/genometools/gt/bin
export PATH=$PATH:/path/to/genometools/gt_libgenometool_lib_header/bin

gt -j 4 ltrdigest -hmms GyDB_collection_n_steinbiss_2009/profiles/*hmm \
                  -encseq EW_assembly_renamed.fa -matchdescstart \
                  < EW_assembly_renamed.ltrh.sorted.gff3 > EW_assembly_renamed.ltrh.sorted.ltrd.gff3
```
Get filter definition
```
wget https://raw.githubusercontent.com/satta/filterama/master/filter_protein_match.lua
```
Filter out elements with no protein domain hits
```
gt select -rule_files filter_protein_match.lua < EW_assembly_renamed.ltrh.sorted.ltrd.gff3 > EW_assembly_renamed.ltrh.sorted.ltrd.filtered.gff3
```
Check how many remained
```
gt stat EW_assembly_renamed.ltrh.sorted.ltrd.filtered.gff3 >> stat.EW_assembly_renamed.gff3.out
```
Extract sequences for each element
```
gt extractfeat -type LTR_retrotransposon -encseq EW_assembly_renamed.fa \
                -matchdescstart < EW_assembly_renamed.ltrh.sorted.ltrd.filtered.gff3 \
                > EW_assembly_ltrh.sorted.ltrd.filtered.sequences.fasta
```




