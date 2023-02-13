# 2. Concatenating, filtering, and classifying repeats


Merging individual repeat library from `RepeatModuler`, `LTRharvest&LTRdigest`, `TransposonPSI`, `SINE database`

```
cat path/to/repeatmoduler/output/consensi.fa.classified \
    path/to/LTRharvest&LTRdigest/pipeline/output/EW_assembly_ltrh.sorted.ltrd.filtered.sequences.fasta \
    path/to/TransposonPSI/output/EW_assembly_renamed.TPSI.allHits.chains.bestPerLocus.fa \
    path/to/SINE/database/SINEs.bnk > Combined.rep.library.fasta
```
## Filtering
We filtered out sequences less than 50bp.
```
module load seqtk/1.3-gimkl-2018b
seqtk seq -L 50 Combined.rep.library.fasta > Combined.rep.library.minlen50.fasta
```
Then we removed redundant sequences using usearch
```
module load USEARCH/11.0.667-i86linux32
usearch -cluster_fast Combined.rep.library.minlen50.fasta -id 0.8 -consout Combined.rep.library.minlen50.usearch.fasta
```
## RepeatClassifier
We than used repeatclassifier to classify the repeat library

`Script used for repeatclassifier`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 16
#SBATCH --job-name RepClas.EW
#SBATCH --mem=10G
#SBATCH --time=48:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load RepeatModeler/2.0.2-Miniconda3

RepeatClassifier -consensi Combined.rep.library.minlen50.usearch.fasta
```
We filtered out all the sequences with Unknown hits after classification and blasted (blastx, evalue <1e-10) against [reviewed insect uniprot database](https://www.uniprot.org/uniprot/?query=reviewed:yes%20taxonomy:50557)

`Filteringout Unknown hits`
```
grep "Unknown" Combined.rep.library.minlen50.usearch.fasta.classified > UnknownIDs.txt
sed -i 's/>//g' UnknownIDs.txt
faSomeRecords Combined.rep.library.minlen50.usearch.fasta.classified UnknownIDs.txt Unknown_repeats.fasta
```
can also use seqtk instead of faSomeRecords `seqtk subseq Combined.rep.library.minlen50.usearch.fasta.classified UnknownIDs.txt > Unknown_repeats.fasta`

Blasting against downloaded Uniprot database

```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 8
#SBATCH --job-name blast.ew
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load BLAST/2.12.0-GCC-9.2.0

# First prepare the database
makeblastdb -in path/to/UniprotSeqs/uniprot-taxonomy__Insecta+[50557]_-filtered-reviewed_yes.fasta -dbtype prot

# Blastx
blastx -query Unknown_repeats.fasta -db path/to/uniprot-taxonomy__Insecta+[50557]_-filtered-reviewed_yes.fasta -evalue 1e-10 \
        -num_threads 10 -max_target_seqs 1 -outfmt '6 qseqid sseqid evalue bitscore sgi sacc stitle' -out Blast_out.txt
```
To filter the repeat library from the blastx report
```
awk -F "\t" '{print $1,$7}' Blast_out.txt  | sort | uniq | grep -i -v "transposon" | grep -i -v "Copia protein" | grep -i -v "mobile element" | \
grep -i -v "transposable"  | grep -i -v "transposase" | awk '{print $1}' > Unknowns_with_Port_hit.txt
```
```
faSomeRecords -exclude Combined.rep.library.minlen50.usearch.fasta.classified Unknowns_with_Port_hit.txt EW_Rep_CombinedLibrary.lib.minlen50.nr.classified.filtered.fa
```
The final Repeat library is `EW_Rep_CombinedLibrary.lib.minlen50.nr.classified.filtered.fa`

PS: `faSomeRecords` script if needed is [here](https://github.com/santiagosnchez/faSomeRecords) 

