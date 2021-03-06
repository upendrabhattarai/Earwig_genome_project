# Maker pipeline
We used maker pipeline for the genome annotation.
Some very useful resources we went through to make our pipeline work are:

[Resource 1](http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/MAKER_Tutorial_for_WGS_Assembly_and_Annotation_Winter_School_2018#Training_ab_initio_Gene_Predictors), [Resource 2](https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2),[Resource 3](https://bioinformaticsworkbook.org/dataAnalysis/GenomeAnnotation/Intro_To_Maker.html#gsc.tab=0),[Resource 4](https://github.com/guillemylla/Crickets_Genome_Annotation/blob/master/G_bimaculatus/Protein_Coding_Genes.md)

## Downloading and formatting publicly available data
We downloaded all the dermapteran protein and mRNA-seq data from NCBI (date: 12-Dec-2021)
We also downloaded all the dermapteran protein data from uniprot database (date: 12-Dec-2021)

Both protein data (NCBI and Uniprot) were merged and placed in a folder `proteins`
We separated `Forficula auricularia` mRNA data from `non-forficula auricularia` using cbdtools as below. 
Reads from `Forficula auricularia` were placed in a folder `est`
and `non-forficula auricularia` reads were placed in `altest` folder.
We also placed the `Trinity.fasta` generated after denovo assembly of the mRNA-seq data in `est` folder.

```
#indexing 
cdbfasta ncbi-dermaptera-mrna.fasta

grep "auricularia" ncbi-dermaptera-mrna.fasta | \
sed 's/>//g'| path/to/installataion/folder/cdbfasta/cdbyank \
ncbi-dermaptera-mrna.fasta.cidx > est/F.auriculariaEST.fasta

grep -v "auricularia" ncbi-dermaptera-mrna.fasta | \
sed 's/>//g'| /path/to/installataion/folder/cdbfasta/cdbyank \
ncbi-dermaptera-mrna.fasta.cidx > altest/Not_F.auriculariaEST.fasta
```
## First round of Maker
Creating Maker config files:
```
module load MAKER/2.31.9-gimkl-2020a
maker -CTL
```
This will create `maker_exe.ctl` `maker_bopt.ctl` `maker_opts.ctl`.
We edited the `maker_opts.ctl` file as below for the first round of maker

`maker_opts.ctl`
```
#-----Genome (these are always required)
genome=/path/to/assembly/file/EW_assembly.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=1 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=1 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=1 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=1 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=1 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=1 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=1 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est=path/to/est/Trinity.fasta,path/to/est/F.auriculariaEST.fasta #set of ESTs or assembled mRNA-seq in fasta format
altest=path/to/altest/Not.F.auriculariaEST.fasta #EST/cDNA sequence file in fasta format from an alternate organism
est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=path/to/proteins/merged_ncbi_uniprot_dermaptera_protein.fasta  #protein sequence file in fasta format (i.e. from mutiple oransisms)
protein_gff=  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org=all #select a model organism for RepBase masking in RepeatMasker
rmlib=path/to/repeat/libraray/EW_Rep_CombinedLibrary.lib.minlen50.nr.classified.filtered.fa #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein=/scale_wlg_persistent/filesets/opt_nesi/CS400_centos7_bdw/MAKER/2.31.9-gimkl-2020a/data/te_proteins.fasta #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff= #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm= #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species= #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=300000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1000 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=1 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP=path/to/tmp #specify a directory other than the system default temporary directory for temporary files
```
Below is the content of `maker_bopts.ctl`

```
#-----BLAST and Exonerate Statistics Thresholds
blast_type=ncbi+ #set to 'ncbi+', 'ncbi' or 'wublast'

pcov_blastn=0.8 #Blastn Percent Coverage Threhold EST-Genome Alignments
pid_blastn=0.85 #Blastn Percent Identity Threshold EST-Genome Aligments
eval_blastn=1e-10 #Blastn eval cutoff
bit_blastn=40 #Blastn bit cutoff
depth_blastn=0 #Blastn depth cutoff (0 to disable cutoff)

pcov_blastx=0.5 #Blastx Percent Coverage Threhold Protein-Genome Alignments
pid_blastx=0.4 #Blastx Percent Identity Threshold Protein-Genome Aligments
eval_blastx=1e-06 #Blastx eval cutoff
bit_blastx=30 #Blastx bit cutoff
depth_blastx=0 #Blastx depth cutoff (0 to disable cutoff)

pcov_tblastx=0.8 #tBlastx Percent Coverage Threhold alt-EST-Genome Alignments
pid_tblastx=0.85 #tBlastx Percent Identity Threshold alt-EST-Genome Aligments
eval_tblastx=1e-10 #tBlastx eval cutoff
bit_tblastx=40 #tBlastx bit cutoff
depth_tblastx=0 #tBlastx depth cutoff (0 to disable cutoff)

pcov_rm_blastx=0.5 #Blastx Percent Coverage Threhold For Transposable Element Masking
pid_rm_blastx=0.4 #Blastx Percent Identity Threshold For Transposbale Element Masking
eval_rm_blastx=1e-06 #Blastx eval cutoff for transposable element masking
bit_rm_blastx=30 #Blastx bit cutoff for transposable element masking

ep_score_limit=20 #Exonerate protein percent of maximal score threshold
en_score_limit=20 #Exonerate nucleotide percent of maximal score threshold
```
`maker_exe.ctl` file has the path for all the executables set in the server
```
#-----Location of Executables Used by MAKER/EVALUATOR
makeblastdb=/opt/nesi/CS400_centos7_bdw/RMBlast/2.10.0-GCC-9.2.0/bin/makeblastdb #location of NCBI+ makeblastdb executable
blastn=/opt/nesi/CS400_centos7_bdw/RMBlast/2.10.0-GCC-9.2.0/bin/blastn #location of NCBI+ blastn executable
blastx=/opt/nesi/CS400_centos7_bdw/RMBlast/2.10.0-GCC-9.2.0/bin/blastx #location of NCBI+ blastx executable
tblastx=/opt/nesi/CS400_centos7_bdw/RMBlast/2.10.0-GCC-9.2.0/bin/tblastx #location of NCBI+ tblastx executable
formatdb= #location of NCBI formatdb executable
blastall= #location of NCBI blastall executable
xdformat= #location of WUBLAST xdformat executable
blasta= #location of WUBLAST blasta executable
RepeatMasker=/opt/nesi/CS400_centos7_bdw/RepeatMasker/4.1.0-gimkl-2020a/RepeatMasker #location of RepeatMasker executable
exonerate=/opt/nesi/CS400_centos7_bdw/exonerate/2.2.0-GCC-9.2.0/bin/exonerate #location of exonerate executable

#-----Ab-initio Gene Prediction Algorithms
snap=/opt/nesi/CS400_centos7_bdw/KorfSNAP/2013-11-29-GCC-9.2.0/bin/snap #location of snap executable
gmhmme3=/scale_wlg_persistent/filesets/opt_nesi/CS400_centos7_bdw/GeneMark-ES/4.62-GCC-9.2.0/gmhmme3 #location of eukaryotic genemark executable
gmhmmp= #location of prokaryotic genemark executable
augustus=/opt/nesi/CS400_centos7_bdw/AUGUSTUS/3.3.3-gimkl-2020a/bin/augustus #location of augustus executable
fgenesh= #location of fgenesh executable
tRNAscan-SE=/scale_wlg_persistent/filesets/opt_nesi/CS400_centos7_bdw/tRNAscan-SE/2.0.7-GCC-9.2.0/bin/ #location of trnascan executable
snoscan= #location of snoscan executable

#-----Other Algorithms
probuild=/scale_wlg_persistent/filesets/opt_nesi/CS400_centos7_bdw/GeneMark-ES/4.62-GCC-9.2.0/probuild #location of probuild executable (required for genemark)
```

We ran maker job with the following script with 400Gb of total memory in 36 cores.
`maker.sl`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks=36
#SBATCH --job-name Maker.EW.round1
#SBATCH --mem=400Gb
#SBATCH --time=5-00:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load MAKER/2.31.9-gimkl-2020a

srun maker -q -c $SLURM_NTASKS -base F.auricularia_Maker_R1 maker_opts.ctl maker_bopts.ctl maker_exe.ctl
```
Togenarate a combined output
```
cd F.auricularia_Maker_R1
fasta_merge -d F.auricularia_Maker_R1_master_datastore_index.log
gff3_merge -d F.auricularia_Maker_R1_master_datastore_index.log
```
## Maker With Ab Initio Gene Predictors

## Training SNAP
We used the output from the maker's first round to train SNAP

```
mkdir snap
cd snap

module load MAKER/2.31.9-gimkl-2020a

cp path/to/F.auricularia_Maker_R1_maker_output/F.auricularia_Maker_R1.all.gff ./
maker2zff F.auricularia_Maker_R1.all.gff
fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna
mkdir params # next step will produce a lots of intermediate files so we will create a subdirectory and perform it there
cd params
forge ../export.ann ../export.dna
cd ../
hmm-assembler.pl F.auricularia_Maker_R1 ./params > F.auricularia_Maker_R1.hmm
```
## Training Augustus
We trained augustus using `Braker`, which we will be using in our maker pipeline.
Busco is another popular program used to train Augustus for genome annotation. We experienced that Running busco with the whole genome assembly to train Augustus took forever. However, after the first round of maker, we can excise the regions containing mRNA annotation with 1000bp flanking regions each side and use that set to train Augustus with Busco. You can find the script [here](augustus_with_busco_from_annotated_subset.md) if you want to do so.

# Maker 2nd round
For the second round of Maker we did ab initio gene prediction using trained model from `SNAP`, and `Augustus`, as well as the braker annotation `braker.gff` as `pred_gff` file.

First lets extract the mapping information from all our input data on the first round of maker as .gff files.
So that we don't have to run blast again. This will save significant run time.

To do so,
```
mkdir maker_R1_outputs
cd maker_R1_outputs
cp path/to/F.auricularia_Maker_R1_maker_output/F.auricularia_Maker_R1.all.gff ./

# Extract protein alignments
awk '{ if ($2 == "protein2genome") print $0 }' F.auricularia_Maker_R1.all.gff > F.auricularia_Maker_all_R1_protein2genome.gff
# Extract transcript alignments
awk '{ if ($2 == "est2genome") print $0 }' F.auricularia_Maker_R1.all.gff > F.auricularia_Maker_all_R1_est2genome.gff
# Extract repeat alignments
awk '{ if ($2 ~ "repeat") print $0 }' F.auricularia_Maker_R1.all.gff > F.auricularia_Maker_all_R1_repeats.gff
```
We then backedup the config files from the first round and edited `maker_opts.ctl` as follows and ran maker job with different base name in the same folder where we ran first round of maker.

`maker_opts.ctl`
```
#-----Genome (these are always required)
genome=/path/to/assembly/file/EW_assembly.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff=path/to/maker_R1_outputs/maker.gene.models.gff #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est= #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff=path/to/maker_R1_outputs/F.auricularia_Maker_all_R1_est2genome.gff #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=  #protein sequence file in fasta format (i.e. from mutiple oransisms)
protein_gff=path/to/maker_R1_outputs/F.auricularia_Maker_all_R1_protein2genome.gff  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org= #select a model organism for RepBase masking in RepeatMasker
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff=path/to/maker_R1_outputs/F.auricularia_Maker_all_R1_repeats.gff #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm=path/to/snap.output/F.auricularia_Maker_R1.hmm #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species=ForAuricularia #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=300000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1000 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=1 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=1 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP=path/to/tmp #specify a directory other than the system default temporary directory for temporary files
```
Submitting the maker job round 2
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks=36
#SBATCH --job-name Maker.EW.R2
#SBATCH --mem=250Gb
#SBATCH --time=1-00:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load MAKER/2.31.9-gimkl-2020a
export AUGUSTUS_CONFIG_PATH=/nesi/nobackup/uoo02752/nematode/busco_downloads/MyAugustusConfig_Ew

srun maker -q -c $SLURM_NTASKS -base F.auricularia_Maker_R2 maker_opts.ctl maker_bopts.ctl maker_exe.ctl
```

## Training SNAP for the third round of maker
We trained SNAP again using the output from the second round of maker.
Scripts as above.

We again need to extract `est2genome` `protein2genome` and `repeats` mapping files.


To do so,
```
mkdir maker_R1_outputs
cd maker_R2_outputs
cp path/to/F.auricularia_Maker_R2_maker_output/F.auricularia_Maker_R2.all.gff ./

# Extract protein alignments
awk '{ if ($2 ~ "protein2genome") print $0 }' F.auricularia_Maker_R2.all.gff > F.auricularia_Maker_all_R2_protein2genome.gff
# Extract transcript alignments
awk '{ if ($2 ~ "est2genome") print $0 }' F.auricularia_Maker_R2.all.gff > F.auricularia_Maker_all_R2_est2genome.gff
# Extract repeat alignments
awk '{ if ($2 ~ "repeat") print $0 }' F.auricularia_Maker_R2.all.gff > F.auricularia_Maker_all_R2_repeats.gff
```
We need to prepare the configuration file `maker_opts.ctl` to run the third round of maker.

`maker_opts.ctl`
```
#-----Genome (these are always required)
genome=path/to/assembly/file/EW_assembly.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff=path/to/maker/gene/models/from/second/round/maker_annotation.gff #MAKER derived GFF3 file
est_pass=1 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=1 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=1 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=1 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=1 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=1 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=1 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est= #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff=path/to/F.auricularia_R2.all.est2genome.gff #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=  #protein sequence file in fasta format (i.e. from mutiple oransisms)
protein_gff=path/to/F.auricularia_R2.all.protein2genome.gff  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org=all #select a model organism for RepBase masking in RepeatMasker
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff=path/to/F.auricularia_R2.all.repeats.gff #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm=path/to/snap/output/after/second/round/of/maker/F.auricularia_R2.Maker.hmm #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species=ForAuricularia #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff=path/to/braker.gff #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=300000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1000 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=1 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP=/nesi/nobackup/uoo02752/tmp_maker/ #specify a directory other than the system default temporary directory for temporary files
```

Togenarate a final output
```
cd F.auricularia_Maker_R2
fasta_merge -d F.auricularia_Maker_R2_master_datastore_index.log
gff3_merge -d F.auricularia_Maker_R2_master_datastore_index.log
```
This gives us the final Maker annotation. We will post process the output to add functional annotation.

To calculate the statistics and different outputs from the maker annotation we used `gaas_maker_merge_outputs_from_datastore.pl` script from [GAAS tools](https://github.com/NBISweden/GAAS). It has got easy conda installation.

```
gaas_maker_merge_outputs_from_datastore.pl -i path/to/maker/output/directory/maker_output/ -o gaas_maker_output
```





