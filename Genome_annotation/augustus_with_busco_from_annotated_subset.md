# Training Augustus with Busco using annotated regions from the Maker output after evidence based annotation.

To do so, we took the regions containing mRNA annotations from our first round of Maker output with 1000bp on each side.
`F.auricularia_Maker_R1.all.gff` is the annotation file from Maker round 1.
To excise the region we need:
```
mkdir augustus
cd augustus
module load BEDTools/2.29.2-GCC-9.2.0

awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' path/to/F.auricularia_Maker_R1_maker_output/F.auricularia_Maker_R1.all.gff | \
awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }' | \
bedtools getfasta -fi path/to/EW_assembly.fasta -bed - -fo F.auricularia_Maker_R1_transcripts1000.fasta
```

We had to rename sequence header to make busco run on these extracted reads
```
awk '/^>/{print ">seq_" ++i; next}{print}' F.auricularia_Maker_R1_transcripts1000.fasta > F.auricularia_Maker_R1_transcripts1000_header_renamed.fasta
```

Run Busco to train augustus 
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 36
#SBATCH --job-name busco_maker.ew
#SBATCH --mem=10G
#SBATCH --time=20:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load BUSCO/5.2.2-gimkl-2020a
cp -r $AUGUSTUS_CONFIG_PATH path/to/Augustus/config/file/MyAugustusConfig_Ew
export AUGUSTUS_CONFIG_PATH=path/to/Augustus/config/file/MyAugustusConfig_Ew

busco --in path/to/F.auricularia_RL1st_maker.transcripts1000_header_renamed.fasta \
--out busco_ew_insecta -c 36 -m genome --long --augustus --augustus_parameters='--progress=true' \
-l insecta_odb10
```
It took ~18 hours to finish augustus training, while trying to train augustus from busco with whole genome assembly as input didn't even finish after 10 days.

We have the augustus trained model in `busco_ew_insecta/run_insecta_odb10/augustus_output/retraining_parameters/BUSCO_busco_ew_insecta/` from busco working directory.
we renamed the files and copied it to the species directory inside `path/to/Augustus/config/file/MyAugustusConfig_Ew`

```
cd busco_ew_insecta/run_insecta_odb10/augustus_output/retraining_parameters/BUSCO_busco_ew_insecta/
rename  "BUSCO_busco_ew_insecta" "ForficulaAuricularia" *
sed -i 's/BUSCO_busco_ew_insecta/ForficulaAuricularia/g' ForficulaAuricularia_parameters.cfg
sed -i 's/BUSCO_busco_ew_insecta/ForficulaAuricularia/g' ForficulaAuricularia_parameters.cfg.orig1
cd ..
mv BUSCO_busco_ew_insecta ForficulaAuricularia
cp -r  ForficulaAuricularia path/to/Augustus/config/file/MyAugustusConfig_Ew/species/
```
After this we can simply pass `ForficulaAuricularia` as augustus species in `maker_opt.ctl` file and export Augustus config path to MyAugustusConfig_Ew while submitting Maker job as we did above.
```
export AUGUSTUS_CONFIG_PATH=path/to/Augustus/config/file/MyAugustusConfig_Ew
```


