# Long read assembly pipeline
---
## 1. Basecalling (Guppy)
We used Guppy version 5.0.7 to basecall the `fast5` files produced by minions. It was ran in `GPU partition` with the following script.

```
#!/bin/bash -e

#SBATCH --job-name=Guppy_EW                 
#SBATCH --account=uoo02752              
#SBATCH --time=10:00:00                
#SBATCH --partition=gpu                 #guppy runs faster in gpu partition in Nesi, than other partition
#SBATCH --gres=gpu:1                    #some configuration for gpu partition, suggested by Nesi support
#SBATCH --mem=6G                                
#SBATCH --ntasks=4                              
#SBATCH --cpus-per-task=1               
#SBATCH --output=%x-%j.out              
#SBATCH --error=%x-%j.err               
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz

module load ont-guppy-gpu/5.0.7
guppy_basecaller -i ../path/to/all/fast5/files/ \
                 -s ./path/to/output/directory/ \
                 --flowcell FLO-MIN106 \
                 --kit SQK-LSK109 \
                 --num_callers 4 -x auto \
                 --recursive \                                  # recursive flag looks for fast5 files in subfolders as well
                 --trim_barcodes --disable_qscore_filtering     # disabled the quality filtering to get all the fastq files produced in one folder instead of pass, and fail. we will do QC later
