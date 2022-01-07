# Quickmerge
We used [Quickmerge](https://github.com/mahulchak/quickmerge) to merge `Supernova assembly` into `Flye assembly`. Flye assembly had better completeness and contiguity so it was used as a primery assembly during merging.
To do so, we used merge_wrapper script from quickmerge with `-l` value equals to the `N50` value of Supernova assembly `i.e 30358` and other parameters as in the script below.

`Script for quickmerge`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --partition=large
#SBATCH --job-name quickmerge.ew
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH=/nesi/nobackup/uoo02752/bin/quickmerge:/nesi/nobackup/uoo02752/bin/quickmerge/MUMmer3.23:$PATH

merge_wrapper.py assembly.fasta EW_10xSN.fasta \
                  -hco 5.0 -c 1.5 -l 30358 -ml 5000 -p EW_Flye_SN_merge 
```
