# 1. Repeat library preparation
## iii. TransposonPSI
To install the software and for more detail information [see here](http://transposonpsi.sourceforge.net/).

`Script to run TransposonPSI`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --job-name transposonPSI
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH="/path/to/conda/environment/:$PATH" for dependencies

/path/to/transposonPSI.pl path/to/assembly/EW_assembly_renamed.fa nuc
```
