[Supernova (v.2.1.1)](https://support.10xgenomics.com/de-novo-assembly/software/pipelines/latest/using/running) was used to assemble the linked reads with the following script. We used 660 million reads as input and all other default parameters.

```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --job-name Supernova_EW
#SBATCH --mem=400G
#SBATCH --time=168:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load Supernova/2.1.1

supernova run --id=EW_10xSN \
              --fastqs=/path/to/linked-read/fastq/files \
              --maxreads=660000000
```

After the completion of the assembly we produced the assembly in fasta file with `pseudohap` style with the following script

```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --job-name Supernova_EW_fasta
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load Supernova/2.1.1

supernova mkoutput \
        --style=pseudohap \
        --asmdir=path/to/working-directory-of-supernova/EW_10xSN/outs/assembly \      # this is a path to the output directory `assembly` created by supernova
        --outprefix=./EW_10xSN
```

After running this script we will get the assembly called `EW_10xSN.fasta.gz`
We checked the assembly statistics with the [Quast](https://github.com/ablab/quast). [Script here](quast.sh)
