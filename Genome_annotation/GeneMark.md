# GeneMark


`Script: GeneMark`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 20
#SBATCH --job-name Genemark.EW
#SBATCH --mem=10G
#SBATCH --time=05:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread


module load GeneMark-ES/4.62-GCC-9.2.0

/scale_wlg_persistent/filesets/opt_nesi/CS400_centos7_bdw/GeneMark-ES/4.62-GCC-9.2.0/gmes_petap.pl \
            --sequence path/to/assembly/file/EW_assembly.fasta \
            --ES --cores $SLURM_NTASKS
```
