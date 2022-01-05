## Longranger-Basic
We need to run Longranger-basic in our linked read data to produced barcoded fastq file to use as an input for `Arks`

`Script to run longranger-basic`
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --job-name longranger_basicsEW
#SBATCH --mem=50G
#SBATCH --time=50:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH="/path/to/installation/directory/of/longranger-2.2.2:$PATH"

longranger basic --id=Earwig \
                  --fastqs=path/to/the/directory/containing/all/linked-read/fastq/files
```

Output file will be at: `./Earwig/outs/barcoded.fastq.gz`
