#!/bin/bash
#SBATCH --time=24:00:00   # walltimea
#SBATCH --account=bradburd1
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G   # memory per CPU
#SBATCH --mail-user=hancockz@umich.edu   # email address
#SBATCH --mail-type=BEGIN,END,FAIL

module load Bioinformatics
module load stacks/2.65

populations -P /nfs/turbo/lsa-bradburd/Zach/haustorius_seqs/stacks/ --popmap new_pop_file.txt -p 10 -r 0.75 --fstats --genepop --structure --write-random-snp --vcf


exit 0