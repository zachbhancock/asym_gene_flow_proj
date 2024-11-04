#!/bin/bash
#SBATCH --time=168:00:00   # walltimea
#SBATCH --account=bradburd1
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G   # memory per CPU
#SBATCH --mail-user=hancockz@umich.edu   # email address
#SBATCH --mail-type=BEGIN,END,FAIL

module load Bioinformatics
module load stacks

denovo_map.pl -T 8 -M 6 -o /nfs/turbo/lsa-bradburd/Zach/haustorius_seqs/stacks/ --samples /nfs/turbo/lsa-bradburd/Zach/haustorius_seqs/redo_de_multiplexed/ --popmap pop_file.txt --paired


exit 0
