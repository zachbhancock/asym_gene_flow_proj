#!/bin/bash

module load Bioinformatics
module load stacks/2.65

populations -P "your_directory"/stacks/ --popmap new_pop_file.txt -p 10 -r 0.75 --fstats --genepop --structure --write-random-snp --vcf


exit 0
