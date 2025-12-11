#!/bin/bash

module load Bioinformatics
module load stacks

denovo_map.pl -T 8 -M 6 -o "your_directory"/stacks/ --samples "your_directory"/ --popmap pop_file.txt --paired

exit 0
