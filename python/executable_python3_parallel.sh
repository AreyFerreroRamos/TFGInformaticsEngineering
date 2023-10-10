#!/bin/bash

# Author: Arey Ferrero Ramos.
# Date: June 9, 2023. Version: 3.
# Description: This script runs a command that executes, applying parallelism, the section of the project made in SEES:lab that performs nestednes assessment.
#	Input:
#		-The type of abundance matrix.
#		-The number of randomized matrices for the assessment.
#	Output:
#		-The value of nestednes of the abundance matrix.
#		-The p-value of the nestedness value of the abundance matrix.

num_matrices=$2
num_procs=$(nproc --all)
num_procs_parallel=$(($num_procs-1))
num_matrices_proc=$(( num_matrices / num_procs_parallel ))

>executable_parallel_tmp.sh
for (( num_proc=1; num_proc<=$num_procs_parallel; num_proc++ )); do
	echo "python3 nestedness_assessment.py ../input_files/count_Genus_all.tsv ../input_files/metadata.csv ../input_files/sp_code.txt \"$1\" \"$num_matrices_proc\"" >> executable_parallel_tmp.sh
done

parallel --eta --bar --jobs $num_procs_parallel :::: executable_parallel_tmp.sh

exit 0
