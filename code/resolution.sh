#!/bin/bash
mkdir -p ../../data/$4/$5_expanded_DNase

sbatch 7_reconstruct_model_DNase.sb $1 $2 $3 ../../data/$4/$5_expaned_frag/ ../../data/$4/$5_expanded_DNase ../../data/$4/Genomic_fragment_expanded_$6
