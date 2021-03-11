#!/bin/bash
mkdir -p ../../data/$6/$4_$6/$5_result
sbatch 7_reconstruct_model_DNase.sb $1 $2 $3 ../../data/$4/$5_frag/ ../../data/$6/$4_$6/$5_result ../../data/$6/Genomic_fragment
#Rscript 8_ensemble_structure.R ../../data/$4/$5_IF.Rdata $2 ../../data/$4_$6/$5_result ../../data/$4/$5_PD.Rdata ensembled_structure
