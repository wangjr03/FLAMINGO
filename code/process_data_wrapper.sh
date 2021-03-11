#!/bin/bash
Rscript 1_preprocess_HiC.R $1 $2 $5
Rscript 2_generate_fragment_data.R $2 $5
Rscript 3_generate_fragment_DNase.R $2 $5
Rscript 4_generate_DNase_profile.R $2/Genomic_fragment $3
Rscript 5_impute_DNase_dist.R $2/Genomic_fragment
Rscript 6_generate_backbone_data.R $4 $2 $5
