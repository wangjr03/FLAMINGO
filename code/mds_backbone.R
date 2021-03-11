args <- commandArgs(T)

# 1. input pairwise distance
# 2. output file

pd <- read.table(args[1])
pd <- as.matrix(pd)
library(MASS)

rm_id <- which(apply(pd,1,min) > 3)

pd <- pd[-rm_id,-rm_id]
pd[which(pd>3)] <- 3
p_t <- isoMDS(pd,k=3)
setwd(args[2])
p_t <- p_t$points
save(p_t,file="MDS_backbone_structure.Rdata"  )