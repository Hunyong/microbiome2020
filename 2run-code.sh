#!/bin/sh
module load r/3.5.2  #In r/3.5.0 and MAST_1.8.2, FromMatrix gives error. Switch r to 3.5.2
Rscript C02.01.simulation-Basic-ijk-batch.R $1 $2 $3 $4 $5 $6

# ijk = scenarios
# $4 = model (3 = zlin)
# $5 = perturb (5 = 50%)
# $6 = n (80, 400)
