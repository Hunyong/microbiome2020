#!/bin/sh
#In r/3.5.0 and MAST_1.8.2, FromMatrix gives error. Switch r to 3.5.2
#In r/3.5.2 scran package does not work. Switch to 4.0.3. Scran is not any more used, and this new version often encounters fatal error. Going back to 3.5.2
module load r/4.0.3  
Rscript C02.01.simulation-Basic-ijk-batch.R $1 $2 $3 $4 $5 $6

# ijk = scenarios
# $4 = model (3 = zlin)
# $5 = perturb (5 = 50%)
# $6 = n (80, 400)
