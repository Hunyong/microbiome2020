#!/bin/sh
module load r/3.5.2
Rscript C01.01.goodness_of_fit.R $1 $2;

# for i in {1..2}; do # ZOE 1, 2
#   for j in {1..3}; do # tpm, rpk, asin
#     sbatch S01.gof.sh $i $j;
#   done;
# done;