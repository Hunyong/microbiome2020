#!/bin/sh
module load r/4.0.3  #In r/3.5.0 and MAST_1.8.2, FromMatrix gives error. Switch r to 3.5.2
Rscript C02.01.simulation-DS2update.R $1 $2 $3 $4 $5 $6
