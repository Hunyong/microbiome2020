#!/bin/bash
#SBATCH --job-name=plot
#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=10g
#SBATCH -n 1
#SBATCH -c 12
#SBATCH -t 03-00:00:00
#SBATCH --output=./Report/slurm-%j.out
#SBATCH --mail-type=end

# cd is really important
cd /proj/cosd_lab/dwuWrite/yixiang/microbiome2020

module load r/4.0.3  
Rscript C02.02.metrics_plot.R
Rscript C02.02.power_plot.R