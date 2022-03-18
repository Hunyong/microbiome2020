# core
for l in 1 2 3; do # ZILN
  for j in 1 3 5; do 
    for k in 7 9 10 12 25 27 28 30 43 45 46 48; do 
      for replica in 1 2 3 4 5 6 7 8 9 10; do
        sbatch --output=./Report/slurm-%j.out --time=64:00:00 --mem 10000 bash/2run-code-metrics.sh 0 $j $k $l 5 80 0.1 $replica;
        sbatch --output=./Report/slurm-%j.out --time=64:00:00 --mem 10000  bash/2run-code-metrics.sh 0 $j $k $l 5 400 0.1 $replica;
        done;
    done; 
  done; 
done;
