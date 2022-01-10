# core
for zoe in 1 2 3; do # ZILN
  for j in 2; do 
    for type in gene genebact bact; do 
      for ngene in 10000 100000; do
        for nsignal in 100 300 1000; do
          for sim in 1..10; do
            sbatch --time=16:00:00 --mem 10000 bash/runR.sh C02.11.simulation-craft.R  $zoe $j $type $sim $nsignal $ngene 1
          done;
        done;
      done;
    done; 
  done; 
done;