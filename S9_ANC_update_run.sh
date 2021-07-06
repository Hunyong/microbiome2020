# ZINB
for l in {1..3}; do # ZINB, ZIG, ZILN
  for j in {0..0}; do 
    for k in {1..54}; do 
      sbatch --time=2:00:00 --mem 10000 S9_DS2_update.sh 0 $j $k $l 5 80;
      sbatch --time=3:00:00 --mem 10000  S9_DS2_update.sh 0 $j $k $l 5 400; 
      if [[ "$l" = 3 ]]; then
        sbatch --time=2:00:00 --mem 10000  S9_DS2_update.sh 0 $j $k $l 2.5 80;
        sbatch --time=3:00:00 --mem 10000  S9_DS2_update.sh 0 $j $k $l 2.5 400;
        sbatch --time=2:00:00 --mem 10000  S9_DS2_update.sh 0 $j $k $l 0 80;
        sbatch --time=3:00:00 --mem 10000  S9_DS2_update.sh 0 $j $k $l 0 400;
      fi;
    done; 
  done; 
done;