# core
for l in 3; do # ZILN
  for j in 1 3 5; do 
    for k in 7 9 10 12  25 27 28 30 43 45 46 48; do 
      if [[ "$j" = 1 ]]; then
      sbatch --time=16:00:00 --mem 10000 2run-code.sh 0 $j $k $l 5 80;   # i ranges 1:16
      sbatch --time=64:00:00 --mem 10000  2run-code.sh 0 $j $k $l 5 400; # i ranges 1:16 
      else 
      sbatch --time=10:00:00 --mem 10000 2run-code.sh 0 $j $k $l 5 80;   # i ranges 1:10
      sbatch --time=40:00:00 --mem 10000  2run-code.sh 0 $j $k $l 5 400; # i ranges 1:10
      fi;
    done; 
  done; 
done;

# extra
for l in {1..3}; do # ZINB, ZIG, ZILN
  for j in {1..5}; do 
    for k in {1..54}; do 
      sbatch --time=10:00:00 --mem 10000 2run-code.sh 0 $j $k $l 5 80;
      sbatch --time=40:00:00 --mem 10000  2run-code.sh 0 $j $k $l 5 400; 
      if [[ "$l" = 3 ]]; then
        sbatch --time=10:00:00 --mem 10000  2run-code.sh 0 $j $k $l 2.5 80;
        sbatch --time=40:00:00 --mem 10000  2run-code.sh 0 $j $k $l 2.5 400;
        sbatch --time=10:00:00 --mem 10000  2run-code.sh 0 $j $k $l 0 80;
        sbatch --time=40:00:00 --mem 10000  2run-code.sh 0 $j $k $l 0 400;
      fi;
    done; 
  done; 
done;

# ijk = scenarios
# i is being looped inside Rscripts.
# $4 = model (1 = zinb, 2 = zig, 3 = zlin)
# $5 = perturb (5 = 50%)
# $6 = n (80, 400)
