# ZINB
for l in {1..3}; do # ZINB, ZIG, ZILN
  for j in {1..5}; do 
#    for k in {1..54}; do 
    for k in {1,7,10,13,16,19,22,25,28,31,34}; do 
      sbatch --time=10:00:00 --mem 10000 2run-code.sh 0 $j $k $l 5 80;
    done; 
  done; 
done;

for l in {2..3}; do # ZIG, ZILN
  for j in {1..5}; do 
    for k in {1,7,13,19,25,31}; do 
      sbatch --time=40:00:00 --mem 10000  2run-code.sh 0 $j $k $l 5 400; 
    done; 
  done; 
done;

for l in {3..3}; do # ZILN
  for j in {1..5}; do 
#    for k in {1..54}; do 
    for k in {19,25,31}; do 
      sbatch --time=10:00:00 --mem 10000  2run-code.sh 0 $j $k $l 2.5 80;
      sbatch --time=40:00:00 --mem 10000  2run-code.sh 0 $j $k $l 2.5 400;
      sbatch --time=10:00:00 --mem 10000  2run-code.sh 0 $j $k $l 0 80;
      sbatch --time=40:00:00 --mem 10000  2run-code.sh 0 $j $k $l 0 400;
    done; 
  done; 
done;