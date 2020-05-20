for i in {1..10}; do for j in {1..5}; do for k in {1..54}; do sbatch --time=25:00:00 2run-code.sh $i $j $k 3 5 400; done; done; done;

# ijk = scenarios
# $4 = model (3 = zlin)
# $5 = perturb (5 = 50%)
# $6 = n (80, 400)
