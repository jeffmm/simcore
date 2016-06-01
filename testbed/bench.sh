ncores=(1 2 4)
afiles=(argon_108.inp argon_2916.inp argon_78732.inp)
ftype="$1"

echo "Running argon benchmarks"

for j in "${afiles[@]}"
do
  echo "FILE $j"
  for k in "${ncores[@]}"
  do
    OMP_NUM_THREADS=$k ./xtime.rb bin/ArNe_sim $j $ftype
  done
done
