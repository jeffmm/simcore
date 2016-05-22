echo "Building with no optimizations, single threaded"

make remove
make CFG=test

echo "Running test on Argon 108 atoms"

echo "Brute Force Scheme"
./xtime.rb bin/ArNe_sim argon_108.inp brute &> brute_test.txt

echo "Cell List Scheme"
./xtime.rb bin/ArNe_sim argon_108.inp cells &> cell_test.txt

echo "Neighbor List All Pairs Scheme"
./xtime.rb bin/ArNe_sim argon_108.inp allpairs &> neighbor_allpairs_test.txt
