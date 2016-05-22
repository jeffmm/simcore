echo "Building with no optimizations, single threaded"

make remove
make CFG=test

echo "Running test on Argon 108 atoms"
echo "Static test first"
echo ""

echo "Brute Force Scheme (static)"
./xtime.rb bin/ArNe_sim argon_108_static.inp brute &> tests/static_brute_test.txt

echo "Cell List Scheme (static)"
./xtime.rb bin/ArNe_sim argon_108_static.inp cells &> tests/static_cell_test.txt

echo "Neighbor List All Pairs Scheme (static)"
./xtime.rb bin/ArNe_sim argon_108_static.inp allpairs &> tests/static_neighbor_allpairs_test.txt


echo ""
echo "Dynamic test"

echo "Brute Force Scheme (dynamic)"
./xtime.rb bin/ArNe_sim argon_108_dynamic.inp brute &> tests/dynamic_brute_test.txt
