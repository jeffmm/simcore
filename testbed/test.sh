make remove
make CFG=test

echo "Running test on Argon 108 atoms"
echo "Building with no optimizations, single threaded"
echo ""
echo "Static tests:"

echo "Brute Force Scheme (static)"
./xtime.rb bin/ArNe_sim argon_108_static.inp brute &> tests/static_brute_test.txt

echo "Cell List Scheme (static)"
./xtime.rb bin/ArNe_sim argon_108_static.inp cells &> tests/static_cell_test.txt

echo "Neighbor List All Pairs Scheme (static)"
./xtime.rb bin/ArNe_sim argon_108_static.inp allpairs &> tests/static_neighbor_allpairs_test.txt

echo "Neighbor List Cells (static)"
./xtime.rb bin/ArNe_sim argon_108_static.inp neighborcells &> tests/static_neighbor_cells_test.txt


echo ""
echo "Dynamic tests:"

echo "Brute Force Scheme (dynamic)"
./xtime.rb bin/ArNe_sim argon_108_dynamic.inp brute &> tests/dynamic_brute_test.txt

echo "Cell List Scheme (dynamic)"
./xtime.rb bin/ArNe_sim argon_108_dynamic.inp cells &> tests/dynamic_cell_test.txt

echo "Neighbor List All Pairs Scheme (dynamic)"
./xtime.rb bin/ArNe_sim argon_108_dynamic.inp allpairs &> tests/dynamic_neighbor_allpairs_test.txt

echo "Neighbor List Cells Scheme (dynamic)"
./xtime.rb bin/ArNe_sim argon_108_dynamic.inp neighborcells &> tests/dynamic_neighbor_cells_test.txt

make remove
make
echo "Building with full optimizations, OpenMP"
echo ""
echo "Static tests (optimized):"

echo "Brute Force Scheme (static, optimized)"
./xtime.rb bin/ArNe_sim argon_108_static.inp brute &> tests/static_brute_opt.txt

echo "Cell List Scheme (static, optimized)"
./xtime.rb bin/ArNe_sim argon_108_static.inp cells &> tests/static_cell_opt.txt

echo "Neighbor List All Pairs Scheme (static, optimized)"
./xtime.rb bin/ArNe_sim argon_108_static.inp allpairs &> tests/static_neighbor_allpairs_opt.txt

echo "Neighbor List Cells  Scheme (static, optimized)"
./xtime.rb bin/ArNe_sim argon_108_static.inp neighborcells &> tests/static_neighbor_cells_opt.txt


echo ""
echo "Dynamic tests (optimized):"

echo "Brute Force Scheme (dynamic, optimized)"
./xtime.rb bin/ArNe_sim argon_108_dynamic.inp brute &> tests/dynamic_brute_opt.txt

echo "Cell List Scheme (dynamic, optimized)"
./xtime.rb bin/ArNe_sim argon_108_dynamic.inp cells &> tests/dynamic_cell_opt.txt

echo "Neighbor List All Pairs Scheme (dynamic, optimized)"
./xtime.rb bin/ArNe_sim argon_108_dynamic.inp allpairs &> tests/dynamic_neighbor_allpairs_opt.txt

echo "Neighbor List Cells Scheme (dynamic, optimized)"
./xtime.rb bin/ArNe_sim argon_108_dynamic.inp neighborcells &> tests/dynamic_neighbor_cells_opt.txt
