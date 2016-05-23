ftype="$1"

echo "argon_sim 108 atoms"
./xtime.rb bin/ArNe_sim argon_108.inp $ftype
echo "argon_sim 2916 atoms"
./xtime.rb bin/ArNe_sim argon_2916_long.inp $ftype
echo "argon_sim 78732 atoms"
./xtime.rb bin/ArNe_sim argon_78732_long.inp $ftype
