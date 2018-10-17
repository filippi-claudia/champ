sed -n '/NATURAL ORBITALS OF MCSCF/,/END/p' $1 > tmp1.$$
sed -n '/OPTIMIZED MCSCF MO/,/END/p'        $1 > tmp2.$$
natorb_numlines=$(wc -l tmp1.$$ | awk '{print ($1)}')
sed '$d' tmp1.$$
sed -n "${natorb_numlines},\$p" tmp2.$$
rm tmp1.$$ tmp2.$$

