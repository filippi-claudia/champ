grep "Iteration "        IterData.log | sed 's/Iteration number //'             > $$_tmp1
grep "current_energy "   IterData.log | sed 's/^cur.* =//'  | sed 's/Har.*$//'  > $$_tmp2
grep "max_displ_bl "     IterData.log | sed 's/^max.* =//'  | sed 's/Boh.*$//'  > $$_tmp3
grep "xyz_rms_displ "    IterData.log | sed 's/^xyz.* =//'  | sed 's/Boh.*$//'  > $$_tmp4
grep "max_grad_zmat_bl " IterData.log | sed 's/^max.* =//'  | sed 's/Har.*$//'  > $$_tmp5
grep "rms_grad_zmat "    IterData.log | sed 's/^rms.* =//'  | sed 's/Har.*$//'  > $$_tmp6
echo "# iter    current_energy                max_displ_bl     xyz_rms_displ   max_grad_zmat_bl                rms_grad_zmat" 
paste $$_tmp1 $$_tmp2 $$_tmp3 $$_tmp4 $$_tmp5 $$_tmp6
rm $$_tmp*

