
csh FTAN_to_traveltime_Eqk.csh 
unlimit
echo "-R/-125/-80/24/50" > region.txt
csh correct_2pi_curvature.csh ./region.txt  .

csh tool_do_surface_slow_lapalce_all.csh ./region.txt .
