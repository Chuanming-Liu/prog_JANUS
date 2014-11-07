unlimit
set pos1 = `pwd`
csh tmp.csh  /media/disk-1/Data_cut /home/jiayi/work/SC/Info/evt.txt /home/jiayi/work/SC/Info/station_g.lst
#run rm_b stations => 355 for 25~40s
csh correct_2pi_curvature.csh ./region.txt .
csh tool_do_surface_slow_lapalce_all.csh ./region.txt .

# run rm_b and GZDJT => 354 for 25-70s
# also run the ani with 1.2deg 25points  => v1_1.2deg
set pos2 = $pos1/rm_sta_GZDJT_390
cd $pos2
csh do_all.csh
cd $pos1

# run rm_b and rm some JS stations(6) => 355-6 for 30~50s
set pos3 = $pos1/rm_sta_JS_390
cd $pos3
csh do_all.csh
cd $pos1







