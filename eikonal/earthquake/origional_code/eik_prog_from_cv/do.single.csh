set sh_dir = ~/PROGS_64/EARTHQUAKES
set fsta1 = event.short.lst

set snr1 = 5
set snr2 = 1
set per = 50

set lon1 = 235
set lon2 = 285
set lat1 = 20
set lat2 = 50
set n1 = 251
set n2 = 151

set event = 20061225200058
#set inf = 20061225200058_3pi.${per}
set inf = $event.${per}

awk '{if ($14>'$snr1' && $15 > '$snr2' && $10>0) print $7,$8,$11,$9,$13,$3,1; else print $7,$8,$12,$10,$13,$3,0;}' $inf > $inf.input
${sh_dir}/correct_2pi_v1_jy_ph.cv.time.earth $inf.input ${per} 20   # make corrections
csh ${sh_dir}/PLOT/plot.correct.csh $inf.input.c.txt

csh ${sh_dir}/C_plot_travel_US $inf.input.c.txt $lon1 $lon2 $lat1 $lat2
csh ${sh_dir}/C_plot_travel_am_US  $inf.input.c.txt $lon1 $lon2 $lat1 $lat2

${sh_dir}/correct_travel_time_curvature_v1_270_c $inf.input.c.txt $per $lon1 $lat1 $n1 $n2

csh ${sh_dir}/C_plot_travel_US $inf.input.c.txt_v2 $lon1 $lon2 $lat1 $lat2
csh ${sh_dir}/C_plot_travel_T0.2_US  $inf.input.c.txt_v2 $lon1 $lon2 $lat1 $lat2
csh ${sh_dir}/C_plot_travel_am_US $inf.input.c.txt_v2 $lon1 $lon2 $lat1 $lat2

set cdist = 200
/home/weisen/PROGS_64/EARTHQUAKES/travel_time_to_slow_map_v4_v1_cv_us_v1_2_c_single_azi_US $fsta1 ${per} 0.2 $lon1 $n1 $lat1 $n2 $cdist $event

ls $inf.input.c.txt > teve
~/PROGS_64/EARTHQUAKES/amp_HD_to_amp_gradient_HD_to_amp_laplace_HD_input_small_region_US teve $per $lon1 $lat1 $n1 $n2

sort -nk2 $inf.input.c.txt_v2_am.HD > t1
sort -nk2 $inf.input.c.txt_am_laplace.txt.HD > t2
sdiff t1 t2 | awk '{if ($3>0.) print $1,$2,$7*('$per'**2)/($3*(2*3.1415926)**2); else print $1,$2,-1.;}' > correction
sort -nk2 slow_azi_$inf.txt.HD.2.v2 > t3
sdiff t3 correction | awk '{if ($3>0. && $11>-0.5) print $1,$2,($3*$3-$11)**0.5,$3,$11; else {print $1,$2,-999,-999,-999}}' > $inf.slow.ac
#echo $event > teve1
#${sh_dir}/FIT/fit_tr_amp_to_second_order_input_small_region teve1 50s_iso_ani_v1.iso $lon1 $lat1 $n1 $n2 $per > log.log
