#!/bin/csh
if ( $#argv != 2 )then
echo "usage: do_laplace [region_inputfile] [sac_path: sac_path/persec_snr/20***.ph.txt***]"
exit
endif

set snr_cri =  10
set distlst = ( 0 15000)
cd $argv[2]
foreach per ( 50 60 70  )
#foreach per ( 30 )
    echo "\n\nworking on "$per"sec period...\n"
    foreach i ( 1   )
	set j = `echo $i | awk '{print $1+1}'`
        set dis1 = ${distlst[$i]}
        set dis2 = ${distlst[$j]}
	echo $dis1 $dis2
	cd $per"sec_"$snr_cri"snr_${dis1}_${dis2}dist"
	ls 20*.ph.txt_v2 | cut -d. -f1 > event_lst_ph
	#ls 20*.am.txt_v2 | cut -d. -f1 > event_lst_am
	#python ../get_samst.py event_lst_ph event_lst_am event_lst_sam
#########################
foreach i ( ) #this plotting part has been done in last step, in script 'correct_2pi_curvature.csh'
wc 20*ph.txt_v2 | grep "0       0       0" | awk '{print "rm",$4}' > rm_0_0_0.csh
echo "removing 0_0_0 files..."
csh rm_0_0_0.csh
ls 20*ph.txt_v2 | awk -v region=$argv[1] '{print "/home/jiayi/Script/GMT/C_plot_travel",$1,region,"\n/home/jiayi/Script/GMT/C_plot_travel_T0.2",$1,region}' > plot_all.csh
echo "converting ph to HD..."
csh plot_all.csh
ls 20*.ph.txt_v2 | cut -d. -f1 > event_lst
#/home/tianye/Programs/Earthquake/travel_time_to_slow_map_v4_v3_v2to10_270_noampHD event_lst $per
#/home/tianye/Programs/Earthquake/slow_maps_to_iso_map_ani_data_v5_n50_270_noamp event_lst 0 360 18 $per"s_iso_ani_v1"
#grep -v " 0 999" $per"s_iso_ani_v1.iso"  | awk '{print $1,$2,$3}' > $per"s_iso_ani_v1.1"

ls 20*ph.txt_v2 | awk -v region=$argv[1] '{print "/home/jiayi/Script/GMT/C_plot_travel_am",$1,region}' > plot_all_am.csh
echo "converting am to HD..."
csh plot_all_am.csh
end
#############################
#set evlst = /home/jiayi/work/SC/run_eikonal/evt.txt
	set evlstph = event_lst_ph
	set evlstam = event_lst_am
	set dist_cri = 500
	echo 'removing laplace.txt.HD.....'
	#foreach event ( `more $evlst` )
		rm -f 20*_am_laplace.txt.HD
	#end
	rm -f slow_azi*
	echo 'begin to do travel_time_to_slow_map.....................'
	/home/jiayi/progs/jy/eikonal/earthquake/travel_time_to_slow_map_v4_v3_v2to10_270_noampHD_dist $evlstph  $per $dist_cri
	echo 'begin to run slow map to iso map..................'
	/home/jiayi/progs/jy/eikonal/earthquake/slow_maps_to_iso_map_ani_data_v5_n50_270_noamp_25p $evlstph  0 360 18 $per"s_iso_ani_v1"
	more ${per}s_iso_ani_v1.iso | awk '{if($3>0)print $0}' > ${per}s_iso_ani_v1.1.iso 
	echo 'begin to run amp to amp grad laplace ....'
	#/home/jiayi/progs/jy/eikonal/earthquake/amp_HD_to_amp_gradient_HD_to_amp_laplace_HD_jy $evlstam $per
	echo 'begin to run slow laplace to iso map ani ........'
	#/home/jiayi/progs/jy/eikonal/earthquake/slow_laplace_maps_to_iso_map_ani_data_v5_n50_270_noamp_scale_25p event_lst_sam 0 360 18 $per"s_laplace_iso_ani_v1" 1
	grep -v '0 9999 0' $per's_iso_ani_v1.iso' | awk '{print $1,$2,$3}' > $per's_iso_ani_v1.1'
	#grep -v " 0 9999" $per"s_laplace_iso_ani_v1_scale_1.iso"  | awk '{print $1,$2,$3}' > $per"s_laplace_iso_ani_v1.1_scale_1.iso"
	cd ..
#exit
    end
end
