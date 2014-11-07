#!/bin/csh
if ( $#argv != 2 )then
echo "usage: do_laplace [region_inputfile] [sac_path: sac_path/evnm.per]"
exit
endif

set exe = /home/jixi7887/progs/jy/eikonal/earthquake/process_ph_amp_seperate
set snr_cri =  10
set distlst = ( 0 20)
set grid = 0.2
set dirph = $argv[2]
cd $argv[2]
foreach per ( 40 )
#foreach per (  8 20 25 30 35 40 45  )
#foreach per ( 8 10 12 14 16 18 20 22 24 25 26 28  30 32 34 35 36 38 40 42 44 45 46 48 )
    echo "\n\nworking on "$per"sec period...\n"
    foreach i ( 1  )
	set j = `echo $i | awk '{print $1+1}'`
        set dis1 = ${distlst[$i]}
        set dis2 = ${distlst[$j]}
	echo $dis1 $dis2
	set dirph = . #../$per"sec_"$snr_cri"snr_${dis1}_${dis2}dist" # the dir that contains *ph.txt* files
	##############change here ###########
	ls $dirph/*.ph.txt_v2 | cut -d/ -f2 | cut -d. -f1,2 > event_lst_ph #evnm.per
	#####################################
	set evlstph = event_lst_ph
	set dist_cri = 300
	#set dtcri = 1
	echo 'removing laplace.txt.HD.....'
	rm -f slow_azi*
	echo 'begin to do travel_time_to_slow_map.....................'
	$exe/travel_time_to_slow_map_v4_v3_v2to10_270_noampHD_jy_g0.4_4direction_sepDir $evlstph  $per $dist_cri $grid  $dirph
	echo 'begin to run slow map to iso map..................'
	$exe/slow_maps_to_iso_map_ani_data_v5_n50_270_noamp_25p_g0.4_sepDir $evlstph  -180 180 18 $per"s_iso_ani_v1" $grid $dirph
	more ${per}s_iso_ani_v1.iso | awk '{if($3>0)print $0}' > ${per}s_iso_ani_v1.1.iso 

	#====clean======
	#rm -f *.ph.txt_v1 *.ph.txt_v2 
	#rm -f *.ph.txt.HD *.ph.txt_v2.HD *.ph.txt_v2.HD_0.2
#exit
    end
end
