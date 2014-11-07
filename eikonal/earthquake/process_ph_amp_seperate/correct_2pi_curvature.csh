#!/bin/csh
#echo on
if ($#argv != 2) then
  echo "USAGE: tool_do_AtoD [region_infile] [sac_path : sac_path/evt.per.ph.txt]"
  exit 1
endif
#/home/tianye/data_Eikonal/SAC_XR/80sec_10snr_960dis/temp
#/home/tianye/data_Eikonal/SAC_TA/80sec_10snr_960dis/#/home/tianye/data_Eikonal/SAC_XR/80sec_10snr_960dis/temp
set here = `pwd`

cd $argv[2]

set exe = /home/jixi7887/progs/jy/eikonal/earthquake/process_ph_amp_seperate
set REG = `more $argv[1]`
#echo $REG
set sta_num_cri = 8
set cur_tm_cri = 0.005
set cur_amp_cri = 0.0625 #1/16   
set cla = 40.5
set clo = 255.0 #-90.0 #-105.0
set deg = 0.2
set snrcri = 10

foreach per (  40   )
#foreach per ( 8 20 25 30 35 40 45 )
#foreach per ( 8 10 12 14 16 18 20 22 24 25 26 28  30 32 34 35 36 38 40 42 44 45 46 48 )
    foreach i(  1  )
	set j = `echo $i | awk '{print $1+1}'`
	echo "=====WORKING ON PER="$per" ================"
##########################################
############################################
	foreach event( ` ls *.$per | cut -d. -f1 ` )
	#foreach event ( 20120518020039 )
	
		echo "-------------DO phase $event"

		$exe/correct_2pi_v1_jy_ph $event.$per $per $sta_num_cri $cla $clo $snrcri #==> evt.per.ph.txt_v1
		echo "2"
		$exe/C_plot_travel_g0.4 $event.$per.ph.txt_v1 $here/$argv[1] $deg #==> evt.per.ph.txt_v1.HD
		echo "---------------$event"
		$exe/correct_travel_time_curvature_v1_270_jy_ph_g0.4 $event.$per $per $sta_num_cri $cur_tm_cri $deg #==> evt.per.ph.txt_v2
		echo "3"
		rm -f `wc $event.$per.'ph.txt_v2' | grep "0       0       0" | awk '{print $4}'`
		$exe/C_plot_travel_g0.4 $event.$per'.ph.txt_v2' $here/$argv[1] $deg #==> evt.per.ph.txt_v2.HD
		$exe/C_plot_travel_T0.2_g0.4 $event.$per'.ph.txt_v2' $here/$argv[1] $deg #==> evt.per.ph.txt_v2.HD_0.2
	end #foreach event
    end #foreach i
end #foreach per
