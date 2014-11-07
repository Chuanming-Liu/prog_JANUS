#!/bin/csh
#echo on
if ($#argv != 2) then
  echo "USAGE: tool_do_AtoD [region_infile] [sac_path : sac_path/sec_snr_dis/evt.ph.txt]"
  exit 1
endif
#/home/tianye/data_Eikonal/SAC_XR/80sec_10snr_960dis/temp
#/home/tianye/data_Eikonal/SAC_TA/80sec_10snr_960dis/#/home/tianye/data_Eikonal/SAC_XR/80sec_10snr_960dis/temp
set here = `pwd`

cd $argv[2]

set REG = `more $argv[1]`
#echo $REG
set sta_num_cri = 8
set cur_tm_cri = 0.005
set cur_amp_cri = 0.0625 #1/16   
set cla = 40.5
set clo = -105.0 #-90.0 #-105.0
set deg = 0.2
set snr = 10

set distlst = ( 0 20 )
foreach per (  40   )
#foreach per ( 8 20 25 30 35 40 45 )
#foreach per ( 8 10 12 14 16 18 20 22 24 25 26 28  30 32 34 35 36 38 40 42 44 45 46 48 )
    foreach i(  1  )
	set j = `echo $i | awk '{print $1+1}'`
	set dis1 = ${distlst[$i]}
	set dis2 = ${distlst[$j]}
	cd $per"sec_${snr}snr_${dis1}_${dis2}dist"
	echo "=====WORKING ON PER="$per" ================"
	echo $per"sec_${snr}snr_${dis1}_${dis2}dist"
##########################################
#rm -f ampwrong_* pha2pi_* phawrong_* *_v1
#rm -f *ph.txt_v1.HD *_am.txt_v1.HD cur_tm_big* cur_amp_big* *.ph.txt_v2* *_am.txt_v2*
############################################
	mkdir trash
	wc *.ph.txt | awk '{if($1<'$sta_num_cri')print "mv", $4 ,"trash/"}' | csh
	foreach event( ` ls *.ph.txt | cut -d. -f1 ` )
	
		echo "-------------DO phase $event"

		/home/jiayi/progs/jy/eikonal/earthquake/correct_2pi_v1_jy_ph $event $per $sta_num_cri $cla $clo
		/home/jiayi/Script/GMT/C_plot_travel_g0.4 $event".ph.txt_v1" $here/$argv[1] $deg
		echo "---------------$event"
		/home/jiayi/progs/jy/eikonal/earthquake/correct_travel_time_curvature_v1_270_jy_ph_g0.4 $event $per $sta_num_cri $cur_tm_cri $deg
		rm -f `wc $event'.ph.txt_v2' | grep "0       0       0" | awk '{print $4}'`
		/home/jiayi/Script/GMT/C_plot_travel_g0.4 $event'.ph.txt_v2' $here/$argv[1] $deg
		/home/jiayi/Script/GMT/C_plot_travel_T0.2_g0.4 $event'.ph.txt_v2' $here/$argv[1] $deg
	end #foreach event
	cd ..
    end #foreach i
end #foreach per
