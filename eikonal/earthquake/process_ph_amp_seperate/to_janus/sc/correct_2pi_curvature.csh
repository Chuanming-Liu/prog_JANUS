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
set sta_num_cri = 50
set cur_tm_cri = 0.005
set cur_amp_cri = 0.0625 #1/16   
set cla = 30.0
set clo = 110.0
set distlst = ( 0 15000 )
foreach per ( 50 60 70   )
    foreach i(  1  )
	set j = `echo $i | awk '{print $1+1}'`
	set dis1 = ${distlst[$i]}
	set dis2 = ${distlst[$j]}
	cd $per"sec_10snr_${dis1}_${dis2}dist"
	echo "=====WORKING ON PER="$per" ================"
##########################################
#rm -f ampwrong_* pha2pi_* phawrong_* *_v1
#rm -f *ph.txt_v1.HD *_am.txt_v1.HD cur_tm_big* cur_amp_big* *.ph.txt_v2* *_am.txt_v2*
############################################
	foreach event(`ls 20*.ph.txt | cut -d. -f1`)
	
		#rm -f ampwrong_*$event*  pha2pi_*$event*  phawrong_*event*
		#echo "get v1"
		echo "-------------DO phase"	
		/home/jiayi/progs/jy/eikonal/earthquake/correct_2pi_v1_jy_ph $event $per $sta_num_cri $cla $clo
		/home/jiayi/Script/GMT/C_plot_travel $event".ph.txt_v1" $here/$argv[1]
		/home/jiayi/progs/jy/eikonal/earthquake/correct_travel_time_curvature_v1_270_jy_ph $event $per $sta_num_cri $cur_tm_cri
		rm -f `wc $event'.ph.txt_v2' | grep "0       0       0" | awk '{print $4}'`
		/home/jiayi/Script/GMT/C_plot_travel $event'.ph.txt_v2' $here/$argv[1]
		/home/jiayi/Script/GMT/C_plot_travel_T0.2 $event'.ph.txt_v2' $here/$argv[1]
		#echo "---------DO amp"
		#/home/jiayi/progs/jy/eikonal/earthquake/correct_2pi_v1_jy_am $event $per $sta_num_cri $cla $clo
		#rm -f  ${event}_ph.txt_v1.HD ${event}_am.txt_v1.HD  cur_tm_big_${event}_v1.txt cur_amp_big_${event}_v1.txt 
		#/home/jiayi/Script/GMT/C_plot_travel_amv1_jy  $event".am.txt_v1" $here/$argv[1]
		#/home/jiayi/progs/jy/eikonal/earthquake/correct_travel_time_curvature_v1_270_jy_am $event $per $sta_num_cri  $cur_amp_cri
		#rm -f `wc $event'.am.txt_v2' | grep "0       0       0" | awk '{print $4}'`
		#/home/jiayi/Script/GMT/C_plot_travel_amv2_jy $event'.am.txt_v2' $here/$argv[1]
	

	end #foreach event
	cd ..
    end #foreach i
end #foreach per
