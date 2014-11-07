#!/bin/csh
# get the sec_snr/evt.ph.txt

if ( $#argv != 3)then
echo "USAGE: xxx.csh 1[SAC_path] (#:SAC_path/evt/evt.sta.BHZ.sac) 2[evt.txt] 3[station.lst]"
exit
endif

set here = `pwd`
python /home/jiayi/progs/jy/eikonal/get_evt_sta.py $argv[2] $argv[3] #output evt_sta.txt

set sacpath = $argv[1]
#cd $sacpath
set snr_cri = 10
set distlst = ( 0 15000 )
foreach per (  50 60 70  )
#foreach per( 40 )
	foreach i ( 1 )
		set  j = `echo $i | awk '{print $1+1}'`
		set dis_cri_small = ${distlst[$i]}
		set dis_cri_big = ${distlst[$j]}
		echo $dis_cri_small $dis_cri_big
	
	echo "per="$per 
	#/home/jiayi/progs/jy/eikonal/earthquake/get_event_period_trvt_amp_dist_cri $per $here/evt_sta.lst $snr_cri $dis_cri_small $dis_cri_big $sacpath
	/home/jiayi/progs/jy/eikonal/earthquake/get_event_period_trvt_amp_dist_cri_2DISP1 $per $here/evt_sta.lst $snr_cri $dis_cri_small $dis_cri_big $sacpath
	end
end

cd $here
