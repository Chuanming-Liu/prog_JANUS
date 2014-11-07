

#rm -f evt_sta.lst 
set here = `pwd`
set pathfile = /data/jiayi/work/WUS/gen_email_earthquake/sacInfo_test/pathfile
#python /home/jiayi/progs/jy/eikonal/earthquake/get_evt_sta_from_pathfile_eqk.py $pathfile evt_sta.lst

set snr_cri = 10
#SAC_path/evt/evt.sta.BHZ.sac
set sacpath = /data/jiayi/work/WUS/do_eik/compare_jg

set distlst = ( 0 20 )
#foreach per (  10 12 14 16 18 20 22 24 25 26 28  30 32 34 35 36 38 40 42 44 45 46 48 )
#foreach per ( 14 20 25 30 35 ) 
#foreach per ( 25 30 35 40 45 )
foreach per (  40 )
        set dis_cri_small = ${distlst[1]}
        set dis_cri_big = ${distlst[2]}
        echo $dis_cri_small $dis_cri_big

        echo "per="$per
	#/home/jiayi/progs/jy/eikonal/earthquake/get_event_period_trvt_amp_dist_cri_2DISP1_love_WUS $per $here/evt_sta.lst $snr_cri $dis_cri_small $dis_cri_big $sacpath 
	/home/jiayi/progs/jy/eikonal/earthquake/get_event_period_trvt_amp_dist_cri_1DISP1_love_WUS $per $here/evt_sta.lst $snr_cri $dis_cri_small $dis_cri_big $sacpath 

#        /home/jiayi/progs/jy/eikonal/noise/get_event_period_trvt_amp_dist_cri_2DISP1 $per $here/evt_sta.lst $snr_cri $dis_cri_small $dis_cri_big $sacpath $cmp
end

cd $here



