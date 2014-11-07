
# this script reads in a station list and a period list, and
# outputs group velocity curves between every pair of stations
# on every period in the period list. 


# the following file names must be set:

# 	format: sta n lat lon

set evt = 'pre_events.lst'
set sta = 'pre_station.lst'
#more $evt
# 	format: per
#set perlist = '/home/weisen/2_PW/GLOBAL_DISP/perlist'
set perlist = './perlist'
# 	root name where gp vel maps reside. file names
#	end in '_R_per' like '_R_125'. perlist must have periods
#	for which grv has group velocity maps.

# phase velocity maps 
set grv = '/home/weisen/2_PW/GLOBAL_DISP/PHASE_VEL/smpkolya_phv'

# group  velocity maps 
# set grv = 'GROUP_VEL/int_CUB2'

#	pathfile file name

set pathfile = 'pathfile'


# mk_pathfile makes file pathfile
#	format: n1 n2 sta1 sta2 xlat1 xlon1 xlat2 xlon2

#/home/weisen/2_PW/GLOBAL_DISP/mk_pathfile_earthquake << !
#./GLOBAL_DISP/mk_pathfile_earthquake << !
#$evt
#$sta
#$pathfile
#!

# grvel_predict reads in perlist & pathfile and outputs PREDICTION_R & PREDICTION_L
# the following file name must be set:
#~levshin/bin/mhr_grvel_predict $pathfile $grv $perlist

./mhr_grvel_predict/mhr_grvel_predict_earth_v3_cv_for_sm $pathfile $grv $perlist
#/home/weisen/2_PW/GLOBAL_DISP/mhr_grvel_predict/mhr_grvel_predict_earth_v2 $pathfile $grv $perlist
