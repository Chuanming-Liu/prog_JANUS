#!/bin/csh
# this is used to run the  run_eik_EqkTravTime_PreProcess.C code
#

set year = $argv[1]
set evtlst = /home/jixi7887/work/US/irisInfo/from_cv/$year.events.lst
set datadir =  /lustre/janus_scratch/jixi7887/US/eikonal_eqk/$year
set codedir = /home/jixi7887/progs/jy/eikonal/earthquake/process_ph_amp_together
set stanumcri = 20
set snrlatecri = 5
set snrpreccri = 10
set lonmin = 235
set lonmax = 285
set latmin = 20
set latmax = 50


awk '{print $1}' $evtlst > temp_${year}.lst
set evtlst = temp_${year}.lst
set perlst = 


$codedir/run_eik_EqkTravTime_PreProcess $evtlst $perlst $stanumcri $snrlatecri $snrpreccri $lonmin $lonmax $latmin $latmax $codedir $datadir

