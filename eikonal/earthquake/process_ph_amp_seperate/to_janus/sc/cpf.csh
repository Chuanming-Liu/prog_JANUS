#!/bin/csh
#link all the _v2 info to this dir


foreach per ( 25 30 40  50 60 70 )
	set dir = ${per}sec_10snr_0_15000dist
	if( ! -e $dir )then
		mkdir $dir
	else
		rm -f -r $dir
		mkdir $dir
	endif
	echo $per
	rm $dir/* -f
	ln -f ../run_eikonal_varT1/$dir/20*ph.txt ./$dir
	
end


