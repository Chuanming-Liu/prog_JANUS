# this is used to get the fit of azimuthal anisotropy
# has the same function as previous get_1psi2psi code, but written in python and much easier to control
set dirout = fit_psi
mkdir -p $dirout
foreach psitype ( 4psi  )
foreach per ( 20 )
	echo "working on period "$per
	#set input =  ${per}sec_20snr_0_4000dist/${per}s_iso_ani_v1.ani
	set input =  ${per}sec_10snr_0_4000dist/${per}s_iso_ani_v1.ani
	#set psitype = 4psi
	set outmap = ${psitype}_${per}_map.txt
	set outpointdir =  ${psitype}_${per}_point # use NO to turn off this output
	python /home/jiayi/progs/jy/python/get_psi.py $input $outmap $outpointdir $psitype
	\mv $outmap $dirout/
	if ( -e  $dirout/$outpointdir )then
		rm -fr $dirout/$outpointdir
	endif
	mv $outpointdir $dirout/
end
end

