echo "input 1] name (for the output disp file name) 2] wave_type (L or R) 3] modelfile"
set name = $argv[1]
set wave = $argv[2]
if ( $wave == L )then
	set nm = 'lov'
else if ( $wave == R )then
	set nm = 'ray'
else
	echo "wrong input for the parameter wave"
	exit
endif
# get the tdisp96.dat file
rm -f tdisp96.dat
#tprep96 -M modelTI.d -d dist.txt -DT 1 -NPTS 100 -NMOD 1 -$wave -PMIN 5 -PMAX 100
tprep96 -M $argv[3] -d dist.txt -DT 1 -NPTS 100 -NMOD 1 -$wave -PMIN 5 -PMAX 100
# get disp file
tdisp96

mv tdisp96.cv.lov tdisp96.cv.lov_$name 
