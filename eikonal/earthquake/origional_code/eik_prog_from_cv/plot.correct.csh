gmtset HEADER_FONT = 4
gmtset HEADER_FONT_SIZE = 14
gmtset HEADER_OFFSET = 0.07c
gmtset ANNOT_FONT_PRIMARY = 4
gmtset ANNOT_FONT_SECONDARY = 4
gmtset LABEL_FONT = 4
gmtset ANNOT_OFFSET_PRIMARY = 0.08c
gmtset ANNOT_OFFSET_SECONDARY = 0.08c
gmtset LABEL_FONT_SIZE = 12p
gmtset ANNOT_FONT_SIZE_PRIMARY = 14p
gmtset LABEL_OFFSET = 0.12c

gmtset PLOT_DEGREE_FORMAT = +ddd:mm:ss
gmtset BASEMAP_TYPE = plain

set input = $1
set outps = $input.ps
set REG1 = -R235/290/25/50
set JXY1 = -Jm0.22
set BON1 = -B5f1/5f1WSne

set tBON = $BON1":."phase":"

#awk '{if ($9>0 && $14 > 10 && $15 > 10) print $0}' $input > t.dat
awk '$7>0' $input > g.dat
awk '{print $1,$2,$3}' g.dat > ph.dat
awk '$7<=0' $input > b.dat
psbasemap $REG1 $JXY1 $tBON -P -K -X3 -Y20.8 > $outps

#awk '{print $7,$8,$11}' t.dat > ph.dat
set bb = `python ~/PROGS_64/python_scripts/compute_mean_std.rs.py ph.dat 3 | awk '{print $1-2*$2}'`
set ee = `python ~/PROGS_64/python_scripts/compute_mean_std.rs.py ph.dat 3 | awk '{print $1+2*$2}'`
makecpt -Crainbow -T$bb/$ee/50 -Z -I > t.cpt
nearneighbor -R -I0.2 -S1.2 -N4/3 ph.dat -G1.grd
#xyz2grd ph.dat -R -I0.2 -G1.grd
grdimage 1.grd -J -R -P -O -K -Ct.cpt >> $outps
surface ph.dat -R -I0.2 -T0 -G2.grd
grdcontour 2.grd -J -R -P -O -K -Ct.cpt >> $outps
psxy ph.dat -R -J -P -O -K -Ct.cpt -St0.15 -W1 >> $outps
pscoast -A5000 -R -J -P -O -K -Na -W3 -S147/224/255 >> $outps
psscale -D13.2/3.6/7/0.3 -P -O -K -Ct.cpt >> $outps

awk '{if ($5>0.) print $1,$2,$5}' g.dat > amp.dat
set tBON = $BON1":."amp":"
psbasemap $REG1 $JXY1 $tBON -P -K -O -Y-9.8 >> $outps
set bb = `python ~/PROGS_64/python_scripts/compute_mean_std.rs.py amp.dat 3 | awk '{print $1-3*$2}' | awk '{if ($1<0) print 0; else print $1;}'`
set ee = `python ~/PROGS_64/python_scripts/compute_mean_std.rs.py amp.dat 3 | awk '{print $1+3*$2}'`
set kk = `echo $bb $ee | awk '{print ($2-$1)/5.}'`
makecpt -Crainbow -T$bb/$ee/${kk} -Z > t.cpt
nearneighbor -R -I0.2 -S1.2 -N4/3 amp.dat -G1.grd
#xyz2grd amp.dat -R -I0.2 -G1.grd
grdimage 1.grd -J -R -P -O -K -Ct.cpt >> $outps
#surface amp.dat -R -I0.2 -T0 -G2.grd
#grdcontour 2.grd -J -R -P -O -K -Ct.cpt >> $outps
psxy amp.dat -R -J -P -O -K -Ct.cpt -St0.15 -W1 >> $outps
pscoast -A5000 -R -J -P -O -K -Na -W3 -S147/224/255 >> $outps
psscale -D13.2/3.6/7/0.3 -P -O -K -Ct.cpt >> $outps

pwd | psxy -R -J -P -O >> $outps

