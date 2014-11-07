# this is used to plot the velocity map for a single event
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

set tBON = $BON1":."vel":"

awk '{if($3>0)print $1+360,$2,1/$3}' $input > ph.dat
awk '$7<=0' $input > b.dat
psbasemap $REG1 $JXY1 $tBON -P -K -X3 -Y20.8 > $outps
set cpt = US.20.Lov.cpt

nearneighbor -R -I0.2 -S1.2 -N4/3 ph.dat -G1.grd
grdimage 1.grd -J -R -P -O -K -C$cpt >> $outps
surface ph.dat -R -I0.2 -T0 -G2.grd
#psxy ph.dat -R -J -P -O -K -C$cpt -St0.15 -W1 >> $outps
pscoast -A5000 -R -J -P -O -K -Na -W3 -S147/224/255 >> $outps
psscale -D13.2/3.6/7/0.3 -P -O -K -C$cpt >> $outps

pwd | psxy -R -J -P -O >> $outps

