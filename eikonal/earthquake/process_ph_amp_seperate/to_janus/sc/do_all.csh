#!/bin/csh
#csh FTAN_to_traveltime.csh /media/disk-2/Data_cut /home/jiayi/work/SC/Info/evt.txt /home/jiayi/work/SC/Info/station_g.lst
csh correct_2pi_curvature.csh ./region.txt .
csh tool_do_surface_slow_lapalce_all.csh ./region.txt .
