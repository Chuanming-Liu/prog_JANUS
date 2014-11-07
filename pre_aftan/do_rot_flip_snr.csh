rm -f process.txt
#/home/jiayi/progs/jy/stack/STACK sta.lst temp_evt.lst ./STACK_ALL_EN
awk '{print $1,$2,$3}' irisInfo/station_combine.lst | sort -d > temp_sta.lst
echo "doing rotation..." >> process.txt
/home/jiayi/progs/jy/rotate_filp/lfrotate_all_64  temp_sta.lst ./stack_cv_jy_fc
mv STACK_R? /media/JAPAN/jiayi/WUS/
mv STACK_?R /media/JAPAN/jiayi/WUS/
echo "doing flip and snr..." >> process.txt
foreach ch (  TT  )
#foreach ch ( )
  echo "=======working on $ch========" >> process.txt
  echo "=======working on $ch========"
  set dir = STACK_${ch}
  set i = 0
  foreach sta ( `ls -d $dir/* | awk '{split($0,a,"'$dir'/");print a[2]}'` )
    echo '===station' $i $sta >> process.txt
    echo '===station' $i $sta 
    echo "doing split" >> process.txt
    rm -f temp_saclst
    ls $dir/$sta/*$ch | awk '{split($0,a,"'$sta'/");print a[2]}' > temp_saclst
    /home/jiayi/progs/jy/rotate_filp/flip_3000 $dir/$sta temp_saclst
#    rm -f $dir/$sta/*_snr.txt
    @ i = $i + 1
  end
  echo "doing snr " >> process.txt
  /home/jiayi/progs/jy/pre_aftan/spectral_snr_f_V2_cor_WChina $dir temp_sta.lst
  echo "end snr" >> process.txt
end
#rm -f temp.lst temp_saclst
echo "doing aftan ..." >> process.txt
python find_eT.py
#foreach ch ( RR TT RT TR )
foreach ch ( TT )
/home/jiayi/progs/jy/aftan/with_jump_no_amp/aftani_c_pgl_cor_love_with_pred_model Infosaclst_eT_snr5_T35_60_${ch}.txt
end
echo "end aftan ..." >> process.txt
