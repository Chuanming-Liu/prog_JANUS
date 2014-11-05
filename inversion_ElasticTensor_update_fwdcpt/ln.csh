

set dir = ../inversion_ElasticTensor
foreach prog ( string_split.C generate_Bs.C gen_random_cpp.C INITstructure_BS.h CALpara_isolay_BS_newV2L_changeEtaSpace.C CALgroup_smooth_BS.C CALmodel_LVZ_ET_BS.C ASC_rw.C BIN_rw_Love.C CALinv_isolay_rf_parallel_saveMEM_BS_updateK_eachjump_parallel.C para_avg_multiple_gp_v4.C  )
	ln -s $dir/$prog .
end
