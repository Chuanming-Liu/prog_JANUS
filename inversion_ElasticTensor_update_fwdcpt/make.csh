#----make the main code--
g++   -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp   Main_do_inversion_parallel_BS_multipleGp_Feb15_updateKeachjump_changeEtaSpace_parallel_cptLkernel.C -std=c++11 -o Main_do_inv_restoreunc_cptLkernel_v2likemis #-O3

#----make the model_avg code--
g++ -I /home/jixi7887/Tool/package/boost_1_54_0 -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4  -fopenmp -O3  CALavg_getposteria_v9_cptLkernel.C -o CALavg_getposteria_v9_restoreunc_cptLkernel_v2likemis -std=c++11

