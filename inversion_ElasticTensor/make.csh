#g++ -I /home/jixi7887/Tool/package/boost_1_54_0 -I /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 Main_do_inversion_parallel.C -o test -fopenmp -O3
#g++ -I /home/jixi7887/Tool/package/boost_1_54_0 -I /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 Main_do_inversion_parallel_BS.C -o Main_do_inversion_parallel_BS -fopenmp -O3
#g++ -I /home/jixi7887/Tool/package/boost_1_54_0 -I /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 Main_do_inversion_parallel_BS_2step.C -o Main_do_inversion_parallel_BS_2step -fopenmp -O3
#g++ -I /home/jixi7887/Tool/package/boost_1_54_0 -I /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 Main_do_inversion_parallel_BS_test.C -o Main_do_inversion_parallel_BS_test -fopenmp -O3
#g++ -I /home/jixi7887/Tool/package/boost_1_54_0  -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp -O3  Main_do_inversion_parallel_BS_multipleGp_Feb15.C -o test
#g++ -I /home/jixi7887/Tool/package/boost_1_54_0  -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp -O3  Main_do_inversion_parallel_BS_multipleGp_Feb15_largec.C -o test_largec
#g++ -I /home/jixi7887/Tool/package/boost_1_54_0  -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp -O3  Main_do_inversion_parallel_BS_multipleGp_Feb15_updateK.C -o test_updateK2_2%
#g++ -I /home/jixi7887/Tool/package/boost_1_54_0  -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp -O3  Main_do_inversion_parallel_BS_multipleGp_Feb15_largec_updateK.C -o test_largec_updateK

#g++ -I /home/jixi7887/Tool/package/boost_1_54_0  -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp -O3  Main_temp.C -o Main_temp
#g++ -I /home/jixi7887/Tool/package/boost_1_54_0  -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp -O3  Main_temp_2.C -o Main_temp_2B
#g++ -I /home/jixi7887/Tool/package/boost_1_54_0  -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp -O3  Main_temp.C -o Main_temp
#g++  -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp -O3  Main_temp_2_compare_MK.C -o Main_temp_2_compare_MK -std=c++11

#g++ -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp -O3  Main_do_inversion_parallel_BS_get_prior.C  -o test_get_prior_changeEtaSpace1.1_writebin -std=c++11
#g++ -I /home/jixi7887/Tool/package/boost_1_54_0  -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp -O3  Main_do_inversion_parallel_BS_get_prior.C  -o test_get_prior_changeEtaSpace1.1_writebin3

#g++ -I /home/jixi7887/Tool/package/boost_1_54_0  -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp -O3 Main_do_inversion_parallel_BS_multipleGp_Feb15_updateK_changeEtaSpace_changeInputc.C -o Main_do_inversion_parallel_BS_multipleGp_Feb15_updateK_changeEtaSpace0.8_changeInputc2
#g++ -I /home/jixi7887/Tool/package/boost_1_54_0  -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp -O3 Main_do_inversion_parallel_BS_multipleGp_Feb15_updateK_changeEtaSpace_changeInputc.C -o test_changeEtaSpace1.1_changeInputc0_constantc

#g++ -I /home/jixi7887/Tool/package/boost_1_54_0  -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp -O3  Main_do_inversion_parallel_BS_multipleGp_Feb15_updateK_changeInputc.C -o test_updateK_changeInputEtaEQ0.9
#g++ -I /home/jixi7887/Tool/package/boost_1_54_0  -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp -O3 Main_do_inversion_parallel_BS_multipleGp_Feb15_updateK_parallel.C -o test_parall
#g++ -I /home/jixi7887/Tool/package/boost_1_54_0  -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp -O3 Main_do_inversion_parallel_BS_multipleGp_Feb15_updateK_changeEtaSpace.C -o test_updateK_changeEtaSpace1.1

#g++ -I /home/jixi7887/Tool/package/boost_1_54_0  -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp -O3 Main_do_inversion_parallel_BS_multipleGp_Feb15_updateKeachjump_changeEtaSpace.C -o test_updateKeachjump_changeEtaSpace1.1_test

#g++ -I /home/jixi7887/Tool/package/boost_1_54_0  -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp -O3 Main_do_inversion_parallel_BS_multipleGp_Feb15_updateKeachjump_changeEtaSpace_parallel.C -o test_updateKeachjump_changeEtaSpace1.1_parall_3_dy
#g++  -O3 -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp   Main_do_inversion_parallel_BS_multipleGp_Feb15_updateKeachjump_changeEtaSpace_parallel.C -o test_updateKeachjump_changeEtaSpace1.1_parall_Apr25_vpvsCM1.5_restoreUnc -std=c++11  # this is the code used to generate Jul paper.
g++  -O3 -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp   Main_do_inversion_parallel_BS_multipleGp_Feb15_updateKeachjump_changeEtaSpace_parallel.C -o test_updateKeachjump_changeEtaSpace1.1_parall_Apr25_vpvsCM1.5_restoreUnc_weitLoveMore4 -std=c++11  
#g++  -O3 -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp   Main_do_inversion_parallel_BS_multipleGp_Feb15_updateKeachjump_changeEtaSpace_parallel.C -o test_updateKeachjump_changeEtaSpace1.1_parall_Apr25_vpvsCM1.5_restoreUnc_RAamp_rmRAPpos -std=c++11  
#g++  -O3 -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp   Main_do_inversion_parallel_BS_multipleGp_Feb15_updateKeachjump_changeEtaSpace_parallel.C -o test3kacc20kitter -std=c++11 # for the Tibet inversion
#g++  -O3 -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp   Main_do_inversion_parallel_BS_multipleGp_Feb15_updateKeachjump_changeEtaSpace_parallel.C -o HVtest -std=c++11 # for the HVtest
#g++  -g  -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4   -fopenmp Main_do_inversion_parallel_BS_multipleGp_Feb15_updateKeachjump_changeEtaSpace_parallel_gdb.C -o test_updateKeachjump_changeEtaSpace1.1_parall_Apr25_vpvsgdb2 -std=c++11

#g++ -I /home/jixi7887/Tool/package/boost_1_54_0  CALavg_getposteria_v3_prior.C -o CALavg_getposteria_v3_prior
#g++ -I /home/jixi7887/Tool/package/boost_1_54_0  CALavg_getposteria_v7.C -o CALavg_getposteria_v7_idphi8
#g++ -I /home/jixi7887/Tool/package/boost_1_54_0 -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4  -fopenmp -O3  CALavg_getposteria_v8.C -o CALavg_getposteria_v8_restoreunc #  this is the code used to generate Jul paper.
#g++ -I /home/jixi7887/Tool/package/boost_1_54_0 -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4  -fopenmp -O3  CALavg_getposteria_v8.C -o CALavg_getposteria_v8_restoreunc_0_100_1km  # get model at every 1km
#g++ -I /home/jixi7887/Tool/package/boost_1_54_0 -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4  -fopenmp -O3  CALavg_getposteria_v8.C -o CALavg_getposteria_v8_test



#g++ -I /home/jixi7887/Tool/package/boost_1_54_0 -O3 -I  /home/jixi7887/Tool/C++Eigen/eigen-3.1.4 -fopenmp Main_computeDisp_BS.C -o Main_computeDisp_BS -std=c++11
