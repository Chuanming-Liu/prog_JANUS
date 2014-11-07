INST_DIR = $(HOME)/bin
BIN = Main_ani_rdKernel_isolay_rf_parallel_single_point
#BIN = test2
FC = gfortran 
CC = g++ 
FFLAGS = -O3 -Wall -ffixed-line-length-none
CFLAGS = -O3 -Wall -fopenmp
#LDLIBS = -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6 -L/home/jixi7887/Tool/package/boost_1_54_0
LDLIBS =  -I /home/jixi7887/Tool/package/boost_1_54_0

FOBJS = Main_ani_rdKernel_isolay_rf_parallel_single_point.o \
	./RF/theo.o ./RF/qlayer.o ./RF/four1.o

$(BIN) : $(FOBJS)
	$(CC) $(FOBJS) -o $(BIN) ${LDLIBS} -fopenmp 
#$(CC) $(FOBJS) -o $(BIN) ${LDLIBS} -fopenmp -lg2c -fno-double-underscore -fno-underscore 

Main_ani_rdKernel_isolay_rf_parallel_single_point.o : Main_ani_rdKernel_isolay_rf_parallel_single_point.C
	$(CC) -fopenmp $(LDLIBS) -O3 -c $^

clean ::
	rm -f $(BIN) core.* $(FOBJS)

		
