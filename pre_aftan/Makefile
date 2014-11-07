#
INST_DIR = $(HOME)/bin
BIN  = spectral_snr_f_V2_WUS_eqk #spectral_snr_f_V2_cor_WChina
#spectral_snr_f_V2_30_150

fflags = -O -ffixed-line-length-none
cflags = -O -m64 -fopenmp

#LDLIBS = -lfftw3 -lgfortran
LDLIBS = -L/home/wshen/lib -lfftw3
FFLAGS = $(DBG) $(fflags) 
CFLAGS = $(DBG) $(cflags)

#FC = g77
FC = gfortran
CC = gcc

DBG =
FOBJS = spectral_snr_f_V2_WUS_eqk.o filter4.o swapn.o
#FOBJS = temp.o filter4.o swapn.o
#FOBJS = spectral_snr_f_V2.o filter4.o swapn.o
all : $(BIN)
$(BIN) : $(FOBJS)
	$(FC) $(FFLAGS) $(CFLAGS) $(FOBJS) -o $(BIN) $(LDLIBS) -lm 

install : $(BIN)
	install -s $(BIN) $(INST_DIR)

clean :
	rm -f $(BIN) core $(FOBJS)
