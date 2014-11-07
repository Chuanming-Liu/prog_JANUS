#
INST_DIR = $(HOME)/bin
BIN  = spectral_snr_f_V2_cv
#BIN = spectral_snr_f_V2

fflags = -O -ffixed-line-length-none
cflags = -O

LDLIBS = -lfftw3

FFLAGS = $(DBG) $(fflags)
CFLAGS = $(DBG) $(cflags)

FC = g77
CC = gcc

DBG =
FOBJS = spectral_snr_f_V2_cv.o filter4_cv.o swapn.o
#FOBJS = spectral_snr_f_V2.o filter4.o swapn.o
#FOBJS = temp.o filter4.o swapn.o

$(BIN) : $(FOBJS)
	$(CC) $(FFLAGS) $(CFLAGS) $(FOBJS) -o $(BIN) $(LDLIBS) -m64 -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6

install :: $(BIN)
	install -s $(BIN) $(INST_DIR)

clean ::
	rm -f $(BIN) core $(FOBJS)
