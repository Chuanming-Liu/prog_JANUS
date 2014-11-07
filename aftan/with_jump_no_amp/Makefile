#
INST_DIR = $(HOME)/bin
BIN  = aftani_c_pgl_cor_love_with_pred_model_WUS_eqk # aftani_c_pgl_cor_love_with_pred_model

fflags = -O2 -Wall -ffixed-line-length-none -m64
cflags = -O2 -Wall -m64

LDLIBS =  -L/home/wshen/lib  -lfftw3

FFLAGS = $(DBG) $(fflags)

CFLAGS = $(DBG) $(cflags)

#FC = g77
FC=gfortran
CC = gcc
DBG =
FOBJS = driver_c_with_pred_model_love_WUS_eqk.o   swapn.o aftanpg.o aftanipg.o ftfilt.o fmax.o taper.o \
	trigger.o tapers.o tgauss.o mspline.o pred_cur.o misc.o phtovel.o
all : $(BIN)
$(BIN) : $(FOBJS)
	$(FC) $(FFLAGS) $(CFLAGS) $(FOBJS) -o $(BIN) $(LDLIBS) 

install : $(BIN)
	install -s $(BIN) $(INST_DIR)

clean :
	rm -f $(BIN) core $(FOBJS)
