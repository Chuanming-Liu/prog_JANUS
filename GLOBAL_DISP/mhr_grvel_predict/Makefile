#
BIN=mhr_grvel_predict_earth_v3_cv_for_sm_eqk
 
#FFLAGS =  -e -Nl30 
FFLAGS =  -e -Nl30 -ffixed-line-length-none
FC = gfortran
CFLAGS = 
 
DBG = 
$(BIN) : mhr_grvel_predict_earth_v3_cv_for_sm_eqk.o azidl.o spline.o  project.o new_path.o integr.o
	
	$(FC) $(FFLAGS) $(CFLAGS)  mhr_grvel_predict_earth_v3_cv_for_sm_eqk.o  azidl.o spline.o project.o new_path.o integr.o -o $(BIN) $(LDFLAGS) $(LDLIBS)

clean ::
	rm -f $(BIN) core.* $(FOBJS) *.o
