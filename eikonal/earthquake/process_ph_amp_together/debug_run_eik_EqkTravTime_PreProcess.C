// this reads in the travel time measurements from get_dist_measurements_parallel.C
// then do the correct_2pi, plot travel time/amp map and correct_time/amp curvature
//
//
//

#define MAIN
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>

#include "/home/jixi7887/progs/jy/eikonal/earthquake/dx.h"
#include "/home/jixi7887/progs/jy/eikonal/earthquake/dy.h"
#include "/home/jixi7887/progs/jy/eikonal/earthquake/process_ph_amp_together/correct_2pi_v1_jy_ph.cv.time.earth_parallel.C"
#include "/home/jixi7887/progs/jy/eikonal/earthquake/process_ph_amp_together/correct_travel_time_curvature_v1_270_c_parallel.C"
#include "/home/jixi7887/progs/jy/eikonal/earthquake/process_ph_amp_together/travel_time_to_slow_map_v4_v1_cv_us_v1_2_c_single_azi_US_parallel.C"
#include "/home/jixi7887/progs/jy/eikonal/earthquake/process_ph_amp_together/amp_HD_to_amp_gradient_HD_to_amp_laplace_HD_input_small_region_US_parallel.C"

#define NSTR 300
#define NEVTMAX 2000
#define NPERMAX 20
using namespace std;

int main(int argc, char *argv[]){
  if(argc!=12){
	printf("Input:\n1]event_list(3 column,evnm,evlo,evla)\n2]perlst (1 column, per)\n3]sta_num_cri\n4]snr_late_cri\n5]snr_precursor_cri\n6]lonmin 7]lonmax 8]latmin 9]latmax\n10]codedir (codedir/C_plot...)\n11] data_dir (datadir/evnm.per)\n");
	exit(-1);
  } 
  
  char evtfile[NSTR],perlstfile[NSTR],codedir[NSTR],datadir[NSTR];
  int stanumcri,snrcri1,snrcri2,npts_x,npts_y;
  float lonmin,lonmax,latmin,latmax ;
  float evlolst[NEVTMAX],evlalst[NEVTMAX];

  sprintf(evtfile,"%s",argv[1]);
  sprintf(perlstfile,"%s",argv[2]);
  stanumcri=atoi(argv[3]);
  snrcri1=atoi(argv[4]);
  snrcri2=atoi(argv[5]);
  lonmin=atof(argv[6]);
  lonmax=atof(argv[7]);
  latmin=atof(argv[8]);
  latmax=atof(argv[9]);
  sprintf(codedir,"%s",argv[10]);
  sprintf(datadir,"%s",argv[11]);
  npts_x=(int)((lonmax-lonmin)/0.2)+1;
  npts_y=(int)((latmax-latmin)/0.2)+1;
  

  int Nevt,i,Nper;
  char evnmlst[NEVTMAX][30];
  FILE *f1;
  int perlst[NPERMAX];

  //--read in evt name ---
  if((f1=fopen(evtfile,"r"))==NULL){
	printf("ERROR, cannot open event name file %s to read\n",evtfile);
	exit(-1);
  }

  for(i=0;;i++){
	if (i>NEVTMAX-1){
		printf("ERROR, the NEVTMAX(%d) is not big enough\n",NEVTMAX);
		exit(-1);
	}
	if(fscanf(f1,"%s %f %f",&evnmlst[i][0],&evlolst[i],&evlalst[i])!=3)break;
  }
  fclose(f1);
  Nevt=i;
  printf("test-- read in %d events\n",Nevt);

  //--read in period list ---
  if((f1=fopen(perlstfile,"r"))==NULL){
	printf("ERROR, cannot open perlst name file %s to read\n",perlstfile);
	exit(-1);
  }

  for(i=0;;i++){
	if(i>NPERMAX-1){
		printf("ERROR, the NER(%d) is not big enough\n",NPERMAX);
		exit(-1);
	}
	if(fscanf(f1,"%d",&perlst[i])!=1)break;
  }
  fclose(f1);     
  Nper=i;
  printf("test-- read in %d T\n",Nper);

  //---
  char str[500],evnm[30],datanm[500],rmstr[1000];
  int per,flaglst[NEVTMAX];
  float evlo,evla;
  for(int iper=0;iper<Nper;iper++){
	per=perlst[iper];
  	printf("test-- per  = %d\n",per);
	int count=0;
  	#pragma omp parallel for default(none) shared(flaglst,count,per,datadir,evnmlst,evlolst,evlalst,npts_x,npts_y,codedir,stanumcri,snrcri1,snrcri2,Nevt,lonmin,lonmax,latmin,latmax) private(str,evnm,datanm,rmstr,evlo,evla) num_threads(2)
	for(int ievt=0;ievt<Nevt;ievt++){
		//#pragma omp critical (countnum)
		//{
			count++;
		//}
		sprintf(evnm,"%s",evnmlst[ievt]);
		evlo=evlolst[ievt];
		evla=evlalst[ievt];
		//int flag=1;
		flaglst[ievt]=1;
		printf("--STEP1 work on %s.%d n=%d\n",evnm,per,count);
		sprintf(datanm,"%s/%s.%d",datadir,evnm,per);
		if(access(datanm,F_OK)!=0){
			printf("data file %s does not exist!\n",datanm);
			//continue;
			flaglst[ievt]=0;
		}
		if(flaglst[ievt]>0){
			flaglst[ievt]=1;
			if((do_correct_2pi(datanm,per,stanumcri,snrcri1,snrcri2))!=1){
				printf("--%s.%d correct 2pi failed\n",evnm,per);
				//continue;
				flaglst[ievt]=0;
			}//==> output datanm.input.c.txt
		}
	}//for ievt

	count=0;
  	#pragma omp parallel for default(none) shared(flaglst,count,per,datadir,evnmlst,evlolst,evlalst,npts_x,npts_y,codedir,stanumcri,snrcri1,snrcri2,Nevt,lonmin,lonmax,latmin,latmax) private(str,evnm,datanm,rmstr,evlo,evla) num_threads(2)
	for(int ievt=0;ievt<Nevt;ievt++){
			count++;
		if(flaglst[ievt]>0){
			sprintf(evnm,"%s",evnmlst[ievt]);
                	evlo=evlolst[ievt];
                	evla=evlalst[ievt];
			flaglst[ievt]=1;
			printf("--STEP2, work on %s.%d n=%d\n",evnm,per,count);
			sprintf(datanm,"%s/%s.%d",datadir,evnm,per);
			sprintf(str,"csh %s/C_plot_travel_US %s.input.c.txt %f %f %f %f\n",codedir,datanm,lonmin,lonmax,latmin,latmax);
			system(str);
			sprintf(str,"csh %s/C_plot_travel_am_US %s.input.c.txt %f %f %f %f\n",codedir,datanm,lonmin,lonmax,latmin,latmax);
			system(str);
			sprintf(rmstr,"rm -f %s.input.c.txt.tomo.grd ",datanm);
			system(rmstr);
		}//if flag
	}//for ievt
	count=0;
  	#pragma omp parallel for default(none) shared(flaglst,count,per,datadir,evnmlst,evlolst,evlalst,npts_x,npts_y,codedir,stanumcri,snrcri1,snrcri2,Nevt,lonmin,lonmax,latmin,latmax) private(str,evnm,datanm,rmstr,evlo,evla) num_threads(2)
	for(int ievt=0;ievt<Nevt;ievt++){
			count++;
		if(flaglst[ievt]>0){
			sprintf(evnm,"%s",evnmlst[ievt]);
                	evlo=evlolst[ievt];
                	evla=evlalst[ievt];
			//flaglst[ievt]=1;
			printf("--STEP3, work on %s.%d n=%d\n",evnm,per,count);
                        sprintf(datanm,"%s/%s.%d.input.c.txt",datadir,evnm,per);
			if((do_correct_curvature(datanm,per,lonmin,latmin,npts_x,npts_y))!=1){
				printf("--%s.%d correct curvature failed\n",evnm,per);
				flaglst[ievt]=0;			
			}//if
		}//if flag
	}//for ievt

	count=0;
  	#pragma omp parallel for default(none) shared(flaglst,count,per,datadir,evnmlst,evlolst,evlalst,npts_x,npts_y,codedir,stanumcri,snrcri1,snrcri2,Nevt,lonmin,lonmax,latmin,latmax) private(str,evnm,datanm,rmstr,evlo,evla) num_threads(2)
	for(int ievt=0;ievt<Nevt;ievt++){
			count++;
		if(flaglst[ievt]>0){
			sprintf(evnm,"%s",evnmlst[ievt]);
                	evlo=evlolst[ievt];
                	evla=evlalst[ievt];
			printf("--STEP4, work on %s.%d n=%d\n",evnm,per,count);
                        sprintf(datanm,"%s/%s.%d.input.c.txt",datadir,evnm,per);
			sprintf(str,"csh %s/C_plot_travel_US %s_v2 %f %f %f %f\n",codedir,datanm,lonmin,lonmax,latmin,latmax);
                printf(str);
                system(str);
                sprintf(str,"csh %s/C_plot_travel_T0.2_US %s_v2 %f %f %f %f\n",codedir,datanm,lonmin,lonmax,latmin,latmax);
                printf(str);
                system(str);
                sprintf(str,"csh %s/C_plot_travel_am_US %s_v2 %f %f %f %f\n",codedir,datanm,lonmin,lonmax,latmin,latmax);
                printf(str);
                system(str);

                sprintf(rmstr,"%s %s_v2.tomo.grd",rmstr,datanm);
		system(rmstr);
		}
	}
/*
	count=0;
  	#pragma omp parallel for default(none) shared(flaglst,count,per,datadir,evnmlst,evlolst,evlalst,npts_x,npts_y,codedir,stanumcri,snrcri1,snrcri2,Nevt,lonmin,lonmax,latmin,latmax) private(str,evnm,datanm,rmstr,evlo,evla) num_threads(2)
	for(int ievt=0;ievt<Nevt;ievt++){
			count++;
		if(flaglst[ievt]>0){
			sprintf(evnm,"%s",evnmlst[ievt]);
                	evlo=evlolst[ievt];
                	evla=evlalst[ievt];
			printf("--STEP5, work on %s.%d n=%d\n",evnm,per,count);
			if((do_travel_time_to_slow_map(evnm,evlo,evla,per,0.2,lonmin,npts_x,lonmax,npts_y,5000.,datadir,datadir))!=1){
                        	printf("--%s.%d do_traveltime_to_slow_map failed\n",evnm,per);
				flaglst[ievt]=0;
			}//==> slow_azi_evnm.per.txt.HD.2.v2
		}
	}

        #pragma omp parallel for default(none) shared(flaglst,count,per,datadir,evnmlst,evlolst,evlalst,npts_x,npts_y,codedir,stanumcri,snrcri1,snrcri2,Nevt,lonmin,lonmax,latmin,latmax) private(str,evnm,datanm,rmstr,evlo,evla) num_threads(2)
        for(int ievt=0;ievt<Nevt;ievt++){
                        count++;
                if(flaglst[ievt]>0){
                        sprintf(evnm,"%s",evnmlst[ievt]);
                        evlo=evlolst[ievt];
                        evla=evlalst[ievt];
                        printf("--STEP6, work on %s.%d n=%d\n",evnm,per,count);
                        sprintf(datanm,"%s/%s.%d.input.c.txt",datadir,evnm,per);
			char tmpstr[50];
                	sprintf(tmpstr,"%s.%d.input.c.txt",evnm,per);
			if((do_amp_to_laplace(tmpstr,per,lonmin,latmin,npts_x,npts_y,datadir,datadir,codedir))!=1){
                        	printf("--%s.%d do_amp_to_laplace failed\n",evnm,per);
				flaglst[ievt]=0;
                	}//==>evnm_am_gradx(y).txt_v2 evnm_am_laplace.txt.HD
			sprintf(rmstr,"%s %s_am_gradx.txt_v2.tomo.grd %s_am_grady.txt_v2.tomo.grd",rmstr, datanm,datanm);
                	system(rmstr);

                }
        }
*/
	
		/*
		if(flag>0){
		//#pragma omp critical (csh1) 
		//{
		sprintf(str,"csh %s/C_plot_travel_US %s.input.c.txt %f %f %f %f\n",codedir,datanm,lonmin,lonmax,latmin,latmax);
		//printf(str);
		system(str);
		sprintf(str,"csh %s/C_plot_travel_am_US %s.input.c.txt %f %f %f %f\n",codedir,datanm,lonmin,lonmax,latmin,latmax);
		//printf(str);
		system(str);
		//}
		sprintf(rmstr,"rm -f %s.input.c.txt.tomo.grd ",datanm);

		sprintf(datanm,"%s.input.c.txt",datanm);
		flag=1;
		if((do_correct_curvature(datanm,per,lonmin,latmin,npts_x,npts_y))!=1){
			printf("--%s.%d correct curvature failed\n",evnm,per);
			system(rmstr);
			//continue;
			flag=0;
		}//==> output datanm,input.c.txt_v2
		if(flag>0){
		//#pragma omp critical (csh2)
		//{
		sprintf(str,"csh %s/C_plot_travel_US %s_v2 %f %f %f %f\n",codedir,datanm,lonmin,lonmax,latmin,latmax);
		printf(str);
		system(str);
		sprintf(str,"csh %s/C_plot_travel_T0.2_US %s_v2 %f %f %f %f\n",codedir,datanm,lonmin,lonmax,latmin,latmax);
		printf(str);
		system(str);
		sprintf(str,"csh %s/C_plot_travel_am_US %s_v2 %f %f %f %f\n",codedir,datanm,lonmin,lonmax,latmin,latmax);
		printf(str);
		system(str);

		sprintf(rmstr,"%s %s_v2.tomo.grd",rmstr,datanm);
		//}
		flag=1;
		if((do_travel_time_to_slow_map(evnm,evlo,evla,per,0.2,lonmin,npts_x,lonmax,npts_y,5000.,datadir,datadir))!=1){
			printf("--%s.%d do_traveltime_to_slow_map failed\n",evnm,per);
			system(rmstr);
			//continue;
			flag=0;
		}//==> slow_azi_evnm.per.txt.HD.2.v2

		if(flag>0){
		char tmpstr[50];
		sprintf(tmpstr,"%s.%d.input.c.txt",evnm,per);
		flag=1;
		if((do_amp_to_laplace(tmpstr,per,lonmin,latmin,npts_x,npts_y,datadir,datadir,codedir))!=1){
			printf("--%s.%d do_amp_to_laplace failed\n",evnm,per);
			system(rmstr);
			//continue;
			flag=0;
		}//==>evnm_am_gradx(y).txt_v2 evnm_am_laplace.txt.HD
		if(flag>0){
		sprintf(rmstr,"%s %s_am_gradx.txt_v2.tomo.grd %s_am_grady.txt_v2.tomo.grd",rmstr, datanm,datanm);
		system(rmstr);
		}}}}}
		//}
	}//for ievt*/
  }//for iper


  return 1;
}



