// input the events.lst station.lst
// output sta_dist
// with array in it.
// output eventname, arayname, staname, dist, stalon, stalat
// input arrya, stanm, stalon, stalat
// fevnet:  eventname, evelont, evelat

#define MAIN
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <fstream>
#include <istream>
#include <time.h>
#include <omp.h>
//#include <ifstream>
using namespace std;

int get_dist2(double lat1, double lon1, double lat2, double lon2, double *dist, double *azi, double *bazi);
bool fexists(const char *filename);
void get_phvel(char filename[300],double per,double *phvel_out,double *amp_temp2);


///////////////////////////////////////////////////////////////////////////////////////////
void get_phvel(char filename[300],double per,double *phvel_out,double *grvel_out)
{
  *phvel_out=0;
  *grvel_out=0;
  FILE *f1;
  double cp,ap_1,ap_2,phvel_1,phvel_2,snr,wd,gv_1,gv_2,tri,gb1,gb2;
  double amp1,amp2;
  int k,i;
  if((f1=fopen(filename,"r"))==NULL)
    {
      cout<<"no file for "<<filename<<endl;
      return;
    }
  ap_1=0;
  int flag=0;
  for(;;)
    {
      if(fscanf(f1,"%d %lf %lf %lf %lf %lf %lf",&k,&cp,&ap_2,&gv_2,&phvel_2,&amp2,&tri)==EOF) {
        //cout<<"end_of_file!!!! "<<filename<<endl;
        break;
        }
      if(ap_2>per)
        {
          if(ap_1==0) 
            {
              fclose(f1);
              return;
            }
          *phvel_out=(phvel_2-phvel_1)/(ap_2-ap_1)*(per-ap_1)+phvel_1;
          *grvel_out=(gv_2-gv_1)/(ap_2-ap_1)*(per-ap_1) + gv_1;
	  flag=1;
          //*amp_temp2=(amp2-amp1)/(ap_2-ap_1)*(per-ap_1)+amp1;
          break;
        }
      ap_1=ap_2;
      phvel_1=phvel_2;

      //amp1=amp2;
      gv_1=gv_2;
      i++;
    }
  fclose(f1);
}


//////////////////////////////////////////////////////////////////////////////////////////////
void get_snr(char filename[300],double per,double *amp,double *snr1, double *snr2)
{
  *amp=0.;
  *snr1 = 0.;
  *snr2 = 0.;
  FILE *f1;
  double cp,ap_1,ap_2,amp_1,amp_2,s1_1,s1_2,s2_1,s2_2;
  double amp1,amp2;
  int k,i;
  if((f1=fopen(filename,"r"))==NULL)
    {
      cout<<"no file for "<<filename<<endl;
      return;
    }
  ap_1=0;
  for(;;)
    {
      if(fscanf(f1,"%lf %lf %lf %lf",&ap_2,&amp_2,&s1_2,&s2_2)==EOF) {
        //cout<<"end_of_file!!!! "<<filename<<endl;
        break;
        }
      if(ap_2>per)
        {
          if(ap_1==0)
            {
              fclose(f1);
              return;
            }
          //*phvel_out=(phvel_2-phvel_1)/(ap_2-ap_1)*(per-ap_1)+phvel_1;
          //*grvel_out=(gv_2-gv_1)/(ap_2-ap_1)*(per-ap_1) + gv_1;
          *amp = (amp_2-amp_1)/(ap_2-ap_1)*(per-ap_1)+amp_1;
          *snr1 = (s1_2-s1_1)/(ap_2-ap_1)*(per-ap_1)+s1_1; 
          *snr2 = (s2_2-s2_1)/(ap_2-ap_1)*(per-ap_1)+s2_1; 

          //*amp_temp2=(amp2-amp1)/(ap_2-ap_1)*(per-ap_1)+amp1;
          break;
        }
      ap_1 = ap_2;
      s1_1 = s1_2;
      s2_1 = s2_2;
      amp_1 = amp_2;
      i++;
    }
  fclose(f1);
}
/////////////////////////////////////////////////////////////////////////////////////////


bool fexists(const char *filename)
{
  fstream ifile(filename);
  return ifile;
}

int get_dist2(double lat1, double lon1, double lat2, double lon2, double *dist, double *azi, double *bazi)
{
double pi;
pi = 4.0*atan(1.0);
double cva = 6378.137;
double cvb = 6356.7523142;
double f = 1/298.257223563;
double L = 0.00;
double jcvA, jcvB;
//if (lon1 < 0.) lon1 = lon1 + 360.;
//if (lon2 < 0.) lon2 = lon2 + 360.;
L = lon1-lon2;
//if (L > 180.000)  L =360.000000 - L;
//if (L < -180.000) L =  360.000 - abs(L);
//L = fabs(L);
double U1 = 0;
U1 = atan((1-f)*tan(lat1/180*pi));
double U2 = 0;
U2 = atan((1-f)*tan(lat2/180*pi));
double cv,cv1,cv2,cv3,cv4,cv5,cvC,numda1;
L = L*pi/180;
double numda = L;
numda1 = numda;
int cc=0;

do {
  numda = numda1;
  cv1 =  sqrt( (cos(U2)*sin(numda))*(cos(U2)*sin(numda))+ (cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(numda))*(cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(numda)) ); // cv1 sin(quan)
  cv2 = sin(U1)*sin(U2)+ cos(U1)*cos(U2)*cos(numda);
//cout<<"cv1 cv2: "<<cv1<<" "<<cv2<<endl;
  cv = atan2(cv1,cv2);
//cout<<"cv: "<<cv<<endl;
  cv3 = cos(U1)*cos(U2)*sin(numda)/sin(cv);
//cout<<cos(U1)*cos(U2)*sin(numda)<<endl;
//cout<<sin(cv)<<endl;
//cout<<"cv3: "<<cv3<<endl;
  cv4 = 1 - cv3*cv3;

if (cv4 == 0)
   cv4 = 0.0000000001;

//cout<<"cv4: "<<cv4<<endl;
  cv5 = cos(cv) - 2*sin(U1)*sin(U2)/cv4;
//cout<<"cv5: "<<cv5<<endl;
  cvC = f/16*cv4*(4 + f*(4 - 3*cv4));
  numda1 = L + (1-cvC)*f*cv3*(cv + cvC*cv1*(cv5 + cvC*cv2*(-1 +2*cv5*cv5)));
//cout<<"numda1 "<<numda1<<endl;

  cc = cc+1;
  if (cc>100000) {
     fprintf(stderr,"cannot compute : %g %g %g %g %g %g\n",lon1,lat1,lon2,lat2,numda,numda1);
     numda = numda1;
     break;
     }
} while (fabs(numda - numda1) > 0.0000000001);
//            //
double mius, cvA, cvB, deltacv,s;
mius = cv4*(cva*cva - cvb*cvb)/(cvb*cvb);
cvA = 1+mius/16384*(4096 + mius*(-768 + mius*(320 - 175*mius)));
//            ////cout<<"cvA "<<cvA<<endl;
cvB = mius/1024*(256+ mius*(-128 + mius*(74 - 47*mius)));
deltacv = cvB*cv1*(cv5 +cvB/4*(cv2*(-1 + 2*cv5*cv5)-cvB/6*cv5*(-3+4*cv1*cv1)*(-3+4*cv5*cv5) ));
//            ////cout<<"delatacv "<<deltacv<<"cvb "<<cvb<<"cv "<<cv<<endl;
s = cvb * cvA *(cv - deltacv);
//fprintf (stderr,"jcvA1: %g jcvA2 %g\n",cos(U2)*sin(numda1),cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(numda1));
//fprintf (stderr,"jcvB1: %g jcvB2 %g\n",cos(U1)*sin(numda1),-cos(U2)*sin(U1)+sin(U2)*cos(U1)*cos(numda1));
//fprintf(stderr,"check atan2: %g\n",atan2(1,-1)*180/pi);
jcvA = atan2( (cos(U2)*sin(numda1)),(cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(numda1)))*180/pi;
jcvB = atan2( (cos(U1)*sin(numda1)),(-sin(U1)*cos(U2)+sin(U2)*cos(U1)*cos(numda1)))*180/pi;
//fprintf(stderr,"jcvA ; %g jcvB: %g\n",jcvA, jcvB);

if (jcvB>=180) jcvB = jcvB-180;
else jcvB = 180 - jcvB;

if (jcvA>=180) jcvA = jcvA-180;
//else if (jcvA<-90 )
//   jcvA = jcvA+180;
 else if (jcvA <0 )
       jcvA = -jcvA;
//     else if (jcvA <=90)
//        jcvA = 180-jcvA ;
     else jcvA = 360 - jcvA;

*dist = s;
*azi = jcvA;
*bazi = jcvB;
return 1;
}

# define NSTMAX 2000
int main(int argn, char *arg[]) {
  clock_t t;
  int i,j,k,nper;
  char stanm[NSTMAX][10], ntnm[NSTMAX][10];
  double stalat[NSTMAX], stalon[NSTMAX],evlat,evlon,tevlat, tevlon;
  double tper,phv,grv,tph,tgr,amp,snr1,snr2;
  double pers[300];
  //char eventnm[10000][20];
  char teve[20],evnm[20],fdir[300],fsac[300],fdisp[300],fsnr[300],pnm[300];
  FILE *ff1, *ff2, *ffout,*fper;
  int nst = 0;
  int nev = 0;
  char chan[6];
//  double la1,lo1,la2,lo2,az,baz,dist;
//  az = 0;
//  baz = 0;
//  dist = 0;
//  double temp[4];
  // read station information
  if (argn != 8) {
    fprintf(stderr,"please input [station.lst] [event.lst] [event] [dir] [perlist] [out.file] [chan]\n");
    printf("[station.lst] stnm,stlo,stla,ntnm\n[event.lst] evnm,evlo,evla\n[perlist] per\nfsac=fdir/evnm.ntnm[i].stanm[i].chan.sac\n[out.file] output_dir=[out.file].per\n");
    return 0;
    }

  if ((ff1 = fopen(arg[1],"r"))==NULL) {
    fprintf(stderr,"cannot open file %s\n",arg[1]);
    fclose(ff1);
    return 0;
    }
  for(i=0;;i++){
	if(i>NSTMAX-1){
		printf("the defined Nstmax(%d) is not large enough! reset it!\n",NSTMAX);
		return 0;
	}
	if (fscanf(ff1,"%s %lf %lf %s",&(stanm[i][0]),&(stalon[i]),&(stalat[i]),&(ntnm[i][0]))!=4) break;
  }
  fclose(ff1);

  strcpy(evnm,arg[3]);
  strcpy(fdir,arg[4]);
  strcpy(chan,arg[7]);
  //fprintf(stderr,"read file %s ok!!\n",arg[1]);
  nst = i;

  if ((ff1 = fopen(arg[2],"r"))==NULL) {
    fprintf(stderr,"cannot open file %s\n",arg[2]);
    fclose(ff1);
    return 0;
    }
  for (i=0;;i++) {
    if (fscanf(ff1,"%s %lf %lf", &(teve[0]),&(tevlon),&(tevlat))!=3) break;
    if (strcmp(teve,evnm) != 0) continue;
    else {evlon = tevlon; evlat = tevlat;break;}
    }
  fclose(ff1);
  //fprintf(stderr,"find event %s information lon=%f lat=%f\n",evnm,evlon,evlat);

  //fprintf(stderr,"station number: %d \n",nst);

  if ((ff1 = fopen(arg[5],"r"))==NULL) {
    fprintf(stderr,"cannot open perlist file %s\n",arg[5]);
    fclose(ff1);
    return 0;
    }
  for (i=0;;i++) {
    if (fscanf(ff1,"%lf", &(pers[i]))!=1) break;
    }
  nper = i;
  //cout<<nper<<" "<<pers[0]<<endl;
  //abort();
  //ffout = fopen(arg[6],"w");
  t = clock();

  double az,baz,dist;
  double tla1,tlo1,tla2,tlo2;
  for (j=0;j<nper;j++) {
    k = 0;
    tper = pers[j];
    sprintf(pnm,"%s.%g\0",arg[6],tper);
    //fprintf(stderr,"%s T=%g\n",pnm,tper);
    ffout = fopen(pnm,"w");
    amp = -1.;snr1 = -1.;snr2 = -1.;   
    phv = -1.; grv = -1.;
    #pragma omp parallel for default(none) shared(tper,nst,fdir,evnm,ntnm,stanm,chan,stalat,stalon,evlat,evlon,ffout,k) private(fsac,fdisp,fsnr,tla1,tlo1,tla2,tlo2,phv,grv,amp,snr1,snr2,dist,az,baz,tph,tgr,i)
    for (i=0;i<nst;i++) {
      sprintf(fsac,"%s/%s.%s.%s.%s.sac",fdir,evnm,ntnm[i],stanm[i],chan);
      if(access(fsac,F_OK)!=0)continue;
      //if(access(fsac,F_OK)!=0){printf("file %s not exist!\n",fsac);continue;}
      sprintf(fdisp,"%s/%s.%s.%s.%s.sac_2_DISP.1",fdir,evnm,ntnm[i],stanm[i],chan);
      if(access(fdisp,F_OK)!=0)continue;
      //if(access(fdisp,F_OK)!=0){printf("file %s not exist!\n",fdisp);continue;}
      sprintf(fsnr,"%s/%s.%s.%s.%s.sac_amp_snr",fdir,evnm,ntnm[i],stanm[i],chan);
      if(access(fsnr,F_OK)!=0)continue;
      //if(access(fsnr,F_OK)!=0){printf("file %s not exist!\n",fsnr);continue;}

      tla1 = stalat[i];
      tlo1 = stalon[i];
      tla2 = evlat;
      tlo2 = evlon;
      //fprintf(stderr,"check %s\n",fsac) ;     

      get_phvel(fdisp,tper,&phv,&grv);
      get_snr(fsnr,tper,&amp,&snr1,&snr2);

      if (get_dist2(tla1, tlo1, tla2, tlo2, &dist, &az, &baz)==0) {
        printf("cannot compute distance!!!");
	continue;
        //return 0;
        }
      if (phv > 0. && grv > 0.) {
        tph = dist/phv;
        tgr = dist/grv;
        }
      else {tph = -1.; tgr = -1.;continue;}
      #pragma omp critical (writefile)
      {
      fprintf(ffout,"%s %s %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", evnm,ntnm[i],stanm[i],dist,az,baz,tlo1,tla1,phv,grv,tph,tgr,amp,snr1,snr2);
      k = k + 1;
      }//pragma
      }//for i<nst
  fclose(ffout);
  printf("%s, %d stations read\n",pnm,k);
  }//for j<nper
  t = clock() - t;
  //printf("the time consumed(CPU time): %g\n",((float)t)/CLOCKS_PER_SEC);
  return 1;
  }//main

