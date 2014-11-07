
#define MAIN
//#include "/home/yingjie/ulisse/progs/NOISE_CODA/HEAD_NOISE/mysac.h"
#include "/home/jixi7887/progs/jy/HEAD/head_noise/64_mysac.h"
//#include "/home/jiayi/progs/jy/HEAD/head_noise/64_sac_db.h"
#include "/home/jixi7887/progs/jy/HEAD/head_noise/64_koftan.h"
#include "/home/jixi7887/progs/jy/HEAD/head_noise/gl_const.h"
#include "/home/jixi7887/progs/jy/HEAD/head_noise/mymacro.h"
#include <omp.h>
#define SLEN 200000

// This version do not need prior minV and manV 
// I just set this two values to 1.5 and 4 becuase the arrival-time window  calculated
// from the velocity range should contain surface wave signals. 
// this is good for paths not too long. If the path is too long,
// the arrival time of slow velocity will be close to the end of 5000 s. 
// In this case, you have to use the orignal version with predicted minV and maxV. 


/* Modified version of spectral_snr.c code to include Misha's filter for faster narrow-band filtering for larger data sets. The program is written for a 5-50s period band. Modify minP/maxP if changes are desired. This version does not write out the narrow-band filtered waveforms (as SAC files) */

/* 
 *  CURRENTLY SET TO USE THE MULTIPLICATION FACTOR FOR CALCULATING NARROW BAND FILTERS
 */


/* FUNCTION PROTOTYPES */
SAC_HD *read_sac (char *fname, float *sig, SAC_HD *SHD, int nmax);
void write_sac (char *fname, float *sig, SAC_HD *SHD);
void filter4_(double *f1, double *f2, double *f3, double *f4, int *npow, double *dt, 
	      int *n, float seis_in[], float seis_out[]);
int get_snr(char *filename);


/*--------------------------------------------------------------------------*/
SAC_HD *read_sac (char *fname, float *sig, SAC_HD *SHD, int nmax)
/*--------------------------------------------------------------------------*/
/* function to read sac files given the name, fname. The function outputs the time signal to the pointer sig, fills the header SHD, if the signal has fewer than nmax points */
{
  FILE *fsac;

  if((fsac = fopen(fname, "rb")) == NULL) {
    printf("could not open sac file to read%s \n", fname);
    exit(1);
  }

  if ( !fsac ) {
    /*fprintf(stderr,"file %s not find\n", fname);*/
    return NULL;
  }

  if ( !SHD ) SHD = &SAC_HEADER;

  fread(SHD,sizeof(SAC_HD),1,fsac);

  if ( SHD->npts > nmax ) {
    fprintf(stderr,"ATTENTION !!! %s npts is limited to %d.\n",fname,nmax);
    SHD->npts = nmax;
  }

  fread(sig,sizeof(float),(int)(SHD->npts),fsac);
  fclose (fsac);

/*-------------  calcule de t0  ----------------*/
   {
     int eh, em ,i;
     float fes;
     char koo[9];

     for ( i = 0; i < 8; i++ ) koo[i] = SHD->ko[i];
     koo[8] = '\0';

     SHD->o = SHD->b + SHD->nzhour*3600. + SHD->nzmin*60 +
     SHD->nzsec + SHD->nzmsec*.001;

     sscanf(koo,"%d%*[^0123456789]%d%*[^.0123456789]%g",&eh,&em,&fes);

     SHD->o  -= (eh*3600. + em*60. + fes);
   /*-------------------------------------------*/}

   return SHD;
}


/*--------------------------------------------------------------------------*/
void write_sac (char *fname, float *sig, SAC_HD *SHD)
/*--------------------------------------------------------------------------*/
{
  FILE *fsac;
  int i;
  if((fsac = fopen(fname, "wb"))==NULL) {
    printf("could not open sac file to write\n");
    exit(1);
  }

  if ( !SHD ) {
    SHD = &SAC_HEADER;
  }

  SHD->iftype = (int)ITIME;
  SHD->leven = (int)TRUE;
  SHD->lovrok = (int)TRUE;
  SHD->internal4 = 6L;
  SHD->depmin = sig[0];
  SHD->depmax = sig[0];

  for ( i = 0; i < SHD->npts ; i++ ) {
    if ( SHD->depmin > sig[i] ) {
      SHD->depmin = sig[i];
    }
    if ( SHD->depmax < sig[i] ) {
      SHD->depmax = sig[i];
    }
   }

  fwrite(SHD,sizeof(SAC_HD),1,fsac);
  fwrite(sig,sizeof(float),(int)(SHD->npts),fsac);

  fclose (fsac);
}


float sig0[SLEN];
char fname[300];
/*--------------------------------------------------------------*/
int get_snr(char *fname )
/*--------------------------------------------------------------*/
{
  FILE *fp1, *fp2, *fd, *fsp;
  int k,i, iswap, cnt, cnt2, *ih, *nsamples;
  int nf = 20;
  int npow;
  float f_cent;
  static int n;
  double dt, delta, multfac, addfac;
  double maxP= 160.0; 	
  double minP= 8.0  ;
  double per[nf],f[nf],snr[nf];
  double b, e, fb, fe, step;
  double dist, minV, maxV, minT, maxT, window, signalmax, noisemax;
  double f1, f2, f3, f4;
  static float seis_in[7000],seis_out[400000];
  float hed[158];
  char name[160], fname1[500], fname2[500];
  
  /* printf("fname is %s\n\n", fname); */

  if ( read_sac (fname, sig0, &SAC_HEADER, SLEN) == NULL )
    {
      fprintf(stderr,"File %s not found\n", fname);
      return 0;
    }

  fb=1.0/maxP;
  fe=1.0/minP;
  step=(log(fb)-log(fe))/(nf-1);

/* Calculate the frequency steps */
  for(k=0;k<nf;k++)
    {
      f[k]=exp(log(fe)+k*step);
      per[k]=1/f[k];
    }

  b=SAC_HEADER.b;
  e=SAC_HEADER.e;
  dist = SAC_HEADER.dist;
  if (dist!=dist) dist=-1;
//  fprintf(stderr,"distance read from head file: %g\n",dist);

  if (dist < 0)
    {
    fprintf(stderr,"dist: %g. There is problem with the dist.\n",dist);
    return 0;
    }
  //  minV = SAC_HEADER.user4;
  // maxV = SAC_HEADER.user5;

  minV = 2.0;
  maxV = 5.0;
  minT = dist/maxV - maxP;
  maxT = dist/minV + 2*maxP;

  if(minT<b)
    minT=b;
  if( maxT> ( e - 1000) )
    maxT= (e - 1000);
  window=maxT-minT;

/* Call to Misha's filter, filter4.f */
/* ==========================================================
 * Parameters for filter4 function:
 * Input parameters:
 * f1,f2   - low corner frequences, f2 > f1, Hz, (double)
 * f3,f4   - high corner frequences, f4 > f3, Hz, (double)
 * npow    - power of cosine tapering,  (int)
 * dt      - sampling rate of seismogram in seconds, (double)
 * n       - number of input samples, (int)
 * seis_in - input array length of n, (float)
 * Output parameters:
 * seis_out - output array length of n, (float)
 * ==========================================================
 */

  cnt2=0;
  for(k=1;k<nf-1;k++)
    {

       f2=f[k+1];
       f3=f[k-1];

//      float fwd = 0.01;
//      f2 = f[k]-fwd/2.;
//      f3 = f[k]+fwd/2.;

      f_cent = 2/(f2+f3); 

/* printf("cnt %d %f\n", cnt2, f_cent); */
//      multfac=2.5;
//      f1=f2/multfac;
//      f4=f3*multfac;

      addfac= .01;
      f1=f2 - addfac;
      f4=f3 + addfac; 

/* read sac header and seismogram */
      iswap = 0;
      if((fd = fopen(fname,"r")) == NULL) {
        printf("Cannot find file %s.\n",fname);
      }
      fread(hed,sizeof(float),158,fd);
      ih = (int *)(&hed[76]);
      if(*ih > 100 || *ih < -100) iswap = 1;
      if(iswap) swapn((unsigned char *)hed,(int)(sizeof(float)),158);
      dt = hed[0];
      delta = hed[50];
      nsamples = (int *)(&hed[79]);
      fread(seis_in,sizeof(float),*nsamples,fd);
      fclose(fd);
      if(iswap) swapn((unsigned char *)seis_in,(int)(sizeof(float)),*nsamples);
      /* printf("Delta= %f, Dt= %f, Nsamples= %d\n",delta,dt,*nsamples); */
      n = *nsamples;

/* call to filter4, which ouputs the filtered seismogram as the seis_out array */
      npow=1;
      filter4_(&f1,&f2,&f3,&f4,&npow,&dt,&n,seis_in,seis_out);
  
      for(cnt=0; cnt<n; cnt++) {
        if(seis_out[cnt] > 0) {
          sig0[cnt]=seis_out[cnt];
        }
        else {
          sig0[cnt]=(-1)*seis_out[cnt];
        }
      }

/* write narrow-band filtered files - commmented out in the nowrt program */
      
//      SAC_HEADER.t1 = minT;
//      SAC_HEADER.t2 = maxT;
//      sprintf(fname2,"%s_snr.txt_%d", fname, cnt2);
//       printf("filtered file: %s\n", fname2);
//      cnt2++;
//      write_sac(fname2, seis_out, &SAC_HEADER); 
//

/* loop over sig0 to find maxima within signal and noise time periods */
      signalmax=0;
      noisemax=0;

      for(i=(int)minT;i<maxT;i++)
	{
	  if(sig0[i]>signalmax) 
	    signalmax=sig0[i];
	}

      float sum = 0.;
      for(i=(int)(maxT+500);i<e;i++)
	{
	   sum = sum+sig0[i]*sig0[i];

	}

      noisemax = sqrt(sum/(e-(int)(maxT+500)));
      snr[k]=signalmax/noisemax;      

    }

  sprintf(fname1,"%s_snr.txt",fname);
  fp2=fopen(fname1,"w");
  for(i=1;i<nf-1;i++)
    {
      fprintf(fp2,"%lf %lf\n",per[i],snr[i]);
    }
  fclose(fp2);
  return 1;
}


//SAC_DB sdb;
/*--------------------------------------------------------------*/
#define Nevtmax 20000 // ATTENTION HERE, the length of list! indicate the maximum number of events
int main (int argc, char *argv[])
/*--------------------------------------------------------------*/
{
  FILE *ff,*ff1;
  char filename[500];
  char evtnmlst[Nevtmax][30],dir[100],str[500],evtnm[30]; 
  int i,Nevtnm;
  float dump;
//check input command line arguments
  if(argc != 3) {
    fprintf(stderr,"USAGE: spectral_snr_f [rotated_dir] [evtnm.lst]\n");
    fprintf(stderr,"Req.s a list of all cut files for spectral SNR values.\n");
    exit(-1);
  }

  if((ff = fopen(argv[2], "r"))==NULL) {
    fprintf(stderr,"Cannot open file list: %s\n", argv[2]);
    exit(1);
  }
  Nevtnm=0;
  for(i=0;;i++){
	if(i>Nevtmax){printf("### Hey, the Nevtmax (%d) is not big enough!! reset it!\n",Nevtmax);exit(0);}
	if(fscanf(ff,"%s",&evtnmlst[i])==EOF)break;
	Nevtnm++;
  }  
  fclose(ff);

  printf("read in %d events\n",Nevtnm);
  sprintf(dir,argv[1]);
  #pragma omp parallel for default(none) shared(evtnmlst,dir,Nevtnm) private(i,evtnm,filename,ff1,str) num_threads(1)
  for(i=0;i<Nevtnm;i++)
    {
      //if(fscanf(ff,"%s %g %g",evtnm[i],&dump,&dump)==EOF)
      //if(fscanf(ff,"%s",&evtnm)==EOF)
	//break;
      sprintf(evtnm,"%s",evtnmlst[i]);
      printf("----working on evt %d %s\n",i,evtnm);
      sprintf(str,"ls %s/%s/*.LHT.sac > %s/%s/temp.lst\n",dir,evtnm,dir,evtnm);
      //printf("%s\n",str);
      system(str);
      sprintf(str,"%s/%s/temp.lst",dir,evtnm);
      if ((ff1=fopen(str,"r"))==NULL)
      {
	printf("Problem! cannot open cor_sac list\n");
	exit(1);
      }
      for (;;){
      
      if(fscanf(ff1,"%s",filename)==EOF) // read evt-sta.sac pair name
	break;
      if(get_snr(filename)==0) continue;
      }
      fclose(ff1);
    }

  return 1;
}
