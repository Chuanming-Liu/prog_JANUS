// this version, mainly change the way of reading/storing data, put data in sub dir, e.g., dir/evnm/... instead of dir/...
#define MAIN
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include "/home/jixi7887/progs/jy/useful/string_split.C"
using namespace std;

#define NEVTMAX 10000
#define NEVTMAX2 6000
#define NPERMAX 20
// RIGHT NOW THIS CODE CANNOT USE PARALLEL, SO DON"T USE -fopenmp WHILE COMPILING

// NEVTMAX is the max # of read in event
// NEVTMAX2 is the max # of events with reasobale slowness for each grip point
int main(int na, char *arg[])
{
  if(na!=14)
    {
      cout<<"usage:travel_time_to_velocity_map \n1]event.lst (1 column evnm) 2]azi_min 3]azi_max 4]N_bin 5]out_name 6]dx 7]lonmin 8]lonmax 9]latmin 10]latmax 11]perlst 12]indir 13]outdir"<<endl;
      return 0;
    }
  FILE *ff,*fin,*fout,*file1,*file_iso,*file_ani;
  int i,j,k;
  int npts_x,npts_y;
  char buff1[300],name_iso[300],name_ani[300],name_ani_n[300],indir[200],outdir[200];
  double lat,lon,t_lat,t_lon,radius,pi;
  int perlst[NPERMAX],per;
  int t_i,t_j,nsta;
  int ii,jj,kk,kkk,min_n;
  int marker_i,marker_j,Nper;
  double min,max,d_bin,dump1,dump2,dump3,dump4,dump5;
  int N_bin;
  fstream mff;
  string line;
  vector<string> v;
  min=atof(arg[2]);
  max=atof(arg[3]);
  N_bin=atoi(arg[4]);
  sprintf(indir,"%s",arg[12]);
  sprintf(outdir,"%s",arg[13]);
  //  cout<<min<<" "<<max<<" "<<N_bin<<endl;
  double hist[N_bin];
  double slow_sum1[N_bin];
  double slow_un[N_bin];
  double azi_w2[N_bin], azi_std[N_bin];
  double azi_weight_sum[N_bin];

  d_bin=(max-min)/N_bin;



  radius=6371.1391285;
  pi=4.0*atan(1.0);
  double dx,dy,x0,y0,x1,y1,temp,lat_temp,temp2,temp3,temp4,trash1,trash2;
//  npts_x=63;
//  npts_y=61;

  dx=atof(arg[6]);//degree
  dy=atof(arg[6]);//degree
  x0=atof(arg[7]);
  x1=atof(arg[8]);
  y0=atof(arg[9]);
  y1=atof(arg[10]);

  if((file1=fopen(arg[11],"r"))==NULL){
	printf("cannot open perlst file %s to read!\n",arg[11]);
	exit(0);
  }
  for(i=0;;i++){
	if(fscanf(file1,"%d",&perlst[i])!=1)break;
	printf("test--i=%d per =%d\n",i,perlst[i]);
  }	
  Nper=i;
  fclose(file1);

  char event_name[300],event_namelst[NEVTMAX][300];
  int Nevt,ievt;

  if((file1=fopen(arg[1],"r"))==NULL){
	printf("Cannot open evlst file %s to read\n",arg[1]);
	exit(0);
  }
  for(i=0;;i++)
  {
	if(i>NEVTMAX-1){
		printf("ERROR,NEVTMAX (%d) is not big enough!\n",NEVTMAX);
		exit(0);
	}
      if(fscanf(file1,"%s",&event_namelst[i][0])==EOF)break;
  }
  Nevt=i;
fclose(file1);


dx=atof(arg[6]);//degree
dy=atof(arg[6]);//degree
npts_x=int((x1-x0)/dx+1);
npts_y=int((y1-y0)/dy+1);

//printf("%f,%f,%f,%f,%d,%d",x0,x1,y0,y1,npts_x,npts_y);
//return 1;

  fprintf(stderr,"Memory check!!\n");
  double slow[npts_x][npts_y][NEVTMAX2],slow_weight[npts_x][npts_y][NEVTMAX2];
  double azi[npts_x][npts_y][NEVTMAX2];
  double weight[NEVTMAX2];
  double weight_sum;
  int n[npts_x][npts_y];
  double slow_sum[npts_x][npts_y],slow_std[npts_x][npts_y];

  
  fprintf(stderr,"Memory enough!!\n");
  
printf("read in %d periods, %d events\n",Nper,Nevt);
for(int iper=0;iper<Nper;iper++){
  per=perlst[iper];
  printf("working on period = %d\n",per);

  //--initialize---
  sprintf(name_iso,"%s/%s.%d.iso",outdir,arg[5],per);
  sprintf(name_ani,"%s/%s.%d.ani",outdir,arg[5],per);
  sprintf(name_ani_n,"%s/%s.%d_ani_n",outdir,arg[5],per);  
  weight_sum=0.;
  for(i=0;i<NEVTMAX2;i++){weight[i]=0.;}

  for(i=0;i<npts_x;i++)
    {
      for(j=0;j<npts_y;j++)
	{
	  slow_sum[i][j]=0;
	  slow_std[i][j]=0;
	  n[i][j]=0;
 	  for(k=0;k<NEVTMAX2;k++){
		slow[i][j][k]=0.;
		slow_weight[i][j][k]=0.;
		azi[i][j][k]=0.;
	  }
	}
    }
   //---

  //nsta=0;

  #pragma omp parallel for default(none) shared(per,indir,Nevt,event_namelst,slow,azi,n,slow_weight,x0,y0,dx,dy,x1,y1,min,max) private(ievt,event_name,fin,buff1,lon,lat,temp,temp2,i,j,mff,line,v)
  for(ievt=0;ievt<Nevt;ievt++)  // 1)for each point i j, count how many times(events),N, it is covered, n[i][j]=N. Record all the info for each point.
    {
	//int flag=1;
  	char year[10];
	sprintf(event_name,"%s",event_namelst[ievt]);
      	sprintf(year,"%c%c%c%c",event_name[0],event_name[1],event_name[2],event_name[3]);
	printf("test--- %d,%s --> %s\n",ievt,event_name,year);
      //if(fscanf(file1,"%s",&event_name)==EOF)
//      if(fscanf(file1,"%s %lf %lf %lf %lf",&event_name,&dump1,&dump2,&dump3,&dump4)==EOF)
        //break;
      //nsta++;
      //################# the name of read in file #############
      //sprintf(buff1,"%s/%s/slow_azi_%s.%d.txt.HD.2.v2",indir,year,event_name,per);
      sprintf(buff1,"%s/%s/%s/slow_azi_%s.%d.txt.HD.2.v2",indir,year,event_name,event_name,per); // changed Sep 23
      //sprintf(buff1,"%s/%d/slow_azi_%s.%d.txt.HD.2.v2",indir,per,event_name,per);
///lustre/janus_scratch/wshen/VV_HELMHOLTZ/TOMO/24/slow_azi_20130527110608.24.txt.HD.2.v2
      //cout<<buff1<<endl;
      mff.open(buff1);
      if ( not mff.is_open()){
		//cout<<buff1<<"does not exist!\n exit!!\n";
		continue;}
      //printf("test--- %d,%s --> %s\n",ievt,event_name,year);
      //printf("passed\n");

      while(getline(mff,line)){
		v.clear();
		Split(line,v," ");
		if(v.size()!=9)continue;
		lon=atof(v[0].c_str());
		lat=atof(v[1].c_str());
		temp=atof(v[2].c_str());
		temp2=atof(v[3].c_str());
		//printf("%f %f %f %f\n",lon,lat,temp,temp2);
	 i=int((lon-x0)/dx+0.1);
	 j=int((lat-y0)/dy+0.1);
	 
	 if(temp<0.5&&temp>0.1)
	   {
	     if(lon>x1+0.01||lon<x0-0.01||lat>y1+0.01||lat<y0-0.01){
		printf("Hey, out of range location, %f %f (%f,%f) (%f,%f)\n",lon,lat,x0,x1,y0,y1);
	     	continue;}
	     if(temp2<min or temp2>max){
		printf("Hey, out of range value! %s: (%f,%f) %f (%f,%f)\n",buff1,lon,lat,temp2,min,max);
		exit(0);
	     }
	     slow[i][j][n[i][j]]=temp;
	     
	     azi[i][j][n[i][j]]=temp2;
	     
	     slow_weight[i][j][n[i][j]]=1;
	     n[i][j]++;
	     if(n[i][j]>NEVTMAX2-1){
		printf("Hey, the NEVTMAX2(%d) is not big enough, there are too many events for each grid point! Reset it!\n",NEVTMAX2);
		exit(0);

	    }
	     //	     slow_sum[i][j]+=temp;
	   }

	}
      mff.close();
    }//for ievt

printf("test-- hey, finish step1\n");
  //fclose(file1);
  if((file_iso=fopen(name_iso,"w"))==NULL){
	printf("Cannot open file %s to write!\n",name_iso);
	exit(0);
  }
  nsta=50;  // PARAMETER : NSTA This nsta is different from that above
  double w2;
  double temp_slow_sum;
  int temp_n;
  #pragma omp parallel for default(none) shared(slow,slow_weight,npts_x,npts_y,n,nsta,azi,weight,slow_sum,slow_std,x0,y0,dx,dy) private(i,j,file_iso,w2,weight_sum,k,kk,temp_slow_sum,temp,temp_n)
  for(i=0;i<npts_x;i++)
    {
      for(j=0;j<npts_y;j++)
	{
	  if(n[i][j]<0.5*nsta) // 2) if the point i,j is coverd by different events less than 50%nsta, write 0 9999
	    {
	      #pragma omp critical (writeiso)
	  	{
	      fprintf(file_iso,"%lf %lf 0 9999 %d\n",x0+i*dx,y0+j*dy,n[i][j]);
	    	}
	      continue;
	  }
	  w2=0;
	  weight_sum=0;
	  for(k=0;k<n[i][j];k++)
	    {
//  among all events that covers point i j, for each event k, count how may events are there within 25 degree,N[k] => weight[k]=1/N[k]*slow_weightijk=1/N[k];
	      weight[k]=0;
	      for(kk=0;kk<n[i][j];kk++) 
		{
		  if(fabs(azi[i][j][kk]-azi[i][j][k])<25)
		    weight[k]++;
		}
	      weight[k]=1/weight[k]*slow_weight[i][j][k];  //???? HOW ABOUT IF WEIGTH[K] i.e. N[K] == 0??
	      weight_sum+=weight[k];

	    }
	  for(k=0;k<n[i][j];k++) // 4) for each points, get the weighted slowness_SUM...=> get weighted std
	    {
	      weight[k]=weight[k]/weight_sum;
	      slow_sum[i][j]+=weight[k]*slow[i][j][k];
	      w2+=weight[k]*weight[k];
	    }

	  temp_slow_sum=slow_sum[i][j];
	  temp=0;
	  for(k=0;k<n[i][j];k++)
	    {
	      temp+=weight[k]*(slow[i][j][k]-slow_sum[i][j])*(slow[i][j][k]-slow_sum[i][j]);
	    }
	  slow_std[i][j]=sqrt(temp/(1-w2)); // 1st std

	  w2=0;
	  weight_sum=0;
	  slow_sum[i][j]=0;
	  temp_n=0;
	  for(k=0;k<n[i][j];k++)
	    {
	      if(fabs(slow[i][j][k]-temp_slow_sum)>2.0*slow_std[i][j]) // 5) if the slowness is outside 2 std, ignore it => get new series of weigt[k], and get weight_sum
		continue;
	      weight_sum+=weight[k];
	      temp_n++; // number of events that have std within 2*1st_std
	    }

	  for(k=0;k<n[i][j];k++)
	    {
	      if(fabs(slow[i][j][k]-temp_slow_sum)>2.0*slow_std[i][j]) // 6) ignore events outside 2*1st_std => based on the new weight_sum got in step 5), normalize the weght[k], => get new slow_sum  => get newly weighted std -- 2ed_std
		continue;
	      weight[k]=weight[k]/weight_sum;
	      slow_sum[i][j]+=weight[k]*slow[i][j][k];
	      w2+=weight[k]*weight[k];
	    }
	  temp=0;
	  for(k=0;k<n[i][j];k++)
	    {
	      if(fabs(slow[i][j][k]-temp_slow_sum)>2.0*slow_std[i][j])
		continue;
	      temp+=weight[k]*(slow[i][j][k]-slow_sum[i][j])*(slow[i][j][k]-slow_sum[i][j]);
	    }
	  slow_std[i][j]=sqrt(temp/(1-w2)); // 2ed std

	  temp=slow_std[i][j]*sqrt(w2)/slow_sum[i][j]/slow_sum[i][j]; // ??? WHAT IS THIS ?????

	  #pragma omp critical (writeiso2)
	  {
	  fprintf(file_iso,"%lf %lf %lf %lf %d\n",x0+i*dx,y0+j*dy,1/slow_sum[i][j],temp,temp_n); // the 1/slow_sum is the velocity with all the events outside 2*1st_std removed, temp is related to 2ed_std
	  }

	}
    }
  fclose(file_iso);
  printf("test--- finish step2\n");

  if((file_ani=fopen(name_ani,"w"))==NULL){
	printf("Cannot open file %s to write!\n",name_ani);
	exit(0);
  }
  if((fout=fopen(name_ani_n,"w"))==NULL){
	printf("Cannot open file %s to write!\n",name_ani_n);
	exit(0);
  }


  double test_lon,test_lat;

  int nstep1, nstep2,npoint;
  nstep1 = 6; // the region for average ani would be a rectangle with length (grid*nstep1*2)
  nstep2 = 1; // the points used to average ani would be (nstep1*2/nstep1+1)**2
  npoint = 100;
  for(i=0;i<npts_x;i++)
    {
      for(j=0;j<npts_y;j++)
        {

	  for(k=0;k<N_bin;k++) // ???? N_bin ????
	    {
	      hist[k]=0;
	      slow_sum1[k]=0;
	      slow_un[k]=0;
	      azi_weight_sum[k]=0;
	      azi_w2[k]=0;
	      azi_std[k]=0;
	    }
	  if(i-nstep1<0||i+nstep1>=npts_x||j-nstep1<0||j+nstep1>=npts_y) // ignore few points on the edge
	    continue;
	  kkk=0;
	  //	  min_n=999999999;
	  for(ii=i-nstep1;ii<=i+nstep1;ii+=nstep2) // get the sum of event number that covering 9 points around point i j => kkk
	    {
	      for(jj=j-nstep1;jj<=j+nstep1;jj+=nstep2)
		{
		  if(n[ii][jj]<nsta*0.5) //nsta was set previously
		    continue;
		  kkk+=n[ii][jj];

		}
	    }
	  fprintf(fout,"%lf %lf %d\n",x0+i*dx,y0+j*dy,kkk);
	  if(kkk<npoint*nsta*0.5||n[i][j]<nsta*0.5) // if point i j itself has event number <50%nsta ,or kkk <9*50%nsta, ignore this point. Otherwise, for the 9 quasi-nearby points, compute the total azimuthal distribution hist[], ???slow_sum1distribution
	    continue;
	  for(ii=i-nstep1;ii<=i+nstep1;ii+=nstep2)
	    {
	      for(jj=j-nstep1;jj<=j+nstep1;jj+=nstep2)
		{
		  if(n[ii][jj]<nsta*0.5)
		    continue;
		  for(k=0;k<n[ii][jj];k++)
		    {
		      if(azi[ii][jj][k]>max||azi[ii][jj][k]<min) // azi outside [azi_min,azi_max]
			{
			  fprintf(stderr,"out of range!! %f (%f,%f)\n",azi[ii][jj][k],min,max);
			  //continue;
			  return 1;
			}
		      //		      fprintf(stderr,"oliver!! %d %d\n",ii,jj);
		      hist[int((azi[ii][jj][k]-min)/d_bin)]++;
		      azi_weight_sum[int((azi[ii][jj][k]-min)/d_bin)]+=slow_weight[ii][jj][k];
		      azi_w2[int((azi[ii][jj][k]-min)/d_bin)]+=slow_weight[ii][jj][k]*slow_weight[ii][jj][k];
		      slow_sum1[int((azi[ii][jj][k]-min)/d_bin)]+=(slow[ii][jj][k]-slow_sum[ii][jj])*slow_weight[ii][jj][k];

		    }
		}
	    }

	  kk=0;
	  
	  for(k=0;k<N_bin;k++) // for those total azi(9 quasi-nearby points), azi is sperated into N_bin parts, with more than 10 events, normalize the azi_w2 and slow_sum1 ???? kk the number of azi-quadrant with more than 10 events
	    {
	      if(hist[k]>=10)
		{
		  kk++;
		  azi_w2[k]=azi_w2[k]/azi_weight_sum[k]/azi_weight_sum[k];
		  slow_sum1[k]=slow_sum1[k]/azi_weight_sum[k];
		}
	    }

	  for(ii=i-nstep1;ii<=i+nstep1;ii+=nstep2) // compute the azimuth-distribution of std got from the 9 quasi-nearby points who has events>50%nsta THE STD OF WHAT??? 
            {
              for(jj=j-nstep1;jj<=j+nstep1;jj+=nstep2)
                {
                  if(n[ii][jj]<nsta*0.5)
                    continue;
                  for(k=0;k<n[ii][jj];k++)
                    {
                      if(azi[ii][jj][k]>max||azi[ii][jj][k]<min)
                        {
                          //fprintf(stderr,"out of range!!");
			  fprintf(stderr,"out of range!! %f (%f,%f)\n",azi[ii][jj][k],min,max);
                          return 1;
                        }
		      azi_std[int((azi[ii][jj][k]-min)/d_bin)]+=slow_weight[ii][jj][k]/azi_weight_sum[int((azi[ii][jj][k]-min)/d_bin)]*(slow[ii][jj][k]-slow_sum[ii][jj]-slow_sum1[int((azi[ii][jj][k]-min)/d_bin)])*(slow[ii][jj][k]-slow_sum[ii][jj]-slow_sum1[int((azi[ii][jj][k]-min)/d_bin)]);
                    }
                }
            }


	  
	  fprintf(file_ani,"%lf %lf %d\n",x0+i*dx,y0+j*dy,kk); // kk of each points, kk--number of azi-quadrant with >10 events
	  for(k=0;k<N_bin;k++) // for the azi-quadrant with events>10, compute the azimuthal distribution of azi_std => get temp ??????
	    {
	      if(hist[k]>=10)
		{
		  
		  azi_std[k]=sqrt(azi_std[k]/(1-azi_w2[k]));
		  temp=azi_std[k]*sqrt(azi_w2[k])/slow_sum[i][j]/slow_sum[i][j];
		 
		  fprintf(file_ani,"%lf %lf %lf\n",min+(0.5+k)*d_bin,1/(slow_sum[i][j]+slow_sum1[k]),temp);
		}
	    }

	}
    }
  
  fclose(file_ani);
  fclose(fout);
  printf("test--, finish Step3\n");
}//for iper
  return 0;
}
