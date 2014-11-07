#define MAIN
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#ifdef _OPENMP
 # include <omp.h>
#else
 #define omp_get_thread_num() 0
#endif
#define NEVTMAX 11000
// this version combine the eqk (slow_azi_stnm*.txt) and the noise (slow_azi_stnm*) together to genrate the iso ani results
// this version, mainly change the way of reading/storing data, put data in sub dir, e.g., dir/evnm/... instead of dir/...
//
using namespace std;

int main(int na, char *arg[])
{
  if(na!=14)
    {
      cout<<"usage:travel_time_to_velocity_map eqk_event.lst min max N_bin out_name x0 y0 npts_x npts_y per noise_event(sta).lst dir_to_eqk_vel dir_to_noise_slow_azi"<<endl;
      return 0;
    }
  FILE *ff,*fin,*fout,*file1,*file_iso,*file_ani;
  int i,j,k;
  int npts_x,npts_y;
  char buff1[300],sta1[10],name_iso[100],name_ani[100],name_ani_n[100],direqk[200],dirnoise[200];
  double lat,lon,t_lat,t_lon,radius,pi,sta1_lon,sta1_lat;
  int t_i,t_j,nsta;
  int ii,jj,kk,kkk,min_n;
  double min,max,d_bin;
  int N_bin;
  int per;
  min=atof(arg[2]);
  max=atof(arg[3]);
  N_bin=atoi(arg[4]);
  sprintf(name_iso,"%s.iso",arg[5]);
  sprintf(name_ani,"%s.ani",arg[5]);
  sprintf(name_ani_n,"%s_ani_n",arg[5]);  
  per = atoi(arg[10]);
  sprintf(direqk,arg[12]);
  sprintf(dirnoise,arg[13]);
  //  cout<<min<<" "<<max<<" "<<N_bin<<endl;
  double hist[N_bin];
  double slow_sum1[N_bin];
  double slow_un[N_bin];
  double azi_w2[N_bin], azi_std[N_bin];
  double azi_weight_sum[N_bin];
  //  for(i=0;i<N_bin;i++)
  //{
  //  hist[i]=0;
  //  slow_sum[i]=0;
  //  slow_un[i]=0;
  //}
  d_bin=(max-min)/N_bin;



  radius=6371.1391285;
  pi=4.0*atan(1.0);
  double dx,dy,x0,y0,x1,y1,temp,lat_temp,temp2,temp3,temp4,trash1,trash2;
//  npts_x=176;
  npts_x=atoi(arg[8]);
//  npts_y=126;
  npts_y=atoi(arg[9]);
  fprintf(stderr,"Memory check!!\n");
  double slow[npts_x][npts_y][NEVTMAX],slow_weight[npts_x][npts_y][NEVTMAX];
  double azi[npts_x][npts_y][NEVTMAX];
  double weight[NEVTMAX];
  double weight_sum;
  int n[npts_x][npts_y],neqk[npts_x][npts_y];
  double slow_sum[npts_x][npts_y],slow_std[npts_x][npts_y],check_sum[npts_x][npts_y];
  double trash3,trash4,trash5;
  char tstr[300];
  // double dx_km[npts_y],dy_km;
  
  fprintf(stderr,"Memory enough!!\n");
  
  dx=0.2;//degree
  dy=0.2;//degree
//  x0=235;
  x0=atof(arg[6]);
  if(x0<0)x0=x0+360.;
//  y0=25;
  y0=atof(arg[7]);
  x1=x0+(npts_x-1)*dx;
  y1=y0+(npts_y-1)*dy;
  //for(j=1;j<npts_y-1;j++)
  //  {
  //  lat_temp=y0+j*dy;
  //  lat_temp=atan(0.993277*tan(lat_temp/180*pi))*180/pi;
  //  dx_km[j]=radius*sin((90-lat_temp)/180*pi)*dx/180*pi;
  //}
  //dy_km=radius*dy/180*pi;
  for(i=0;i<npts_x;i++)
    {
      for(j=0;j<npts_y;j++)
	{
	  slow_sum[i][j]=0;
	  n[i][j]=0;
	  neqk[i][j]=0;
	}
    }
  char event_name[300];
  char events_lst[NEVTMAX][20];
  int marker_i,marker_j;
  marker_i=int((245.0-x0+0.001)/dx);
  marker_j=int((40.0-y0+0.001)/dy);
  if((file1=fopen(arg[1],"r"))==NULL){//read in event name list
	printf("Cannot open file %s to read\n",arg[1]);
	exit(0);
  }
  nsta=0;
  int cv1,cv2,neveqk,nevnoise;
  printf("ok here!\n");
  cv2 = 0;  
  for (;;) {
    	if(fscanf(file1,"%s",&(events_lst[cv2][0]))==EOF) break;
    	cv2++;
  	if(cv2>NEVTMAX){
		printf("Hey, event number is larger than the NEVTMAX %d>%d\n",cv2,NEVTMAX);
		exit(0);
  	}
    }
  fclose(file1);
  neveqk = cv2;
  //read in the noise event(station)
  if((file1=fopen(arg[11],"r"))==NULL){
	printf("Cannot open file %s to read\n",arg[11]);
	exit(0);
  }
  cv2=0;
  for(;;){
	if(fscanf(file1,"%s",&(events_lst[cv2+neveqk][0]))==EOF)break;
	cv2++;
  	if(neveqk+cv2>NEVTMAX){
		printf("Hey, event+station number is larger than the NEVTMAX %d>%d\n",neveqk+cv2,NEVTMAX);
		exit(0);
  	}
  }
  nevnoise=cv2;
  fclose(file1);
  printf("ok here!\n");
 
  cv2 = 0;
  char year[10];
  // read in the vel file of earthquake
  for(ii=0;ii<neveqk;ii++){
	cv2=cv2+1;
	sprintf(event_name,"%s\0",events_lst[ii]); //#2005/slow_azi_20051230182646.50.txt.HD.2.v2
	sprintf(year,"%c%c%c%c",event_name[0],event_name[1],event_name[2],event_name[3]);
	sprintf(buff1,"%s/%s/%s/slow_azi_%s.%d.txt.HD.2.v2",direqk,year,event_name,event_name,per); // changed Sep25, 2014
	
      if((fin=fopen(buff1,"r"))==NULL)
       {
         //cout<<buff1<<" not exist!!"<<endl;
         //return 1;
         continue;
       }
         //cout<<buff1<<" do exist!!"<<endl; //===test
       nsta++;
       cv1=0;
       if (fmod(cv2,500)==1) printf("now read %d th %s\n",cv2,buff1);
       for(;;)
       {
         if (fgets(tstr,300,fin) == NULL) break;
	 if(sscanf(tstr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&lon,&lat,&temp,&temp2,&trash1,&trash2,&trash3,&trash4,&trash5)!=9)continue;
	  temp2=temp2+180; // change from (-180,180) to (0,360) //===NOT SURE IF THIS STEP IS RIGHT OR NOT
	 if(lon<0)lon=lon+360.;
         if(lon>x1+0.01||lon<x0-0.01|lat>y1+0.01||lat<y0-0.01)
           continue;
         i=int((lon-x0)/dx+0.1);
         j=int((lat-y0)/dy+0.1);
         if(temp<0.5&&temp>0.1)
           {
             slow[i][j][n[i][j]]=temp;
             azi[i][j][n[i][j]]=temp2;
             slow_weight[i][j][n[i][j]]=1;
             n[i][j]++;
	     neqk[i][j]++;
             cv1++;
           }
       }
      fclose(fin);
  }
  printf("\n\n=======\nafter read in eqk only: n[30][73]=%d\n====\n",n[30][73]);
  //read in the slow_azi file of noise
  for(ii=neveqk;ii<neveqk+nevnoise;ii++)
    {
      cv2 = cv2+1;
      sprintf(event_name,"%s\0",events_lst[ii]); //#slow_azi_TA.P42A.txt.HD
      //nsta++;
      sprintf(buff1,"%s/slow_azi_%s.txt.HD",dirnoise,event_name);
      if((fin=fopen(buff1,"r"))==NULL)
       {
	 //cout<<buff1<<" not exist!!"<<endl;
	 //return 1;
	 continue;
       }
       cv1=0;
       nsta++;
       if (fmod(cv2,500)==1) printf("now read %d th %s\n",cv2,buff1);
       for(;;)
       {
	 if (fgets(tstr,300,fin) == NULL) break;
         if (sscanf(tstr,"%lf %lf %lf %lf",&lon,&lat,&temp,&temp2) != 4) continue;
	 if(lon<0)lon=lon+360.;
	 if(lon>x1+0.01||lon<x0-0.01|lat>y1+0.01||lat<y0-0.01)
	   continue;
	 i=int((lon-x0)/dx+0.1);
	 j=int((lat-y0)/dy+0.1);
	 if(temp<0.5&&temp>0.1)
	   {
	     slow[i][j][n[i][j]]=temp;
	     azi[i][j][n[i][j]]=temp2;
	     slow_weight[i][j][n[i][j]]=1;
	     n[i][j]++;
 	     cv1++;
	     //printf("i=%d j=%d, lon=%f lat=%f, slowness=%f azi=%f,n=%d\n",i,j,lon,lat,temp,temp2,n[i][j]);//---test----
	     //exit(0);//----test----
	     //if(i==30 and j==73){printf("n[%d][%d]=%d\n",i,j,n[i][j]);}
	   }
       }
      //printf("ok here! %s %d %lf %lf\n",buff1,cv1,x1,x0);
      fclose(fin);
    }
  printf("\n\n=======\nafter read in eqk&noise only: n[30][73]=%d neqk[30][73]=%d\n====\n",n[30][73],neqk[30][73]);
  file_iso=fopen(name_iso,"w");
  printf("read ok! read in %d events(%d)+station(%d), %d of them has vel+slow_azi file\n",neveqk+nevnoise,neveqk,nevnoise,nsta);
  nsta=50;
  double w2;
  double temp_slow_sum;
  int temp_n;
  for(i=0;i<npts_x;i++)
    {
      
      #pragma omp parallel for default(none) private (j,k,temp,w2,weight_sum,kk,weight,temp_slow_sum,temp_n,slow_std) shared (file_iso,i,slow_weight,n,azi,slow,x0,y0,npts_x,npts_y,dx,dy,nsta,slow_sum,neqk) 
      for(j=0;j<npts_y;j++)
	{
	  if(n[i][j]<0.5*nsta)
	    {
	      #pragma omp critical(writei1)
	      {
	      fprintf(file_iso,"%lf %lf 0 9999 %d\n",x0+i*dx,y0+j*dy,n[i][j]);
	      }
	      continue;
	    }
	  w2=0;
	  weight_sum=0;
	  for(k=0;k<n[i][j];k++)
	    {
	      weight[k]=0;
	      for(kk=0;kk<n[i][j];kk++)
		{
		  if(fabs(azi[i][j][kk]-azi[i][j][k])<25)
		    weight[k]++;
		}
	      weight[k]=1/weight[k]*slow_weight[i][j][k];
	      weight_sum+=weight[k];
	      //if(i==marker_i&&j==marker_j)
	      //fprintf(stderr,"%lf %lf %lf\n",1/slow[i][j][k],weight[k],slow_weight[i][j][k]);
	    }
	  for(k=0;k<n[i][j];k++)
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
	  slow_std[i][j]=sqrt(temp/(1-w2));
	  w2=0;
	  weight_sum=0;
	  slow_sum[i][j]=0;
	  temp_n=0;
	  for(k=0;k<n[i][j];k++)
	    {
	      if(fabs(slow[i][j][k]-temp_slow_sum)>2.0*slow_std[i][j])
		continue;
	      weight_sum+=weight[k];
	      temp_n++;
	    }
	  for(k=0;k<n[i][j];k++)
	    {
	      if(fabs(slow[i][j][k]-temp_slow_sum)>2.0*slow_std[i][j])
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
	  slow_std[i][j]=sqrt(temp/(1-w2));

	  temp=slow_std[i][j]*sqrt(w2)/slow_sum[i][j]/slow_sum[i][j];
	  #pragma omp critical(writei2)
	  {
	  //fprintf(file_iso,"%lf %lf %lf %lf %d\n",x0+i*dx,y0+j*dy,1/slow_sum[i][j],temp,temp_n);
	  fprintf(file_iso,"%lf %lf %lf %lf %d  Nall %d Neqk %d Nin2stdall %d\n",x0+i*dx,y0+j*dy,1/slow_sum[i][j],temp,temp_n,n[i][j],neqk[i][j],temp_n);
	  }
	}
    }
  fclose(file_iso);
  file_ani=fopen(name_ani,"w");
  fout=fopen(name_ani_n,"w");  
  double test_lon,test_lat;
  //test_lon=240.6;
  //test_lat=36.6;
  printf("iso ok!\n");
  for(i=0;i<npts_x;i++)
    {
      #pragma omp parallel for default(none) private (j,k,ii,jj,kk,slow_sum1,slow_un,azi_weight_sum,azi_w2,azi_std,kkk,hist,temp) shared (file_ani,fout,i,n,nsta,azi,slow_weight,max,min,slow_sum,x0,y0,npts_y,npts_x,N_bin,dx,dy,d_bin,slow,check_sum) 
      for(j=0;j<npts_y;j++)
        {

	  for(k=0;k<N_bin;k++)
	    {
	      hist[k]=0;
	      slow_sum1[k]=0;
	      slow_un[k]=0;
	      azi_weight_sum[k]=0;
	      azi_w2[k]=0;
	      azi_std[k]=0;
	    }
	  if(i-3<0||i+3>=npts_x||j-3<0||j+3>=npts_y)
	    continue;
	  kkk=0;
	  for(ii=i-3;ii<=i+3;ii+=3)
	    {
	      for(jj=j-3;jj<=j+3;jj+=3)
		{
		  if(n[ii][jj]<nsta*0.5)
		    continue;
		  kkk+=n[ii][jj];
		}
	    }
 	  #pragma omp critical(writeo1)
 	  {
	  fprintf(fout,"%lf %lf %d\n",x0+i*dx,y0+j*dy,kkk);
	  }
	  if(kkk<9*nsta*0.5||n[i][j]<nsta*0.5)
	    continue;
	  for(ii=i-3;ii<=i+3;ii+=3)
	    {
	      for(jj=j-3;jj<=j+3;jj+=3)
		{
		  if(n[ii][jj]<nsta*0.5)
		    continue;
		  for(k=0;k<n[ii][jj];k++)
		    {
		      if(azi[ii][jj][k]>max||azi[ii][jj][k]<min)
			{
			  printf("azi(%.1f) out of range!![%.1f~%.1f]",azi[ii][jj][k],min,max);
			  exit(0);
			  //return 1;
			  //continue;
			}
		      //		      fprintf(stderr,"oliver!! %d %d\n",ii,jj);
		      hist[int((azi[ii][jj][k]-min)/d_bin)]++;
		      azi_weight_sum[int((azi[ii][jj][k]-min)/d_bin)]+=slow_weight[ii][jj][k];
		      azi_w2[int((azi[ii][jj][k]-min)/d_bin)]+=slow_weight[ii][jj][k]*slow_weight[ii][jj][k];
		      slow_sum1[int((azi[ii][jj][k]-min)/d_bin)]+=(slow[ii][jj][k]-slow_sum[ii][jj])*slow_weight[ii][jj][k];
		      //slow_un[int((azi[ii][jj][k]-min)/d_bin)]+=(slow[ii][jj][k]-slow_sum[ii][jj])*(slow[ii][jj][k]-slow_sum[ii][jj]);
		    }
		}
	    }
	  //	  fprintf(stderr,"tobal!! %d %d\n",i,j);
	  kk=0;
	  
	  for(k=0;k<N_bin;k++)
	    {
	      if(hist[k]>=10)
		{
		  kk++;
		  azi_w2[k]=azi_w2[k]/azi_weight_sum[k]/azi_weight_sum[k];
		  slow_sum1[k]=slow_sum1[k]/azi_weight_sum[k];
		}
	    }
	  //	  fprintf(stderr,"tobal_1!! %d %d\n",i,j);

	  for(ii=i-3;ii<=i+3;ii+=3)
            {
              for(jj=j-3;jj<=j+3;jj+=3)
                {
                  if(n[ii][jj]<nsta*0.5)
                    continue;
                  for(k=0;k<n[ii][jj];k++)
                    {
                      if(azi[ii][jj][k]>max||azi[ii][jj][k]<min)
                        {
			  printf("azi(%.1f) out of range!![%.1f~%.1f]",azi[ii][jj][k],min,max);
                          exit(0);
                        }
		      azi_std[int((azi[ii][jj][k]-min)/d_bin)]+=slow_weight[ii][jj][k]/azi_weight_sum[int((azi[ii][jj][k]-min)/d_bin)]*(slow[ii][jj][k]-slow_sum[ii][jj]-slow_sum1[int((azi[ii][jj][k]-min)/d_bin)])*(slow[ii][jj][k]-slow_sum[ii][jj]-slow_sum1[int((azi[ii][jj][k]-min)/d_bin)]);
                    }
                }
            }

	  #pragma omp critical(write_azi)
          {	  
	  fprintf(file_ani,"%lf %lf %d\n",x0+i*dx,y0+j*dy,kk);
          
	  for(k=0;k<N_bin;k++)
	    {
	      if(hist[k]>=10)
		{
		  azi_std[k]=sqrt(azi_std[k]/(1-azi_w2[k]));
		  temp=azi_std[k]*sqrt(azi_w2[k])/slow_sum[i][j]/slow_sum[i][j];
		  fprintf(file_ani,"%lf %lf %lf %lf %lf\n",min+(0.5+k)*d_bin,1/(slow_sum[i][j]+slow_sum1[k]),temp,azi_std[k],sqrt(azi_w2[k]));
		}
	    }
          }
	}
    }
  fclose(file_ani);
  fclose(fout);
  return 0;
}
