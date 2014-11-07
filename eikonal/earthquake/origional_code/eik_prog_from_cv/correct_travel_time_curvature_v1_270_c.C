#define MAIN
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
using namespace std;

double get_dist(double lat1,double lon1,double lat2,double lon2)
{
  double theta,pi,temp;
  double radius=6371;
  pi=4.0*atan(1.0);

  lat1=atan(0.993277*tan(lat1/180*pi))*180/pi;
  lat2=atan(0.993277*tan(lat2/180*pi))*180/pi;

  temp=sin((90-lat1)/180*pi)*cos(lon1/180*pi)*sin((90-lat2)/180*pi)*cos(lon2/180*pi)+sin((90-lat1)/180*pi)*sin(lon1/180*pi)*sin((90-lat2)/180*pi)*sin(lon2/180*pi)+cos((90-lat1)/180*pi)*cos((90-lat2)/180*pi);
  if(temp>1)
    {
      cout<<"warning cos(theta)>1 and correct to 1!!"<<temp<<endl;
      temp=1;
    }
  if(temp<-1)
    {
      cout<<"warning cos(theta)<-1 and correct to -1!!"<<temp<<endl;
      temp=-1;
    }
  theta=fabs(acos(temp));
  return theta*radius;
}


int main(int na, char *arg[])
{
  if(na!=7)
    {
      cout<<"usage:correct_tr_t_curvature event_name period lon1 lat1 n1 n2"<<endl;
      return 0;
    }
  FILE *ff,*fin,*fin_amp,*fin2,*fin3,*fout,*file1;
  int i,j;
  int npts_x,npts_y;
  char buff1[300],sta1[10],path[300];
  double lat,lon,lat2,lon2,t_lat,t_lon,radius,pi,sta1_lon,sta1_lat;
  int t_i,t_j;
  int marker_nn,marker_EN[2][2],marker_E,marker_N;
  double period,dist;
  period=atof(arg[2]);

  radius=6371.1391285;
  pi=4.0*atan(1.0);
  double dx,dy,x0,y0,x1,y1,temp,temp1,temp2,temp3,temp4,temp5,lat_temp;
//  npts_x=176;
//  npts_y=126;
  npts_x = atoi(arg[5]);
  npts_y = atoi(arg[6]);
  fprintf(stderr,"Memory check!!\n");
  double tr_t[npts_x][npts_y],amp[npts_x][npts_y];
  double dx_km[npts_y],dy_km;
  fprintf(stderr,"Memory enough!!\n");
  
  dx=0.2;//degree
  dy=0.2;//degree
//  x0=235;
//  y0=25;
  x0 = atof(arg[3]);
  y0 = atof(arg[4]);
  x1=x0+(npts_x-1)*dx;
  y1=y0+(npts_y-1)*dy;
  for(j=1;j<npts_y-1;j++)
    {
      lat_temp=y0+j*dy;
      lat_temp=atan(0.993277*tan(lat_temp/180*pi))*180/pi;
      dx_km[j]=radius*sin((90-lat_temp)/180*pi)*dx/180*pi;
    }
  dy_km=radius*dy/180*pi;
  sprintf(buff1,"%s.HD",arg[1]);
  if((fin=fopen(buff1,"r"))==NULL)
    {
      cout<<buff1<<" not exist!!"<<endl;
      return 1;
    }
  sprintf(buff1,"%s_am.HD",arg[1]);
  if((fin_amp=fopen(buff1,"r"))==NULL)
    {
      cout<<buff1<<" not exist!!"<<endl;
      return 1;
    }
  
  for(i=0;i<npts_x;i++)
    {
      for(j=0;j<npts_y;j++)
	{
	 tr_t[i][j]=0;
	 amp[i][j]=0;
	}
    }
  
  for(;;)
    {
      if(fscanf(fin,"%lf %lf %lf",&lon,&lat,&temp)==EOF) break;
      if(lon>x1+0.01||lon<x0-0.01|lat>y1+0.01||lat<y0-0.01)
	continue;
      i=int((lon-x0)/dx+0.1);
      j=int((lat-y0)/dy+0.1);
      
      tr_t[i][j]=temp;
    }
  fclose(fin);
  for(;;)
    {
      if(fscanf(fin_amp,"%lf %lf %lf",&lon,&lat,&temp)==EOF) break;
      if(lon>x1+0.01||lon<x0-0.01|lat>y1+0.01||lat<y0-0.01)
        continue;
      i=int((lon-x0)/dx+0.1);
      j=int((lat-y0)/dy+0.1);

      amp[i][j]=temp;
    }
  fclose(fin_amp);


  //      double temp1,temp2;
  marker_nn=0;
  int marker;
  sprintf(buff1,"%s",arg[1]);
  if((fin3=fopen(buff1,"r"))==NULL)
    {
      cout<<buff1<<" not exist!!"<<endl;
      return 1;
    } 
  //  sprintf(buff1,"%s.ph.txt_v2",arg[1]);
  //fout=fopen(buff1,"w");
  double arr[3000][6];
  char names[3000][10];
  char tname[10];
  int ii,jj,ist,ist_old,tflag;
  ist=0;
  ist_old=0;
  sprintf(buff1,"%s_v2",arg[1]);
  fout=fopen(buff1,"w");

  for(;;)
    {
      if(fscanf(fin3,"%lf %lf %lf %lf %lf %s %lf",&lon2,&lat2,&temp2,&temp3,&temp4,&(tname[0]),&temp5)==EOF) 
	{
	  fclose(fin3);
	  break;
	}
      ist_old++;
      if(temp2<2*period || temp5 < 0)
	{
	  arr[ist][0]=lon2;
	  arr[ist][1]=lat2;
	  arr[ist][2]=temp2;
	  arr[ist][3]=temp3;
	  arr[ist][4]=temp4;
          arr[ist][5]=(double)temp5;
          sprintf(names[ist],"%s\0",tname);
	  //  fprintf(fout,"%lf %lf %lf %lf %lf\n",lon2,lat2,temp2,temp3,temp4);
	  ist++;
	  continue;
	}
      i=int((lon2-x0)/dx+0.5);
      j=int((lat2-y0)/dy+0.5);
      if(i<0||j<0||i>npts_x||j>npts_y)
	{
	  ist_old--;
	  cout<<lon2<<" "<<lat2<<" "<<temp<<" out of range"<<endl;
	  continue;
	}
      temp=(tr_t[i+2][j]/-12.0+tr_t[i+1][j]*4.0/3.0+tr_t[i][j]*-5.0/2.0+tr_t[i-1][j]*4.0/3.0+tr_t[i-2][j]/-12.0)/dx_km[j]/dx_km[j];
      temp1=(tr_t[i][j+2]/-12.0+tr_t[i][j+1]*4.0/3.0+tr_t[i][j]*-5.0/2.0+tr_t[i][j-1]*4.0/3.0+tr_t[i][j-2]/-12.0)/dy_km/dy_km;
      temp=temp+temp1;
      if(temp>0.005||temp<-0.005)
	{
	  cout<<lon2<<" "<<lat2<<" curvature "<<temp<<endl;
          tflag = -2;
          fprintf(fout,"%lf %lf %lf %lf %lf %s %d\n",lon2,lat2,temp2,temp3,temp4,tname,tflag);
	  continue;
	}
      temp=(tr_t[i+1][j]-tr_t[i-1][j])/2.0/dx_km[j];
      temp1=(tr_t[i][j+1]-tr_t[i][j-1])/2.0/dy_km;
      if(temp1==0)
	{
	  temp1=0.00001;
	}
      temp=sqrt(temp1*temp1+temp*temp);
      if(temp>0.5 || temp<1/6.5)
	{
	  cout<<lon2<<" "<<lat2<<" vel "<<temp<<endl;
          tflag = -3;
          fprintf(fout,"%lf %lf %lf %lf %lf %s %d\n",lon2,lat2,temp2,temp3,temp4,tname,tflag);
	  continue;
	}
      temp=(amp[i+2][j]/-12.0+amp[i+1][j]*4.0/3.0+amp[i][j]*-5.0/2.0+amp[i-1][j]*4.0/3.0+amp[i-2][j]/-12.0)/dx_km[j]/dx_km[j];
      temp1=(amp[i][j+2]/-12.0+amp[i][j+1]*4.0/3.0+amp[i][j]*-5.0/2.0+amp[i][j-1]*4.0/3.0+amp[i][j-2]/-12.0)/dy_km/dy_km;
      temp=fabs((temp+temp1)/amp[i][j]*period*period/4/pi/pi);
      if(temp>1.0/16.0)
	{
	  cout<<"amp_curvature "<<lon2<<" "<<lat2<<" "<<temp<<endl;
	  continue;
	}
      //marker=0;
      //for(ii=-3;ii<=3;ii+=3)
      //{
      //  for(jj=-3;jj<=3;jj+=3)
      //    {
      //      temp=(amp[i+ii+1][j+jj]-amp[i+ii-1][j+jj])/2.0/dx_km[j+jj];
      //      temp1=(amp[i+ii][j+jj+1]-amp[i+ii][j+jj-1])/2.0/dy_km;
      //      temp=sqrt(temp1*temp1+temp*temp)/amp[i+ii][j+jj]*period;
      //      if(temp>0.15)
      //	{
      //	  marker=1;
      //	  break;
      //	}
      //    }
      //  if(marker==1)
      //    break;
      //}
      //if(marker==1)
      //{
      //  cout<<lon2<<" "<<lat2<<" amp gradient > 0.15 "<<temp<<endl;
      //  continue;
      //}
      //      
      //cout<<"amp_curvature "<<lon2<<" "<<lat2<<" "<<temp<<endl;
      arr[ist][0]=lon2;
      arr[ist][1]=lat2;
      arr[ist][2]=temp2;
      arr[ist][3]=temp3;
      arr[ist][4]=temp4;
      arr[ist][5]=(double)temp5;
      sprintf(names[ist],"%s\0",tname);

      //	  fprintf(fout,"%lf %lf %lf %lf %lf\n",lon2,lat2,temp2,temp3,temp4);
      ist++;
    }
  


  if(ist*1.0/ist_old<0.5)
    {
      cout<<"too many stations removed!! "<<ist<<" and "<<ist_old<<endl;
    }
  else
    if(ist<4)
      {
	cout<<"too less station left!!"<<endl;
      }
    else
      {
	//sprintf(buff1,"%s.ph.txt_v2",arg[1]);
	//fout=fopen(buff1,"w");
	for(i=0;i<ist;i++)
	  {
	    fprintf(fout,"%lf %lf %lf %lf %lf %s %d\n",arr[i][0],arr[i][1],arr[i][2],arr[i][3],arr[i][4],names[i],(int)arr[i][5]); 
	  }
	//fclose(fout);
      }
  fclose(fout);
  return 0;
}
