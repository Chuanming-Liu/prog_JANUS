// use HD and HD_0.8

// espetially for US
// change the marker_nn to 4
// change the distance cri to 150. 
// this version can work on either all events(sta) in the event(sta) list (argc==9) or on one single event(sta) (argc==10)
// I feel the final output of az baz is strange, need double check

/*
#define MAIN
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include "/home/jixi7887/progs/jy/eikonal/earthquake/dx.h"
#include "/home/jixi7887/progs/jy/eikonal/earthquake/dy.h"
*/
using namespace std;
/*
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
*/
int get_dist2(double lat1, double lon1, double lat2, double lon2, double *dist, double *azi, double *bazi)
{
double pi;
pi = 4.0*atan(1.0);
double cva = 6378.137;
double cvb = 6356.7523142;
double f = 1/298.257223563;
double L = 0.00;
double jcvA, jcvB;
L = lon1-lon2;
double U1 = 0;
U1 = atan((1-f)*tan(lat1/180*pi));
double U2 = 0;
U2 = atan((1-f)*tan(lat2/180*pi));
double cv,cv1,cv2,cv3,cv4,cv5,cvC,numda1;
L = L*pi/180;
double numda = L;
numda1 = numda;
do {
  numda = numda1;
  cv1 =  sqrt( (cos(U2)*sin(numda))*(cos(U2)*sin(numda))+ (cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(numda))*(cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(numda)) ); // cv1 sin(quan)
  cv2 = sin(U1)*sin(U2)+ cos(U1)*cos(U2)*cos(numda);
  cv = atan2(cv1,cv2);
  cv3 = cos(U1)*cos(U2)*sin(numda)/sin(cv);
  cv4 = 1 - cv3*cv3;

if (cv4 == 0)
   cv4 = 0.0000000001;

  cv5 = cos(cv) - 2*sin(U1)*sin(U2)/cv4;
  cvC = f/16*cv4*(4 + f*(4 - 3*cv4));
  numda1 = L + (1-cvC)*f*cv3*(cv + cvC*cv1*(cv5 + cvC*cv2*(-1 +2*cv5*cv5)));
} while (fabs(numda - numda1) > 0.0000000001);
double mius, cvA, cvB, deltacv,s;
mius = cv4*(cva*cva - cvb*cvb)/(cvb*cvb);
cvA = 1+mius/16384*(4096 + mius*(-768 + mius*(320 - 175*mius)));
cvB = mius/1024*(256+ mius*(-128 + mius*(74 - 47*mius)));
deltacv = cvB*cv1*(cv5 +cvB/4*(cv2*(-1 + 2*cv5*cv5)-cvB/6*cv5*(-3+4*cv1*cv1)*(-3+4*cv5*cv5) ));
s = cvb * cvA *(cv - deltacv);
jcvA = atan2( (cos(U2)*sin(numda1)),(cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(numda1)))*180/pi;
jcvB = atan2( (cos(U1)*sin(numda1)),(-sin(U1)*cos(U2)+sin(U2)*cos(U1)*cos(numda1)))*180/pi;


if (jcvB>180) jcvB = jcvB-180;
else jcvB = 180 - jcvB;

if (jcvA>180) jcvA = jcvA-180;
 else if (jcvA <0 )
       jcvA = -jcvA;
     else jcvA = 360 - jcvA;

*dist = s;
*azi = jcvA;
*bazi = jcvB;
return 1;
}

//int main(int na, char *arg[])
//evtlst: evnm,evlo,evla
int do_travel_time_to_slow_map( const char *sta1, double sta1_lon, double sta1_lat, double period,double dxdy, double x0, int xn,double y0,int yn,double cdist1,const char *indir, const char *outdir)
{// the effect of cdist: require the grid point has a close station
/*  if(na<9)
    {
      cout<<"usage:travel_time_to_velocity_map [station_morgan.lst] [period] [dx/dy] [x0] [xn] [y0] [yn] [cdist] [optional: sta]"<<endl;
      return 0;
    }
*/
  double ag1,ag2,diffa,diffb;
  FILE *ff,*fin,*fin2,*fin3,*fout,*file1;
  int i,j,cc,tflag;
  int npts_x,npts_y;
  char buff1[300],pflag[5],cvstr[300];
  double lat,lon,lat2,lon2,t_lat,t_lon,radius,pi;
  double lonin[2000],latin[2000],ttin[2000],ain[2000];
  int nn,cvi,cvii;
  int t_i,t_j;
  int marker_nn,marker_EN[2][2],marker_E,marker_N;
  double dist;
  double plat,plon,tdist,cdist;
  double az,baz;
  double dst;
  double xi,xj,xk;
  double yi,yj,yk;
  double loni,lati;
  double ang_n,slow_n;
  //period=atof(arg[2]);
  cdist = period * 4. * 3.;
  /*dxdy = atof(arg[3]);
  x0 = atof(arg[4]);
  xn = atoi(arg[5]);
  y0 = atof(arg[6]);
  yn = atoi(arg[7]);
  cdist1 = atof(arg[8]);
  */
  /*
  if (atoi(arg[8]) == 1) {
      sprintf(pflag,"phase.c");
      }
  else {
      sprintf(pflag,"group");
      }
  */
//  cout<<pflag<<endl;


  radius=6371.1391285;
  pi=4.0*atan(1.0);
  double dx,dy,x1,y1,temp,temp1,temp2,lat_temp;
  npts_x=xn;
  npts_y=yn;
  //fprintf(stderr,"Memory check!!\n");
  double tr_t[npts_x][npts_y];
  double dx_km[npts_y],dy_km[npts_y];
  double mdist1[npts_x][npts_y];
  double mdist2;
  
  
  int reason_n[npts_x][npts_y],cvn_lat;
  //fprintf(stderr,"Memory enough!!\n");
  
  dx=dxdy;//degree
  dy=dxdy;//degree
//  x0=235;
//  y0=26;
  x1=x0+(npts_x-1)*dx;
  y1=y0+(npts_y-1)*dy;
  for(j=1;j<npts_y-1;j++)
    {
      lat_temp=y0+j*dy;
      cvn_lat = int(lat_temp/0.2+0.1);
      if (cvn_lat>449) cvn_lat = 449;
      //lat_temp=atan(0.993277*tan(lat_temp/180*pi))*180/pi;
      //dx_km[j]=radius*sin((90-lat_temp)/180*pi)*dx/180*pi;
      dx_km[j]=dx_km1[cvn_lat]*dx/0.2;
      dy_km[j]=dy_km1[cvn_lat]*dy/0.2;
    }
// dy_km=radius*dy/180*pi;
//  cout<<dy_km<<endl;
//  abort();  

  //file1=fopen(arg[1],"r");
  //clock_t t,t1;
  //t = clock();
  //char AA[] = "334A\0";
  double trash0,trash1,trash2,trash3,trash6;
  int trash4;
  char trash5[6];
  

  //for(;;)
    //{
      
  //if(fscanf(file1,"%s %lf %lf",&sta1,&sta1_lon,&sta1_lat)==EOF)
     //break;
      
      //if (na == 10) {
      //if (strcmp (sta1,arg[9]) != 0) continue; }

      if (sta1_lon < 0)  
	 { sta1_lon = sta1_lon + 360.; }
      //sprintf(buff1,"travel_time_%s.%s.txt_v1.HD",sta1,pflag);
      sprintf(buff1,"%s/%s.%g.input.c.txt_v2.HD",indir,sta1,period);
      if((fin=fopen(buff1,"r"))==NULL)
       {
	 cout<<buff1<<" not exist!!"<<endl;
	 return -1;
       }
      //sprintf(buff1,"travel_time_%s.%s.txt_v1.HD_0.2",sta1,pflag);
      sprintf(buff1,"%s/%s.%g.input.c.txt_v2.HD_0.2",indir,sta1,period);
      if((fin2=fopen(buff1,"r"))==NULL)
	{
	  cout<<buff1<<" not exist!!"<<endl;
	  return -1;
	}
      //      sprintf(buff1,"travel_time_%s.txt_v1",sta1);
      //if((fin3=fopen(buff1,"r"))==NULL)
      //{
      //  cout<<buff1<<" not exist!!"<<endl;
      //  return 1;
      //}
//      sprintf(buff1,"slow_azi_%s.%s.txt.HD.2.v2",sta1,pflag);
      sprintf(buff1,"%s/slow_azi_%s.%g.txt.HD.2.v2",outdir,sta1,period);
      fout=fopen(buff1,"w");
      for(i=0;i<npts_x;i++)
	{
	  for(j=0;j<npts_y;j++)
	    {
	      tr_t[i][j]=0;
              reason_n[i][j] = 0;
	    }
	}
//      sprintf(buff1,"travel_time_%s.%s.txt_v1",sta1,pflag);
      sprintf(buff1,"%s/%s.%g.input.c.txt_v2",indir,sta1,period);
      if((fin3=fopen(buff1,"r"))==NULL)
	{
	  cout<<buff1<<" not exist!!"<<endl;
	  return -1;
	}      
      fclose(fin3);

      nn = 0;
      fin3=fopen(buff1,"r");
      for(cvi=0;cvi<2000;cvi++) {lonin[cvi]=0.;latin[cvi]=0.;ttin[cvi]=-1.;}
      cvi = 0;
      for (;;) {
         if (fgets(cvstr,100,fin3) == NULL) {break;}
         //sscanf(cvstr,"%lf %lf %lf %lf %s %d",&(lonin[cvi]),&(latin[cvi]),&temp2);
         sscanf(cvstr,"%lf %lf %lf %lf %lf %s %d",&(trash0),&(trash1),&trash2, &trash3, &trash6, &(trash5[0]),&trash4);
         if (trash4>0) {
            lonin[cvi] = trash0;
            latin[cvi] = trash1;
            ttin[cvi] = trash2;
            ain[cvi] = trash6;
            cvi = cvi + 1;
            }
         else {continue;}
         //if (fabs(lon2-lon) > cdist1/110. || fabs(lat2-lat) > cdist1/110. ) continue;
         //dist=get_dist(lat,lon,lat2,lon2);
         //if (mdist1[i][j] > dist) {mdist1[i][j] = dist;  }
         }
      nn = cvi;
      fclose(fin3);
      //cout<<"read in file ok! "<<nn<<endl;

      //cout <<"now read "<<buff1<<endl;

      //t1 = clock();
      //fprintf(stderr,"now time is : %g\n",(float)(t1-t)/CLOCKS_PER_SEC);
      //mdist1 = 1000.;

      for(;;)
	{
	  //mdist1 = 1000.;
	  if(fscanf(fin,"%lf %lf %lf",&lon,&lat,&temp)==EOF) break;
	  if(fscanf(fin2,"%lf %lf %lf",&lon2,&lat2,&temp2)==EOF) break;
	  if(lon!=lon2||lat!=lat2)
	    {
	      fprintf(stderr,"HD and HD_0.2 files not compatiable!!\n");
	      return -1;
	    }
	  if(lon>x1+0.01||lon<x0-0.01|lat>y1+0.01||lat<y0-0.01)
            {
	    continue;
            }

	  i=int((lon-x0)/dx+0.1);
	  j=int((lat-y0)/dy+0.1);
	  if(temp<temp2-2||temp>temp2+2||temp<2*period)   // 2 period criterior
	    {tr_t[i][j]=0;
	     reason_n[i][j] = 1;
	     continue;}
          mdist1[i][j] = 1000.;
	  marker_nn = 4;
	  marker_EN[0][0]=0;
	  marker_EN[0][1]=0;
	  marker_EN[1][0]=0;
          marker_EN[1][1]=0;

          tflag = 0;
	  for(cvi=0;cvi<nn;cvi++)    // looking for stations close to the point.
	    {
              lon2 = lonin[cvi]; lat2 = latin[cvi]; 
              // require the grid point has a close station
              if (fabs(lon2-lon) > cdist1/110. || fabs(lat2-lat) > cdist1/110. ) continue;  // too far, continue
	      if (lon2 < 0) lon2 = lon2+360.;

	      if(lon2-lon<0)
		marker_E=0;
	      else
		marker_E=1;

	      if(lat2-lat<0)
                marker_N=0;
              else
                marker_N=1;

              dist=get_dist(lat,lon,lat2,lon2);
	      if(marker_EN[marker_E][marker_N]!=0)
	      	continue;
              // require the grid point has a close station
	      if(dist<cdist1 && dist >= 1)
	      //if(dist>=1)
		{
		  marker_nn--;
		  if(marker_nn==0){tflag = 1;break;}
		  marker_EN[marker_E][marker_N]++;
		}
	    }
     	    if (tflag < 1) {
                temp = 0.;
                temp2 = 0.;
		reason_n[i][j] = 2;
		}
	  tr_t[i][j]=temp;
	}
      fclose(fin);
      fclose(fin2);

      cvi = 0;
      for(i=1;i<npts_x-1;i++) for(j=1;j<npts_y-1;j++){if (tr_t[i][j]>0) cvi++; }
      ////cout<<cvi<<endl;

      //      double temp1,temp2;
      //fclose(fout);
      //t1 = clock();
      //fprintf(stderr,"now time is : %g\n",(float)(t1-t)/CLOCKS_PER_SEC);

      //cout<<"now  wright!"<<endl;

      for(i=0;i<npts_x;i++) {
	  for(j=0;j<npts_y;j++) {
              if (i==0 || i==npts_x-1 || j==0 || j==npts_y-1) {
                reason_n[j][j]=8;
                fprintf(fout,"%lf %lf 0 999 %d\n",x0+i*dx,y0+j*dy,reason_n[i][j]);
                continue;
                }
              if (reason_n[i][j] == 2 || reason_n[i][j] == 1) {
                fprintf(fout,"%lf %lf 0 999 %d\n",x0+i*dx,y0+j*dy,reason_n[i][j]);
                continue;
                }

	      temp1=(tr_t[i+1][j]-tr_t[i-1][j])/2.0/dx_km[j];
	      temp2=(tr_t[i][j+1]-tr_t[i][j-1])/2.0/dy_km[j];
	      if(temp2==0)
		{
		  temp2=0.00001;
		}
	      temp=sqrt(temp1*temp1+temp2*temp2);
	      if(temp>0.6||temp<0.2)
		{
                  reason_n[i][j] = 3;
		  fprintf(fout,"%lf %lf 0 999 %d\n",x0+i*dx,y0+j*dy,reason_n[i][j]);
		}
	      else if ( tr_t[i+1][j]==0||tr_t[i-1][j]==0||tr_t[i][j+1]==0||tr_t[i][j-1]==0 ) 
		{
		  reason_n[i][j] = 4;
                  fprintf(fout,"%lf %lf 0 999 %d\n",x0+i*dx,y0+j*dy,reason_n[i][j]);
  		}
	      else
		{
                   // get mdist2 ////////////////////////////////////////////////////////////////////////////
                   mdist2 = 999.;
                   //cout<<i<<" "<<j<<endl;
                   lon = x0+i*dx; lat = y0+j*dy;
                   //cout<<lon<<" "<<lat<<endl;
                   for (cvi=0;cvi<nn;cvi++) {
                     lon2 = lonin[cvi];lat2 = latin[cvi];
                     //cout<<lon2<<" "<<lat2<<" "<<cvi<<" "<<nn<<" "<<mdist2<<" "<<cvii<<endl;
              	// require the grid point has a close station
                     if (fabs(lon2-lon) > cdist1/111. || fabs(lat2-lat) > cdist1/111. ) continue;
                     dist = 112.*pow((lon2-lon)*(lon2-lon) + (lat2-lat)*(lat2-lat),0.5);
                     if (mdist2 > dist) {mdist2 = dist;cvii = cvi;} 
                     }
                    mdist2 = get_dist(lat,lon,latin[cvii],lonin[cvii]);
                    //cout<<mdist2<<endl;
                    ////////////////////////////////////////////////////////////////////////////////////////////////

		    plat = y0+j*dy;
                    plon = x0+i*dx;
			
                    /////////////////////////// get azimuth ///////////////////////////////////////////
                    baz = 0.;az = 0.;
                    get_dist2(plat,plon,sta1_lat,sta1_lon,&dst,&az,&baz);     
                    //cout<<plat<<" "<<plon<<" "<<sta1_lat<<" "<<sta1_lon<<endl;
                    //cout<<dst<<" "<<az<<" "<<baz<<endl;
		    //abort();
                    ///////////////////////////////////////////////////////////////////////////////////
                    az = az + 180.;
                    az = 90.-az;
                    baz = 90.-baz;
                    if (az > 180.) az = az - 360.;
                    if (az < -180.) az = az + 360.;
                    if (baz > 180.) baz = baz - 360.;
                    if (baz < -180.) baz = baz + 360.;
                    ag1 = az;
                    //ang_n = az;
                    ag2 = atan2(temp2,temp1)/pi*180.;
                    diffa = ag2 - ag1;
                    //cout<<ag2<<" "<<ag1<<" "<<endl;
                    //abort();
                    if (diffa < -180.) diffa = diffa + 360.;
                    if (diffa > 180.) diffa = diffa - 360.;
                    //cout<<diffa<<endl;
                    //abort();

                    if (fabs(plat - sta1_lat)>5. || fabs(plon - sta1_lon)>5.) {
                         fprintf(fout,"%lf %lf %lf %lf %g %g %g %g %g\n",x0+i*dx,y0+j*dy,temp,ag2,mdist2,dst,az,baz,diffa);
                         //fprintf(stderr,"%lf %lf %lf %lf %g %g %g %g\n",x0+i*dx,y0+j*dy,temp,atan2(temp2,temp1)/pi*180,mdist2,dst,az,baz);
                         //cout<<sta1_lat<<sta1_lon<<endl;
                         //abort();
                         continue;
                         }
                    tdist = get_dist(plat,plon,sta1_lat,sta1_lon);
                    if (tdist < cdist-50.)  {
                       fprintf(fout,"%lf %lf 0 999 5\n",plat,plon); // too close to the central station 
                       }
                    else {
                        fprintf(fout,"%lf %lf %lf %lf %g %g %g %g %g\n",x0+i*dx,y0+j*dy,temp,ag2,mdist2,dst,az,baz,diffa);
                       }  
		}
	    }
	}
      fclose(fout);
      //t1 = clock();
      //fprintf(stderr,"now time is : %g\n",(float)(t1-t)/CLOCKS_PER_SEC);
    //}
  //fclose(file1);
  return 1;
}
