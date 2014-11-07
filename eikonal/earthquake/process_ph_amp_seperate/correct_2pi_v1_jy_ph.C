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
      cout<<"usage:correct_2pi filename period least_sta_# centra_lat centra_lon snrcri"<<endl;
      return 0;
    }
  FILE *ff,*f1,*f2,*famp,*fpha1,*fpha2;
  char file_name[400],tempnm[100];
  double dis,dis_min;
  double per,clon,clat;
  per=atof(arg[2]);
  sprintf(tempnm,"%s",arg[1]);
  if((f1 = fopen(tempnm, "r"))==NULL) {
    printf("cannot open %s\n",tempnm);
    exit(1);
  }
  clat=atof(arg[4]);
  clon = atof(arg[5]);
  double snrcri = atof(arg[6]);
  int N=3000; // should be larger than then maximum number of lines in each evt.ph.txt
  double lon[N],lat[N],time[N],vel[N],amp[N],snr;
  int YesNo[N];
  char stnm[N][100];
  int i,j,k,ii,jj,iev,ist;
  char tevnm[20],tntnm[10],tstanm[10];
  double tdist,taz,tbaz,tphv,tgpv,ttph,ttgp,tsnr1,tsnr2,tamp,tlon,tlat;

   ist=0;
   cout<<"working on "<<arg[1]<<endl;
  for(i=0;;i++)
    {
	if(ist>N-1){
		printf("N(%d) is too small!\n",N);
		exit(0);
	}
	if(fscanf(f1,"%s %s %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&tevnm[0],&tntnm[0],&tstanm[0],&tdist,&taz,&tbaz,&tlon,&tlat,&tphv,&tgpv,&ttph,&ttgp,&tamp,&tsnr1,&tsnr2)!=15)break;
	if(tsnr1>snrcri and tsnr2>snrcri and tgpv>0.){
		if(tlon<0)tlon=tlon+360; // lon--[0,360]
		lon[ist]=tlon;lat[ist]=tlat;time[ist]=ttph;vel[ist]=tphv;amp[ist]=tamp;
		sprintf(stnm[ist],"%s.%s",tntnm,tstanm);
     	        YesNo[ist]=0;
		ist++;
	}
//      if(fscanf(f1,"%lf %lf %lf %lf %lf %s %g",&lon[i],&lat[i],&time[i],&vel[i],&amp[i],&stnm[i],&snr)==EOF) 
//	break;
/*	if (time[i]<1)
		{printf("wrong time stnm=%s time=%f\n",stnm[i],time[i]);
		i=i-1;}
*/ 
    }
  ist=i;

  cout<<ist<<endl;
  fclose(f1);
  if(ist<atof(arg[3]))
    {
      printf("Only %d stations\n",ist);
      exit(1);
    }
  int order[ist];
  double fact[ist];
  for(i=0;i<ist;i++)  // 
    {
      //      fact[i]=lat[i]-lon[i];
//      fact[i]=(lat[i]-40)*(lat[i]-40)+(lon[i]-245)*(lon[i]-245);  // PARAMETER : center position of the array
      fact[i]=(lat[i]-clat)*(lat[i]-clat)+(lon[i]-clon)*(lon[i]-clon);  // PARAMETER : center position of the array
    }
  order[0]=0;
  for(i=1;i<ist;i++) 
// Order the stations according to fact[i], in list order[k], there stores the serial number of station, and the serial number is ordered by the value of fact[i]. 
//e.g. if sta i: 1 2 3, fact[i]: 100 500 200, 
//     order[i]: 1 3 2
    {
      for(j=0;j<i;j++)
	{
	  //	  if(fact[i]>fact[order[j]])
	  if(fact[i]<fact[order[j]])
	  break;
	}
      for(k=i;k>j;k--)
	{
	  order[k]=order[k-1];
	}
      order[j]=i;      
    }
 // abort();
  int mark;
  double temp;
  char buff[300];
  sprintf(buff,"%s.ph.txt_v1",arg[1]);
  //ff=fopen(buff,"w");
  // fprintf(ff,"%lf %lf %lf %lf %lf\n",lon[order[0]],lat[order[0]],time[order[0]],vel[order[0]],amp[order[0]]);
  int ino=0,flag;
  //  for(i=0;i<ist;i++)
 // cout<<i<<" "<<lon[order[i]]<<" "<<lat[order[i]]<<endl;
  char tempchar[100];
/*  sprintf(tempchar,"ampwrong_%s",arg[1]);
  famp = fopen(tempchar,"w");
  sprintf(tempchar,"pha2pi_%s",arg[1]);
  fpha1 = fopen(tempchar,"w");
  sprintf(tempchar,"phawrong_%s",arg[1]);
  fpha2 = fopen(tempchar,"w");
*/
  for(i=1;i<ist;i++)
/* ..1) for each station order[i] find its nearest station order[mark], if the amp of these two station different too much, erase the record of station[order[i]] with that of station[order[mark]], Otherwise
*/

    {
	//if(i%50==0)printf("test-- i=%d out of%d\n",i,ist);
      dis_min=999999999;
      for(j=0;j<i;j++)
	{
	  dis=get_dist(lat[order[i]],lon[order[i]],lat[order[j]],lon[order[j]]);
	  if(dis<dis_min)
	    {
	      dis_min=dis;
	      mark=j;
	      // cout<<mark<<endl;
	    }
	}
	
/*      if(fabs(amp[order[i]]-amp[order[mark]])>10*amp[order[i]]||fabs(amp[order[i]]-amp[order[mark]])>10*amp[order[mark]])
	{
	  //	  if(i==137)
	//   fprintf(stderr,"%d %d\n",order[i],order[mark]);
//	  cout<<dis<<" "<<amp[order[i]]<<" "<<amp[order[mark]]<<endl;
	//  fprintf(famp,"%10s %10g %10g  %10s %10g %10g\n",lon[order[i]],lat[order[i]],stnm[order[i]],lon[order[mark]],lat[order[mark]],stnm[order[mark]]);
	  YesNo[order[i]]=999;
	  ino++;
	  lon[order[i]]=lon[order[mark]];
	  lat[order[i]]=lat[order[mark]];
	  time[order[i]]=time[order[mark]];
	  vel[order[i]]=vel[order[mark]];
	  amp[order[i]]=amp[order[mark]];
	  strcpy(stnm[order[i]],stnm[order[mark]]);
	  continue;
	}
*/
      dis=vel[order[i]]*time[order[i]];
      temp=dis/vel[order[mark]];

	//flag=1;
      for(;;)// Otherwise,(amp doesn't differ too much), compare the arrival times, if dt > n*T/2: t_new=t+m*pi. Also, erase this station if the modified dt > 6s. vel_new=dist/t_new
	{
/// THIS PART MAY NEED TO BE CHANGED! IF WE MODIFY THE ARRIVAL TIME WITHOUT DOING ANYTHING TO THE AMP, THE AMP ISN"T THE AMPLITUDE AT THE MODIFIED ARRIVAL TIME
	  if(time[order[i]]-temp>per/2)
	   {
	   // fprintf(fpha1,"%10s %10g %10g %10.4f %10g  %10s %10g %10g %10.4f %10g\n",stnm[order[i]],lon[order[i]],lat[order[i]],vel[order[i]],time[order[i]],stnm[order[mark]],lon[order[mark]],lat[order[mark]],vel[order[mark]],temp);
	    time[order[i]]-=per;
	   // flag=0;
	   }
	  else
	    if(temp-time[order[i]]>per/2)
	   { 
	 // if(flag)
	 // fprintf(fpha1,"%10s %10g %10g %10.4f %10g  %10s %10g %10g %10.4f %10g\n",stnm[order[i]],lon[order[i]],lat[order[i]],vel[order[i]],time[order[i]],stnm[order[mark]],lon[order[mark]],lat[order[mark]],vel[order[mark]],temp);
	   time[order[i]]+=per;
	  // flag=0;
  	   }
	    else {
	      break;
		}
	}
      if(fabs(time[order[i]]-temp)>6)
	{
	//fprintf(fpha2,"%10s %10g %10g %10.4f %10g  %10s %10g %10g %10.4f %10g\n",stnm[order[i]],lon[order[i]],lat[order[i]],vel[order[i]],time[order[i]],stnm[order[mark]],lon[order[mark]],lat[order[mark]],vel[order[mark]],temp);
	  //	  if(i==137)
	  //fprintf(stderr,"%d %d %lf %lf\n",i,mark,lon[order[mark]],lat[order[mark]]);
//	  cout<<"misft dt > 6 sec"<<" "<<time[order[i]]<<" "<<temp<<" "<<fabs(time[order[i]]-temp)<< endl;
	  YesNo[order[i]]=999;
	  ino++;
	  lon[order[i]]=lon[order[mark]];
	  lat[order[i]]=lat[order[mark]];
	  time[order[i]]=time[order[mark]];
	  vel[order[i]]=vel[order[mark]];
	  amp[order[i]]=amp[order[mark]];
	  strcpy(stnm[order[i]],stnm[order[mark]]);
	  continue;
	}
      vel[order[i]]=dis/time[order[i]];
      //fprintf(ff,"%lf %lf %lf %lf %lf\n",lon[order[i]],lat[order[i]],time[order[i]],vel[order[i]],amp[order[i]]);
      //      fprintf(ff,"%lf %lf %lf %lf %lf %lf %lf %lf\n",lon[order[i]],lat[order[i]],time[order[i]],vel[order[i]],amp[order[i]],lon[order[mark]],lat[order[mark]],vel[order[mark]]);
      //cout<<lon[order[i]]<<" "<<lat[order[i]]<<" "<<fact[order[i]]<<endl;
    }
  /* fclose(famp);
   fclose(fpha1);
   fclose(fpha2);
   */
  //  cout<<ino<<" "<<ist<<" "<<ino/ist-0.9<<endl;
  if(ino*1.0/ist>0.5||ist-ino<atof(arg[3]))
    {
      printf("too many stations removed!! %d %d\n",ist,ino);
      exit(1);
    }
  ff=fopen(buff,"w");
  
  for(i=0;i<ist;i++)
    {
      if(YesNo[order[i]]==0)
	fprintf(ff,"%lf %lf %lf %lf %lf\n",lon[order[i]],lat[order[i]],time[order[i]],vel[order[i]],amp[order[i]]);	
    }
  
  fclose(ff);
  return 0;
}
