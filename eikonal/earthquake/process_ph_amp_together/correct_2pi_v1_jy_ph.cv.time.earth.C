#define MAIN
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
using namespace std;
// every measurement is marked with cflag: 1--good data 2--do 2pi correction -1--bad data (even after 2pi correction)
// this corrects the travel time and amp together

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
  if(na!=4)
    {
      cout<<"usage:correct_2pi event period least_sta_#"<<endl;
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
  //clat=atof(arg[4]);
  //clon = atof(arg[5]);
  int N=3000; // should be larger than then maximum number of lines in each evt.ph.txt
  double lon[N],lat[N],time[N],vel[N],amp[N],snr;
  char stas[N][6];
  int YesNo[N];
  char stnm[N][100];
  int i,j,k,ii,jj,iev,ist,fflag[N];
  int tf,tflag,nst;
  double lon1[N],lat1[N],time1[N],vel1[N],amp1[N];
  double tlon,tlat,ttime,tvel,s1,s2,tamp,tempa;
  char stas1[N][6],tsta[6];
  int fc[N]; // good,no correction, 1; corrected; 2; rejected -1;
   cout<<"working on "<<arg[1]<<endl;
  j = 0;k = 0;
  s1 = 0; s2 = 0;
  for(i=0;i<N;i++)
    {
      //if(fscanf(f1,"%lf %lf %lf %lf %lf %s %g",&lon[i],&lat[i],&time[i],&vel[i],&amp[i],&stnm[i],&snr)==EOF) 
      //if(fscanf(f1,"%lf %lf %lf %lf %lf",&lon[i],&lat[i],&time[i],&vel[i],&amp[i])==EOF) 
      if(fscanf(f1,"%lf %lf %lf %lf %lf %s %d",&tlon,&tlat,&ttime,&tvel,&tamp,&(tsta),&tflag)==EOF) 
	break;
/*	if (time[i]<1)
		{printf("wrong time stnm=%s time=%f\n",stnm[i],time[i]);
		i=i-1;}
*/
      if (tflag > 0.) {  // store everything useful in lon,lat,time,vel,and stas;
        lon[j]=tlon; lat[j]=tlat; time[j]=ttime; vel[j]=tvel; amp[j]=tamp;strcpy(stas[j],tsta);        
        YesNo[j]=0;
        s1 = s1 + tlon;
        s2 = s2 + tlat;
        j = j + 1; }
      else {
        lon1[k]=tlon; lat1[k]=tlat; time1[k]=ttime; vel1[k]=tvel;amp1[k]=tamp; strcpy(stas1[k],tsta);
        k = k + 1;
        }
    }
  ist=j; // useful ones;
  nst=k; // nonuseful ones;
  //if (ist )
  //clon = s1/ist; clat = s2/ist;
  cout<<ist<<endl;
  fclose(f1);
  if(ist<atof(arg[3]))
    {
      printf("Only %d stations\n",ist);
      exit(1);
    }
  clon = s1/ist; clat = s2/ist;
  cout<<clon<<" "<<clat<<endl;
  //abort();
  int order[ist];
  double fact[ist];
  for(i=0;i<ist;i++)  // 
    {
      //      fact[i]=lat[i]-lon[i];
//      fact[i]=(lat[i]-40)*(lat[i]-40)+(lon[i]-245)*(lon[i]-245);  // PARAMETER : center position of the array
      //if (fflag[i]<0.5) {fact[i]=999.;continue;}
      fact[i]=(lat[i]-clat)*(lat[i]-clat)+(lon[i]-clon)*(lon[i]-clon);  // PARAMETER : center position of the array
    }
  order[0]=0;
  for(i=1;i<ist;i++) 
// Order the stations according to fact[i], in list order[k], there stores the serial number of station, and the serial number is ordered by the value of fact[i]. 
//e.g. if sta i: 1 2 3, fact[i]: 100 500 200, 
//     order[i]: 1 3 2
    {
      //if (fflag[i]<0.5) continue;
      for(j=0;j<i;j++)
	{
	  //	  if(fact[i]>fact[order[j]])
          //if (fflag[j]<0.5) continue;
	  if(fact[i]<fact[order[j]])
	  break;
	}
      for(k=i;k>j;k--)
	{
          //if (fflag[j]<0.5) continue;
	  order[k]=order[k-1];
	}
      order[j]=i;      
    }
 // abort();
  int mark;
  double temp;
  char buff[300];
  sprintf(buff,"%s.c.txt",arg[1]);
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
	
      dis=vel[order[i]]*time[order[i]];
      temp=dis/vel[order[mark]];
      tempa = amp[order[mark]];

/// THIS PART MAY NEED TO BE CHANGED! IF WE MODIFY THE ARRIVAL TIME WITHOUT DOING ANYTHING TO THE AMP, THE AMP ISN"T THE AMPLITUDE AT THE MODIFIED ARRIVAL TIME
      if (fabs(time[order[i]]-temp) > per/2.) {
        fc[order[i]]=2;
	while(time[order[i]]-temp>per/2) {
	   // fprintf(fpha1,"%10s %10g %10g %10.4f %10g  %10s %10g %10g %10.4f %10g\n",stnm[order[i]],lon[order[i]],lat[order[i]],vel[order[i]],time[order[i]],stnm[order[mark]],lon[order[mark]],lat[order[mark]],vel[order[mark]],temp);
	    time[order[i]]-=per;
	   // flag=0;
	   }
	while(temp-time[order[i]]>per/2)
	   { 
	 // if(flag)
	 // fprintf(fpha1,"%10s %10g %10g %10.4f %10g  %10s %10g %10g %10.4f %10g\n",stnm[order[i]],lon[order[i]],lat[order[i]],vel[order[i]],time[order[i]],stnm[order[mark]],lon[order[mark]],lat[order[mark]],vel[order[mark]],temp);
	   time[order[i]]+=per;
	  // flag=0;
  	   }
        }
      else {
        fc[order[i]]=1;
	      //break;
        }
	//}
      if(fabs(time[order[i]]-temp)>10 || fabs(amp[order[i]]-tempa)>20*amp[order[i]]||fabs(amp[order[i]]-tempa)>20*tempa)
	{
	//fprintf(fpha2,"%10s %10g %10g %10.4f %10g  %10s %10g %10g %10.4f %10g\n",stnm[order[i]],lon[order[i]],lat[order[i]],vel[order[i]],time[order[i]],stnm[order[mark]],lon[order[mark]],lat[order[mark]],vel[order[mark]],temp);
          fc[order[i]]=-1;
	  YesNo[order[i]]=999;
	  ino++;
	  lon[order[i]]=lon[order[mark]];
	  lat[order[i]]=lat[order[mark]];
	  time[order[i]]=time[order[mark]];
	  vel[order[i]]=vel[order[mark]];
	  amp[order[i]]=amp[order[mark]];
	  //strcpy(stnm[order[i]],stnm[order[mark]]);
	  continue;
	}
      vel[order[i]]=dis/time[order[i]];
    }
  if(ino*1.0/ist>0.9||ist-ino<atof(arg[3]))
    {
      printf("too many stations removed!! %d %d\n",ist,ino);
      exit(1);
    }
  ff=fopen(buff,"w");
  
  for(i=0;i<ist;i++)
    {
      //if(YesNo[order[i]]==0)
	//fprintf(ff,"%lf %lf %lf %lf %lf\n",lon[order[i]],lat[order[i]],time[order[i]],vel[order[i]],amp[order[i]]);	
	fprintf(ff,"%lf %lf %lf %lf %lf %s %d\n",lon[order[i]],lat[order[i]],time[order[i]],vel[order[i]],amp[order[i]],stas[order[i]],fc[order[i]]);	
    }
  for (i=0;i<nst;i++) {
    fprintf(ff,"%lf %lf %lf %lf %lf %s 0\n",lon1[i],lat1[i],time1[i],vel1[i],amp[i],stas1[i]);
    }
  fclose(ff);
  return 0;
}
