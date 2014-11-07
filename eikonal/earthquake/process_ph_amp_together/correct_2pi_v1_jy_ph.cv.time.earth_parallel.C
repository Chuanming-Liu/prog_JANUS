/*
#define MAIN
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
*/
using namespace std;
// every measurement is marked with cflag: 1--good data 2--do 2pi correction -1--bad data (even after 2pi correction)
// this corrects the travel time and amp together
// this code can hardly be parallized, since the the 2pi correction data depend on each other.

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


//int main(int na, char *arg[])
int do_correct_2pi(char *infile,double per,int stnumcri, double snrcri1, double snrcri2)
{ 
  /*if(na!=6)
    {
      cout<<"usage:correct_2pi 1]input_travel_time_file(the one from get_dist_measurements_parallel)\n2]period \n3]least_sta_#\n4]snr_criteri for end_sig_SNR\n5]snr_criteria for precurcer_sig_SNR\n"<<endl;
      return 0;
    }*/
  FILE *ff,*f1,*f2,*famp,*fpha1,*fpha2;
  char file_name[400]; //,tempnm[100];
  double dis,dis_min;///,snrcri1,snrcri2;
  double clon,clat;
  //per=atof(arg[2]);
  //snrcri1=atof(arg[4]);
  //snrcri2=atof(arg[5]);
  //sprintf(tempnm,"%s",arg[1]);
  if((f1 = fopen(infile, "r"))==NULL) {
    printf("cannot open %s\n",infile);
    return -1;
  }
  //clat=atof(arg[4]);
  //clon = atof(arg[5]);
  int N=3000; // should be larger than then maximum number of lines in each evt.ph.txt
  double lon[N],lat[N],time[N],vel[N],amp[N],snr;
  char stas[N][10];
  int YesNo[N];
  char stnm[N][100];
  int i,j,k,ii,jj,iev,ist,fflag[N];
  int tf,tflag,nst;
  double lon1[N],lat1[N],time1[N],vel1[N],amp1[N];
  double tlon,tlat,ttime,tvel,s1,s2,tamp,tempa;
  char stas1[N][10],tsta[10];
  int fc[N]; // good,no correction, 1; corrected; 2; rejected -1;
  char tevnm[20],tntnm[10],tstanm[10];
  double tdist,taz,tbaz,tphv,tgpv,ttph,ttgp,tsnr1,tsnr2;
  // cout<<"working on "<<tempnm<<endl;
  j = 0;k = 0;
  s1 = 0; s2 = 0;
  for(i=0;;i++)
    {
      if(i>N-1){
	printf("ERROR, the N(%d) is not big enough!\n",N);
        return -1;
	}

      if(fscanf(f1,"%s %s %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&tevnm[0],&tntnm[0],&tstanm[0],&tdist,&taz,&tbaz,&tlon,&tlat,&tphv,&tgpv,&ttph,&ttgp,&tamp,&tsnr1,&tsnr2)!=15)break;
      if(tsnr1>snrcri1 and tsnr2>snrcri2 and tgpv>0.){
		vel[j]=tphv;lon[j]=tlon;lat[j]=tlat;time[j]=ttph;amp[j]=tamp;strcpy(stas[j],tstanm);
		YesNo[j]=0;
		s1=s1+tlon;s2=s2+tlat;
		j=j+1;
	}
      else{
		lon1[k]=tlon;lat1[k]=tlat;time1[k]=ttgp;vel1[k]=tgpv;amp1[k]=tamp;strcpy(stas1[k],tstanm);
		k=k+1;
      }	
    }//for i
  ist=j; // useful ones;
  nst=k; // nonuseful ones;
  fclose(f1);
  if(ist<stnumcri)
    {
      printf("Only %d stations\n",ist);
      return -1;
    }
  clon = s1/ist; clat = s2/ist;
  //cout<<clon<<" "<<clat<<endl;
  //abort();
  int order[ist];
  double fact[ist];
  for(i=0;i<ist;i++)  // 
    {
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
  //sprintf(buff,"%s.c.txt",arg[1]);
  sprintf(buff,"%s.input.c.txt",infile);
  int ino=0,flag;
  char tempchar[100];
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
	    time[order[i]]-=per;
	   // flag=0;
	   }
	while(temp-time[order[i]]>per/2)
	   { 
	 // if(flag)
	   time[order[i]]+=per;
	  // flag=0;
  	   }
        }
      else {
        fc[order[i]]=1;
	      //break;
        }
      if(fabs(time[order[i]]-temp)>10 || fabs(amp[order[i]]-tempa)>20*amp[order[i]]||fabs(amp[order[i]]-tempa)>20*tempa)
	{
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
    }//for
  if(ino*1.0/ist>0.9||ist-ino<stnumcri)
    {
      printf("too many stations removed!! %d %d\n",ist,ino);
      return -1;
    }
  ff=fopen(buff,"w");
  
  for(i=0;i<ist;i++)
    {
	fprintf(ff,"%lf %lf %lf %lf %lf %s %d\n",lon[order[i]],lat[order[i]],time[order[i]],vel[order[i]],amp[order[i]],stas[order[i]],fc[order[i]]);	
    }
  for (i=0;i<nst;i++) {
    fprintf(ff,"%lf %lf %lf %lf %lf %s 0\n",lon1[i],lat1[i],time1[i],vel1[i],amp[i],stas1[i]);
    }
  fclose(ff);
  return 1;
}
