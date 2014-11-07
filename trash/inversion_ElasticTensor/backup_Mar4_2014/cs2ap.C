#include <iostream>
#include <algorithm>

#include <vector>
#include <fstream>
#include <cmath>
#define _USE_MATH_DEFINES

using namespace std;


int cs2ap(double Ac,double As, double uncAc, double uncAs, double &amp,double &phi,double &uncamp, double &uncphi, int phiflag){

  double T=M_PI*2/phiflag;
  amp=sqrt(Ac*Ac+As*As);

  phi=atan2(Ac,As);
  phi=phi/phiflag;
	
  phi=T/4.-phi; //fast axis direction
  while(phi>T)
          phi=phi-T;
  while(phi<0)
          phi=phi+T;

  phi=phi*180./M_PI; //rad2deg

  //compute the unc of amp and phi
  uncamp = sqrt(pow(Ac*uncAc,2)+pow(As*uncAs,2))/amp; //(Ac*uncAc+As*uncAs)/amp;
  if(fabs(As)>0 and fabs(Ac)>0){uncphi = 180./M_PI*1/(1+pow(Ac/As,2))*sqrt(pow(uncAc/As,2)+pow(Ac*uncAs/As/As,2));}//(1/(1+pow(Ac/As,2))*(uncAc/As-Ac*uncAs/As/As))*180/M_PI;}
  else {uncphi=0.;}

  return 1;
}//cs2ap

int main(int argc, char *argv[]){
  if(argc!=5){
	printf("input 1)Ac 2)As 3)uncAc 4)uncAs\n");
	exit(0);
  }
  double As,Ac,amp,phi,uncAs,uncAc,uncamp,uncphi;
  

Ac=atof(argv[1]);
As=atof(argv[2]);
uncAc=atof(argv[3]);
uncAs=atof(argv[4]);

cs2ap(Ac,As,uncAc,uncAs,amp,phi,uncamp,uncphi,2);
  printf("%.4f_%.4f_%.4f_%.4f\n",amp,phi,uncamp,uncphi);  

  return 1;
}

