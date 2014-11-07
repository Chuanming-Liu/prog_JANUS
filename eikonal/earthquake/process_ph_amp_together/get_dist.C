#define MAIN
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
using namespace std;
int get_dist(double lat1, double lon1, double lat2, double lon2, double *dist, double *azi, double *bazi)
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

int main ()
{
	printf("268.000000 31.000000 0.278913 104.328912 52.1675 8479.88 98.6531 100.464 5.67581                 |  268.000000 31.000000 0.278367 104.460828 52.1675 9832.66 51.1006 44.536 53\n");
	//
	//double stlo=-94.9829,stla=43.8565;
	double stlo=268.0,stla=31.0;
	double evlo=-80.159,evla=-44.806;
	double dist,az,baz;
	get_dist(evla,evlo,stla,stlo,&dist,&az,&baz);
	printf("hey azi=%f baz=%f dist=%f\n",az,baz,dist);
	az = az + 180.;
                    az = 90.-az;
                    baz = 90.-baz;
                    if (az > 180.) az = az - 360.;
                    if (az < -180.) az = az + 360.;
                    if (baz > 180.) baz = baz - 360.;
                    if (baz < -180.) baz = baz + 360.;
	printf("hey azi=%f baz=%f dist=%f\n",az,baz,dist);
	get_dist(stla,stlo,evla,evlo,&dist,&az,&baz);
	printf("hey azi=%f baz=%f dist=%f\n",az,baz,dist);
return 1;
}

