#include<iostream>
#include<algorithm>
#include<vector>
#include<cmath>
#include<fstream>
#include"./string_split.C"
#include"./generate_Bs.C"
#include"./gen_random.C"
#include"./INITstructure.h"



using namespace std;
int main(int argc, char *argv[]){
 vector<int> foo (3,0),bar,hah;
 int i;

 
 for(i=0;i<3;i++){
   bar.push_back(i);
   hah.push_back(2*i);
 }
 hah.push_back(10);
 foo=hah;
 cout<<foo.size()<<endl;
 foo=bar;
 for(i=0;i<foo.size();i++){
   printf(" foo[%d]=%d\n",i,foo[i]);
 
 }

 cout<<"hihihihi\n";
 for(i=8;i<8;i++){
   printf("hihi %d\n",i);
 }


 int n=3;
 vector<int> aa(n,100),bb;
 for(i=0;i<n;i++){
   cout<<i<<" "<<aa[i]<<endl;  
 }
 cout<<aa.size()<<endl;
 
 vector<int> x1(3,0);
 i=x1.size();
 n=atoi(argv[1]);
 vector<vector<int> > x2(3,vector<int>(3,1));//,x3(atoi(argv[1]),x1);

 double x3[n][3];
 //cout<<"x3.size()="<<x3.size();
 
 for(i=0;i<x2.size();i++){
   for(int j=0;j<x2[i].size();j++){
     printf("%d ",x2[i][j]);
   }
   printf("\n");
 }
  
 vector<int> a1(2,0),a2(2,1),a3;
 //double a1[2]={0.},a2[2];
 a2=a1;
 a1[0]=100;
 a3=a2*100;
 
 for(i=0;i<a3.size();i++){
   printf("a2=%d, a3=%d\n",a2[i],a3[i]);
 }

 /*for(i=0;i<x3.size();i++){
   for(int j=0;j<x3[i].size();j++){
     printf("%d ",x3[i][j]);
   }
   printf("\n");
 }
 for(i=0;i<n;i++){
   for(int j=0;j<3;j++){
     printf("%g ",x3[i][j]);
   }
   printf("\n");
 }
 */

 return 1;
}
