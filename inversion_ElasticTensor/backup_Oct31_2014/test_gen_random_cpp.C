// this is used to test if the gen_random code generates good random numbers
//double gen_random_normal(double mean, double sigma)
//
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <chrono>
#include <random>
#include "gen_random_cpp.C"
//#include "gen_random.C"

using namespace std;

int main(){
  double mean,sigma,rand1,rand2;
  mean=0,sigma=1;
  FILE *out;
  if((out=fopen("rand_output.txt","w"))==NULL){
	printf("Cannot open file to write\n");
	exit(0);
  }
  for (int i=0;i<1e6;i++){
	rand1=gen_random_normal(mean,sigma);
	rand2=gen_random_unif01();
  	fprintf(out,"%d %g %g\n",i,rand1,rand2);
  }
  fclose(out);
  return 1;
}

