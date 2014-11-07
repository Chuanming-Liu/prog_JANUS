// this generates random number with the C++ <random> without using <boost>
// need gcc/4.8.0 or higher
// **** need a global generator defined in front of the main program*** 
// if producing generator inside the gen_random function, the resulting distribution won't be right.

// put these line in front of the main program
//#include <chrono>
//#include <random>
//unsigned seed=chrono::system_clock::now().time_since_epoch().count();
//default_random_engine generator (seed);

double gen_random_unif01(){
  double rand;
  uniform_real_distribution<double> distribution(0.,1.);
  rand=distribution(generator);
  return rand;
}

double gen_random_normal( double mean, double sigma){
  double rand;
  normal_distribution<double> distribution(mean,sigma);
  rand=distribution(generator);
  return rand;
}



 


