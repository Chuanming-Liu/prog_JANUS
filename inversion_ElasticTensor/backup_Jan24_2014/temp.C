// test the set_intersection function
/*template <class InputIterator1, class InputIterator2, class OutputIterator>
  OutputIterator set_intersection (InputIterator1 first1, InputIterator1 last1,
                                   InputIterator2 first2, InputIterator2 last2,
                                   OutputIterator result)
{
  while (first1!=last1 && first2!=last2)
  {
    if (*first1<*first2) ++first1;
    else if (*first2<*first1) ++first2;
    else {
      *result = *first1;
      ++result; ++first1; ++first2;
    }
  }
  return result;
}

*/

#include <iostream>
#include <vector>
#include <fstream>
using namespace std;

vector<int> get_intersection_index(vector<double> lst1, vector<double> lst2){
  vector<int> index1,index2;
  vector<double> result;
  vector<int>::iterator first1,last1,first2,last2;
  int i1=0,i2=0;
  first1=lst1.begin();
  last1=lst1.end();
  while(first1!=last1 and first2!=last2){
	if (*first1<*first2){++first1;++i1;}
	else if (*first2<*first1){++first2;++i2}
	else{
		*result=*first;
		++result;++first1;++first2;
		index1.push_back(i1);index2.push_back(i2);
		++i1;++i2;
	}

  }
  for(int i=0;i<result.size();i++){
	printf("result%d = %g\n",i,result[i]);
  }
  return index1;
}


int main (){
//	int first[]={10,20,30}
//	int second[]={20,30,40,50}
vector<double> first(3);
first.push_back(10);first.push_back(20.);first.push_back(30.);
vector<double> second(4);
second.push_back(20.);second.push_back(30.);second.push_back(40.);second.push_back(50);

vector<double> v(10);
vector<int> index;

index=get_intersection(first,second);
for(i=0;i<index.size();i++){
	printf("index%d = %d\n",i,index[i]);
}

return 0;
}

