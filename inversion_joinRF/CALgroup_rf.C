/*#include<fstream>
#include<iostream>
#include<vector>
//#include"INITstructure.h"
#include<algorithm>
*/
//#include"/home/jiayi/progs/jy/HEAD/head_c++/generate_Bs.C"
// CAUTION, since the forward RF code cannot handle model with more than 100 layers, so the crust part cannot have too many layers!
//
// this is the updategrouprf prgram which layerize model. This is different from the smooth version, the ouput points indicate layer model , not smooth model.
// In smooth model, the model is in fact not layerized, but just use some vel AT some depth to represent the points ON the model line
////    *|
////     |__
////        *|
////         |
////         ----*
// make a line going through these points, and you get the smooth model
//
// Here, It's different. It is layerized model, you CANNOT make a line going through these points to represnt the model. 
// // *|
// //  |___
// //      |
// //     *|
// //      ----|
// //          |
// //          *
//
using namespace std;
/*========CONTENT===========
double get_AvgVs(double vsv, double vsh)
int updategrouprf1(groupdef &group)
int updategrouprf2(groupdef &group)
int updategrouprf3(groupdef &group)
int updategrouprfBs(groupdef &group) //calculate the Bspline function
int updategrouprf(groupdef &group)
//==========================
*/
//class groupcal{
//public:
/* 
	double get_AvgVs(double vsv, double vsh)
	{
	  double vs;
	  vs=sqrt((2*pow(vsv,2)+pow(vsh,2))/3.0);
	  return vs;
	}
*/
//----------------------------------------------------- 
        int updategrouprf1(groupdef &group)//layer the model from layered model
        {
	  int i;
          if(group.flag!=1){cout<<"wrong flag, here updating layered model, flag=1, group.flag="<<group.flag<<endl;return 0;}
	  //if (group.thick1.size()+group.value1.size()!=0){cout<<"problem! the thick1 or value1 is not empty!"<<endl;return 0;}
	  group.thick1.clear();
	  group.Rvalue1.clear();
	  group.Lvalue1.clear();
	  group.Avalue1.clear();
          group.nlay=group.np;
          for (i=0;i<group.np;i++)
                { group.thick1.push_back(group.thick*group.ratio[i]);
                  group.Lvalue1.push_back(group.Lvalue[i]);
		  group.Rvalue1.push_back(group.Rvalue[i]); 
		  group.Avalue1.push_back(get_AvgVs(group.Rvalue1[i],group.Lvalue1[i]));
		}
        }//updategrouprf1
//----------------------------------------------------- 
        int updategrouprfBs(groupdef &group) //calculate the Bspline function     
        { 
          int nBs,order;
          double factor;
          if(group.thick>=150)group.nlay=30;//60;
	  else if(group.thick<10)group.nlay=3;//5;
	  else if(group.thick<20)group.nlay=5;//10;
	 // else if(group.thick<70)group.nlay=40;// changed on June 7, 2012
	  else group.nlay=30;//30;
          nBs=group.np;
          if(nBs<4) order=3;
          else order=4; //cubic Bspline
          factor = 2.;
	  gen_B_spline(nBs,order,0.,group.thick,factor,group.nlay,group.Bsplines);
          group.flagBs=1;
	  return 1;
      }//updategrouprfBs
//----------------------------------------------------- 
        int updategrouprf2(groupdef &group)//layer the model from Bspline model   
        {
	  int i,j,nnlay,nBs;
	  double tmpR,tmpL;
          if(group.flag!=2){cout<<"wrong flag, here updating Bspline model, flag=2, group.flag="<<group.flag<<endl;return 0;}
	  //if (group.thick1.size()+group.value1.size()!=0){cout<<"problem! the thick1 or value1 is not empty!"<<endl;return 0;}
	  group.thick1.clear();group.Rvalue1.clear();group.Lvalue1.clear();group.Avalue1.clear();
          if(group.flagBs!=1)updategrouprfBs(group);
	  nnlay=group.nlay;
	  nBs=group.np;
	  for (i=0;i<nnlay;i++)
		{
		  tmpR=0.;tmpL=0.;
		  for (j=0;j<nBs;j++){tmpR=tmpR+group.Bsplines[j*nnlay+i]*group.Rvalue[j];tmpL=tmpL+group.Bsplines[j*nnlay+i]*group.Lvalue[j];} //???
		  group.Rvalue1.push_back(tmpR);
		  group.Lvalue1.push_back(tmpL);
		  group.Avalue1.push_back(get_AvgVs(group.Rvalue1[i],group.Lvalue1[i]));
		  group.thick1.push_back(group.thick/nnlay);
		}
	  return 1;		
        }//updategrouprf2

//----------------------------------------------------- 
        int updategrouprf3(groupdef &group)//layer the model from gradient model
        {
	  int tn,j;
	  double dh,dvL,dvR;
	  if(group.flag!=4){cout<<"wrong flag, here updating gradient model, flag=4, group.flag="<<group.flag<<endl;return 0;}
	//  if (group.thick1.size()+group.value1.size()!=0){cout<<"problem! the thick1 or value1 is not empty!"<<endl;return 0;}
	  group.thick1.clear();group.Rvalue1.clear();group.Lvalue1.clear();group.Avalue1.clear();
	  if(group.thick>2.)tn=max(int(group.thick/0.5),10);
	  else if(group.thick<0.5)tn=2;
	  else tn=4;
	  dh=group.thick/double(tn);
	    dvR=(group.Rvalue[1]-group.Rvalue[0])/(tn-1.);//by default, only two values,[0] is smaller than [1]
	    dvL=(group.Lvalue[1]-group.Lvalue[0])/(tn-1.);
	    for(j=0;j<tn;j++)
		{
		  group.Rvalue1.push_back(group.Rvalue[0]+j*dvR);
		  group.Lvalue1.push_back(group.Lvalue[0]+j*dvL);
		  group.Avalue1.push_back(get_AvgVs(group.Rvalue1[j],group.Lvalue1[j]));
		  group.thick1.push_back(dh);
		}
	  group.nlay=tn;
	  return 1;
        }//updategrouprf3
//----------------------------------------------------- 
	int updategrouprf4(groupdef &group)//layer the water layer
	{ // have some problems here *************** what's the input format for water layer???? ***** 
	  if(group.flag!=5){cout<<"wrong flag, here updating water model, flag=5, group.flag="<<group.flag<<endl;return 0;}
	  group.thick1.clear();group.Rvalue1.clear();group.Lvalue1.clear();
	  group.Rvalue1.push_back(group.Rvalue[0]-1);
	  group.Lvalue1.push_back(-1);
	  group.Avalue1.push_back(-1);
	  group.thick1.push_back(group.thick-0.);
	  group.nlay=1;
	  return 1;
	}//updategrouprf4
//----------------------------------------------------- 
        int updategrouprf(groupdef &group)
        {
	  if(group.flag==1)updategrouprf1(group);//lay
	  else if(group.flag==2)updategrouprf2(group);//Bs
	  else if(group.flag==4)updategrouprf3(group);//grad
	  else if (group.flag==5)updategrouprf4(group);//water layer
	  else cout<<"wrong value in group.flag : "<<group.flag<<endl;
        }//updategrouprf
//};//groupcal
