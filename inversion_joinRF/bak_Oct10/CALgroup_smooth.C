/*#include<fstream>
#include<iostream>
#include<vector>
//#include"INITstructure.h"
#include<algorithm>
*/
//#include"/home/jiayi/progs/jy/HEAD/head_c++/generate_Bs.C"
using namespace std;

/*========CONTENT===========
this is for the smooth model instead of layered model. 
Here,for the BSpline model: there are some 0-thickness layers 
and for the layered model: the nlay of the layered model is doubled

double get_AvgVs(double vsv, double vsh)
int updategroup1(groupdef &group)
int updategroup2(groupdef &group)
int updategroup3(groupdef &group)
int updategroupBs(groupdef &group) //calculate the Bspline function
int updategroup(groupdef &group)
//==========================
*/
//class groupcal{
//public:

	double get_AvgVs(double vsv, double vsh)
	{
	  double vs;
	  vs=sqrt((2*pow(vsv,2)+pow(vsh,2))/3.0);
	  return vs;
	}
//----------------------------------------------------- 
        int updategroup1(groupdef &group)//layer the model from layered model
        {
	  int i;double tVs;
          if(group.flag!=1){cout<<"wrong flag, here updating layered model, flag=1, group.flag="<<group.flag<<endl;return 0;}
	  //if (group.thick1.size()+group.value1.size()!=0){cout<<"problem! the thick1 or value1 is not empty!"<<endl;return 0;}
	  group.thick1.clear();
	  group.Rvalue1.clear();
	  group.Lvalue1.clear();
	  group.Avalue1.clear();
          group.nlay=group.np*2;
          for (i=0;i<group.np;i++)
                {
		  tVs = get_AvgVs(group.Rvalue[i],group.Lvalue[i]);
		  group.thick1.push_back(group.thick*group.ratio[i]);
                  group.Lvalue1.push_back(group.Lvalue[i]);
		  group.Rvalue1.push_back(group.Rvalue[i]); 
		  group.Avalue1.push_back(tVs);
		  
		  group.thick1.push_back(0.);
                  group.Lvalue1.push_back(group.Lvalue[i]);
		  group.Rvalue1.push_back(group.Rvalue[i]); 
		  group.Avalue1.push_back(tVs);
		}
        }//updategroup1
//----------------------------------------------------- 
        int updategroupBs(groupdef &group) //calculate the Bspline function     
        {
          int nBs,order;
          double factor;
          if(group.thick>=150)group.nlay=30;//60;
	  else if(group.thick<10)group.nlay=3;//5;
	  else if(group.thick<20)group.nlay=5;//10;
          else group.nlay=30;//30;
          nBs=group.np;
          if(nBs<4) order=3;
          else order=4; //cubic Bspline
          factor = 2.;
	  gen_B_spline(nBs,order,0.,group.thick,factor,group.nlay,group.Bsplines);
          group.flagBs=1;
	  return 1;
      }//updategroupBs
//----------------------------------------------------- 
//    *|
//     |__
//	  *|
//	   |
//	   ----*  
// layer the model from Bspline model. The last layer has thick1==0, the picture above shows how we parameterize Bspline into "3-node" model, line the nodes and you can obtain the model you want.
//
        int updategroup2(groupdef &group)//layer the model from Bspline model  (smooth version) 
        {
	  int i,j,nnlay,nBs;
	  double tmpR,tmpL;
          if(group.flag!=2){cout<<"wrong flag, here updating Bspline model, flag=2, group.flag="<<group.flag<<endl;return 0;}
	  //if (group.thick1.size()+group.value1.size()!=0){cout<<"problem! the thick1 or value1 is not empty!"<<endl;return 0;}
	  group.thick1.clear();group.Rvalue1.clear();group.Lvalue1.clear();group.Avalue1.clear();
          if(group.flagBs!=1)updategroupBs(group);
	  nnlay=group.nlay;
	  nBs=group.np;
	  for (i=0;i<nnlay-1;i++)
		{
		  tmpR=0.;tmpL=0.;
		  for (j=0;j<nBs;j++){tmpR=tmpR+group.Bsplines[j*nnlay+i]*group.Rvalue[j];tmpL=tmpL+group.Bsplines[j*nnlay+i]*group.Lvalue[j];} //???
		  group.Rvalue1.push_back(tmpR);
		  group.Lvalue1.push_back(tmpL);
		  group.Avalue1.push_back(get_AvgVs(group.Rvalue1[i],group.Lvalue1[i]));
		  group.thick1.push_back(group.thick/(nnlay-1));
		}
	  i=nnlay-1;tmpR=0.;tmpL=0.;
	  for (j=0;j<nBs;j++){tmpR=tmpR+group.Bsplines[j*nnlay+i]*group.Rvalue[j];tmpL=tmpL+group.Bsplines[j*nnlay+i]*group.Lvalue[j];}
	  group.Rvalue1.push_back(tmpR);
	  group.Lvalue1.push_back(tmpL);
	  group.Avalue1.push_back(get_AvgVs(group.Rvalue1[i],group.Lvalue1[i]));
	  group.thick1.push_back(0.);
	  return 1;		
        }//updategroup2

//----------------------------------------------------- 
        int updategroup3(groupdef &group)//layer the model from gradient model
        {
	  int tn,j;
	  double dh,dvL,dvR;
	  if(group.flag!=4){cout<<"wrong flag, here updating gradient model, flag=4, group.flag="<<group.flag<<endl;return 0;}
	//  if (group.thick1.size()+group.value1.size()!=0){cout<<"problem! the thick1 or value1 is not empty!"<<endl;return 0;}
	  group.thick1.clear();group.Rvalue1.clear();group.Lvalue1.clear();group.Avalue1.clear();
	  if(group.thick>2.)tn=max(int(group.thick/0.5),10);
	  else if(group.thick<0.5)tn=2;
	  else tn=4;
	  dh=group.thick/double(tn-1);
	    dvR=(group.Rvalue[1]-group.Rvalue[0])/(tn-1.);//by default, only two values,[0] is smaller than [1]
	    dvL=(group.Lvalue[1]-group.Lvalue[0])/(tn-1.);
	    for(j=0;j<tn-1;j++)
		{
		  group.Rvalue1.push_back(group.Rvalue[0]+j*dvR);
		  group.Lvalue1.push_back(group.Lvalue[0]+j*dvL);
		  group.Avalue1.push_back(get_AvgVs(group.Rvalue1[j],group.Lvalue1[j]));
		  group.thick1.push_back(dh);
		}
	   j=tn-1;
	   group.Rvalue1.push_back(group.Rvalue[0]+j*dvR);
	   group.Lvalue1.push_back(group.Lvalue[0]+j*dvL);
	   group.Avalue1.push_back(get_AvgVs(group.Rvalue1[j],group.Lvalue1[j]));
	   group.thick1.push_back(0.);	
	  group.nlay=tn;
	  return 1;
        }//updategroup3
//----------------------------------------------------- 
	int updategroup4(groupdef &group)//layer the water layer
	{ // have some problems here *************** what's the input format for water layer???? ***** 
	  if(group.flag!=5){cout<<"wrong flag, here updating water model, flag=5, group.flag="<<group.flag<<endl;return 0;}
	  group.thick1.clear();group.Rvalue1.clear();group.Lvalue1.clear();
	  group.Rvalue1.push_back(group.Rvalue[0]-1);
	  group.Lvalue1.push_back(-1);
	  group.Avalue1.push_back(-1);
	  group.thick1.push_back(group.thick-0.);
	
	  group.Rvalue1.push_back(group.Rvalue[0]-1);
	  group.Lvalue1.push_back(-1);
          group.Avalue1.push_back(-1);
          group.thick1.push_back(0.);

	  group.nlay=2;
	  return 1;
	}//updategroup4
//----------------------------------------------------- 
        int updategroup(groupdef &group)
        {
	  if(group.flag==1)updategroup1(group);//lay
	  else if(group.flag==2)updategroup2(group);//Bs
	  else if(group.flag==4)updategroup3(group);//grad
	  else if (group.flag==5)updategroup4(group);//water layer
	  else cout<<"wrong value in group.flag : "<<group.flag<<endl;
        }//updategroup
//};//groupcal
