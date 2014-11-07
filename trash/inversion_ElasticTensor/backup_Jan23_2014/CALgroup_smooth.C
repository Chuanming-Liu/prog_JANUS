/*#include<fstream>
#include<iostream>
#include<vector>
//#include"INITstructure.h"
#include<algorithm>
*/
//#include"/home/jiayi/progs/jy/HEAD/head_c++/generate_Bs.C"
using namespace std;

/*========CONTENT===========
this is for the smooth model instead of layered model (for MINEOS instead of Hermman's code)
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
//	     *|(v1)
//	      |(h1)
//(v1 ,h1'=0) *----*(v2)
//	     	   |(h2)
//	           |
//(v2 ,h2'=0)      *----	
//
//
        int updategroup1(groupdef &group)//layer the model from layered model
        {
	  int i;
	  double tVs;
          if(group.flag!=1){cout<<"wrong flag, here updating layered model, flag=1, group.flag="<<group.flag<<endl;return 0;}
	  //if (group.thick1.size()+group.value1.size()!=0){cout<<"problem! the thick1 or value1 is not empty!"<<endl;return 0;}
	  group.thick1.clear();
	  group.vsvvalue1.clear();
	  group.vshvalue1.clear();
	  group.vpvvalue1.clear();
	  group.vphvalue1.clear();
	  group.etavalue1.clear();
	  group.thetavalue1.clear();
	  group.phivalue1.clear();
	  group.rhovalue1.clear();
          group.nlay=group.np*2;
          for (i=0;i<group.np;i++)
                {
		  group.thick1.push_back(group.thick*group.ratio[i]);
                  group.vsvvalue1.push_back(group.vsvvalue[i]);
		  group.vshvalue1.push_back(group.vshvalue[i]); 
                  group.vpvvalue1.push_back(group.vpvvalue[i]);
		  group.vphvalue1.push_back(group.vphvalue[i]); 
		  group.etavalue1.push_back(group.etavalue[i]);
		  group.thetavalue1.push_back(group.thetavalue[i]);
		  group.phivalue1.push_back(group.phivalue[i]);
		  group.rhovalue1.push_back(group.rhovalue[i]);

		  group.thick1.push_back(0.);
                  group.vsvvalue1.push_back(group.vsvvalue[i]);
		  group.vshvalue1.push_back(group.vshvalue[i]); 
                  group.vpvvalue1.push_back(group.vpvvalue[i]);
		  group.vphvalue1.push_back(group.vphvalue[i]); 
		  group.etavalue1.push_back(group.etavalue[i]);
		  group.thetavalue1.push_back(group.thetavalue[i]);
		  group.phivalue1.push_back(group.phivalue[i]);
		  group.rhovalue1.push_back(group.rhovalue[i]);
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
	  printf("nBs=%d order=%d group.thick=%g factor=%g nlay=%d\n",nBs,order,group.thick,factor,group.nlay);
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
        int updategroup2(groupdef &group)//layer the model from Bspline model   
        {
	  int i,j,nnlay,nBs;
	  double tmpvsv,tmpvsh,tmpvpv,tmpvph,tmpeta,tmptheta,tmpphi,tmprho;
          if(group.flag!=2){cout<<"wrong flag, here updating Bspline model, flag=2, group.flag="<<group.flag<<endl;return 0;}
	  //if (group.thick1.size()+group.value1.size()!=0){cout<<"problem! the thick1 or value1 is not empty!"<<endl;return 0;}
	  group.thick1.clear();
	  group.vsvvalue1.clear();
	  group.vshvalue1.clear();
	  group.vpvvalue1.clear();
	  group.vphvalue1.clear();
	  group.etavalue1.clear();
	  group.thetavalue1.clear();
	  group.phivalue1.clear();
	  group.rhovalue1.clear();

          if(group.flagBs!=1)updategroupBs(group);
	  nnlay=group.nlay;
	  nBs=group.np;
	  for (i=0;i<nnlay-1;i++)
		{
		  //tmpR=0.;tmpL=0.;
		  tmpvsv=tmpvsh=tmpvpv=tmpvph=tmpeta=tmptheta=tmpphi=tmprho=0.;
		  for (j=0;j<nBs;j++){
		    tmpvsv=tmpvsv+group.Bsplines[j*nnlay+i]*group.vsvvalue[j];
		    tmpvsh=tmpvsh+group.Bsplines[j*nnlay+i]*group.vshvalue[j];
		    tmpvpv=tmpvpv+group.Bsplines[j*nnlay+i]*group.vpvvalue[j];
		    tmpvph=tmpvph+group.Bsplines[j*nnlay+i]*group.vphvalue[j];
		    tmpeta=tmpeta+group.Bsplines[j*nnlay+i]*group.etavalue[j];
		    tmptheta=tmptheta+group.Bsplines[j*nnlay+i]*group.thetavalue[j];
		    tmpphi=tmpphi+group.Bsplines[j*nnlay+i]*group.phivalue[j];
		    tmprho=tmprho+group.Bsplines[j*nnlay+i]*group.rhovalue[j];
		  } //for j
		  group.vsvvalue1.push_back(tmpvsv);
		  group.vshvalue1.push_back(tmpvsh);
		  group.vpvvalue1.push_back(tmpvpv);
		  group.vphvalue1.push_back(tmpvph);
		  group.etavalue1.push_back(tmpeta);
		  group.thetavalue1.push_back(tmptheta);
		  group.phivalue1.push_back(tmpphi);
		  group.rhovalue1.push_back(tmprho);
		  
		  group.thick1.push_back(group.thick/(nnlay-1));
		}//for i
	  i=nnlay-1;tmpvsv=tmpvsh=tmpvpv=tmpvph=tmpeta=tmptheta=tmpphi=tmprho=0.;
	  for (j=0;j<nBs;j++){
		    tmpvsv=tmpvsv+group.Bsplines[j*nnlay+i]*group.vsvvalue[j];
		    tmpvsh=tmpvsh+group.Bsplines[j*nnlay+i]*group.vshvalue[j];
		    tmpvpv=tmpvpv+group.Bsplines[j*nnlay+i]*group.vpvvalue[j];
		    tmpvph=tmpvph+group.Bsplines[j*nnlay+i]*group.vphvalue[j];
		    tmpeta=tmpeta+group.Bsplines[j*nnlay+i]*group.etavalue[j];
		    tmptheta=tmptheta+group.Bsplines[j*nnlay+i]*group.thetavalue[j];
		    tmpphi=tmpphi+group.Bsplines[j*nnlay+i]*group.phivalue[j];
		    tmprho=tmprho+group.Bsplines[j*nnlay+i]*group.rhovalue[j];
	  }//for j
	  group.vsvvalue1.push_back(tmpvsv);
	  group.vshvalue1.push_back(tmpvsh);
	  group.vpvvalue1.push_back(tmpvpv);
	  group.vphvalue1.push_back(tmpvph);
	  group.etavalue1.push_back(tmpeta);
	  group.thetavalue1.push_back(tmptheta);
	  group.phivalue1.push_back(tmpphi);
	  group.rhovalue1.push_back(tmprho);

	  group.thick1.push_back(0.);
	  return 1;		
        }//updategroup2

//----------------------------------------------------- 
        int updategroup3(groupdef &group)//layer the model from gradient model
        {
	  int tn,j;
	  double dh,dvsv,dvsh,dvpv,dvph,deta,dtheta,dphi,drho;
	  if(group.flag!=4){cout<<"wrong flag, here updating gradient model, flag=4, group.flag="<<group.flag<<endl;return 0;}
	//  if (group.thick1.size()+group.value1.size()!=0){cout<<"problem! the thick1 or value1 is not empty!"<<endl;return 0;}
	  group.thick1.clear();
	  group.vsvvalue1.clear();
	  group.vshvalue1.clear();
	  group.vpvvalue1.clear();
	  group.vphvalue1.clear();
	  group.etavalue1.clear();
	  group.thetavalue1.clear();
	  group.phivalue1.clear();
	  group.rhovalue1.clear();

	  if(group.thick>2.)tn=max(int(group.thick/0.5),10);
	  else if(group.thick<0.5)tn=2;
	  else tn=4;
	  dh=group.thick/double(tn-1);
	  dvsv=(group.vsvvalue[1]-group.vsvvalue[0])/(tn-1.);//by default, only two values,[0] is smaller than [1]
	  dvsh=(group.vshvalue[1]-group.vshvalue[0])/(tn-1.);
	  dvpv=(group.vpvvalue[1]-group.vpvvalue[0])/(tn-1.);
	  dvph=(group.vphvalue[1]-group.vphvalue[0])/(tn-1.);
	  deta=(group.etavalue[1]-group.etavalue[0])/(tn-1.);
	  dtheta=(group.thetavalue[1]-group.thetavalue[0])/(tn-1.);
	  dphi=(group.phivalue[1]-group.phivalue[0])/(tn-1.);
	  drho=(group.rhovalue[1]-group.rhovalue[0])/(tn-1.);

	  
	  for(j=0;j<tn-1;j++)
		{
		  group.vsvvalue1.push_back(group.vsvvalue[0]+j*dvsv);
		  group.vshvalue1.push_back(group.vshvalue[0]+j*dvsh);
		  group.vpvvalue1.push_back(group.vpvvalue[0]+j*dvpv);
		  group.vphvalue1.push_back(group.vphvalue[0]+j*dvph);
		  group.etavalue1.push_back(group.etavalue[0]+j*deta);
		  group.thetavalue1.push_back(group.thetavalue[0]+j*dtheta);
		  group.phivalue1.push_back(group.phivalue[0]+j*dphi);
		  group.rhovalue1.push_back(group.rhovalue[0]+j*drho);
		 
		  group.thick1.push_back(dh);
		}
	   j=tn-1;
	   group.vsvvalue1.push_back(group.vsvvalue[0]+j*dvsv);
   	   group.vshvalue1.push_back(group.vshvalue[0]+j*dvsh);
           group.vpvvalue1.push_back(group.vpvvalue[0]+j*dvpv);
  	   group.vphvalue1.push_back(group.vphvalue[0]+j*dvph);
	   group.etavalue1.push_back(group.etavalue[0]+j*deta);
	   group.thetavalue1.push_back(group.thetavalue[0]+j*dtheta);
	   group.phivalue1.push_back(group.phivalue[0]+j*dphi);
	   group.rhovalue1.push_back(group.rhovalue[0]+j*drho);

	   group.thick1.push_back(0.);	
	  group.nlay=tn;
	  return 1;
        }//updategroup3
//----------------------------------------------------- 
	int updategroup4(groupdef &group)//layer the water layer
	{ // have some problems here *************** what's the input format for water layer???? ***** 
	  if(group.flag!=5){cout<<"wrong flag, here updating water model, flag=5, group.flag="<<group.flag<<endl;return 0;}
	  group.thick1.clear();
	  group.vsvvalue1.clear();
	  group.vshvalue1.clear();
	  group.vpvvalue1.clear();
	  group.vphvalue1.clear();
	  group.etavalue1.clear();
	  group.thetavalue1.clear();
	  group.phivalue1.clear();
	  group.rhovalue1.clear();

	  group.vsvvalue1.push_back(0.);
	  group.vshvalue1.push_back(0.);
	  group.vpvvalue1.push_back(group.vpvvalue[0]);
	  group.vphvalue1.push_back(group.vphvalue[0]);
	  group.etavalue1.push_back(1.0);
	  group.thetavalue1.push_back(0.);	  	  
	  group.phivalue1.push_back(0.);
	  group.rhovalue1.push_back(1.02);

	  group.thick1.push_back(group.thick-0.);
	
	  group.vsvvalue1.push_back(0.);
	  group.vshvalue1.push_back(0.);
	  group.vpvvalue1.push_back(group.vpvvalue[0]);
	  group.vphvalue1.push_back(group.vphvalue[0]);
	  group.etavalue1.push_back(1.0);
	  group.thetavalue1.push_back(0.);	  	  
	  group.phivalue1.push_back(0.);
	  group.rhovalue1.push_back(1.02);

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
	  else {cout<<"### ERROR updategroup: wrong value in group.flag : "<<group.flag<<endl;exit(0);}
        }//updategroup
//};//groupcal
