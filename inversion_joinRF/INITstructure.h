#include<iostream>
#include<vector>
#include<cstdlib>
using namespace std;

struct groupdef {
int np,nlay;
int flagBs; // flag for Bspline, it will change to 1 once the Bspling fuctions are stored in
int flag;   //1--layered model 2--Bspline 3--gradient 5-- water layer
double thick,vpvs; //thickness of each group,, vpvs ratio
vector<double> ratio,Rvalue,Rvalue1,thick1,Bsplines;
vector<double> Lvalue,Lvalue1,Avalue1;
//ratio, from readin file of layered model;value from readin file,vel;
//value1, layered vel; thick1, layered h;
//Bsplines, Bspline function the this group;
};

struct layermoddef {
int nlayer;
vector<double> vsv,vsh,vsavg,vpvs,vp,rho,qs,qp,thick;
//parameters for layered model
};

struct rfdef {
int nrfo,rt;
double L,misfit;
vector<double> to,rfo,unrfo,tn,rfn,tnn,rfnn;
//parameters for receiver function
//to , rfo and unrfo are readin file
//tn, rfn are computed ones.
};

struct dispdef{
int npper,ngper,fphase,fgroup;
//fphase,fgroup indicate if group/phase disp is read(1) or not(0)
double pL,pmisfit,gL,gmisfit,L,misfit;
//period, liability,misfit
double L2,misfit2,pL2,gL2,pmisfit2,gmisfit2;// this is for temperary use
vector<double> pper,pvelo,pvel,unpvelo,gper,gvelo,gvel,ungvelo,period1;
//readin period, readin vel, model-derived vel, readin unc;
//period1 stores the period that will be computed 
};

struct datadef{
rfdef rf;
dispdef Rdisp,Ldisp;
double p,L,misfit;
//weight of rf;
};
///////////////////////////////////////////////////////////////////////////////////////////////
struct paradef{
int Rnpara,Lnpara,flag;
//number of parameter, vel and h; onece updated para, falg =1, after 1st model_derived para;
double L,misfit;
vector<double> Rparameter,Lparameter;
vector< vector<string> > Rpara0,Lpara0;
vector<vector<double> > Rspace1,Lspace1;
//??,??,range of para
};

struct modeldef{
int ngroup,flag,cc;
double tthick;
//#of readin groups;if layered falg=1;
vector<groupdef> groups; // a sequence of groups
layermoddef laym0;
datadef data;
};
