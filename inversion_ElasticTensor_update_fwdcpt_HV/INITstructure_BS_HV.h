//#include<iostream>
//#include<vector>
// this version added something to enable smooth vel profile (Bspline)
// this version, modified dispdef, added the h/v ratio info
#include<cstdlib>
//using namespace std;
/*
 * before using the compute_kernel, remember to update the disp inside model (disp should be cpt from the model) with" updatemodel(model,flagupdaterho) + compute_dispMineos(model,PREM,Nprem,Rflag,Lflag)" 
 * before using Vkernel2Lkernel, should keep the mod and para updated by using para2mod; usually Vpara2Lpara, which has para2mod inside, is performed before Vkernel2Lkernel, so normally, don't need another para2mod before Vkernel2Lkerenl
 *

*/
struct groupdef {
int np,nlay; //number of parameters for each para in each group (e.g., #of para describing vsv in group 1), filled in readmod; # of layers in each group, initialized in readmod;
int flagBs; // flag for Bspline, it will change to 1 once the Bspling fuctions are stored in
int flag;   //1--layered model 2--Bspline 3--Bspline(need to be changed to point) 4--gradient 5-- water layer 6--grid(line the grids and then u get the model, similar to MINEOS's input)
int pflag; // indicate the type of input model, if not specified, vsv=vsh,vpv=vph,eta=1, and vp=vs*vpvs; 1--vs; 2--vsv,vsh; 3--vsv,vsh,vp; 4--vsv,vsh,vpv,vph; 5--vsv,vsh,vpv,vph,eta;
int flagcpttype; // indicate the type of forward computation that will be used for this group; 1--use Vkernel to do all cpt (model is TI or iso) 2--use Vkernel(for RA) and Lovekernel(for AZ) 3--use Vkernel(for RA) and Azikernel(for AZ) 4--use Lovekernel (for both RA and AZ) 
// Vkernel-- dC/dV, dC/dh; 
// Lovekernel-- dC/dX, (dC/dh);
// Azikernel-- dC/dVcos, dC/dVsin;

double thick,vpvs; //thickness of each group,, vpvs ratio
//double anitheta,aniphi; 
vector<double> ratio,thick1,Bsplines;
vector<double> vsvvalue,vsvvalue1,vshvalue,vshvalue1,vpvvalue,vpvvalue1,vphvalue,vphvalue1,etavalue,etavalue1,rhovalue,rhovalue1,thetavalue1,thetavalue,phivalue,phivalue1;
//the rho is computed from vp_avg during readmodAniso

//vector<double> ratio,Rvalue,Rvalue1,thick1,Bsplines;
//vector<double> Lvalue,Lvalue1,Avalue1;
//vector<double> AAvalue,AAvalue1,CCvalue,CCvalue1,LLvalue,LLvalue1,FFvalue,FFvalue1,NNvalue,NNvalue1;
//vector<double> Bsvalue,Bcvalue,Bsvalue1,Bcvalue1,Gcvalue,Gsvalue,Gcvalue1,Gsvalue1,Hcvalue,Hsvalue,Hcvalue1,Hsvalue1,Ecvalue,Esvalue,Ecvalue1,Esvalue1;
//ratio, from readin file of layered model;value from readin file,vel;
//value1, layered vel; thick1, layered h;
//Bsplines, Bspline function the this group;
};

struct layermoddef {
int nlayer;
vector<double>vpvs,rho,qs,qp,thick;
vector<double> vsv,vsh,vpv,vph,eta,theta,phi;

//vector<double> vsv,vsh,vsavg,vpvs,vp,rho,qs,qp,thick;
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
int fhv,nhvper;
double pL,pmisfit,gL,gmisfit,L,misfit,hvmisfit,hvL;
//period, liability,misfit
double L2,misfit2,pL2,gL2,pmisfit2,gmisfit2;// this is for temperary use
vector<double> pper,pvelo,pvel,unpvelo,gper,gvelo,gvel,ungvelo,period1;
vector<double> hvper,hvratioo,hvratio,unhvratioo; // the H/V ratio dispersion curve.
//vector<double> paziampo,paziphio,paziamp,paziphi;
//readin period, readin vel, model-derived vel, readin unc;
//period1 stores the period that will be computed 
};

struct datadef{
rfdef rf;
dispdef Rdisp,Ldisp;
dispdef AziampRdisp,AziampLdisp,AziphiRdisp,AziphiLdisp;
double p,L,misfit;
//weight of rf;
};
///////////////////////////////////////////////////////////////////////////////////////////////
struct paradef{
int npara,flag; 
//number of parameter, vel and h; onece updated para, falg =1, after 1st model_derived para;
// the length of Lovepara is the same as V para
double L,misfit;
vector<double> parameter,LoveRAparameter;
vector< vector<double> > para0,LoveAZparameter;
// in para0[i]: 0)type_flag1 1)dv_type 2)dv 3)sigma 4)groupid 5)valueid 6)type_flag2 7)LVflag 8)RayleighWave_flag 9)Lovewave_flag 10)AZ_flag
// type_flag1:(filled in readpara) 0--value/Bcoeff; 1--gp thickness; -1--vpvs
// type_flag2:(filled in readpara) 1-vsv; 2-vsh; 3-vpv; 4-vph; 5-eta; 6-theta; 7-phi; 8-rho; 9-h; 10-dVsvcos; 11-dVsvsin; 12-vpvs;
// the rest flags are fill in mod2para when the para.flag<1
// LVflag:(filled in mod2para) 1--use Vpara&Vkernel to do all forward computation; -1--use Lovepara&Lovekernel to do all/part of the forward computation
// RayleighWave_flag(filled in mod2para): 1- this para affects Rayleighwave_RA/AZ_dispersion; 0-does not affect, 2-- affects AZ disp, not RA disp
// LoveWave_flag(filled in mod2para): same as above
// AZ_flag(filled in mod2para):2--this is a 2phi parameter; 4--4phi parameter; 0--not a AZ parameter

vector<vector<double> > space1;

//vector<double> Rparameter,Lparameter;
//vector<double> AAparameter,CCparameter,FFparameter,LLparameter,NNparameter;
////vector< vector<string> > Rpara0,Lpara0;
//vector<vector<double> > Rspace1,Lspace1;
//vector< vector<string> > AApara0,CCpara0,FFpara0,LLpara0,NNpara0;
//vector<vector<double> > AAspace1,CCspace1,FFspace1,LLspace1,NNspace1;
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
