//read and write binary files
/*#include<iostream>
#include<fstream>
#include<vector>
*/
// this version, also write the LoveRAparameter[] and LoveAZparameter[][]
#include <fstream>
using namespace std;

//int write_bin(modeldef &model,char *fbname,paradef &para,int sign, int iiter, int iaccp)
int write_bin(modeldef &model,ofstream &out,paradef &para,int sign, int iiter, int iaccp)
{
//ofstream out(fbname,ios::in|ios::binary);
//if(!out)
//{cout<<"####cannot open file to write binary!!\n";exit (0);}
//-----parameter---------
int *ids,*n_len,*size,i;
double *Lmis,*param,*paramLRA,*paramLAZ1,*paramLAZ2,*gthick,*para_laym0;
ids=new int[3];
n_len=new int[3];
Lmis=new double[9];
param=new double[para.npara];
paramLRA=new double[para.npara];
paramLAZ1=new double[para.npara];
paramLAZ2=new double[para.npara];
gthick=new double[model.ngroup];
para_laym0=new double[11];
size = new int[6];
//------value passing---------
//--ids
ids[0]=sign; 
ids[1]=iiter; 
ids[2]=iaccp;
size[0]=sizeof(int)*3;
//--Lmis
Lmis[0]=model.data.L; 
Lmis[1]=model.data.misfit; 
Lmis[2]=model.data.rf.misfit; 
Lmis[3]=model.data.Rdisp.misfit;
Lmis[4]=model.data.Ldisp.misfit;
Lmis[5]=model.data.AziampRdisp.misfit;
Lmis[6]=model.data.AziampLdisp.misfit;
Lmis[7]=model.data.AziphiRdisp.misfit;
Lmis[8]=model.data.AziphiLdisp.misfit;
size[1]=sizeof(double)*9;
//--n_len
n_len[0]=para.npara; 
n_len[1]=model.ngroup;
n_len[2]=model.laym0.nlayer;
size[2]=sizeof(int)*3;
//--param
for(i=0;i<para.npara;i++){param[i]=para.parameter[i];}
size[3]=sizeof(double)*para.npara;
//--param LoveRAparameter
for(i=0;i<para.npara;i++){paramLRA[i]=para.LoveRAparameter[i];}
//--param LoveAZparameter1
for(i=0;i<para.npara;i++){paramLAZ1[i]=para.LoveAZparameter[i][0];}
//--param LoveAZparameter2
for(i=0;i<para.npara;i++){paramLAZ2[i]=para.LoveAZparameter[i][1];}
//--gthick
for(i=0;i<model.ngroup;i++)gthick[i]=model.groups[i].thick;
size[4]=sizeof(double)*model.ngroup;
//--para_laym0
size[5]=sizeof(double)*11;
//------write binary file----------------
/*
sign iiter iaccp
model.data.L   model.data.misfit    ...rf.misfit    ...Rdisp.misfit Ldisp.misfit
para.Rnpara+para.Lnpara   model.ngroup    model.laym0.nlayer
para.R/Lparameter[]
model.groups[].thick
model.laym0.thick[i] ...vsv[i] ...vsh[i] ...vpv[i] ..vph[i] ..eta[i] ...rho[i] ...qp[i] ...qs[i] ...theta[i] ...phi[i]
*/
out.write((char *)ids,size[0]);
out.write((char *)Lmis,size[1]);
out.write((char *)n_len,size[2]);
out.write((char *)param,size[3]);
out.write((char *)paramLRA,size[3]);
out.write((char *)paramLAZ1,size[3]);
out.write((char *)paramLAZ2,size[3]);
out.write((char *)gthick,size[4]);
for(i=0;i<model.laym0.nlayer;i++)
{
   para_laym0[0]=model.laym0.thick[i];
   para_laym0[1]=model.laym0.vsv[i];
   para_laym0[2]=model.laym0.vsh[i];
   para_laym0[3]=model.laym0.vpv[i];
   para_laym0[4]=model.laym0.vph[i];
   para_laym0[5]=model.laym0.eta[i];
   para_laym0[6]=model.laym0.rho[i];
   para_laym0[7]=model.laym0.qp[i];
   para_laym0[8]=model.laym0.qs[i];
   para_laym0[9]=model.laym0.theta[i];
   para_laym0[10]=model.laym0.phi[i];
   out.write((char *)para_laym0,size[5]);
}
//-----END of writing------------------
//out.close();
delete [] ids;
delete [] n_len;
delete [] Lmis;
delete [] param;
delete [] gthick;
delete [] para_laym0;
delete [] size;
return 1;
}


int read_bin(vector<modeldef> &modelall, vector<paradef> &paraall,char *fbname, vector<int> &signall, vector<int> &iiterall, vector<int> &iaccpall)
{
ifstream in(fbname,ios::in| ios::binary);
if(!in)
{cout<<"cannot open binary file to read! return\n"<<fbname<<endl;return 0;}
/*
sign iiter iaccp
model.data.L   model.data.misfit    ...rf.misfit    ...disp.misfit
para.npara   model.ngroup    model.laym0.nlayer
para.parameter[]
model.groups[].thick
model.laym0.thick[i] ...vs[i] ...vp[i] ...rho[i]
*/
/*
sign iiter iaccp
model.data.L   model.data.misfit    ...rf.misfit    ...Rdisp.misfit Ldisp.misfit
para.Rnpara+para.Lnpara   model.ngroup    model.laym0.nlayer
para.R/Lparameter[]
model.groups[].thick
model.laym0.thick[i] ...vsv[i] ...vsh[i] ...vp[i] ...rho[i] ...qp[i] ...qs[i]
*/
in.seekg(0,ios::end);
int SSIZE = (int)in.tellg();
in.seekg(0,ios::beg);
int N=0;
groupdef tgp;
modeldef model;
paradef para;
initmodel(model);
model.flag=1;//the readin model is a layered model.

while(in.tellg()<SSIZE){
N=N+1;
//cout<<"N="<<N<<endl;
int *size,*ids,*n_len,i;
double *misL,*gthick,*param,*paramLRA,*paramLAZ1,*paramLAZ2,*para_laym0;
double a=1.0,depth=0.,tqp,tqs;
int sign,iiter,iaccp;
vector<double> tv;
//--------initialization of the vectors -----------
model.groups.clear();
model.laym0.thick.clear();
model.laym0.vsv.clear();
model.laym0.vsh.clear();
model.laym0.vpv.clear();
model.laym0.vph.clear();
model.laym0.eta.clear();
model.laym0.rho.clear();
model.laym0.qs.clear();
model.laym0.qp.clear();
model.laym0.theta.clear();
model.laym0.phi.clear();
para.parameter.clear();
para.LoveRAparameter.clear();
para.LoveAZparameter.clear();

//-------------------------------------------------
size = new int[6];
ids = new int[3];
n_len = new int[3];
misL=new double[9];
size[0]=sizeof(int)*3;
size[1]=sizeof(double)*9;
size[2]=sizeof(int)*3;
size[5]=sizeof(double)*11;
in.read((char *)ids,size[0]);//sign iiter iaccp
in.read((char *)misL,size[1]);//model.data.L   model.data.misfit    ...rf.misfit    ...disp.Rmisfit ...disp.Lmisfit
in.read((char *)n_len,size[2]);//para.npara   model.ngroup    model.laym0.nlayer
size[3]=sizeof(double)*n_len[0];//double*para.npara
size[4]=sizeof(double)*n_len[1];//double*model.ngroup
param=new double[size[3]];
paramLRA=new double[size[3]];
paramLAZ1=new double[size[3]];
paramLAZ2=new double[size[3]];
gthick=new double[size[4]];
para_laym0=new double[size[5]];
in.read((char *)param,size[3]);//para.parameter[]
in.read((char *)paramLRA,size[3]);//para.LoveRAparameter[]
in.read((char *)paramLAZ1,size[3]);//para.LoveAZparameter[][0]
in.read((char *)paramLAZ2,size[3]);//para.LoveAZparameter[][1]
in.read((char *)gthick,size[4]);//model.groups[].thick
//----passing value to input structures -------------------------
for(i=0;i<n_len[2];i++)//model.laym0.nlayer
{
  in.read((char *)para_laym0,size[5]);
  model.laym0.thick.push_back(para_laym0[0]);
  model.laym0.vsv.push_back(para_laym0[1]);
  model.laym0.vsh.push_back(para_laym0[2]);
  model.laym0.vpv.push_back(para_laym0[3]);
  model.laym0.vph.push_back(para_laym0[4]);
  model.laym0.eta.push_back(para_laym0[5]);
  model.laym0.rho.push_back(para_laym0[6]);
  model.laym0.qp.push_back(para_laym0[7]);
  model.laym0.qs.push_back(para_laym0[8]);
  model.laym0.theta.push_back(para_laym0[9]);
  model.laym0.phi.push_back(para_laym0[10]);
  depth=depth+para_laym0[0];
  //if(depth<18.0) //here we used the Q in crust and mantle. haven't take the water,sediment into account.
  //{tqp=1400.;tqs=600;}
  //else
  //{tqp=200.;tqs=80.0;}
  //model.laym0.qs.push_back(tqs);
  //model.laym0.qp.push_back(tqp);
   
}
sign=ids[0];
iiter=ids[1];
iaccp=ids[2];
model.data.L=misL[0];
model.data.misfit=misL[1];
model.data.rf.misfit=misL[2];
model.data.Rdisp.misfit=misL[3];
model.data.Ldisp.misfit=misL[4];
model.data.AziampRdisp.misfit=misL[5];
model.data.AziampLdisp.misfit=misL[6];
model.data.AziphiRdisp.misfit=misL[7];
model.data.AziphiLdisp.misfit=misL[8];

para.misfit=model.data.misfit; //----there was a bug here, fixed on May 17, 2012
para.npara=n_len[0];
model.ngroup=n_len[1];
model.laym0.nlayer=n_len[2];
for(i=0;i<para.npara;i++){para.parameter.push_back(param[i]);}
for(i=0;i<para.npara;i++){para.LoveRAparameter.push_back(paramLRA[i]);}
for(i=0;i<para.npara;i++){tv.clear();tv.push_back(paramLAZ1[i]);tv.push_back(paramLAZ2[i]);para.LoveAZparameter.push_back(tv);}
for(i=0;i<model.ngroup;i++){model.groups.push_back(tgp);model.groups[i].thick=gthick[i];}
//int read_bin(vector<modeldef> &modelall, vector<paradef> &paraall,char *fbname, vector<int> &signall, vecort<int> &iiterall, vector<int> &iaccpall)

//if (sign==1 and model.data.misfit<1.5 )
if (sign ==1 )
{
modelall.push_back(model);
paraall.push_back(para);
signall.push_back(sign);
iiterall.push_back(iiter);
iaccpall.push_back(iaccp);
}//only stores acceptable models
//----------------------------------------------------------------
delete [] size;
delete [] ids;
delete [] n_len;
delete [] misL;
delete [] gthick;
delete [] param;
delete [] para_laym0;
}
in.close();
return 1;
}








