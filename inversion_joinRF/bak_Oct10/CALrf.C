// this is used to do the forward calculation for receiver function based on the model parameters (or layerized model?)
//
//
using namespace std;
// oh wait, the rf requires layerized model while Mineos use smooth model. Will this cause any inconsistency? 
// And, be careful about translating smooth model into layerized model.
// so, the Mineos and RF have totally different way of layerizing model. The updatemodelTibet and updategroup program should be different.
// Then how can we differentiate these two models? Only use the RF layer model while computing RF, which means, this RF layer model only exist here in this program?

extern "C"
{
void theo_ (int *n,float *fbeta,float *h,float *vps,float *qa,float *qb,float *fs,float *din,float *a0,float *c0,float *t0,int *nd,float *rx);
}
/*===========CONTENT=============
int compute_rf(modeldef &model,float depcri1, float depcri2,float qpcri, float qscri)
=================================
*/
//--------------------------------
// do para2mod before this subroutine

int compute_rf(modeldef &model, float depcri1, float depcri2,float qpcri, float qscri)
{
  int i,newnlayer,nl,nn,N=100,NN=1000; // N denote the maximum layers in the velocity model
  float slow,pi,din,rt;
//  double tvs[100]={0.},tvpvs[100]={0.},tqs[100]={0.},tqp[100]={0.},trho[100]={0.},tthick[100]={0.},rx[1000]={0.};
  float *tvs,*tvpvs,*tqs,*tqp,*trho,*tthick,*rx;
  modeldef RFmodel; // this is only used for storing the layerized model. The computed rf will store in the input model

  tvs = new float[N]; 
  tvpvs=new float[N]; 
  tqs=new float[N]; 
  tqp=new float[N]; 
  trho=new float[N]; 
  tthick=new float[N]; 
  rx=new float[NN];

  for(i=0;i<N;i++)
	{tvs[i]=tvpvs[i]=tqs[i]=tqp[i]=trho[i]=tthick[i]=0.;}
  for (i=0;i<NN;i++) rx[i]=0.;

  RFmodel=model;
  updatemodelRF(RFmodel,depcri1,depcri2,qpcri,qscri);
 
  newnlayer=RFmodel.laym0.nlayer<100?(RFmodel.laym0.nlayer+1):100; //The forward RF code can not handle model with too many layers, # of layer should be < 100
  nl=0;
  for(i=0;i<newnlayer-1;i++)
	{
	  if(RFmodel.laym0.vsv[i]>0.)
		{
		  tvs[nl]=RFmodel.laym0.vsv[i]; // RF is sensitive to Vsv not Vsh
	  	  tvpvs[nl]=RFmodel.laym0.vpvs[i];
		  tqs[nl]=RFmodel.laym0.qs[i];
		  tqp[nl]=RFmodel.laym0.qp[i];
		  tthick[nl]=RFmodel.laym0.thick[nl];
		  nl=nl+1;
		}//if
	}//for i
  tvs[nl]=tvs[nl-1];
  tvpvs[nl]=tvpvs[nl-1];
  tqs[nl]=tqs[nl-1];
  tqp[nl]=tqp[nl-1];
  tthick[nl]=0.; 
  nl=nl+1;

  nn=1000;
  slow=0.06;
  pi=atan(1)*4;
//  cout<<"===== pi="<<pi<<endl;
  din=180.*asin(tvs[nl-1]*tvpvs[nl-1]*slow)/pi;
  rt=model.data.rf.rt;

  float gau = 2.5;
  float dt = 0.005;
  float t0 = 0.;
  int nn1 = (int)newnlayer;

  /*FILE *ftemp; 
  ftemp=fopen("./temp_mod_rf.txt","w"); 
  for (i=0;i<nn1;i++) fprintf(ftemp,"%g %g %g %g \n",tvs[i],tthick[i],tvpvs[i],tqp[i]);
  fclose(ftemp);
  */
  theo_(&nn1,tvs,tthick,tvpvs,tqp,tqs,&rt,&din,&gau,&dt,&t0,&nn,rx);
  model.data.rf.tn.clear();
  model.data.rf.rfn.clear();
  //ftemp = fopen("check_rf.txt","w");
  for(i=0;i<nn;i++)
	{
	  //fprintf(ftemp,"%g %g\n",i*1./model.data.rf.rt,rx[i]);
	  model.data.rf.tn.push_back(i*1./model.data.rf.rt);
  	  model.data.rf.rfn.push_back(rx[i]);
	}//for i
  //fclose(ftemp);

  delete[] tvs;
  delete[] tvpvs;
  delete[] tqs;
  delete[] tqp;
  delete[] trho;
  delete[] tthick;
  delete[] rx;
  return 1;

}//compute rf


