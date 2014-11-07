/* content
int readPREM(const char *PREMnm, vector<vector<double> >  &PREM,int &Nprem)
int write_modMineos(modeldef &model, const char *outname,vector<vector<double> > PREM,int Nprem,int &Nmod)
int interpolate(dispdef indisp,dispdef &outdisp)
int interpolate_model(layermoddef inlay, layermoddef &outlay, double dh)
int read_dispMineos(dispdef &indisp,const char* Moutputnm, int Nmod)
int compute_dispMineos(modeldef &model,vector<vector<double> > PREM,int Nprem, int Rsurflag,int Lsurflag)
int compute_diff(dispdef disp1, dispdef disp2, vector<vector<double> > &pveldiff,vector<vector<double> > &gveldiff)
int compute_Vkernel_single_para(paradef para, int ip,modeldef model, vector<vector<double> > PREM,int Nprem, int Rflag, int Lflag,int flagupdaterho, vector<double> &trkp1, vector<double> &trkg1, vector<double> &tlkp1, vector<double> &tlkg1)
int compute_Vkernel(paradef para,modeldef &model,vector<vector<vector<double> > > &kernel,vector<vector<double> > PREM,int Nprem, int Rflag,int Lflag, int flagupdaterho)
int computeLovekernel(vector<double> Vkernel1,double rho,double vel,vector<double> &Lkernel1, int flagupdaterho)
int computeLovekerneleta(vector<double> Vkernel1,double c, vector<double> &Lkernel1)
int Vkernel2Lkernel(paradef para,modeldef model,vector<vector<vector<double> > > Vkernel,vector<vector<vector<double> > > &Lkernel, int flagupdaterho)
int Vpara2ET2LoveCoeff(vector<double> Vparameter, double[6][6] &ET,double[8] &RAcoeff, double[8][2] &AZcoeff)
int getGroupVPara(groupdef group, vector<vector<double> > &Vparameter2, int flagupdaterho)
int Vpara2Lovepara(paradef para,modeldef model,int flagupdaterho)
int compute_RAdisp(modeldef &model, paradef para, modeldef refmod, paradef refpara, vector<vector<vector<double> > > Vkernel, vector<vector<vector<double> > > Lkernel)
int cs2ap(double Ac,double As, double &amp,double &phi, int phiflag, int RLflag)
int compute_AZdisp(modeldef &model, paradef para, modeldef refmod, paradef refpara, vector<vector<vector<double> > > Vkernel, vector<vector<vector<double> > > Lkernel)
int compute_dispKernel(modeldef &model,paradef para, modeldef refmodel, paradef refpara, vector<vector<vector<double> > > Vkernel, vector<vector<vector<double> > >  Lkernel, int Rflag, int Lflag,int Razflag, int Lazflag)
int get_misfitKernel(modeldef &model,paradef &para,modeldef refmodel, paradef refpara, vector<vector<vector<double> > > Vkernel, vector<vector<vector<double> > >  Lkernel,int Rflag, int Lflag, int AziampRflag, int AziampLflag, int AziphiRflag, int AziphiLflag, float inpamp, float inpphi, int flagupdaterho)
int get_misfitMineosRA(modeldef &model,paradef &para, vector<vector<double> > PREM,int Nprem, int Rflag,int Lflag, int flagupdaterho)
int get_index(int var,vector<int> vv)
int Bsp2Point(modeldef BSmodel,paradef BSpara,modeldef &Pmodel,paradef &Ppara,int flagupdaterho)
*/
using namespace std;
#define _USE_MATH_DEFINES

int readPREM(const char *PREMnm, vector<vector<double> >  &PREM,int &Nprem)
{
  int i,j;
  vector<string> v;
  vector<double> tmod;
  string line;
  fstream Mf;
  
  PREM.clear();
  //read in PREM model
  Mf.open(PREMnm);
  Nprem=0;
  i=0;
  if(not Mf.is_open()){cout<<"### in compute_idspMineos, cannot open PREMfile "<<PREMnm<<endl;exit(0);}
  while(getline(Mf,line)){
    i++;
    v.clear();
    tmod.clear();
    if (i<4)continue;//the first three line are titles
    Split(line,v," ");
    if (v.size()!=9){cout<<"### check, problem with column # in PREM file,here size= "<<v.size()<<" i="<<i<<endl;exit(0);} //this line can be removed!!!!!!!
    for(j=0;j<9;j++)// by default, the column number should be 9
        tmod.push_back(atof(v[j].c_str()));
    PREM.push_back(tmod);
    Nprem++;    
  }//while
  Mf.close();

}

//---------------------------------------

int write_modMineos(modeldef &model, const char *outname,vector<vector<double> > PREM,int Nprem,int &Nmod)
{//combine the inversion model and the prem model together, write out into a file which will be used as Mineos input model
  FILE *outf;
  int i,N;
  double indepth,inrad,trad;
  vector<double> tmod;
  vector<vector<double> > Minput;
  if (model.flag!=1)
  {
	cout<<"### in write_modMineos, Molde hasn't been interperated!\n";
	exit(0);
  }
  //depth of the input model
  indepth=model.tthick;
  inrad=(6371.0000-indepth)*1000; //radius, meter
  N=0;
//  cout<<"*** test inrad="<<inrad<<" indepth="<<indepth<<endl;
  //before reaching the bottom of my inversion model. fill with PREM model
  for(i=0;i<Nprem;i++)//by default, the prem model has 9 columns and Nprem lines
  {
	trad=PREM[i][0];
	if (trad<inrad)
	 {Minput.push_back(PREM[i]);N=N+1;}
	else
	  break;
  }//for i
  //then, begin fill with my inversion model .
  int nlayer=model.laym0.nlayer;

  trad=inrad-model.laym0.thick[nlayer-1];
  double tmpth=0.;
  for (i=nlayer-1;i>-1;i--){ //revised on Sep 11, 2012
//  for (i=nlayer-1;i>0;i--){ //####################this seems to be wrong to me, shouldn't i reach 0 ??? in this way, the top layer is a single-velocity layer
	tmod.clear();
	tmpth=tmpth+model.laym0.thick[i]*1000.;
	tmod.push_back(trad+model.laym0.thick[i]*1000.);//0
	tmod.push_back(model.laym0.rho[i]*1000.);//1
  	tmod.push_back(model.laym0.vpv[i]*1000.);//2
 	tmod.push_back(model.laym0.vsv[i]*1000.);//3
	tmod.push_back(model.laym0.qp[i]);//4
	tmod.push_back(model.laym0.qs[i]);//5
  	tmod.push_back(model.laym0.vph[i]*1000.);//6
  	tmod.push_back(model.laym0.vsh[i]*1000.);//7
  	tmod.push_back(model.laym0.eta[i]);//8
 	Minput.push_back(tmod);
	

        trad=tmod[0];
	N=N+1;
	if (tmod[0]>6371000.0001){
	  printf("###### in write_modMineos, radius exceeds 6371000 m!!! %g tmpth=%g model.tthick=%g\n",tmod[0],tmpth,model.tthick);
	  exit(0);
	}
  }//for i
  tmod[0]=6371000;
//  Minput.push_back(tmod);//################# in this way, the top layer is a single-velocity layer,this works for present 1lay-sediment model. But in general, this seems wrong to me ?????? Aug23 2012. Need revision!!!
//  N=N+1;
  if((outf=fopen(outname,"w"))==NULL){
	printf("#### cannot open file to write!%s\n",outname);
	exit(0);
  }//if outf
  fprintf(outf,"radius  rho      vpv      vsv      qkappa   qshear   vph      vsh      eta\n");
  fprintf(outf,"1 1.00000 1\n");
  //fprintf(outf,"%d 33 66\n",N);
  fprintf(outf,"%d 25 71\n",N); //AK135 model
  for(i=0;i<N;i++)
  {
//   	cout<<i<<" "<<Minput[i].size()<<endl;
	fprintf(outf,"%8.0f%9.2f%9.2f%9.2f%9.1f%9.1f%9.2f%9.2f%9.5f\n",Minput[i][0],Minput[i][1],Minput[i][2],Minput[i][3],Minput[i][4],Minput[i][5],Minput[i][6],Minput[i][7],Minput[i][8]);
  }
  fclose(outf); 
  Nmod=N; //number of model layers
  return 1;
}//write_modMineos
//--------------------------------------

int interpolate(dispdef indisp,dispdef &outdisp)
{ // based on given period (the input period) and Mineos_T and Mineos_vel, interpolat Mineos_vel into values at the given periods
  // here, require that both pper and gper are sorted list(Increasing). Otherwise this code need to be changed (use id=find(period1.begin(),period1.end(),model.data.disp.gper[i]);like what we did for the Hermman inversion code) 
  vector<double> inper,refperiod;
  double Tref,tpvel,tgvel;
  int i,j,k,n,m;

  inper=indisp.pper;//same as indisp.gper, the per from Mineos is decreasing. While the refperiod is increasing
  outdisp.pvel.clear();outdisp.gvel.clear();

  if(outdisp.fphase>0){
    refperiod=outdisp.pper;
    n=outdisp.npper;
    k=indisp.npper;
    if (inper[0]<refperiod[n-1] or inper[k-1]>refperiod[0])
    {  
	printf("#### in interpolate, the range of refperiod is outside that of Mineos period!\nReset para for Mineos\n");
	printf("#### refT_min=%g refT_max=%g inT[0]=%g inT[k-1]=%g\n",refperiod[0],refperiod[n-1],inper[0],inper[k-1]);
	exit(0);
    }
    m=0;
    for(i=0;i<n;i++){
      Tref=refperiod[i];
      for(j=k-1;j>-1;j--){
	if(inper[j]<Tref)continue;
	tpvel=(indisp.pvel[j]-indisp.pvel[j-1])/(inper[j]-inper[j-1])*(Tref-inper[j-1])+indisp.pvel[j-1];
	outdisp.pvel.push_back(tpvel);
	k=j-1;//k is changing!
	m++;
	break;
      }//for j
    }//for i
    if(m!=n){printf("### in interpolate phase, inconsistancy between output disp.size() and inputT.size()!\n### outdisp.size()=%d inT.size()=%d\n",m,n);exit(0);}
  }//if fphase

  refperiod.clear();
  if(outdisp.fgroup>0){
    refperiod=outdisp.gper;
    n=outdisp.ngper;
    k=indisp.ngper;
    if (inper[0]<refperiod[n-1] or inper[k-1]>refperiod[0])
    {
        printf("#### in interpolate group, the range of refperiod is outside that of Mineos period!\nReset para for Mineos\n");
        printf("#### refT_min=%g refT_max=%g\n",refperiod[0],refperiod[n-1]);
        exit(0);
    }
    m=0;
    for(i=0;i<n;i++){
      Tref=refperiod[i];
      for(j=k-1;j>-1;j--){
        if(inper[j]<Tref)continue;
        tgvel=(indisp.gvel[j]-indisp.gvel[j-1])/(inper[j]-inper[j-1])*(Tref-inper[j-1])+indisp.gvel[j-1];
        outdisp.gvel.push_back(tgvel);
        k=j-1;//k is changing!
        m++;
	break;
      }//for j
    }//for i
    if(m!=n){printf("### in interpolate, inconsistancy between output disp.size() and inputT.size()!\n### outdisp.size()=%d inT.size()=%d\n",m,n);exit(0);}
  }// if fgroup

  return 1;
}//interpolate
//--------------------------------------
/*int interpolate_model(layermoddef inlay, layermoddef &outlay, double dh)
{// this is used to interpolate model to a constant interval model (dh between consecutive points is constant)
    double tthick,dep1,dep2,tdep2;
    double vsv,vsh,vpvs,vp,rho,qs,qp;
    int i,j,N,flag;

    tthick=0.;
    for(i=0;i<inlay.nlayer;i++){tthick=tthick+inlay.thick[i];}
    printf("===total thickness=%g\n",tthick);//--test---

    //    *|
    //     |__
    //        *|
    //         * 
    N=floor(tthick/dh)+2; // number of points ==N, number of layer==N,but the last 2 layers has different thickness(x, and 0) ..
    outlay.nlayer=N;   
    outlay.vsv.push_back(inlay.vsv[0]);
    outlay.vsh.push_back(inlay.vsh[0]);
    outlay.vpvs.push_back(inlay.vpvs[0]);
    outlay.vp.push_back(inlay.vp[0]);
    outlay.rho.push_back(inlay.rho[0]);
    outlay.qp.push_back(inlay.qp[0]);
    outlay.qs.push_back(inlay.qs[0]);
    outlay.thick.push_back(dh);

    dep1=0.;	    
    for(i=1;i<N-1;i++){
      dep1=dep1+dh;
      dep2=0.;flag=-1;
      tdep2=0.;
      if(dep1>tthick){printf("##### problem, the output depth is larger than the Max(input_depth)\n");exit(0);}
      for(j=1;j<inlay.nlayer;j++){//dismis the last layer in inlay, since its thickness is 0
	      dep2=dep2+inlay.thick[j-1];
	      if(dep2>dep1){vsv=(inlay.vsv[j]-inlay.vsv[j-1])/(inlay.thick[j-1])*(dep1-tdep2)+inlay.vsv[j-1];
	      	vsh=(inlay.vsh[j]-inlay.vsh[j-1])/inlay.thick[j-1]*(dep1-tdep2)+inlay.vsh[j-1];
		vpvs=(inlay.vpvs[j]-inlay.vpvs[j-1])/inlay.thick[j-1]*(dep1-tdep2)+inlay.vpvs[j-1];
		vp=(inlay.vp[j]-inlay.vp[j-1])/inlay.thick[j-1]*(dep1-tdep2)+inlay.vp[j-1];
		rho=(inlay.rho[j]-inlay.rho[j-1])/inlay.thick[j-1]*(dep1-tdep2)+inlay.rho[j-1];
	        qs=(inlay.qs[j]-inlay.qs[j-1])/inlay.thick[j-1]*(dep1-tdep2)+inlay.qs[j-1];
		qp=(inlay.qp[j]-inlay.qp[j-1])/inlay.thick[j-1]*(dep1-tdep2)+inlay.qp[j-1];
		flag=1;
		break;
	      }//if
	      else if (fabs(dep1-dep2)<0.001){
	      	vsv=inlay.vsv[j];vsh=inlay.vsh[j];vpvs=inlay.vpvs[j];vp=inlay.vp[j];
		rho=inlay.rho[j];qs=inlay.qs[j];qp=inlay.qp[j];	      
	        flag=1;
		break;
	      }// else if
	      tdep2=dep2;
      }//if j
      if(flag<0){printf("##### cannot interpolate for dep=%g, model tthick=%g dep2=%g\n",dep1,tthick,dep2);exit(0);}

      outlay.vsv.push_back(vsv);
      outlay.vsh.push_back(vsh);
      outlay.vpvs.push_back(vpvs);
      outlay.vp.push_back(vp);
      outlay.rho.push_back(rho);
      outlay.qs.push_back(qs);
      outlay.qp.push_back(qp);
      outlay.thick.push_back(dh);      
    }//if i
    outlay.thick[N-2]=tthick-(N-2)*dh;    
      outlay.vsv.push_back(vsv);
      outlay.vsh.push_back(vsh);
      outlay.vpvs.push_back(vpvs);
      outlay.vp.push_back(vp);
      outlay.rho.push_back(rho);
      outlay.qs.push_back(qs);
      outlay.qp.push_back(qp);
      outlay.thick.push_back(0.);
   return 1;
} // interpolate_model
*/
//--------------------------------------
int read_dispMineos(dispdef &indisp,const char* Moutputnm, int Nmod)
{//from Minoes output(Moutputnm with Nmod lines of model), read in the dispersion information to disp(indisp)
  fstream Mf;
  int i;
  double tper,tpvel,tgvel;
  string line;
  vector<string> v;
  dispdef tdisp;
  vector<double> period1;

  Mf.open(Moutputnm);
  if(not Mf.is_open()){cout<<"### in read_dispMineos, cannot open Moutput file "<<Moutputnm<<endl;exit(0);}
  i=0;
  tdisp.pper.clear();tdisp.pvel.clear();tdisp.gper.clear();tdisp.gvel.clear();
  tdisp.npper=0;tdisp.ngper=0;
  while(getline(Mf,line))
  {
    i++;
    v.clear();
    //line 1~(N+10), start from 1
    if (i<Nmod+12)continue;
    Split(line,v," ");
    tper=1000.0/atof(v[4].c_str());
    tpvel=atof(v[3].c_str());
    tgvel=atof(v[6].c_str());
    tdisp.pper.push_back(tper);
    tdisp.gper.push_back(tper);
    tdisp.pvel.push_back(tpvel);
    tdisp.gvel.push_back(tgvel);
    tdisp.npper++;
    tdisp.ngper++;
  }//while  
  Mf.close();

  //interpolate the tdisp into specified period list, and store the vel in indisp
  interpolate(tdisp,indisp);
/*  **** CHECK HERE ********
  FILE *tempf;
  if((tempf=fopen("temp_disp1.txt","w"))==NULL)cout<<"cannot open file to write\n";
  printf("N pper=%d pvel=%d pvelo=%d gper=%d gvel=%d gvelo=%d\n",indisp.pper.size(),indisp.pvel.size(),indisp.pvelo.size(),indisp.gper.size(),indisp.gvel.size(),indisp.gvelo.size());
  for(i=0;i<indisp.npper;i++)
    {
     if(indisp.fgroup>0)
    fprintf(tempf,"%g %g %g %g %g %g\n",indisp.pper[i],indisp.pvel[i],indisp.pvelo[i],indisp.gper[i],indisp.gvel[i],indisp.gvelo[i]);
     else
    fprintf(tempf,"%g %g %g \n",indisp.pper[i],indisp.pvel[i],indisp.pvelo[i]);
    }
  fclose(tempf);   
  if((tempf=fopen("temp_disp2.txt","w"))==NULL)cout<<"cannot open file to write\n";
  for(i=0;i<tdisp.npper;i++)
    fprintf(tempf,"%g %g %g %g\n",tdisp.pper[i],tdisp.pvel[i],tdisp.gper[i],tdisp.gvel[i]);
  fclose(tempf);   
  ***********************
*/
  return 1;
}//read_dispMineos

//--------------------------------------
int compute_dispMineos(modeldef &model,vector<vector<double> > PREM,int Nprem, int Rsurflag,int Lsurflag, int ipara)
{ //execute mineos_bran, give a model, compute dispersion curves for L and/or R, store value in model.d.R/Ldisp.pvel/gvel
  // !!! REMEMBER TO USE UPDATEMODEL BEFORE USING THIS SUBROUTINE !!!
  int Nmod,jcmp,i;
  string line;
  float wmin,wmax;
  vector<string> v;
  vector<double> tmod;
  char Minmodnm[200],Moutputnm[200];
  char str[500],moddir[100],modnm[100];
  dispdef Rdisp,Ldisp;
  //******parameters********
  wmin=1000./70.;//1000./70.;
  wmax=1000./10.;
  //************************

  //write Mineos input model
  sprintf(Minmodnm,"MineosInputMod_%d.txt",ipara);
  write_modMineos(model,Minmodnm,PREM,Nprem,Nmod);
//  cout<<"**** test Nmod="<<Nmod<<endl;// Nmod record the # of model layers

  //execute Mineos and get disp_calc
  sprintf(moddir,".");
  sprintf(modnm,"MineosInputMod_%d",ipara);
  if(Lsurflag>0){
     jcmp=2;
     sprintf(str,"csh run_Mineos_bran.csh %d %s %s %.1f %.1f\n",jcmp,moddir,modnm,wmin,wmax);
     system(str);
     sprintf(Moutputnm,"%s_T",modnm);
     read_dispMineos(model.data.Ldisp,Moutputnm,Nmod);  	
     model.data.AziampLdisp.pvel.clear();model.data.AziampLdisp.gvel.clear();model.data.AziphiLdisp.pvel.clear();model.data.AziphiLdisp.gvel.clear();

     // bug fixed on Oct 30 2013. modification
     for( i=0;i<model.data.AziampLdisp.npper;i++){
     	model.data.AziampLdisp.pvel.push_back(0.);
     }
     for(i=0;i<model.data.AziphiLdisp.npper;i++)
     	model.data.AziphiLdisp.pvel.push_back(0.);
     for( i=0;i<model.data.AziampLdisp.ngper;i++){
     	model.data.AziampLdisp.gvel.push_back(0.);
     }
     for( i=0;i<model.data.AziphiLdisp.ngper;i++)
     	model.data.AziphiLdisp.gvel.push_back(0.);
  }//if Lsurflag
  if(Rsurflag>0){
     jcmp=3;
     sprintf(str,"csh run_Mineos_bran.csh %d %s %s %.1f %.1f\n",jcmp,moddir,modnm,wmin,wmax);
     //printf("@@@ check, compute_dispMineos:  %s  #################\n",str);
     system(str);
     //get disp info from Mineos output 
     sprintf(Moutputnm,"%s_S",modnm);
     //cout<<"ok000\n";
     read_dispMineos(model.data.Rdisp,Moutputnm,Nmod);//according to Moutputnm, store pvel and gvel into m.d.R/Ldisp.pvel/gvel
     //cout<<"ok00\n";
     model.data.AziampRdisp.pvel.clear();model.data.AziphiRdisp.pvel.clear();model.data.AziampRdisp.gvel.clear();model.data.AziphiRdisp.gvel.clear();
     cout<<"ok0\n";
     // bug fixed on Oct 30 2013. modification
     for( i=0;i<model.data.AziampRdisp.npper;i++){
     	model.data.AziampRdisp.pvel.push_back(0.);
     }
     for(i=0;i<model.data.AziphiRdisp.npper;i++)
     	model.data.AziphiRdisp.pvel.push_back(0.);
     for( i=0;i<model.data.AziampRdisp.ngper;i++){
     	model.data.AziampRdisp.gvel.push_back(0.);
     }
     for(i=0;i<model.data.AziphiRdisp.ngper;i++)
     	model.data.AziphiRdisp.gvel.push_back(0.);
  }//if Rsurflag

  
  return 1; 
}//compute_dispMineos

//--------------------------------------
int compute_diff(dispdef disp1, dispdef disp2, vector<vector<double> > &pveldiff,vector<vector<double> > &gveldiff)
{//veldiff: [T,v1,v2,v1-v2]
  int i;
  vector<double> tdiff;
  //para2mod

  //compute_disp

  //get diff
  pveldiff.clear();gveldiff.clear();
  //********check*******does npper!=0 means fphase>0??***********
  //printf("npper=%d fphase=%d ngper=%d fgroup=%d\n",disp1.npper,disp1.fphase,disp1.ngper,disp1.fgroup);
  //*****************
  if (disp1.npper!=disp2.npper or disp1.ngper!=disp2.ngper){cout<<"### in compute_diff, nT doesn't match!\n";printf("disp1.npper=%d disp2.npper=%d disp1.ngper=%d disp2.ngper=%g\n",disp1.npper,disp2.npper,disp1.ngper,disp2.ngper);exit(0);}
  for(i=0;i<disp1.npper;i++)
  {
	//***********check*********
	if(pow(disp1.pper[i]-disp2.pper[i],2)>0.1){printf("## in compute_diff, period lists don't match!\n###T1=%g T2=%g\n",disp1.pper[i],disp2.pper[i]);exit(0);}
	//************************
	tdiff.clear();
	tdiff.push_back(disp1.pper[i]);
	tdiff.push_back(disp1.pvel[i]);
	tdiff.push_back(disp2.pvel[i]);
	tdiff.push_back(disp1.pvel[i]-disp2.pvel[i]);
	pveldiff.push_back(tdiff);
  }//for i npper

  for(i=0;i<disp1.ngper;i++){
	//***********check*********
	if(pow(disp1.gper[i]-disp2.gper[i],2)>0.1){printf("## in compute_diff, period lists don't match!\n###T1=%g T2=%g\n",disp1.gper[i],disp2.gper[i]);exit(0);}
	//************************
	tdiff.clear();
	tdiff.push_back(disp1.gper[i]);
	tdiff.push_back(disp1.gvel[i]);
	tdiff.push_back(disp2.gvel[i]);
	tdiff.push_back(disp1.gvel[i]-disp2.gvel[i]);
	gveldiff.push_back(tdiff);
 	
  }//for i ngper

  return 1;
}//compute_diff


//---------------------------------------
int compute_Vkernel_single_para(paradef para, int i,modeldef model, vector<vector<double> > PREM,int Nprem, int Rflag, int Lflag,int flagupdaterho, vector<double> &trkp1, vector<double> &trkg1, vector<double> &tlkp1, vector<double> &tlkg1, int ipara){
 //comute the kernel for each paramter;
 // clear the vectors 
 // compute_dispMineos for the input_old_model before using this code, keep the disp in model updated

 modeldef newmodel;
 paradef newpara;
 float dp=0.01,ddp; //dp*100% perturbation
 int j;
 vector<vector<double> > DRpvel,DRgvel,DLpvel,DLgvel;
 
 newmodel=model;
 newpara=para;

 trkp1.clear();trkg1.clear();tlkp1.clear();tlkg1.clear();
 if(Rflag>0 and Lflag==0){
	ddp=para.parameter[i]*dp;
	newpara.parameter[i]=para.parameter[i]+ddp;
	para2mod_static(newpara,model,newmodel);
	updatemodel(newmodel,flagupdaterho);
 	compute_dispMineos(newmodel,PREM,Nprem,1,0,ipara);
	compute_diff(newmodel.data.Rdisp,model.data.Rdisp,DRpvel,DRgvel);
	for(j=0;j<model.data.Rdisp.npper;j++)
		trkp1.push_back(DRpvel[j][3]/ddp);
	for(j=0;j<model.data.Rdisp.ngper;)
		trkg1.push_back(DRgvel[j][3]/ddp);
 	for(j=0;j<model.data.Ldisp.npper;j++)
		tlkp1.push_back(0.);
	for(j=0;j<model.data.Ldisp.ngper;j++)
		tlkg1.push_back(0.);
  }
  else if (Rflag*Lflag>0){
	ddp=para.parameter[i]*dp;
	newpara.parameter[i]=para.parameter[i]+ddp;
	para2mod_static(newpara,model,newmodel);
	updatemodel(newmodel,flagupdaterho);
        printf("@@@ check, compute_Vkernel_single_para, para %g -->%g\n",para.parameter[i],newpara.parameter[i]);
 	compute_dispMineos(newmodel,PREM,Nprem,1,1,ipara);
	compute_diff(newmodel.data.Rdisp,model.data.Rdisp,DRpvel,DRgvel);
	compute_diff(newmodel.data.Ldisp,model.data.Ldisp,DLpvel,DLgvel);	
	//---check--
	//for(j=0;j<model.data.Rdisp.npper;j++)printf("@@@ check, compute_Vkernel_single_para, T=%g,vin=%g vold=%g vnew=%g DRpvel[%d][3]=%g\n",newmodel.data.Rdisp.pper[j],model.data.Rdisp.pvelo[j],model.data.Rdisp.pvel[j],newmodel.data.Rdisp.pvel[j],j,DRpvel[j][3]);
	//
	for(j=0;j<model.data.Rdisp.npper;j++)
                trkp1.push_back(DRpvel[j][3]/ddp);
        for(j=0;j<model.data.Rdisp.ngper;j++)
                trkg1.push_back(DRgvel[j][3]/ddp);
	for(j=0;j<model.data.Ldisp.npper;j++)
		tlkp1.push_back(DLpvel[j][3]/ddp);
	for(j=0;j<model.data.Ldisp.ngper;j++)
		tlkg1.push_back(DLgvel[j][3]/ddp);
  } 
  else if(Rflag==0 and Lflag>0){
	ddp=para.parameter[i]*dp;
	newpara.parameter[i]=para.parameter[i]+ddp;
	para2mod_static(newpara,model,newmodel);
	updatemodel(newmodel,flagupdaterho);
 	compute_dispMineos(newmodel,PREM,Nprem,0,1,ipara);
	compute_diff(newmodel.data.Ldisp,model.data.Ldisp,DLpvel,DLgvel);
	for(j=0;j<model.data.Rdisp.npper;j++)
		trkp1.push_back(0.);
	for(j=0;j<model.data.Rdisp.ngper;)
		trkg1.push_back(0.);
 	for(j=0;j<model.data.Ldisp.npper;j++)
		tlkp1.push_back(DLpvel[j][3]/ddp);
	for(j=0;j<model.data.Ldisp.ngper;j++)
		tlkg1.push_back(DLgvel[j][3]/ddp);	

  }
  //printf("@@@ check, compute_Vkernel_single_para, para %g -->%g\n",para.parameter[i],newpara.parameter[i]);
  return 1;
}//compute_Vkernel_single_para


//---------------------------------------
int compute_Vkernel(paradef para,modeldef model,vector<vector<vector<double> > > &kernel,vector<vector<double> > PREM,int Nprem, int Rflag,int Lflag, int flagupdaterho){ // ; BS
  // purterb the input para, then compute kernel
  //kernel: kernel[kRp[nP][nT],kRg[][],kLp[][],kLg[][]]
  //should keep the model(&its disp) updated before using this subroutine

  int j,ng,nv,ppflag,LVflag;
  vector<vector<double> > trkp2,trkg2,tlkp2,tlkg2;
  vector<double> trkp1,trkg1,tlkp1,tlkg1;
  vector<double> temp1,temp2,temp3,temp4;
  vector<double> Rp0(model.data.Rdisp.npper,0.),Rg0(model.data.Rdisp.ngper,0.),Lp0(model.data.Ldisp.npper,0.),Lg0(model.data.Ldisp.ngper,0.);

  for(j=0;j<para.npara;j++){
	trkp2.push_back(Rp0);
	trkg2.push_back(Rg0);
	tlkp2.push_back(Lp0);
	tlkg2.push_back(Lg0);
  }

  //ttk1.push_back(-999.0);
  //ttk2.push_back(ttk1);

  kernel.clear();

  //if(model.flag==0){updatemodel(model,flagupdaterho);}
  //compute_dispMineos(model,PREM,Nprem,Rflag,Lflag);// this has been done in the main program
  //printf("@@@ check, compute_Vkernel, from initial model\n");
  //for(j=0;j<model.data.Rdisp.npper;j++)printf("  @@@ check, T=%g,vin=%g vold=%g\n",model.data.Rdisp.pper[j],model.data.Rdisp.pvelo[j],model.data.Rdisp.pvel[j]);
//num_threads(1)
  #pragma omp parallel for default(none) shared(model,para,Rp0,Rg0,Lp0,Lg0,PREM,Nprem,Rflag,Lflag,flagupdaterho,trkp2,trkg2,tlkp2,tlkg2) private(ng,ppflag,LVflag,trkp1,trkg1,tlkp1,tlkg1) 
  for(int i=0;i<para.npara;i++){
        /*trkp1=Rp0;
	trkg1=Rg0;
	tlkp1=Lp0;
	tlkg1=Lg0;
	*/
	ng=(int)para.para0[i][4];
	ppflag=(int)para.para0[i][6];
	LVflag=(int)para.para0[i][7];

  	//printf("@@@ check, compute_Vkernel, begin3\n");
	if(LVflag==1 and para.space1[i][2]*(para.space1[i][2]-0.001)<0.){// if it won't be used to cpt LoveKernel and para don't need to be perturbed or scaled; then no need to cpt its partial deriv
		trkp1=Rp0;trkg1=Rg0;tlkp1=Lp0;tlkg1=Lg0;
		printf("@@@ check, compute_Vkernel, dealing with the %dth para, NO Vkernel cpt=====\n",i);
	}
	else if( model.groups[ng].flagcpttype==3 ){// cpt with V+AZ kernel 
		if(ppflag==10 or ppflag==11){// by default, ppflag is increasing, so ppflag==1 for this group has kernel computed already , BUT ATTENTION, THIS PART NEED CHANGE IN PARALLEL CASE!
			continue;
			printf("@@@ check, compute_Vkernel, dealing with the %dth para, NO Vkernel cpt=====\n",i);
			/*
			trkp1=temp1;
			trkg1=temp2;
			tlkp1=temp3;
			tlkg1=temp4;
			*/
		}
		else if(ppflag>5 and ppflag!=9)// not vsv~eta, and not h
			{
			printf("@@@ check, compute_Vkernel, dealing with the %dth para, NO Vkernel cpt=====\n",i);
			trkp1=Rp0;trkg1=Rg0;tlkp1=Lp0;tlkg1=Lg0;
			}
		else {
			printf("@@@ check, compute_Vkernel, dealing with the %dth para, DOOOOO Vkernel cpt=====\n",i);
			compute_Vkernel_single_para(para,i,model,PREM,Nprem, Rflag,Lflag,flagupdaterho,trkp1,trkg1,tlkp1,tlkg1,i);
			//if(ppflag==1){temp1=trkp1;temp2=trkg1;temp3=tlkp1;temp4=tlkg1;}
		}	

	}//if 3
       else if( model.groups[ng].flagcpttype==1 or model.groups[ng].flag==1 or model.groups[ng].flag==6 ){// cpt with Vkernel, or layered model, or point model; BS
		if(ppflag>5 and ppflag!=9)
			{
			printf("@@@ check, compute_Vkernel, dealing with the %dth para, NO Vkernel cpt=====\n",i);
			trkp1=Rp0;trkg1=Rg0;tlkp1=Lp0;tlkg1=Lg0;
			}
		else 
			{
			printf("@@@ check, compute_Vkernel, dealing with the %dth para, DOOOOO Vkernel cpt=====\n",i);
			compute_Vkernel_single_para(para,i,model,PREM,Nprem, Rflag,Lflag,flagupdaterho,trkp1,trkg1,tlkp1,tlkg1,i);
			}
	}//else if 1
	else if( model.groups[ng].flag==2 ){// cpt part/all with Lkernel, and the group is Bspline
		printf("####, compute_Vkernel, computation for this situation is still under construction! Sorry!\n");
		exit(0);
	}
	else{
		printf("####, compute_Vkernel, this is an unconsidered situation, consider it NOW!\nOR add limits to the checkParaModel function to prevent this situation from happening\n");
		exit(0);
	}//else

	# pragma critical (printout)
	{
  	printf("@@@ check, compute_Vkernel, end of the %d th para\n",i);
	trkp2[i]=trkp1;
	trkg2[i]=trkg1;
	tlkp2[i]=tlkp1;
	tlkg2[i]=tlkg1;
	}
	/* 
	trkp2.push_back(trkp1);
	trkg2.push_back(trkg1);
	tlkp2.push_back(tlkp1);
	tlkg2.push_back(tlkg1);
	*/
  }// for i<npara

  for(int i=0;i<para.npara;i++){// this part is especially for parallel running case; get Vkernel for AZcos, AZsin after all Vkernel have been computed.
	ng=(int)para.para0[i][4];
	ppflag=(int)para.para0[i][6];
	if(model.groups[ng].flagcpttype==3){
		if(ppflag==1){//Vsv
			nv=(int)para.para0[i][5];
			temp1=trkp2[i];
			temp2=trkg2[i];
			temp3=tlkp2[i];
			temp4=tlkg2[i];
		}
		else if (ppflag==10 or ppflag==11){//AZcos or AZsin
			if ((int)para.para0[i][5]!=nv){printf("@@@ check compute_Vkernel; the order of parameter has problem!");exit(0);}
			trkp2[i]=temp1;
			trkg2[i]=temp2;
			tlkp2[i]=temp3;
			tlkg2[i]=temp4;
		}				
	}//if flag
  }//for i<npara

  //---check---
  if(trkp2.size()!=para.npara){
	printf("@@@ check compute_Vkernel, the size of kernel(%d) isn't right(%d)!! \n",trkp2.size(),para.npara);
	exit(0);
  }
 /* for(int k=0;k<para.npara;k++){
	printf("ipara=%d size of Rp Rg Lp Lg= %d %d %d %d\n",k,trkp2[k].size(),trkg2[k].size(),tlkp2[k].size(),tlkg2[k].size());
  }*/
  kernel.push_back(trkp2);
  kernel.push_back(trkg2);
  kernel.push_back(tlkp2); 
  kernel.push_back(tlkg2);


  return 1;
}// compute_Vkernel



//#####################################################################################

//---------------------------------------
int computeLovekernel(vector<double> Vkernel1,double rho,double vel,vector<double> &Lkernel1, int flagupdaterho){ 
  //compute the Lkernel (dc/dX) from the Vkernel (dc/dV)

  //for ACLN: dc/dX=dc/d(rho*V^2)=dc/dV/(d(rho*V^2)/dV)=(dc/dV)/(V^2*drho/dV+2*rho*V)
  //SO, if rho is constant, won't change with V --> drho/dV=0 ==> dc/dX=(dc/dV)/(2*rho*V)
  //if not, then drho/dV=d(0.3601*vp+0.541)/dV=... and dc/dX=...

  //for F: dc/dF=(dc/dF)|A,L=(dc/d((A-2*L)*eta))|A,L=1/(A-2*L)*dc/deta
  int i;
  double c;
  Lkernel1.clear();

  if(flagupdaterho==1){//under construction
	printf("###, computeLovekernel, the flagupdaterho==1 situation is still under construction!\n ");
	exit(0);
  }
  else
	c=2*rho*vel;

  for(i=0;i<Vkernel1.size();i++){
	Lkernel1.push_back(Vkernel1[i]/c);
  }//for i
  
  return 1;
}//convertVkernel2Lovekernel
//---------------------------------------
int computeLovekerneleta(vector<double> Vkernel1,double c, vector<double> &Lkernel1){
  Lkernel1.clear();
  for(int i=0;i<Vkernel1.size();i++){
	//Lkernel1.push_back(Vkernel1[i]*c);
	Lkernel1.push_back(Vkernel1[i]);//check kernel F and eta
	//printf("@@@ check, computeLovekerneleta,ip=%d, kernel= %g c=%g\n",i,Vkernel1[i]*c,c); 
  }
  return 1;
}
//---------------------------------------

int Vkernel2Lkernel(paradef para,modeldef model,vector<vector<vector<double> > > Vkernel,vector<vector<vector<double> > > &Lkernel, int flagupdaterho){ // ; BS
// should keep the mod and para updated, use para2mod before using this subroutine;
// para2mod is done in the Vpara2Lovepara, by default, i think the model here is updated, so para2mod doesn't show up here
// this is for iso or TI model, due to the computation capability of MINEOS 
  int i,LVflag,ng,nv,ppflag;  
  vector<vector<double> > trkp2,trkg2,tlkp2,tlkg2;
  vector<double> trkp1,trkg1,tlkp1,tlkg1;
  vector<double> Rp0(model.data.Rdisp.npper,0.),Rg0(model.data.Rdisp.ngper,0.),Lp0(model.data.Ldisp.npper,0.),Lg0(model.data.Ldisp.ngper,0.);
  double rho,A,L,tA,tL,Vv;

  Lkernel.clear();

  for(i=0;i<para.npara;i++){
	trkp1.clear(); trkg1.clear();tlkp1.clear();tlkg1.clear();

	LVflag=(int)para.para0[i][7];
	if(LVflag==1){trkp1=Rp0;trkg1=Rg0;tlkp1=Lp0;tlkg1=Lg0;}// only need Vkernel,so Lkernel-->0
	else if(LVflag==-1){//need Lkernel
		ng=(int)para.para0[i][4];
		nv=(int)para.para0[i][5];
		ppflag=(int)para.para0[i][6];
                if(model.groups[ng].flag==1 or model.groups[ng].flag==6){//layered or point ; BS
			if(ppflag==9){//h
				trkp1=Vkernel[0][i];
				trkg1=Vkernel[1][i];
				tlkp1=Vkernel[2][i];
				tlkg1=Vkernel[2][i];
			}//if ppflag==9
			else if (ppflag>5){
				trkp1=Rp0;trkg1=Rg0;tlkp1=Lp0;tlkg1=Lg0;
			}
			else if (ppflag==5){//eta
				//the real computation of the dv/dF is inside the compute_RAdisp code
				/*
				rho=model.groups[ng].rhovalue[nv];
				A=rho*pow(model.groups[ng].vphvalue[nv],2.0);//value in the group.Xvalue
				L=rho*pow(model.groups[ng].vsvvalue[nv],2.0);
				//---check---
				printf("@@@ Vkernel2Lkernel, check A&L from m.g(%g & %g) vs that from para(%g & %g)\n",A,L,tA,tL);
				if(pow(A-tA,2)>0.01 or pow(L-tL,2)>0.01){
					printf("### Vkernel2Lkernel, A&tA or L&tL differ, bug in code!\n");
					exit(0);
				}
				A=1./(A-2*L);
				*/
				A=1.;
				computeLovekerneleta(Vkernel[0][i],A,trkp1);
				computeLovekerneleta(Vkernel[1][i],A,trkg1);
				computeLovekerneleta(Vkernel[2][i],A,tlkp1);
				computeLovekerneleta(Vkernel[3][i],A,tlkg1);
				
			}
			else{
				rho=model.groups[ng].rhovalue[nv];
				Vv=para.parameter[i];//value in the Vparameter (velocity)
				//Lv;//value in the Loveparameter (elastic coefficient)
				computeLovekernel(Vkernel[0][i],rho,Vv,trkp1,flagupdaterho);
				computeLovekernel(Vkernel[1][i],rho,Vv,trkg1,flagupdaterho);
				computeLovekernel(Vkernel[2][i],rho,Vv,tlkp1,flagupdaterho);
				computeLovekernel(Vkernel[3][i],rho,Vv,tlkg1,flagupdaterho);
				//---check---
				//if(ppflag==4){tA=rho*Vv*Vv;}
				//else if(ppflag==1){tL=rho*Vv*Vv;}
				//---
			}

			//printf("@@@ check, Vkernel2Lkernel, ppflag=%d============\n",ppflag);
			
		}
		else if (model.groups[ng].flag==2){//Bspline
			//under construction
			printf("####, Vkernel2Lkernel, computation for this situation is still under construction! Sorry!\n");
			exit(0);
		}
		else{
			printf("####, Vkernel2Lkernel, this is an unconsidered situation, consider it NOW!\nOR add limits to the checkParaModel function to prevent this situation from happening\n");
			exit(0);
		}	
	}// if LVflag==-1
	else{
		printf("### Vkernel2Lkernel, wrong LVflag!\n");
		exit(0);
	}
 	trkp2.push_back(trkp1);
	trkg2.push_back(trkg1);
	tlkp2.push_back(tlkp1);
	tlkg2.push_back(tlkg1);		  	

  }// for i<npara  
  Lkernel.push_back(trkp2);
  Lkernel.push_back(trkg2);
  Lkernel.push_back(tlkp2); 
  Lkernel.push_back(tlkg2);
  return 1;
}//Vkernel2Lkernel
//---------------------------------------------------------------
int read_kernel(paradef &para, modeldef &model,vector<vector<vector<double> > > &kernel, const char *fRker, const char *fLker,int &Rflag, int &Lflag,vector<vector<double> > PREM,int Nprem)
{
  //kernel: kernel[kRp[nP][nT],kRg[][],kLp[][],kLg[][]]
  // read in kernel from outside file, now it only read in Rph and Lph kernel.
  fstream mff;
  string line;
  vector<string> v;
  vector<vector<double> > tkp2,tkg2,ttk2;
  vector<double> tkp1,Rg0,Lg0,ttk1;
  int nT,nP,i,j;

 
  //compute_dispMineos(model,PREM,Nprem,Rflag,Lflag);

  for(i=0;i<model.data.Rdisp.ngper;i++)Rg0.push_back(0.);
  for(i=0;i<model.data.Ldisp.ngper;i++)Lg0.push_back(0.);
  //ttk1.push_back(-999.);ttk2.push_back(ttk1);  

  kernel.clear();
  if(Rflag>0){
	nT=model.data.Rdisp.npper;
	nP=para.npara;
//	kRp=new double [nP][nT];
	double** kRp= new double*[nP];
	for(i=0;i<nP;++i)kRp[i]=new double[nT];
	
	mff.open(fRker);
	if(not mff.is_open() ){
		printf("#### in read_kernel, cannot open file %s to read\n",fRker);
		return 0;	}
	i=0;
	while(getline(mff,line)){
	  v.clear();
	  Split(line,v," ");//T=v[0], kernel[][inpara][nT]=v[1~Npara]
	  if(v.size()!=nP+1){
		printf("### in read_kernel, the input kernel file has incorrect nP=%d, should be %d\n",v.size()-1,nP);
		return 0;
	  }// if lenv
          for(j=0;j<nP;j++){
		kRp[j][i]=atof(v[j+1].c_str());
	  }// for j
	  i=i+1;
	}//while 
	if (i!=nT){printf("#### in read_kernel, nT is inconsistent, i=%d nT=%d\n",i,nT);return 0;}
	mff.close();

	tkp2.clear();
	tkg2.clear();
	for(i=0;i<nP;i++){
	  tkp1.clear();
	  for(j=0;j<nT;j++)
		tkp1.push_back(kRp[i][j]);
	  tkp2.push_back(tkp1);
	  tkg2.push_back(Rg0);
	}//for i nP	
	kernel.push_back(tkp2);
	kernel.push_back(tkg2);

//	delete [][] kRp;	
	for(i=0;i<nP;++i)delete [] kRp[i];
	delete [] kRp;
  }//Rflag>0
  else {
	kernel.push_back(ttk2);kernel.push_back(ttk2);}
 
  if(Lflag>0){
        nT=model.data.Ldisp.npper;
        nP=para.npara;
//      kRp=new double [nP][nT];
        double** kLp= new double*[nP];
        for(i=0;i<nP;++i)kLp[i]=new double[nT];

        mff.open(fLker);
        if(not mff.is_open() ){
                printf("#### in read_kernel, cannot open file %s to read\n",fLker);
                return 0;       }
        i=0;
        while(getline(mff,line)){
          v.clear();
          Split(line,v," ");//T=v[0], kernel[][inpara][nT]=v[1~Npara]
          if(v.size()!=nP+1){
                printf("### in read_kernel, the input kernel file has incorrect nP=%d, should be %d\n",v.size()-1,nP);
                return 0;
          }// if lenv
          for(j=0;j<nP;j++){
                kLp[j][i]=atof(v[j+1].c_str());
          }// for j
          i=i+1;
        }//while 
        if (i!=nT){printf("#### in read_kernel, nT is inconsistent, i=%d nT=%d\n",i,nT);return 0;}
        mff.close();

        tkp2.clear();
        tkg2.clear();
        for(i=0;i<nP;i++){
          tkp1.clear();
          for(j=0;j<nT;j++)
                tkp1.push_back(kLp[i][j]);
          tkp2.push_back(tkp1);
          tkg2.push_back(Lg0);
        }//for i nP     
        kernel.push_back(tkp2);
        kernel.push_back(tkg2);

        for(i=0;i<nP;++i)delete [] kLp[i];
        delete [] kLp;

  }
  else{kernel.push_back(ttk2);kernel.push_back(ttk2);}

  vector<string>().swap(v);
  return 1;
}//read_kernel

//---------------------------------------
//write_kernel(Vkernel,model0,para1,dirlay,nodeid,lon,lat,Rsurflag,Lsurflag);
int write_kernel(vector<vector<vector<double> > > Vkernel,modeldef model0,paradef para1,char *kernelnmR, char *kernelnmL, int Rsurflag, int Lsurflag){
	//===write kernel===only write the kernel for phvel kernel[0] and kernel[2]
	FILE *fkernel;
	int i,j;
	//printf("test--- Vkernel.size=%d\n",Vkernel.size());
	if(Rsurflag>0){
       		//sprintf(tmpstr,"%s/kernelRp1ani_%s_%.1f_%.1f.txt",dirlay,nodeid,lon,lat); 
       		if((fkernel=fopen(kernelnmR,"w"))==NULL){
			printf("### write_kernel, cannot open file %s to write!!\n",kernelnmR);
			exit(0);
		}
       		
		/*Kernel[4][np][nT] test--
		for(i=0;i<para1.npara;i++){
			printf("for para %d, k[0][%d].size()=%d\n",i,i,Vkernel[0][i].size());
			for(j=0;j<model0.data.Rdisp.npper;j++){
				printf("   T%d, k[0][%d][%d]=%g\n",j,i,j,Vkernel[0][i][j]);
			}
		}
		*/
       		for(i=0;i<model0.data.Rdisp.npper;i++){
			fprintf(fkernel,"%f ",model0.data.Rdisp.pper[i]);
			for(j=0;j<para1.npara;j++){
				fprintf(fkernel," %g",Vkernel[0][j][i]);
			}
			fprintf(fkernel,"\n");
			
       		}
		//cout<<" check finish writting kernel file\n";
       		fclose(fkernel);
   	 }
    	if(Lsurflag>0){
       		//sprintf(tmpstr,"%s/kernelLp1ani_%s_%.1f_%.1f.txt",dirlay,nodeid,lon,lat); 
       		if((fkernel=fopen(kernelnmL,"w"))==NULL){
			printf("### write_kernel, cannot open file %s to write!!\n",kernelnmL);
			exit(0);
		}
       		for(i=0;i<model0.data.Ldisp.npper;i++){
			//printf("test--- write Lkernel %d\n",i);
			fprintf(fkernel,"%f ",model0.data.Ldisp.pper[i]);
			for(j=0;j<para1.npara;j++){
				fprintf(fkernel," %g",Vkernel[2][j][i]);
			}
			fprintf(fkernel,"\n");
        	}
       		fclose(fkernel);
    	}
  return 1;
}//write_kernel
//#####################################################################################
//---------------------------------------
int rotate_ET(Matrix<double,6,6> ETin, double theta, double phi, Matrix<double,6,6> &ETout){
//the -phi is for some interpretation's convenience. 
// in the x-n,y-e,z-down system, -phi ==> CW rotation from N
// in the x-s,y-e,z-up system, -phi==> CCW rotation from S
  Matrix3f a1m,a2m,a3m,am;
  Matrix<double, 6,6> Mm;
  double ox,oy,oz;
  double axx,axy,axz,ayx,ayy,ayz,azx,azy,azz;
  //--rotating angles---
  ox=theta*M_PI/180.;
  oy=0.;
  oz=(-1*phi)*M_PI/180.; 
  //--coordinate transformation matrix ---
  
  a3m<< cos(oz),sin(oz),0.,
	-sin(oz),cos(oz),0.,
	0.,0.,1.;
  a2m<< cos(oy),0.,-sin(oy),
	0.,1.,0.,
	sin(oy),0.,cos(oy);
  a1m<< 1.,0.,0.,
	0.,cos(ox),sin(ox),	
	0.,-sin(ox),cos(ox);
  //
  //am=a1m*a2m*a3m;//rotate Z, then Y, then X ==> due to the symmetry of T.I system, initial rotation around Z does nothing!
  am=a3m*a2m*a1m;//rotate the tensor CW around X, then CW around Y, then CW around Z == rotate the coordinate CCW around X, then CCW around the TRANSFORMED Y', then CCW around the TRANSFORMED Z''

  axx=am(0,0);axy=am(0,1);axz=am(0,2);
  ayx=am(1,0);ayy=am(1,1);ayz=am(1,2);
  azx=am(2,0);azy=am(2,1);azz=am(2,2);
 
  Mm<< axx*axx,       axy*axy,        axz*axz,        2*axy*axz,         2*axz*axx,              2*axx*axy,
        ayx*ayx,       ayy*ayy,        ayz*ayz,        2*ayy*ayz,         2*ayz*ayx,              2*ayx*ayy,
        azx*azx,       azy*azy,        azz*azz,        2*azy*azz,         2*azz*azx,              2*azx*azy,
        ayx*azx,       ayy*azy,        ayz*azz,        ayy*azz+ayz*azy,   ayx*azz+ayz*azx,        ayy*azx+ayx*azy,
        azx*axx,       azy*axy,        azz*axz,        axy*azz+axz*azy,   axz*azx+axx*azz,        axx*azy+axy*azx,
        axx*ayx,       axy*ayy,        axz*ayz,        axy*ayz+axz*ayy,   axz*ayx+axx*ayz,        axx*ayy+axy*ayx;
/*  cout<<"@@@a3m=\n"<<a3m<<"\na1m=\n"<<a1m<<"\nam=\n"<<am<<endl;
  cout<<"@@@ check, rotate_ET, rotating angle theta="<<theta<<"  phi="<<phi<<"\n"<<"the rotating matrix=\n"<<Mm<<endl;
*/
  ETout=(Mm*ETin)*Mm.transpose();//dot(dot(Mm,ETin),Mm.T);

  /*cout<<"\ntest--@@@ ETin=\n"<<ETin<<endl;
  printf("theta=%g phi=%g\n",theta,phi);
  cout<<"\ntset--@@@ ETout=\n"<<ETout<<endl;
  */
  /*
  Eigen::FullPivLU<Matrix2f> lua1(a1);
  Eigen::PartialPivLU<Matrix2f> lua1_2(a1);
  cout<<"a1 inverse LU=\n"<<lua1.inverse()<<endl;
  cout<<"a1 inverse partial LU=\n"<<lua1_2.inverse()<<endl;
  */

  return 1;
}//rotate_ET

//---------------------------------------
int ET2LoveCoeff(Matrix<double,6,6> ET,double RAcoeff[5],double AZcoeff[8][2]){
// compute the ACFLN,BHEG from Cij
	//double ET_eff[13]; // [A C L N F] [Bc Bs Ec Es Gc Gs Hc Hs]


	RAcoeff[0]=(ET(3,3)+ET(4,4))*.5;//L
	RAcoeff[1]=(ET(0,0)+ET(1,1))*.125-ET(0,1)*.25+ET(5,5)*.5;//N
	RAcoeff[2]=ET(2,2);//C
	RAcoeff[3]=.375*(ET(0,0)+ET(1,1))+ET(0,1)*.25+ET(5,5)*.5;//A
	RAcoeff[4]=(ET(0,2)+ET(1,2))*.5;//F


	int c=1;
	AZcoeff[0][0]=(ET(4,4)-ET(3,3))*.5;//Gc
	AZcoeff[0][1]=ET(4,3)*c;//Gs
	AZcoeff[1][0]=(ET(0,0)+ET(1,1))*.125-ET(0,1)*.25-ET(5,5)*.5;//Ec
	AZcoeff[1][1]=(ET(0,5)-ET(1,5))*.5*c;//Es
	AZcoeff[2][0]=AZcoeff[2][1]=0.;
	AZcoeff[3][0]=(ET(0,0)-ET(1,1))*.5;//Bc
	AZcoeff[3][1]=ET(0,5)+ET(1,5)*c;//Bs
	AZcoeff[4][0]=(ET(0,2)-ET(1,2))*.5;//Hc
	AZcoeff[4][1]=ET(2,5)*c;//Hs
	
	/*===check--test---
	printf("@@@ check ET2LoveCoeff\n");
	for(int i=0;i<5;i++){
		printf(" the %dth RApara=%g, AZpara cos=%g sin=%g\n",i,RAcoeff[i],AZcoeff[i][0],AZcoeff[i][1]);
	}*/
	
	return 1;
}// ET2LoveCoeff

//---------------------------------------
//int Vpara2ET(vector<double> Vparameter, double[6][6] &ET, vector<double> eff_Vparameter){
int Vpara2ET2LoveCoeff(vector<double> Vparameter, Matrix<double,6,6> &ET,double RAcoeff[8], double AZcoeff[8][2]){
  // from Vparameter --> get TI ET, --> rotate TI ET --> ET --> RAcoeff&AZcoeff;
  //Vparameter: (vsv,vsh,vpv,vph,eta,theta,phi,rho)
  //		 0   1   2   3    4   5    6   7   
  int i;
  //double** ETTI = new double*[6];
  //for(i=0;i<6;i++)ETTI=new double[6];
  Matrix<double,6,6> ETTI;   

  double A,C,F,L,N,rho;
//  cout<<"@@@ check, Vpara2ET, begin to fill ACFLN\n";
  L=Vparameter[7]*pow(Vparameter[0],2.0);
  N=Vparameter[7]*pow(Vparameter[1],2.0);
  C=Vparameter[7]*pow(Vparameter[2],2.0);
  A=Vparameter[7]*pow(Vparameter[3],2.0);
  F=Vparameter[4]*(A-2*L);
  
  //printf("check--@@@ Vpara2ET2LoveCeoff: F=%g=%g*(%g-2*%g)\n",F,Vparameter[4],A,L);

  /*ETTI={
	{A,A-2N,F,0.,0.,0.},
	{A-2N,A,F,0.,0.,0.},
	{F,F,C,0.,0.,0.},
	{0.,0.,0.,L,0.,0.,0.},
	{0.,0.,0.,0.,L,0.,0.},
	{0.,0.,0.,0.,0.,L,0.},
	{0.,0.,0.,0.,0.,0.,N}	}
  */
  ETTI<< A,  A-2*N,  F,  0.,  0.,  0.,
	A-2*N,   A,  F,  0.,  0.,  0.,
	 F,      F,  C,  0.,  0.,  0.,
	0.,     0.,  0.,  L,  0.,  0.,
	0.,     0.,  0.,  0.,  L,  0.,
	0.,     0.,  0.,  0.,  0., N;

  //printf("\ntest--@@@ check Vpara2ET2LoveCoeff, ETTI initialization OK!\nrho=%g  vsv=%g  vsh=%g  vpv=%g  vph=%g  eta=%g\n",Vparameter[7],Vparameter[0],Vparameter[1],Vparameter[2],Vparameter[3],Vparameter[4]);
  rotate_ET(ETTI,Vparameter[5],Vparameter[6],ET);
  ET2LoveCoeff(ET,RAcoeff,AZcoeff);
	/*===check
	printf("@@@ check Vpara2ET2LoveCoeff, ET2LoveCoeff\n");
	for(int i=0;i<5;i++){
		printf(" the %dth RApara=%g, AZpara cos=%g sin=%g\n",i,RAcoeff[i],AZcoeff[i][0],AZcoeff[i][1]);
	}
	*/
  for(i=5;i<8;i++){
	RAcoeff[i]=Vparameter[i];
	AZcoeff[i][0]=AZcoeff[i][1]=Vparameter[i];
  }  
 
  /*eff_Vparameter.clear();
  rho=Vparameter[7];
  for(i=0;i<4;i++){//vsv~vph
	eff_Vparameter.push_back(pow(RAcoeff[i]/rho,0.5));
  }//for i
  eff_Vparameter.push_back(RAcoeff[4]/(RAcoeff[0]-2*RAcoeff[2])); //eta
  eff_Vparameter.push_back(0.);//theta
  eff_Vparameter.push_back(0.);//phi
  eff_Vparameter.push_back(rho);//rho
  eff_Vparameter.push_back(Vparameter[8]);//h
  */
  //for(i =0;i<6;i++)delete[] ETTI[i];
  //delete [] ETTI;
  return 1;
}// Vpara2ET2LoveCoeff

//---------------------------------------
int getGroupVPara(groupdef group, vector<vector<double> > &Vparameter2, int flagupdaterho){
// get the parameters(vsv,vsh,vpv,vph,eta,theta,phi,rho(,h)) that belong to group i; Vparameter contains nlayer(==group.np;) lines and 9 columns.
// this code should probably use the para instead of group to fill VVparameter; if using group, then the need to do para2mod everytime before using this code, in order to guarantee that the model is updated;

  int i;
  double tvp;
  vector<double> tparameter;
  if ((group.flag-1)*(group.flag-6)!=0){// not layered model or point model ; BS
	printf("#### getGroupVpara, wrong group type, not layerd model!\n");
	exit(0);
  }
  Vparameter2.clear();
  if(flagupdaterho){
    for (i=0;i<group.np;i++){
	tparameter.clear();
	tparameter.push_back(group.vsvvalue[i]);
	tparameter.push_back(group.vshvalue[i]);
	tparameter.push_back(group.vpvvalue[i]);
	tparameter.push_back(group.vphvalue[i]);
	tparameter.push_back(group.etavalue[i]);
	tparameter.push_back(group.thetavalue[i]);
	tparameter.push_back(group.phivalue[i]);
	if(group.flag==5){//water layer
		tparameter.push_back(1.02);
 	}
	else {
		tvp=0.5*(group.vpvvalue[i]+group.vphvalue[i]);
		if(tvp<7.5)
			tparameter.push_back(0.541+0.3601*tvp);
		else
			tparameter.push_back(3.35);
	}
	//tparameter.push_back(group.ratio[i]*group.thick);
	Vparameter2.push_back(tparameter);
    }//for  
  }//if
  else{
    for (i=0;i<group.np;i++){
	tparameter.clear();
	tparameter.push_back(group.vsvvalue[i]);
	tparameter.push_back(group.vshvalue[i]);
	tparameter.push_back(group.vpvvalue[i]);
	tparameter.push_back(group.vphvalue[i]);
	tparameter.push_back(group.etavalue[i]);
	tparameter.push_back(group.thetavalue[i]);
	tparameter.push_back(group.phivalue[i]);
	tparameter.push_back(group.rhovalue[i]);
	//tparameter.push_back(group.ratio[i]*group.thick);
	Vparameter2.push_back(tparameter);
    }  
  }//else
  return 1; 
}// getGroupVPara
//---------------------------------------

int Vpara2Lovepara(paradef &para,modeldef &model,int flagupdaterho){
  //from the parameter to LoveRA/AZparameter
  //
  int i,j,k,Nlovelayer=0,nvmax=0,NLg=0,id,LVflag,ppflag,ng,nv;
  vector<double> LAZp0(2,0.),Vparameter1(8,0.);
  vector<vector<double> > Vparameter2;
  vector<int> laygpid(model.ngroup,-1),Bspgpid(model.ngroup,-1);

  //printf("@@@ check1, Vpara2Lovepara, theta=%g\n",model.groups[1].thetavalue[0]);
  modeldef temp;
  temp=model;
  
  para2mod(para,temp,model);//if p2m was done in the main function, this step would be redundant
  //printf("@@@ check2, Vpara2Lovepara, theta=%g\n",model.groups[1].thetavalue[0]);
  para.LoveRAparameter.clear();
  para.LoveAZparameter.clear();
  for(i=0;i<para.npara;i++){
  	para.LoveRAparameter.push_back(0.);
  	para.LoveAZparameter.push_back(LAZp0);
  }//for i	

  //--obtain parameter that will be used to define the RA/AZparameter3
  for(i=0;i<model.ngroup;i++){
	//---get the info for groups that need Lpara and is layered ---  
        if((model.groups[i].flagcpttype-2)*(model.groups[i].flagcpttype-4)==0 and (model.groups[i].flag-1)*(model.groups[i].flag-6)==0){// group need Love para to do forward cpt, and group is layered or point ; BS
		if(model.groups[i].np>nvmax)nvmax=model.groups[i].np;
		laygpid[i]=NLg;
		NLg++;
	}

	//---groups that need Lpara and is Bspline---
	//-- under construction;
  }//for i
  //printf("@@@ check3, Vpara2Lovepara, theta=%g\n",model.groups[1].thetavalue[0]);
  //===check
  //printf("@@@check, Vpara2Lovepara, # of groups that need Lpara is %d, and maximum layers(gp.np) among these groups is %d\n",NLg,nvmax);
  //

  if(NLg<1){
  	//printf("----Vpara2Lovepara, no group needs Lovepara---\n");
	return 1;
  }
  //--get RA/AZ parameters for layered or piont groups: for every group(m.ngroup) and layer/point(g.np) ; BS
  double RAparameter3[NLg][nvmax][8],AZparameter3[NLg][nvmax][8][2];
  double RAparameter1[8],AZparameter1[8][2];
  Matrix<double,6,6> ET;

  for(i=0;i<model.ngroup;i++){
	id=laygpid[i];
	if(id>=0){
		getGroupVPara(model.groups[i],Vparameter2,flagupdaterho);
		/*===check
		printf("@@@ check, Vpara2Lovepara, getGroupVPara\n");
		for(j=0;j<Vparameter2.size();j++){
			for(k=0;k<Vparameter2[j].size();k++){
				printf("  In group %d, layer %d, the %dth parameter is %g\n",i,j,k,Vparameter2[j][k]);
			}
		}
		*/
		for(j=0;j<model.groups[i].np;j++){
			Vpara2ET2LoveCoeff(Vparameter2[j],ET,RAparameter1,AZparameter1);
			for(k=0;k<8;k++){
				RAparameter3[id][j][k]=RAparameter1[k];
				AZparameter3[id][j][k][0]=AZparameter1[k][0];
				AZparameter3[id][j][k][1]=AZparameter1[k][1];
			}
		
			/*---check---
			printf("\n@@@ check, Vpara2Lovepara, group%d nv%d:\n",i,j);
			for(k=0;k<8;k++){
			  printf(" k=%d RA: %5g, AZ:cos %5g sin %5g\n ",k,RAparameter1[k],AZparameter1[k][0],AZparameter1[k][1]);
			}*/
			
		}//for j<np
	}//if laygpid>0
  }//for i
  //printf("@@@ check4, Vpara2Lovepara, theta=%g\n",model.groups[1].thetavalue[0]);

  //---get RA/AZ parameters for Bspline groups
  //---under construction

  //-- fill the LoveRA/AZparameter vector;
  for(i=0;i<para.npara;i++){
  	LVflag=(int)para.para0[i][7];
  	if (LVflag==1)continue; // leave the parameter i as 0 (filled at the beginning of this subroutine)
	ng=(int)para.para0[i][4];
	nv=(int)para.para0[i][5];
	ppflag=(int)para.para0[i][6];
	if(ppflag==9){//thickness
		para.LoveRAparameter[i]=para.parameter[i];
		para.LoveAZparameter[i][0]=para.parameter[i];
		para.LoveAZparameter[i][1]=para.parameter[i];
	}//if
	else if(ppflag>9)continue;
	else{//ppflag 1~8
		id=laygpid[ng];
                if(id>=0){// this is layered or point group and has Lpara ; BS
			para.LoveRAparameter[i]=RAparameter3[id][nv][ppflag-1];
			para.LoveAZparameter[i][0]=AZparameter3[id][nv][ppflag-1][0];
			para.LoveAZparameter[i][1]=AZparameter3[id][nv][ppflag-1][1];
		}//if laygpid
		//else if(Bspgpid[ng]>0){// this is B-spline group and has Lpara
		//	//under construction
		//}//else if
		else{
			printf("###, Vpara2Lpara, unconsidered situation! gp %d nv %d ppflag %d LVflag=%d model.gp.flagcpttype=%d\n",ng,nv,ppflag,LVflag,model.groups[ng].flagcpttype);
			exit(0);
		}
	}//else
  }//for i
  //printf("@@@ check5, Vpara2Lovepara, theta=%g\n",model.groups[1].thetavalue[0]);
  //printf("@@@ check Eta, p.LRApara=%g, Vpara2=%g\n",para.LoveRAparameter[5],Vparameter2[0][4]);

  // 
  return 1;
}//Vpara2Lovepara

//---------------------------------------
int Lovepara2Vpara(paradef &para, modeldef model)
{ //from the RA part of the Love para, get vsv~vph,eta & theta=pho=0;
  int i,ng,nv,ppflag;
  double rho,L=-1.,A=-1.;
  //cout<<"check N_Vpara vs N_Lpara"<<para.npara<<" "<<para.parameter.size()<<" "<<para.LoveRAparameter.size()<<endl;
  for(i=0;i<para.npara;i++){
	ng=(int)para.para0[i][4];
	if(model.groups[ng].flagcpttype==2 or model.groups[ng].flagcpttype==4){//this group needs to cpt RA disp using the Vparameters, so this Lpara2Vpara is needed, (for the cpttype==4 case, the computation of F_partial_deriv needs information of dL, dA, dEta, so need the effective Vpara!)
		ppflag=(int)para.para0[i][6];
		if((ppflag-7)*(ppflag-6)==0){//ppflag=6 or 7, theta or phi
			para.parameter[i]=0.;
		 }
		else if(ppflag>7)continue;//not related to Lpara, doesn't need to be changed
		else if(ppflag==5){//eta; by default, should meet L and A before meet F;
			if(L<0 or A<0){
				printf("### Lovepara2Vpara, the L and A haven't been defined!!\n");	
				exit(0);
			}

			para.parameter[i]=para.LoveRAparameter[i]/(A-2*L);
			L=-1;A=-1;
		}
		else{
			nv=(int)para.para0[i][5];
			rho=model.groups[ng].rhovalue[nv];
			//printf("@@@ check, Lovepara2Vpara: para%d %g==> %g sqrt(%g/%g)\n",i,para.parameter[i],sqrt(para.LoveRAparameter[i]/rho),para.LoveRAparameter[i],rho);
			para.parameter[i]=sqrt(para.LoveRAparameter[i]/rho);
			if(ppflag==1)
				L=para.LoveRAparameter[i];
			else if(ppflag==4)
				A=para.LoveRAparameter[i];
		}//else
	}//if flagcpttype==2
  }//for i<npara

  return 1;
}
//---------------------------------------
int compute_RAdisp(modeldef &model, paradef para, modeldef refmod, paradef refpara, vector<vector<vector<double> > > Vkernel, vector<vector<vector<double> > > Lkernel, int Rsurflag, int Lsurflag){
  int i,j,k,nT,nP,LVflag,ng;
  vector<double> paradiff;
  vector<double> tdisp, *tdisplist[4];  
  vector<vector<double> > tker,trefdisp;
  vector<vector<vector<double> > > kernel;
  vector<int> Rflag(para.npara,0),Lflag(para.npara,0);
  //vector<double> cflag(para.npara,1.);
  int nTlist[4];
  double tvel;

  kernel=Vkernel;  

  double Fref,F,Aref, Lref,c,Etaref,dA,dL,dEta,Eta,L,A;//check
  int ppflag; //check
  FILE *ftemp;//check

  //printf("@@@ check, compute_RAdisp\n");
  //-- prepare para and kernel list ---
  for(i=0;i<para.npara;i++){
	LVflag=(int)para.para0[i][7];
	ng=(int)para.para0[i][4];
	ppflag=(int)para.para0[i][6];
	if((int)para.para0[i][8]==1)Rflag[i]=1*Rsurflag; //parai8 could be 2
	if((int)para.para0[i][9]!=0)Lflag[i]=1*Lsurflag; //parai9 could be -1
	
	
	if(LVflag==1 or model.groups[ng].flagcpttype==2){
		if((ppflag-10)*(ppflag-11)==0) // a bug, modified on Nov 5, 2013
			paradiff.push_back(0.);
		else
			paradiff.push_back(para.parameter[i]-refpara.parameter[i]);
	}
	else{
		paradiff.push_back(para.LoveRAparameter[i]-refpara.LoveRAparameter[i]);
		//---check---kernel of F and eta
		//ppflag=(int)para.para0[i][6];
		if(ppflag==1){
			Lref=refpara.LoveRAparameter[i];
			L=para.LoveRAparameter[i];
			dL=L-Lref;
			}
		else if(ppflag==4){
			Aref=refpara.LoveRAparameter[i];
			A=para.LoveRAparameter[i];
			dA=A-Aref;
			}
		else if(ppflag==5){
			Fref=refpara.LoveRAparameter[i];
			F=para.LoveRAparameter[i];
			Etaref=refpara.parameter[i];
			Eta=F/(A-2*L);
			dEta=Eta-Etaref;
			}
		/*
		if((int)para.para0[i][6]==5){c=1./(A-2*L);}
		else {c=1;}
		cflag[i]=c;
		cout<<"para"<<i<<"  c="<<c<<endl;
		*/
		if(ppflag==5 and fabs(dA+dL+dEta)>1E-5){
			c=dEta/((Aref-2*Lref)*dEta+Etaref*(dA-2*dL));
			//printf("@@@ check, !!i=%d the c=%g  dEta=%g dA=%g dL=%g, Aref=%g Lref=%g Etaref=%g Eta=%g\n",i,c,dEta,dA,dL,Aref,Lref,Etaref,Eta);
			//printf("@@@ F=%g Fref=%g, dF/F=%g\n",F,Fref,(F-Fref)/Fref*100);
			//ftemp=fopen("temp_c.txt","a");
			//fprintf(ftemp,"%g\n",c);
			//fclose(ftemp);	
			for(j=0;j<4;j++){
				for(k=0;k<kernel[j][i].size();k++){
					kernel[j][i][k]=kernel[j][i][k]*c;
				}
			}//for j
		}//if ppflag==5
		else{
			for(j=0;j<4;j++)
				kernel[j][i]=Lkernel[j][i];		
		}
		//for(j=0;j<4;j++)
		//	kernel[j][i]=Lkernel[j][i];		
	}//else
  }//for	

  /*---check---
  printf("@@@ check, compute_RAdisp, RAparadiff:\n");
  for(i=0;i<para.npara;i++){
	  if((int)para.para0[i][7]==1 or  model.groups[(int)para.para0[i][4]].flagcpttype==2)
	  	printf("\trefpara%d=%8g, paradiff=%8g, diff(%)=%g\n",i,refpara.parameter[i],paradiff[i],paradiff[i]/refpara.parameter[i]*100);
	  else
	  	printf("\trefpara%d=%8g, paradiff=%8g, diff(%)=%g L\n",i,refpara.LoveRAparameter[i],paradiff[i],paradiff[i]/refpara.LoveRAparameter[i]*100);
  }
  */

  tdisplist[0]=&model.data.Rdisp.pvel;
  tdisplist[1]=&model.data.Rdisp.gvel;
  tdisplist[2]=&model.data.Ldisp.pvel;
  tdisplist[3]=&model.data.Ldisp.gvel;
  trefdisp.push_back(refmod.data.Rdisp.pvel);
  trefdisp.push_back(refmod.data.Rdisp.gvel);
  trefdisp.push_back(refmod.data.Ldisp.pvel);
  trefdisp.push_back(refmod.data.Ldisp.gvel);

  nP=para.npara;
  nTlist[0]=model.data.Rdisp.npper;nTlist[1]=model.data.Rdisp.ngper;nTlist[2]=model.data.Ldisp.npper;nTlist[3]=model.data.Ldisp.ngper;

  for(i=0;i<2;i++){
	tdisp.clear();

	tker=kernel[i];
	nT=nTlist[i];
	
	//--check--
	//if(tker[0].size()!=nTlist[i]){cout<<"### in cpt_dispKernel, wrong dimension in Kernel!\n"<<"nT_ker nT_disp:"<<tker[0].size()<<" "<<nTlist[i]<<endl;exit(0);}
	//if(nP!=paradiff.size()){cout<<"### in cpt_dispKernel, wrong dimension in Kernel!\n"<<"nP_ker nP_para:"<<nP<<" "<<paradiff.size()<<endl;exit(0);}
	//


	for(j=0;j<nT;j++){
		tvel=(trefdisp[i][j]);
		for(k=0;k<nP;k++){
			tvel=tvel+paradiff[k]*tker[k][j]*Rflag[k];//*cflag[k];			
		}
		tdisp.push_back(tvel);
	}//for j

	*(tdisplist[i])=tdisp;
  }//for i

  for(i=2;i<4;i++){
	tdisp.clear();

	tker=kernel[i];
	nT=nTlist[i];
	

	for(j=0;j<nT;j++){
		tvel=(trefdisp[i][j]);
		for(k=0;k<nP;k++){
			tvel=tvel+paradiff[k]*tker[k][j]*Lflag[k];			
		}
		tdisp.push_back(tvel);
	}//for j

	*(tdisplist[i])=tdisp;
  }//for i

  return 1;	
}//compute_RAdisp
//---------------------------------------
//---------------------------------------
int cs2ap(double Ac,double As, double &amp,double &phi, int phiflag){
// convert cos and sin to the amplitue and phi where phi is the fast axis (I think it's CW from N, not sure...)
// Ac*cos(Xt)+As*sin(Xt)=sqrt(A^2+B^2)*sin(Xt+phi)=sqrt(A^2+B^2)*sin(X(t+phi')); phi=atan(Ac/As), if As/cos(phi)<0 then phi=phi+pi;  phi'=phi/X 
//# Xpsi: fitting_curve=A*sin(X*(t+Xpsi))+C ==> fast_axis = pi/2/X-Xpsi
//# this can be understood if you draw a picture; +Xpsi moves the sin curve left by Xpsi, and 1/4 period of the sin curve is pi/2/X; the positive peak was at pi/2/X,then is moved left by Xpsi --> the positive peak is at pi/2/X=Xpsi
//

  // Rayleigh can also have 4-phi azi!!	
  //if(RLflag<2)phiflag=2;//for Rayleigh, only 2-phi azi; otherwise(RLflag>=2), Love wave can be either 2 or 4-phi azi, indicated by the input phiflag;

  double T=M_PI*2/phiflag;

  amp=sqrt(Ac*Ac+As*As);
/*  if(pow(Ac+Ac,2)<1E-10){//Ac=As=0
	//printf("  @@@ check, cs2ap, Ac=As=0\n");
  	phi=0.; // the phi here has no meaning since amp=0.
  }
  else if(As*As<1E-10){//As=0 --> Ac*cos(Xt)=Ac*sin(X(t+T/4))
	//printf("  @@@ check, cs2ap, As=0\n");
  	phi=T/4;
	if(Ac<0)
		phi=phi+T/2;
  }
  else if(Ac*Ac<1E-10){//Ac=0 --> As*sin(Xt)
	//printf("  @@@ check, cs2ap, Ac=0\n");
  	phi=0.;
	if(As<0)
		phi=T/2;
  }
  else{
	//printf("  @@@ check, cs2ap, normal\n");
  	phi=atan(Ac/As);
  	if(As/cos(phi)<0.)
		phi=phi+M_PI;
	phi=phi/phiflag;
  }
*/
  phi=atan2(Ac,As);
  phi=phi/phiflag;
  phi=T/4.-phi; //fast axis direction

  while(phi>T)
	  phi=phi-T;
  while(phi<0)
	  phi=phi+T;	  
	  
  phi=phi*180./M_PI; //rad2deg

  return 1;
}//cs2ap


//---------------------------------------
template <typename T> int sgn(T val) {
	    return (T(0) < val) - (val < T(0));
}
//---------------------------------------
int compute_AZdisp(modeldef &model, paradef para, modeldef refmod, paradef refpara, vector<vector<vector<double> > > Vkernel, vector<vector<vector<double> > > Lkernel, int Razflag, int Lazflag){

  //right now, this subroutine only consider phvel AZimuthal (amp&phi) disp curves; gpvel is not taken into consideration right now
  //Also, the code considers either the 2-phi azimuthal component, or the 4-phi component for Love waves based on the choice of Lphi
  //ATTENTION, the Razflag should be max(RAZampflag, RAZphiflag) and so as the Lazflag.
  //since I think the angle in Montagner&Nataf(1986) paper is CCW from S, (not CW from N, as they stated in the paper), I need to modify the phi computed from Loveparameters, but no need to do modifications to the phi computed from AZcos, AZsin. 

  int i,j,k,LVflag,ppflag,nT,nP,Lphi,Rphi;
  vector<double> Cosparadiff,Sinparadiff;
  vector<double> tampdisp,tphidisp, *tdisplist[4]; 
  vector<vector<double> > tker,trefdisp;
  vector<vector<vector<double> > > kernel;
  vector<int> Rflag(para.npara,0),Lflag(para.npara,0),phiflag(para.npara,0);
  double sin,cos,tamp,tphi,amp,phi;
  double Aref,Lref,Etaref,dA,dL,dEta,F,Fref,c,A,L,Eta;
  int nTlist[4];

  //printf("@@@ check, compute_AZdisp\n");
  //########## IMPORTANT PARAMETER ########### CONTROL WHICH RAYLEIGH-WAVE LOVE-WAVE AZI TO COMPUTE 
  Lphi=4;
  Rphi=2;
  //also, important is the period setup in the compute_misfitDISP_single_phi function.
  //########################
  kernel=Vkernel;
  //-- prepare para and kernel list ---
  for(i=0;i<para.npara;i++){
	LVflag=(int)para.para0[i][7];
	ppflag=(int)para.para0[i][6];
	//if(para.para0[i][8]*para.para0[i][10]>0)Rflag[i]=1*Razflag;
	if(para.para0[i][8]>0 and (int)para.para0[i][10]==Rphi)Rflag[i]=1*Razflag;
	if((int)(fabs(para.para0[i][9]*para.para0[i][10]))==Lphi)Lflag[i]=1*Lazflag*sgn(para.para0[i][9]);/// modification here?
	//printf("@@@ check, the %dth Rflag=%d Lflag=%d,p6=%g,para10=%g, sng(p9)=%d !!!!!!!!!\n",i,Rflag[i],Lflag[i],para.para0[i][6],para.para0[i][10],sgn(para.para0[i][9]));

	//printf("@@@ check, compute_AZdisp, %dth para, LVflag=%d, ppflag=%d\n",i,LVflag,ppflag);
	if(LVflag==1){
		if(ppflag==10){
			Cosparadiff.push_back(para.parameter[i]-refpara.parameter[i]);
			Sinparadiff.push_back(0.);
			//printf("@@@ check, compute_AZdisp, %dth para, cos diff=%g, sin diff=%g\n",i,para.parameter[i]-refpara.parameter[i],0);
		}
		else if (ppflag==11){
			Cosparadiff.push_back(0.);
			Sinparadiff.push_back(para.parameter[i]-refpara.parameter[i]);
			//printf("@@@ check, compute_AZdisp, %dth para, cos diff=%g, sin diff=%g\n",i,0,para.parameter[i]-refpara.parameter[i]);
			}
		else{Cosparadiff.push_back(0.);Sinparadiff.push_back(0.);
			//printf("@@@ check, compute_AZdisp, %dth para, cos diff=0, sin diff=0\n",i);
		}    	
	}//if LVflag==1

	else{
		Cosparadiff.push_back(para.LoveAZparameter[i][0]-refpara.LoveAZparameter[i][0]);
		Sinparadiff.push_back(para.LoveAZparameter[i][1]-refpara.LoveAZparameter[i][1]);
		//printf("@@@ check, compute_AZdisp, %dth para,Rflag=%d,Lflag=%d ppflag=%d, cos diff=%g-%g=%g, sin diff=%g-%g=%g (LVflag=-1 %g %g)\n",i,Rflag[i],Lflag[i],ppflag,para.LoveAZparameter[i][0],refpara.LoveAZparameter[i][0],para.LoveAZparameter[i][0]-refpara.LoveAZparameter[i][0],para.LoveAZparameter[i][1],refpara.LoveAZparameter[i][1],para.LoveAZparameter[i][1]-refpara.LoveAZparameter[i][1],para.LoveAZparameter[i][1],refpara.LoveAZparameter[i][1]);
		//---compute kernel of F ---
		if(ppflag==1){
			Lref=refpara.LoveRAparameter[i];
			L=para.LoveRAparameter[i];
			dL=L-Lref;
			}//F=-2*para.LoveRAparameter[i];}
		else if (ppflag==4){
			Aref=refpara.LoveRAparameter[i];
			A=para.LoveRAparameter[i];
			dA=A-Aref;
			}//F=para.LoveRAparameter[i]+F;}
		//modification!!!! may need to be modified further!
		else if (ppflag==5){
			Fref=refpara.LoveRAparameter[i];
			F=para.LoveRAparameter[i];
			Etaref=refpara.parameter[i];
			Eta=F/(A-2*L);
			dEta=Eta-Etaref; //F=para.parameter[i]*F;
		//printf("Fref=%g=%g*(%g-2*%g)\n",Fref,Etaref,Aref,Lref);
		//printf("@@@ check, Eta=%g Etaref=%g\n",Eta,Etaref);
		//printf("@@@ check, FfromLRApara=%g FreffromRApara=%g\n",para.LoveRAparameter[i],refpara.LoveRAparameter[i]);
		
		}
	
		if(ppflag==5 and fabs(dA+dL+dEta)>1E-5){
			c=dEta/((Aref-2*Lref)*dEta+Etaref*(dA-2*dL));
                        //printf("@@@ check, !!i=%d the c=%g  dEta=%g dA=%g dL=%g,Aref=%g Lref=%g Etaref=%g Eta=%g \n",i,c,dEta,dA,dL,Aref,Lref,Etaref,Eta);
                        //printf("@@@ F=%g Fref=%g,dF/F=%g\n",F,Fref,(F-Fref)/Fref*100);
                        for(j=0;j<4;j++){
                                for(k=0;k<kernel[j][i].size();k++){
                                        kernel[j][i][k]=kernel[j][i][k]*c;
                                }
                        }//for j
			
		}//if ppflag==5
		else{
			for(j=0;j<4;j++)
				kernel[j][i]=Lkernel[j][i];
		}//else ppflag==5
	}//else
  }//for i<npara 
  /*---check---
  printf("\n@@@ check, compute_AZdisp, AZparadiff:\n");
  for(i=0;i<para.npara;i++){
	  if((int)para.para0[i][7]==1)
	  	printf("\trefpara%d=(%8.3f,%8.3f), paradiff=(%8.3f,%8.3f)\n",i,refpara.parameter[i],refpara.LoveAZparameter[i][1],Cosparadiff[i],Sinparadiff[i]);
	  else
	  	printf("\trefpara%d=(%8.3f,%8.3f), paradiff=(%8.3f,%8.3f) Lovepara\n",i,refpara.LoveAZparameter[i][0],refpara.LoveAZparameter[i][1],Cosparadiff[i],Sinparadiff[i]);
  }
  */

  tdisplist[0]=&model.data.AziampRdisp.pvel;
  tdisplist[1]=&model.data.AziphiRdisp.pvel;
  tdisplist[2]=&model.data.AziampLdisp.pvel;
  tdisplist[3]=&model.data.AziphiLdisp.pvel;
  trefdisp.push_back(refmod.data.AziampRdisp.pvel);
  trefdisp.push_back(refmod.data.AziphiRdisp.pvel);
  trefdisp.push_back(refmod.data.AziampLdisp.pvel);
  trefdisp.push_back(refmod.data.AziphiLdisp.pvel);

  nP=para.npara;
  nTlist[0]=model.data.AziampRdisp.npper;nTlist[1]=model.data.AziphiRdisp.npper;
  nTlist[2]=model.data.AziampLdisp.npper;nTlist[3]=model.data.AziphiLdisp.npper;

   //test----
   char strsin[300],strcos[300];
  //---compute R AZdisp ---
  if(nTlist[0]==nTlist[1]){
	i=0;
	tampdisp.clear();
	tphidisp.clear();
	
	tker=kernel[i];
	nT=nTlist[i];

	for(j=0;j<nT;j++){
		cos=sin=0.;
		tamp=trefdisp[i][j];
		tphi=trefdisp[i+1][j];
		//sprintf(strsin,"nT=%d sin=  ",j);//test--
		//sprintf(strcos,"nT=%d cos=  ",j);//test--
		for(k=0;k<nP;k++){
			cos=cos+Cosparadiff[k]*tker[k][j]*Rflag[k];
			sin=sin+Sinparadiff[k]*tker[k][j]*Rflag[k];
			/*if(fabs(Sinparadiff[k]*tker[k][j]*Rflag[k])>1E-4){
				sprintf(strsin,"%s  \n\t(nP=%d) sin=%7.3f=%.3f*%.3f*%d  ",strsin,k,Sinparadiff[k]*tker[k][j]*Rflag[k],Sinparadiff[k],tker[k][j],Rflag[k]);
			}//---test---
			if(fabs(Cosparadiff[k]*tker[k][j]*Rflag[k])>1E-4){
				sprintf(strcos,"%s  \n\t(nP=%d) cos=%7.3f=%.3f*%.3f*%d  ",strcos,k,Sinparadiff[k]*tker[k][j]*Rflag[k],Cosparadiff[k],tker[k][j],Rflag[k]);
			}//---test---
			*/
		}//for k<nP
		cs2ap(cos,sin,amp,phi,Rphi);//convert cos,sin to amp,phi
;
		//---test---
		//amp=cos;
		//phi=sin;
		//
		tamp=tamp+amp;
		tphi=tphi+phi;
		tampdisp.push_back(tamp);
		tphidisp.push_back(tphi);
		//printf("\n\n%s\n%s\n",strsin,strcos);
		//printf("\nRayleigh: nT=%d cos=%7g sin=%7f ==> amp=%7f(tamp=%7f) phi=%7f\n",j,cos,sin,amp,tamp,phi);//--test--
	}//for j<nT


	*(tdisplist[i])=tampdisp;
	*(tdisplist[i+1])=tphidisp;
  }//if nT0==nT1
  else{
	i=0;
	tampdisp.clear();
	tker=kernel[i];
	nT=nTlist[i];
	for(j=0;j<nT;j++){
		cos=sin=0.;
		tamp=trefdisp[i][j];
		for(k=0;k<nP;k++){
			cos=cos+Cosparadiff[k]*tker[k][j]*Rflag[k];
			sin=sin+Sinparadiff[k]*tker[k][j]*Rflag[k];
		}
		cs2ap(cos,sin,amp,phi,Rphi);
		//printf("Rayleigh: nT=%d cos=%7g sin=%7f ==> amp=%7f phi=%7f\n",j,cos,sin,amp,phi);//--test--
		tamp=tamp+amp;
		tampdisp.push_back(tamp);
	}//for j<nT
	*(tdisplist[i])=tampdisp;

	i=1;
	tphidisp.clear();
	tker=kernel[i-1];
	nT=nTlist[i];
	for(j=0;j<nT;j++){
		cos=sin=0.;
		tphi=trefdisp[i][j];
		for(k=0;k<nP;k++){
			cos=cos+Cosparadiff[k]*tker[k][j]*Rflag[k];
			sin=sin+Sinparadiff[k]*tker[k][j]*Rflag[k];
		}
		cs2ap(cos,sin,amp,phi,Rphi);
		tphi=tphi+phi;
		tphidisp.push_back(tphi);
	}//for j<nT
	*(tdisplist[i])=tphidisp;

  }//else nT0==nT1

  //--compute L AZdisp--
  //--test--
  //for(i=0;i<nP;i++)printf("Lflag[%d]=%d Rflag=%d\n",i,Lflag[i],Rflag[i]);
  //
  if(nTlist[2]==nTlist[3]){
	i=2;
	tampdisp.clear();
	tphidisp.clear();
	
	tker=kernel[i];
	nT=nTlist[i];

	for(j=0;j<nT;j++){
		//printf("\n@@@ check, iT=%d\n",j);
		//sprintf(strsin,"nT=%d sin=  ",j);//test--
		//sprintf(strcos,"nT=%d sin=  ",j);//test--
		cos=sin=0.;
		tamp=trefdisp[i][j];
		tphi=trefdisp[i+1][j];
		for(k=0;k<nP;k++){
			cos=cos+Cosparadiff[k]*tker[k][j]*Lflag[k];
			sin=sin+Sinparadiff[k]*tker[k][j]*Lflag[k];
			/*--test---
			if(fabs(Lflag[k])>1E-4){//if(fabs(Sinparadiff[k]*tker[k][j]*Lflag[k])>1E-4){
				sprintf(strsin,"%s  \n\t(nP=%d) sin=%7.3f=%.3f*%.3f*%d  ",strsin,k,Sinparadiff[k]*tker[k][j]*Lflag[k],Sinparadiff[k],tker[k][j],Lflag[k]);
			}//---test---
			if(fabs(Lflag[k])>1E-4){//if(fabs(Cosparadiff[k]*tker[k][j]*Lflag[k])>1E-4){
				sprintf(strcos,"%s  \n\t(nP=%d) cos=%7.3f=%.3f*%.3f*%d  ",strcos,k,Sinparadiff[k]*tker[k][j]*Lflag[k],Cosparadiff[k],tker[k][j],Lflag[k]);
			}//---test---
			*/
		}
		cs2ap(cos,sin,amp,phi,Lphi);//convert cos,sin to amp,phi
		
		tamp=tamp+amp;
		tphi=tphi+phi;
		tampdisp.push_back(tamp);
		tphidisp.push_back(tphi);
		//printf("\n\n%s\n%s\n",strsin,strcos);//test---
		//printf("\nLove:  nT=%d cos=%7g sin=%7f ==> amp=%7f(tamp=%7f) phi=%7f\n",j,cos,sin,amp,tamp,phi);//--test--
	}//for j<nT


	*(tdisplist[i])=tampdisp;
	*(tdisplist[i+1])=tphidisp;
  }//if nT2==nT3
  else{
	i=2;
	nT=nTlist[i];
	if(nT>0){
	tampdisp.clear();
	tker=kernel[i];
	for(j=0;j<nT;j++){
		cos=sin=0.;
		tamp=trefdisp[i][j];
		for(k=0;k<nP;k++){
			cos=cos+Cosparadiff[k]*tker[k][j]*Lflag[k];
			sin=sin+Sinparadiff[k]*tker[k][j]*Lflag[k];
		}
		cs2ap(cos,sin,amp,phi,Lphi);
		tamp=tamp+amp;
		tampdisp.push_back(tamp);
	}//for j<nT
	*(tdisplist[i])=tampdisp;
	}

	i=3;
	nT=nTlist[i];
	if(nT>0){
	tphidisp.clear();
	tker=kernel[i-1];
	for(j=0;j<nT;j++){
		cos=sin=0.;
		tphi=trefdisp[i][j];
		for(k=0;k<nP;k++){
			cos=cos+Cosparadiff[k]*tker[k][j]*Lflag[k];
			sin=sin+Sinparadiff[k]*tker[k][j]*Lflag[k];
		}

		cs2ap(cos,sin,amp,phi,Lphi);
		tphi=tphi+phi;
		tphidisp.push_back(tphi);
	}//for j<nT
	*(tdisplist[i])=tphidisp;
	}

  }//else nT2==nT3
  return 1;
}//compute_AZdisp
//---------------------------------------

//---------------------------------------
int compute_dispKernel(modeldef &model,paradef para, modeldef refmodel, paradef refpara, vector<vector<vector<double> > > Vkernel, vector<vector<vector<double> > >  Lkernel, int Rflag, int Lflag,int Razflag, int Lazflag){
  
  compute_RAdisp(model,para,refmodel,refpara,Vkernel,Lkernel,Rflag,Lflag);
  compute_AZdisp(model,para,refmodel,refpara,Vkernel,Lkernel,Razflag,Lazflag);
  
  //compute_misfitDISP();	

  return 1;
}//compute_dispKernel 

//---------------------------------------

int get_misfitKernel(modeldef &model,paradef &para,modeldef refmodel, paradef refpara, vector<vector<vector<double> > > Vkernel, vector<vector<vector<double> > >  Lkernel,int Rflag, int Lflag, int AziampRflag, int AziampLflag, int AziphiRflag, int AziphiLflag, float inpamp, float inpphi, int flagupdaterho){
  Vpara2Lovepara(para,model,flagupdaterho);
  compute_dispKernel(model,para,refmodel,refpara,Vkernel,Lkernel,Rflag,Lflag,max(AziampRflag,AziphiRflag),max(AziampLflag,AziphiLflag));
  compute_misfitDISP(model,Rflag,Lflag,AziampRflag,AziampLflag,AziphiRflag,AziphiLflag,inpamp,inpphi);
  para.L=model.data.L;para.misfit=model.data.misfit;
}

//int compute_dispMineos(modeldef &model,vector<vector<double> > PREM,int Nprem, int Rsurflag,int Lsurflag)
//---------------------------------------
int get_misfitMineosRA(modeldef &model,paradef &para, vector<vector<double> > PREM,int Nprem, int Rflag,int Lflag, int flagupdaterho){
//only compute the misfit for RA dispersion curves. the azimuthal part cannot be computed with Mineos, need to use partial derivative
  modeldef tmodel;
  tmodel=model;
  Lovepara2Vpara(para,tmodel);
  para2mod(para,tmodel,model);
  updatemodel(model,flagupdaterho);
  compute_dispMineos(model,PREM,Nprem,Rflag,Lflag,0);
  compute_misfitDISP(model,Rflag,Lflag,0,0,0,0,0,0);
  para.L=model.data.L;para.misfit=model.data.misfit;
}

//-----------------------------------------------------
        int get_index(int var,vector<int> vv){
        //similar to the index function in python
          int i, id;
          id=-1;
          for(i=0;i<vv.size();i++){
                if(var==vv[i]){
                        id=i;
                        break;
                }
          }

          return id;
        }

//-----------------------------------------------------
	int Bsp2Point(modeldef BSmodel,paradef &BSpara,modeldef &Pmodel,paradef &Ppara,int flagupdaterho){//; BS
	// this is used to transfer the Bspline model/para to Point model/para; and this process will make the Vpara2Lpara, Vkernel2Lkernel ... process much easier.
	// represent the Bspline smooth model with Point (vel(z))    
	// cannot add & in front of BSpara(then the BSpara may change in the para2mod process), something strange will happen, not sure why.
	    modeldef ttmodel;
	    paradef tpara;
	    int i,count,ng,j,ig,nv,ppflag,id,tid,N,k;
	    vector< vector<int> > BSppflag2;
	    vector<int> BSng,BSppflag1;
	    vector<double> tpara0;

    	    para2mod(BSpara,BSmodel,ttmodel);
	    updatemodel(ttmodel,flagupdaterho);
	    Pmodel=ttmodel;
	    Pmodel.flag=0;//indicate Pmodel in not updated (the value1[] is not filled)
        /*--check--
        printf("\n\nBS2P, after p2m\n");
                for(int m=0;m<BSpara.npara;m++){
                                int tng=(int)BSpara.para0[m][4];
                                                if(tng==2 and ((int)BSpara.para0[m][6]-1)*((int)BSpara.para0[m][6]-2)==0){printf("para%d ppflag=%d, nv=%d %g\n",m,(int)BSpara.para0[m][6],(int)BSpara.para0[m][5],BSpara.parameter[m]);}
                }
        */

	    for (i=0;i<ttmodel.ngroup;i++){
		if(ttmodel.groups[i].flag!=3)continue;   //flag==3-->Bs need to be changed to point model 
		
	    	Pmodel.groups[i].vsvvalue=ttmodel.groups[i].vsvvalue1;
	    	Pmodel.groups[i].vshvalue=ttmodel.groups[i].vshvalue1;
	    	Pmodel.groups[i].vpvvalue=ttmodel.groups[i].vpvvalue1;
	    	Pmodel.groups[i].vphvalue=ttmodel.groups[i].vphvalue1;
	    	Pmodel.groups[i].etavalue=ttmodel.groups[i].etavalue1;
	    	Pmodel.groups[i].rhovalue=ttmodel.groups[i].rhovalue1;
	    	Pmodel.groups[i].thetavalue=ttmodel.groups[i].thetavalue1;
	    	Pmodel.groups[i].phivalue=ttmodel.groups[i].phivalue1;
		Pmodel.groups[i].thick1=ttmodel.groups[i].thick1;
		
	    	Pmodel.groups[i].np=ttmodel.groups[i].nlay; // # of parameters
		Pmodel.groups[i].flag=6;
		Pmodel.groups[i].flagBs=ttmodel.groups[i].flagBs;;
		Pmodel.groups[i].vpvs=ttmodel.groups[i].vpvs;
	    }//for i

	    //in para0[i]: 0)type_flag1 1)dv_type 2)dv 3)sigma 4)groupid 5)valueid 6)type_flag2 7)LVflag 8)RayleighWave_flag 9)Lovewave_flag 10)AZ_flag
	    //get the groups that need BS->point, and also the ppflag list of that group
	    for(i=0;i<BSmodel.ngroup;i++){
	    	if(BSmodel.groups[i].flag==3){
			BSng.push_back(i);
			BSppflag1.clear();
			for(j=0;j<BSpara.npara;j++){
				ng=(int)BSpara.para0[j][4];
				nv=(int)BSpara.para0[j][5];
				ppflag=(int)BSpara.para0[j][6];
				if(ng==i and nv==0)BSppflag1.push_back(ppflag);
			}//for j
			BSppflag2.push_back(BSppflag1);
	    	}//if flag==3
	    }//for i
	
	    //fill the Ppara.para0; order by group, then order by layer, then order by ppflag
	    count=0;
	    Ppara.para0.clear();
	    tid=-1;
	    for(i=0;i<BSpara.npara;i++){
	    	ng=(int)BSpara.para0[i][4];
		id=get_index(ng,BSng);
		if(id<0){//doesn't belong to a group that needs BS->point
			Ppara.para0.push_back(BSpara.para0[i]);
			count++;
		}
		else{//there will be nlay*N parameters in this group
		// 0)0 1)NM(no matter) 2)NM 3)matters,=that in BSpara 4)groupid 5)ilay_in_that_group 6)type_flag2 (indicate vsv,vsh... for detail,check INITstructure_Bs.h) 7)will be filled later?(WFL) 8)WFL? 9)WFL? 10)WFL?
			if(id==tid)continue;
			BSppflag1=BSppflag2[id];
			N=BSppflag1.size();
			//the dv, sigma of the para0 here is probably not useful, since we will only perturb the Bspling parameter, not these grip parameters
			for(j=0;j<ttmodel.groups[ng].nlay;j++){
				for(k=0;k<N;k++){
					if(j>0 and BSmodel.groups[ng].np>1 )
						tpara0=BSpara.para0[i+k+N];
					else
						tpara0=BSpara.para0[i+k];
					tpara0[0]=0;
					tpara0[5]=j;
					tpara0[6]=BSppflag1[k];
					Ppara.para0.push_back(tpara0);
					count++;
				}//for k
			}//for j
		}//else
	    	tid=id;
	    }//i<BSpara.npara
	    Ppara.npara=count;
		
	    /*---check --
	    printf("number of para in Ppara=%d\n",Ppara.npara);
	    for(i=0;i<Ppara.npara;i++)
	    	printf("parameter %d\n p1=%g  dv=%g sigma=%g ng=%g nv=%g pflag=%g\n",i,Ppara.para0[i][1],Ppara.para0[i][2],Ppara.para0[i][3],Ppara.para0[i][4],Ppara.para0[i][5],Ppara.para0[i][6]);
	    */
	    tpara=Ppara;
	    mod2para(Pmodel,tpara,Ppara);
	    return 1;
	}//Bsp2Point
