//##########this code may need revison in writting Mineos_model
//##########in this version, the compute kernel is modified on Sep20, 2012. previous perturbation is too small for moho thickness
// this parallel is unfinished, i've only parallelized the compute_dispMineos--Rflag*Lflag>0 part
// this rf version the misfit take the misfit of RF into account, the get_misfitKernel() have been changed

using namespace std;
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
  if(not Mf.is_open()){cout<<"### in compute_dispMineos, cannot open PREMfile "<<PREMnm<<endl;exit(0);}
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
/*  tmod.clear();
  tmod.push_back(inrad);
  tmod.push_back(model.laym0.rho[nlayer-1]*1000);
  tmod.push_back(model.laym0.vp[nlayer-1]*1000);
  tmod.push_back(model.laym0.vsv[nlayer-1]*1000);
  tmod.push_back(model.laym0.qp[nlayer-1]);
  tmod.push_back(model.laym0.qs[nlayer-1]);
  tmod.push_back(model.laym0.vp[nlayer-1]*1000);
  tmod.push_back(model.laym0.vsh[nlayer-1]*1000);
  tmod.push_back(1.0);
  Minput.push_back(tmod);
  N=N+1;
*/
  trad=inrad-model.laym0.thick[nlayer-1];
  double tmpth=0.;
  for (i=nlayer-1;i>-1;i--){ //revised on Sep 11, 2012
//  for (i=nlayer-1;i>0;i--){ //####################this seems to be wrong to me, shouldn't i reach 0 ??? in this way, the top layer is a single-velocity layer
	tmod.clear();
	tmpth=tmpth+model.laym0.thick[i]*1000.;
	tmod.push_back(trad+model.laym0.thick[i]*1000.);//0
	tmod.push_back(model.laym0.rho[i]*1000.);//1
  	tmod.push_back(model.laym0.vp[i]*1000.);//2
 	tmod.push_back(model.laym0.vsv[i]*1000.);//3
	tmod.push_back(model.laym0.qp[i]);//4
	tmod.push_back(model.laym0.qs[i]);//5
  	tmod.push_back(model.laym0.vp[i]*1000.);//6
  	tmod.push_back(model.laym0.vsh[i]*1000.);//7
  	tmod.push_back(1.0);//8
 	Minput.push_back(tmod);
	
/*	tmod[1]=model.laym0.rho[i-1]*1000.;
	tmod[2]=model.laym0.vp[i-1]*1000.;
	tmod[3]=model.laym0.vsv[i-1]*1000.;
	tmod[4]=model.laym0.qp[i-1];
	tmod[5]=model.laym0.qs[i-1];//test---------
	tmod[6]=model.laym0.vp[i-1]*1000.;
	tmod[7]=model.laym0.vsh[i-1]*1000.;
	Minput.push_back(tmod); 
*/
        trad=tmod[0];
	N=N+1;
	if (tmod[0]>6371000.0001){
	  printf("###### in write_modMineos, radius exceeds 6371000 m!!! %g tmpth=%g\n",tmod[0],tmpth);
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
int interpolate_model(layermoddef inlay, layermoddef &outlay, double dh)
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
int compute_dispMineos(modeldef &model,vector<vector<double> > PREM,int Nprem, int Rsurflag,int Lsurflag,int ipara)
{ //execute mineos_bran, give a model, compute dispersion curves for L and/or R, store value in model.d.R/Ldisp.pvel/gvel
  int Nmod,jcmp;
  string line;
  float wmin,wmax;
  vector<string> v;
  vector<double> tmod;
  char Minmodnm[200],Moutputnm[200];
  char str[500],moddir[100],modnm[100];
  dispdef Rdisp,Ldisp;
  //******parameters********
  wmin=1000./100.;//1000./70.;
  wmax=1000./6.;
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
  }//if Lsurflag
  if(Rsurflag>0){
     jcmp=3;
     sprintf(str,"csh run_Mineos_bran.csh %d %s %s %.1f %.1f\n",jcmp,moddir,modnm,wmin,wmax);
     system(str);
     //get disp info from Mineos output 
     sprintf(Moutputnm,"%s_S",modnm);
     read_dispMineos(model.data.Rdisp,Moutputnm,Nmod);//according to Moutputnm, store pvel and gvel into m.d.R/Ldisp.pvel/gvel
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
int compute_kernel(paradef para,modeldef &model,vector<vector<vector<double> > > &kernel,vector<vector<double> > PREM,int Nprem, int Rsurflag,int Lsurflag,int isoflag,float depcri1,float depcri2,float qpcri,float qscri)
{// purterb the input para, then compute kernel
  //kernel: kernel[kRp[nP][nT],kRg[][],kLp[][],kLg[][]]
  //***** may need to add one more case: isoflag==1.
  //the perturbation have been changed. modified on Sep 20, 2012
  int i,j;
  paradef newpara;
  modeldef newmodel;
  double dp,ddp;
  vector<vector<double> > ttk2,DRpvel,DRgvel,DLpvel,DLgvel,tkp2,tkg2,tlkp2,tlkg2;
  vector<double> ttk1,tkp1,tkg1,tlkp1,tlkg1;
  vector<double> Rp0,Rg0,Lp0,Lg0;
  for(i=0;i<model.data.Rdisp.npper;i++)Rp0.push_back(0.);
  for(i=0;i<model.data.Rdisp.ngper;i++)Rg0.push_back(0.);
  for(i=0;i<model.data.Ldisp.npper;i++)Lp0.push_back(0.);
  for(i=0;i<model.data.Ldisp.ngper;i++)Lg0.push_back(0.);
  ttk1.push_back(-999.0);
  ttk2.push_back(ttk1);
  dp=0.02; // 2% change

  kernel.clear();
  newmodel=model;
  
  if(para.Rnpara!=para.Lnpara and Rsurflag*Lsurflag!=0){printf("### in compute_kernel, Rnpara!=Lnpara!!\n");exit(0);}
    
  if(Rsurflag>0 and Lsurflag==0){
	printf(" compute kernel case Rsurflag>0, Lsurflag==0\n");  
	tkp2.clear();tkg2.clear();
	compute_dispMineos(model,PREM,Nprem,1,0,0);
	for(i=0;i<para.Rnpara;i++){
		tkp2.push_back(Rp0);tkg2.push_back(Rg0);
	}

	#pragma omp parallel for default(none) shared(model,para,PREM,Nprem,depcri1,depcri2,qpcri,qscri,dp,tkp2,tkg2,tlkp2,tlkg2)  private (i,j,newpara,newmodel,tkp1,tkg1,tlkp1,tlkg1,DRpvel,DRgvel,DLpvel,DLgvel,ddp)
	for(i=0;i<para.Rnpara;i++){
	  if(para.Rspace1[i][2]<0.001){// if the para don't need to be perturbed, then there is no need to compute its paratial deriv
	    continue;
	  }
	  newpara=para;
	  newmodel=model;
	  tkp1.clear();tkg1.clear();	 
	  //newpara.Rparameter[i]=para.Rparameter[i]+dp;
	  newpara.Rparameter[i]=para.Rparameter[i]*(1+dp);
	  ddp=para.Rparameter[i]*dp;
	  para2mod(newpara,model,newmodel);
	  updatemodelTibet(newmodel,depcri1,depcri2,qpcri,qscri);
	  compute_dispMineos(newmodel,PREM,Nprem,1,0,i);
	  compute_diff(newmodel.data.Rdisp,model.data.Rdisp,DRpvel,DRgvel);
	  for(j=0;j<model.data.Rdisp.npper;j++)
		tkp1.push_back(DRpvel[j][3]/ddp);
	  tkp2[i]=tkp1;
	  for(j=0;j<model.data.Rdisp.ngper;j++)
		tkg1.push_back(DRgvel[j][3]/ddp);
	  tkg2[i]=tkg1;
	}//for i Rnpara

	for(i=0;i<para.Lnpara;i++){
	   tkp2.push_back(Rp0);
	   tkg2.push_back(Rg0);	}//for i Lnpara
 
	kernel.push_back(tkp2);
	kernel.push_back(tkg2);
	kernel.push_back(ttk2);kernel.push_back(ttk2);
  }//Rsurflag
  else if(Rsurflag*Lsurflag>0){
	tkp2.clear();tkg2.clear();
	compute_dispMineos(model,PREM,Nprem,1,1,0);
	for(i=0;i<para.Rnpara+para.Lnpara;i++){
	      tkp2.push_back(Rp0);tkg2.push_back(Rg0);
	      tlkp2.push_back(Lp0);tlkg2.push_back(Lg0);
	}
	#pragma omp parallel for default(none) shared(model,para,PREM,Nprem,depcri1,depcri2,qpcri,qscri,dp,tkp2,tkg2,tlkp2,tlkg2) private (i,j,newpara,newmodel,tkp1,tkg1,tlkp1,tlkg1,DRpvel,DRgvel,DLpvel,DLgvel,ddp)
	for(i=0;i<para.Rnpara;i++){
	  if(para.Rspace1[i][2]<0.001){
		continue;
	  }
	  newpara=para;
	  newmodel=model;
	  tkp1.clear();tkg1.clear();tlkp1.clear();tlkg1.clear();
	  newpara.Rparameter[i]=para.Rparameter[i]*(1+dp);
	  ddp=para.Rparameter[i]*dp;
	  para2mod(newpara,model,newmodel);
	  updatemodelTibet(newmodel,depcri1,depcri2,qpcri,qscri);
	  compute_dispMineos(newmodel,PREM,Nprem,1,1,i);
          compute_diff(newmodel.data.Rdisp,model.data.Rdisp,DRpvel,DRgvel);
	  compute_diff(newmodel.data.Ldisp,model.data.Ldisp,DLpvel,DLgvel);
          for(j=0;j<model.data.Rdisp.npper;j++)
                {tkp1.push_back(DRpvel[j][3]/ddp);}
          tkp2[i]=tkp1;
          for(j=0;j<model.data.Rdisp.ngper;j++)
                tkg1.push_back(DRgvel[j][3]/ddp);
          tkg2[i]=tkg1;
	  for(j=0;j<model.data.Ldisp.npper;j++)
		tlkp1.push_back(DLpvel[j][3]/ddp);
	  tlkp2[i]=tlkp1;
	  for(j=0;j<model.data.Ldisp.ngper;j++)
		tlkg1.push_back(DLgvel[j][3]/ddp);
	  tlkg2[i]=tlkg1;
        }//for i Rnpara
	
	#pragma omp parallel for default(none) shared(model,para,PREM,Nprem,depcri1,depcri2,qpcri,qscri,dp,tkp2,tkg2,tlkp2,tlkg2) private (i,j,newpara,newmodel,tkp1,tkg1,tlkp1,tlkg1,DRpvel,DRgvel,DLpvel,DLgvel,ddp)
	for(i=0;i<para.Lnpara;i++){
	  if(para.Lspace1[i][2]<0.001){
		continue;
	  }	
	  if(atoi(para.Lpara0[i][0].c_str())==1){//group thickness, the effect of thickness have been taken into account in Rpara. 
	    	continue;	
	  }
	  newpara=para;
	  newmodel=model;
	  tkp1.clear();tkg1.clear();tlkp1.clear();tlkg1.clear();
	  newpara.Lparameter[i]=para.Lparameter[i]*(1+dp);
	  ddp=para.Lparameter[i]*dp;
	  para2mod(newpara,model,newmodel);
          updatemodelTibet(newmodel,depcri1,depcri2,qpcri,qscri);
	  compute_dispMineos(newmodel,PREM,Nprem,1,1,i+para.Rnpara);
	  compute_diff(newmodel.data.Rdisp,model.data.Rdisp,DRpvel,DRgvel);
	  compute_diff(newmodel.data.Ldisp,model.data.Ldisp,DLpvel,DLgvel);
	  for(j=0;j<model.data.Rdisp.npper;j++){tkp1.push_back(DRpvel[j][3]/ddp);}
          tkp2[i+para.Rnpara]=tkp1;
          for(j=0;j<model.data.Rdisp.ngper;j++)tkg1.push_back(DRgvel[j][3]/ddp);
          tkg2[i+para.Rnpara]=tkg1;
          for(j=0;j<model.data.Ldisp.npper;j++)tlkp1.push_back(DLpvel[j][3]/ddp);
          tlkp2[i+para.Rnpara]=tlkp1;
          for(j=0;j<model.data.Ldisp.ngper;j++)tlkg1.push_back(DLgvel[j][3]/ddp);
          tlkg2[i+para.Rnpara]=tlkg1;
	}//for i Lnpara

	kernel.push_back(tkp2);kernel.push_back(tkg2);
	kernel.push_back(tlkp2);kernel.push_back(tlkg2);
  }//Rflag and Lflag>0

  else if(Rsurflag==0 and Lsurflag>0){
	printf("Hey pragma hasn't been added here!\n");
	exit(0);
  	tkp2.clear();tkg2.clear();
	compute_dispMineos(model,PREM,Nprem,0,1,0);
	for(i=0;i<para.Rnpara;i++){
	  tkp2.push_back(Lp0);
	  tkg2.push_back(Lg0);
	}// for i Rnpara
	for(i=0;i<para.Lnpara;i++){
	  if(para.Lspace1[i][2]<0.001){//**** ATTENTION, here in this case, if you want to perturb the thickness, in Lpara.in remember to make the thickness's space to be non-zero. even though the thickness is controled by that of Rayleigh only.
		tkp2.push_back(Lp0);tkg2.push_back(Lg0);	      	
		continue;
	  }
	  newpara=para;
	  newmodel=model;
	  tkp1.clear();tkg1.clear();
	  if(atoi(para.Lpara0[i][0].c_str())==1)// the group thickness. para2mod cannot transfer thickness from Lgroup, but only from Rgroup. THIS PART MAY NEED CHANGE ******
	  {newpara.Rparameter[i]=para.Rparameter[i]*(1+dp);
	    ddp = para.Rparameter[i]*dp;}
	  else
	  {newpara.Lparameter[i]=para.Lparameter[i]*(1+dp);
	    ddp=para.Lparameter[i]*dp;}
	  para2mod(newpara,model,newmodel);
	  updatemodelTibet(newmodel,depcri1,depcri2,qpcri,qscri);
	  compute_dispMineos(newmodel,PREM,Nprem,0,1,i+para.Rnpara);
	  compute_diff(newmodel.data.Ldisp,model.data.Ldisp,DLpvel,DLgvel);
	  for(j=0;j<model.data.Ldisp.npper;j++)tkp1.push_back(DLpvel[j][3]/ddp);
	  tkp2.push_back(tkp1);
	  for(j=0;j<model.data.Ldisp.ngper;j++)tkg1.push_back(DLgvel[j][3]/ddp);
	  tkg2.push_back(tkg1);
	}//for i Lnpara
	kernel.push_back(ttk2);kernel.push_back(ttk2);
	kernel.push_back(tkp2);kernel.push_back(tkg2);
  }// if Lflag
  return 1; 
}//compute_kernel
//---------------------------------------------------------------
int read_kernel(paradef &para, modeldef &model,vector<vector<vector<double> > > &kernel, const char *fRker, const char *fLker,int &Rflag, int &Lflag,vector<vector<double> > PREM,int Nprem, int isoflag)
{
  //kernel: kernel[kRp[nP][nT],kRg[][],kLp[][],kLg[][]]
  // read in kernel from outside file, now it only read in Rph and Lph kernel.
  fstream mff;
  string line;
  vector<string> v;
  vector<vector<double> > tkp2,tkg2,ttk2;
  vector<double> tkp1,Rg0,Lg0,ttk1;
  int nT,nP,i,j;

  if(para.Rnpara!=para.Lnpara and isoflag!=1){printf("### in read_kernel, Rnpara!=Lnpara!!\n");exit(0);}
  compute_dispMineos(model,PREM,Nprem,Rflag,Lflag,0);

  for(i=0;i<model.data.Rdisp.ngper;i++)Rg0.push_back(0.);
  for(i=0;i<model.data.Ldisp.ngper;i++)Lg0.push_back(0.);
  ttk1.push_back(-999.);ttk2.push_back(ttk1);  

  kernel.clear();
  if(Rflag>0){
	nT=model.data.Rdisp.npper;
	nP=para.Rnpara+para.Lnpara;
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
        nP=para.Rnpara+para.Lnpara;
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

//---------------------------------------------------------------
int compute_dispKernel(modeldef &model,paradef para,modeldef refmod, paradef refpara, vector<vector<vector<double> > > kernel)
{
  //according to kernel(nP*nT) and para(1*nP), compute disp(1*nT) , Ti=sum(Pj*Kji)
  //***the & before 'model' may need to be deleted***;
  int i,j,k,nT,nP;
  vector<vector<double> > tker,tpara,trefdisp;
  vector<double> tdisp, *tdisplist[4];
  double tvel,Rdp,Ldp;
  vector<double> paradiff;
  
  for(i=0;i<para.Rnpara;i++){
	Rdp=para.Rparameter[i]-refpara.Rparameter[i];
	paradiff.push_back(Rdp);
  }
  for(i=0;i<para.Lnpara;i++){
	Ldp=para.Lparameter[i]-refpara.Lparameter[i];
	paradiff.push_back(Ldp);
  }

  tdisplist[0]=&model.data.Rdisp.pvel;
  tdisplist[1]=&model.data.Rdisp.gvel;
  tdisplist[2]=&model.data.Ldisp.pvel;
  tdisplist[3]=&model.data.Ldisp.gvel;

  trefdisp.push_back(refmod.data.Rdisp.pvel);
  trefdisp.push_back(refmod.data.Rdisp.gvel);
  trefdisp.push_back(refmod.data.Ldisp.pvel);
  trefdisp.push_back(refmod.data.Ldisp.gvel);
  //***check****
  int nTlist[4];
  nTlist[0]=model.data.Rdisp.npper;nTlist[1]=model.data.Rdisp.ngper;nTlist[2]=model.data.Ldisp.npper;nTlist[3]=model.data.Ldisp.ngper;
  //*******  
  for(i=0;i<4;i++){// 4 kernel, Rp Rg Lp Lg
    tdisp.clear();
    if(kernel[i].size()<2)continue;
    tker=kernel[i];
    nT=tker[0].size();
    nP=tker.size();
    //*** check***
    if(nT!=nTlist[i]){cout<<"### in cpt_dispKernel, wrong dimension in Kernel!\n"<<"nT_ker nT_disp:"<<nT<<" "<<nTlist[i]<<endl;exit(0);}
    if(nP!=paradiff.size()){cout<<"### in cpt_dispKernel, wrong dimension in Kernel!\n"<<"nP_ker nP_para:"<<nP<<" "<<paradiff.size()<<endl;exit(0);}
    //************
    for(j=0;j<nT;j++){
	tvel=(trefdisp[i][j]);
	for(k=0;k<nP;k++)tvel=tvel+paradiff[k]*tker[k][j];
	tdisp.push_back(tvel);
    }//for j
    *(tdisplist[i])=tdisp;
  }//for i 
  return 1;
}// compute_dispKernel
//---------------------------------------------------------

int get_misfitKernel(modeldef &model,paradef &para,modeldef refmodel,paradef refpara,vector<vector<vector<double> > > kernel,int Rflag,int Lflag,float depcri1,float depcri2, float qpcri, float qscri,float inp)
{//computing calc_disp using kernel, then get misfit between obs and calc disp
  compute_dispKernel(model,para,refmodel,refpara,kernel);
  compute_rf(model,depcri1,depcri2,qpcri,qscri);
  compute_misfit(model,Rflag,Lflag,inp);
  para.L=model.data.L;para.misfit=model.data.misfit;
  return 1;
}
//------------------------------------------------------
//int get_misfitMineos(modeldef inmodel, modeldef &outmodel,paradef para,int Rflag, int Lflag,float depcri1,float depcri2,float qpcri,float qscri,vector<vector<double> > PREM,int Nprem)
int get_misfitMineos(modeldef &model,paradef &para,int Rflag, int Lflag,float depcri1,float depcri2,float qpcri,float qscri,vector<vector<double> > PREM,int Nprem,float inp)
{//computing calc_disp using Mineos based on input model. then get misfit between obs and clac disp
  //para2mod(para,inmodel,outmodel);
  updatemodelTibet(model,depcri1,depcri2,qpcri,qscri);
  printf("############vsv=%g vsh=%g vp=%g vpvs=%g\n",model.laym0.vsv[0],model.laym0.vsh[0],model.laym0.vp[0],model.groups[0].vpvs);//---test---
  compute_dispMineos(model,PREM,Nprem,Rflag,Lflag,0);
  compute_rf(model,depcri1,depcri2,qpcri,qscri);
  compute_misfit(model,Rflag,Lflag,inp);
  para.L=model.data.L;para.misfit=model.data.misfit;
  return 1;
}
//------------------------------------------------------
