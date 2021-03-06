     int write_model(modeldef &model,char *name, char *outdir)
        {
          char namerf[300],namedispp[300],namedispg[300],namemod[300];
          FILE *ff;
          int i;
          double dep=0.;
          sprintf(namerf,"%s/%s.rf",outdir,name);
          sprintf(namedispp,"%s/%s.p.disp",outdir,name);
          sprintf(namedispg,"%s/%s.g.disp",outdir,name);
          sprintf(namemod,"%s/%s.mod",outdir,name);
          ff=fopen(namerf,"w");
          for (i=0;i<model.data.rf.rfn.size();i++)
                {
                  if(model.data.rf.tn.size()>0 and model.data.rf.tn[i]<=10)
                        fprintf(ff,"%g %g %g %g %g\n",model.data.rf.tn[i],model.data.rf.rfn[i],model.data.rf.to[i],model.data.rf.rfo[i],model.data.rf.unrfo[i]);
                }//for i
         fclose(ff);

         if(model.data.disp.npper>0)
         {
         ff=fopen(namedispp,"w");
         for(i=0;i<model.data.disp.npper;i++)
                {
                  if(model.data.disp.pvel.size()>0)
                        fprintf(ff,"%g %g %g %g\n",model.data.disp.pper[i],model.data.disp.pvel[i],model.data.disp.pvelo[i],model.data.disp.unpvelo[i]);
                }//fori
          fclose(ff);
          }//if


         if(model.data.disp.ngper>0)
         {
         ff=fopen(namedispg,"w");
         for(i=0;i<model.data.disp.ngper;i++)
                {
                  if(model.data.disp.gvel.size()>0)
                        fprintf(ff,"%g %g %g %g\n",model.data.disp.gper[i],model.data.disp.gvel[i],model.data.disp.gvelo[i],model.data.disp.ungvelo[i]);
                }//fori
          fclose(ff);
          }//if

          ff=fopen(namemod,"w");
          for(i=0;i<model.laym0.nlayer;i++) //write layered model in stair-step
                {
                  fprintf(ff,"%g %g\n",dep,model.laym0.vs[i]);
                  dep=dep+model.laym0.thick[i];
                  fprintf(ff,"%g %g\n",dep,model.laym0.vs[i]);
                }
          fclose(ff);
          return 1;
        }//write model


//-----------------------------------------------------	 
	int gen_newpara(paradef inpara, paradef &outpara, int pflag )
	{//change the p.parameter[]
	 /*//#### the default relationship between RAvp&RAvs, eta&RAvs is set here within this function ######	
	 float fRAvp,feta,ceta;
	 fRAvp=0.5; //RAvp=fRAvp*RAvs
	 ceta=1.0; feta=-4.2; //eta=ceta-feta*RAvs
	 //#########################################
	 */
	 outpara=inpara;
	 gen_new_para_single(inpara.space1,outpara.parameter,inpara.npara,pflag);

	  return 1;
	}




//---------------------------------------
int rotate_ET(Matrix<double,6,6> ETin, double theta, double phi, Matrix<double,6,6> &ETout){
  Matrix3f a1m,a2m,a3m,am;
  Matrix<double, 6,6> Mm;
  double ox,oy,oz;
  double axx,axy,axz,ayx,ayy,ayz,azx,azy,azz;
  //--rotating angles---
  ox=theta*M_PI/180.;
  oy=0.;
  oz=phi*M_PI/180.;
  //--coordinate transformation matrix ---
  
  a3m<< cos(oz),sin(oz),0.,
	-sin(oz),cos(oz),0.,
	0.,0.,1.;
  cout<<a3m;
  a2m<< cos(oy),0.,-sin(oy),
	0.,1.,0.,
	sin(oy),0.,cos(oy);
  a1m<< 1.,0.,0.,
	0.,cos(ox),sin(ox),	
	0.,-sin(ox),cos(ox);
  //
  am=a3m*a2m*a1m;
  axx=am(0,0);axy=am(0,1);axz=am(0,2);
  ayx=am(1,0);ayy=am(1,1);ayz=am(1,2);
  azx=am(2,0);azy=am(2,1);azz=am(2,2);
 
  Mm<< axx*axx,       axy*axy,        axz*axz,        2*axy*axz,         2*axz*axx,              2*axx*axy,
        ayx*ayx,       ayy*ayy,        ayz*ayz,        2*ayy*ayz,         2*ayz*ayx,              2*ayx*ayy,
        azx*azx,       azy*azy,        azz*azz,        2*azy*azz,         2*azz*azx,              2*azx*azy,
        ayx*azx,       ayy*azy,        ayz*azz,        ayy*azz+ayz*azy,   ayx*azz+ayz*azx,        ayy*azx+ayx*azy,
        azx*axx,       azy*axy,        azz*axz,        axy*azz+axz*azy,   axz*azx+axx*azz,        axx*azy+axy*azx,
        axx*ayx,       axy*ayy,        axz*ayz,        axy*ayz+axz*ayy,   axz*ayx+axx*ayz,        axx*ayy+axy*ayx;
  ETout=(Mm*ETin)*Mm.transpose();//dot(dot(Mm,ETin),Mm.T);


  // 
  /*
  Eigen::FullPivLU<Matrix2f> lua1(a1);
  Eigen::PartialPivLU<Matrix2f> lua1_2(a1);
  cout<<"a1 inverse LU=\n"<<lua1.inverse()<<endl;
  cout<<"a1 inverse partial LU=\n"<<lua1_2.inverse()<<endl;
  */

  return 1;
}

int rotate_ET(double[6][6] ETin,double theta,double phi,double[6][6] &ETout){
	double ox,oy,oz;
	double[3][3] a1m,a2m,a3m,am,M1m,M2m,M3m,M4m;
	double Mm[6][6]; //,ETout[6][6];

	//--rotating angles---
	ox=theta*M_PI/180.;
	oz=phi*M_PI/180.;
	//--coordinate transformation matrix ---
	a3m={
		{cos(oz),sin(oz),0},
		{-sin(0z),cos(oz),0},
		{0,0,1}
	};
	a2m={
		{cos(oy),0,-sin(oy)},
		{0,1,0},
		{sin(oy),0,cos(oy)}
	};
	a1m={
		{1,0,0},
		{0,cos(o1),sin(o1)},
		{0,-sin(o1),cos(o1)}
	};
	/*
	!am=a3m*a2m*a1m;
	//--define coefficient transformation matrix --
	axx=am[0,0];axy=am[0,1];axz=am[0,2];
	ayx=am[1,0];ayy=am[1,1];ayz=am[1,2];
	azx=am[2,0];azy=am[2,1];azz=am[2,2];

	Mm={
		{axx*axx,	axy*axy,	axz*axz,	2*axy*axz,         2*axz*axx,              2*axx*axy},
		{ayx*ayx,	ayy*ayy,	ayz*ayz,	2*ayy*ayz,         2*ayz*ayx,              2*ayx*ayy},
		{azx*azx,	azy*azy,	azz*azz,	2*azy*azz,         2*azz*azx,              2*azx*azy},
		{ayx*azx,       ayy*azy,    	ayz*azz,	ayy*azz+ayz*azy,   ayx*azz+ayz*azx,        ayy*azx+ayx*azy},
		{azx*axx,       azy*axy,    	azz*axz,	axy*azz+axz*azy,   axz*azx+axx*azz,        axx*azy+axy*azx},
		{axx*ayx,      	axy*ayy,   	axz*ayz,	axy*ayz+axz*ayy,   axz*ayx+axx*ayz,        axx*ayy+axy*ayx}
	};
	!ETout=dot(dot(Mm,ETin),Mm.T);
	*/
      return 1;	
} // rotate_ET









//---------------------------------------

//---------------------------------------

//---------------------------------------
int Vkernel2Lovekernel(Vkernel,para,group, int inng, &LoveRAkernel,&LoveAZkernel,flagupdaterho){
  // Vkernel: kernel[kRp[nP][nT],kRg[][],kLp[][],kLg[][]]
  // the order of P is the same as that in the para
  // this Love kernel should match the Loveparameter (only for one group of parameter, but could contain many layers) from Vpara2Lovepara; so, only
/*
 the size of kernel[0-3] that is related to para that belongs to group_ng should be np*6 (vsv~eta+h)
 find them out, and, for each Rp~Lg, sort them [[layer1:Kvsv[nT],Kvsh[nT],...,Keta[nT]]  [layer2]  [layer3]...] //OR maybe do this step in the checkParaModel step, requiring the ordering of para input.

 then, compute the corresponding Love kernel, KA[nT]~KF[nT],and Kh[nT]
 then, order then in terms of the LoveRA and LoveAZ parameters
 in fact, some of the kernel is almost zero, and no need to compute at all (e.g., dLove/dA; Love~LNGE Rayleigh~ALCFGBH)
*/

  vector<double> tLkernel; //kernel for one Love parameter, has size of nT

  int i;
  for(i=0;i<4;i++){
	tk.pusb_back([][])
  }
  for(i=0;i<para.npara;i++){
    	flag=para.para0[i][0];//flag, explainations are in mod2para;
    	ng=para.para0[i][4];
    	nv=para.para0[i][5];
    	pflag=para.para0[i][6];

   	if(ng!=inng)continue;
   	if(flag!=0)continue;
 	if(pflag>5)continue;


    	vel=para.parameter[i];
    	//----check----
	if(pflag==1){
	  printf("@@@ check, Vkernel2Lovekernel vel_from_para=%g, vel_from_model=%g\n",para.parameter[i],group.vsvvalue[nv]);
 	}
	//---

	if(flagupdaterho==1){
	  rho=group.rhovalue[nv];
	  for(nk=0;nk<4;nk++){
	    computeLovekernel(Vkernel[nk][i],rho,vel,tLkernel,pflag);
	    tk[nk].push_back(tLkernel); //Vkernel[nk][i] has size of nT
	  }//for nk
	}//if flagupdaterho
	else{
	  for(nk=0;nk<4;nk++){
	    computeLovekernel(Vkernel[nk][i],-1.,vel,tLkernel,pflag);
	    tk[nk].push_back(tLkernel);
	  }
	}//else flagupdaterho


  }//for i

  //tk[][][] ==> [[kRp: (layer1:A[nT],C[nT],..N[nT]) (layer2 ...) ...]  [kRg ...] [kLp]  [kLg]] == LoveRAparameter 
  //LoveAZparametr=[[kRp: (layer1: BcBsEcEsGcGsHcHs)] ... ] = [[kRp:(layer1: AANNLLFF)]]
	  
  printf("@@@check, Vkernel2Lovekernel,length of Love_kRp=%d\n ",tk[0].size());
  return 1;
}//Vkernel2Lovekernel

//---------------------------------------
int readPREM(const char *PREMnm, vector<vector<double> >  &PREM,int &Nprem)
int write_modMineos(modeldef &model, const char *outname,vector<vector<double> > PREM,int Nprem,int &Nmod)
int interpolate(dispdef indisp,dispdef &outdisp)
int interpolate_model(layermoddef inlay, layermoddef &outlay, double dh)
int read_dispMineos(dispdef &indisp,const char* Moutputnm, int Nmod)
int compute_dispMineos(modeldef &model,vector<vector<double> > PREM,int Nprem, int Rsurflag,int Lsurflag)
int compute_diff(dispdef disp1, dispdef disp2, vector<vector<double> > &pveldiff,vector<vector<double> > &gveldiff)
int compute_Vkernel_single_para(paradef para, int ip,modeldef model, vector<vector<double> > PREM,int Nprem, int Rflag, int Lflag,int flagupdaterho, vector<double> &trkp1, vector<double> &trkg1, vector<double> &tlkp1, vector<double> &tlkg1)
int compute_Vkernel(paradef para,modeldef &model,vector<vector<vector<double> > > &kernel,vector<vector<double> > PREM,int Nprem, int Rflag,int Lflag, int flagupdaterho)
int computeLovekernel(vector<double> Vkernel1,double rho,double vel,double vpvs,vector<double> &Lkernel1, int flagupdaterho)
int computeLovekerneleta(vector<double> Vkernel1,double c, vector<double> &Lkernel1)
int Vkernel2Lkernel(paradef para,modeldef model,vector<vector<vector<double> > > Vkernel,vector<vector<vector<double> > > &Lkernel, int flagupdaterho)
int Vpara2ET2LoveCoeff(vector<double> Vparameter, double[6][6] &ET,double[8] &RAcoeff, double[8][2] &AZcoeff)
int getGroupVPara(groupdef group, vector<vector<double> > &Vparameter2, int flagupdaterho)
int Vpara2Lovepara(paradef para,modeldef model,int flagupdaterho)
int compute_RAdisp(modeldef &model, paradef para, modeldef refmod, paradef refpara, vector<vector<vector<double> > > Vkernel, vector<vector<vector<double> > > Lkernel)
int cs2ap(double Ac,double As, double &amp,double &phi, int phiflag, int RLflag)
int compute_AZdisp(modeldef &model, paradef para, modeldef refmod, paradef refpara, vector<vector<vector<double> > > Vkernel, vector<vector<vector<double> > > Lkernel)




