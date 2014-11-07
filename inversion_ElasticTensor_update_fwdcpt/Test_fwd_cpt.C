// this code contains 3ways to compute the RAdisp
// 1) Mineos, 2) using Lkernel, 3) using Vkernel
// hihihi

int get_RAmodpara(paradef inpara, modeldef inmodel, paradef &RApara, modeldef &RAmodel,int flagupdaterho ){
  //should be Bspline para and model 

  RApara=inpara;
  Lovepara2Vpara(RApara,inmodel);//from the RA part of the Love para, get vsv~vsh, eta, &theta=phi=0
  para2mod(RApara,inmodel,RAmodel);//para.parameter[]-> m.Xvalue[] -> para.parameter[]
  updatemodel(RAmodel,flagupdaterho);// m.g[].Xvalue[] -> m.laym0.X[]

  return 1;
}
//-----------------------------------
int write_RAdisp(char *foutR, char *foutL, modeldef RAmodel){
  FILE *Routdisp,*Loutdisp;
  if((Routdisp=fopen(foutR,"w"))==NULL){
        printf("cannot open file to write %s\n",foutR);
        exit(0);
  }
  if((Loutdisp=fopen(foutL,"w"))==NULL){
        printf("cannot open file to write %s\n",foutL);
        exit(0);
  }
  write_ASCdisp_single(RAmodel.data.Rdisp,Routdisp);
  write_ASCdisp_single(RAmodel.data.Ldisp,Loutdisp);
  return 1;
}
//-----------------------------------
int computeRAdisp_Mineos(modeldef RAmodel,vector<vector<double> >  PREM, int Nprem){
  //input Bspline para and model
  // compute the RA disp curve with Mineos
  char foutnmR[100],foutnmL[100];
 
  sprintf(foutnmR,"RAdisp_R_Mineos.txt");
  sprintf(foutnmL,"RAdisp_L_Mineos.txt");
  compute_dispMineos(RAmodel,PREM,Nprem,1,1,0);
  write_RAdisp(foutnmR,foutnmL,RAmodel);  

  return 1;
}

//-----------------------------------
int computeRAdisp_Lkernel(paradef RApara, modeldef RAmodel,paradef refpara, modeldef refmodel, vector<vector<vector<double> > > Vkernel, vector<vector<vector<double> > > Lkernel){
  //input Bspline para and model
  // compute the RA disp curve with Love kernel
  char foutnmR[100],foutnmL[100];
  int flagupdaterho=0;

  Vpara2Lovepara(RApara,RAmodel,flagupdaterho);
  for (int i=0;i<RAmodel.ngroup;i++){
	RAmodel.groups[i].flagcpttype=4;//compute RAdisp with Lkernel, AZdisp with Lkernel
  }
  compute_RAdisp(RAmodel,RApara,refmodel,refpara,Vkernel,Lkernel,1,1);
  sprintf(foutnmR,"RAdisp_R_Lkernel.txt");
  sprintf(foutnmL,"RAdisp_L_Lkernel.txt");
  write_RAdisp(foutnmR,foutnmL,RAmodel);  

  return 1;
}
//-----------------------------------
int computeRAdisp_Vkernel(paradef RApara, modeldef RAmodel,paradef refpara, modeldef refmodel, vector<vector<vector<double> > > Vkernel, vector<vector<vector<double> > > Lkernel){
  //input Bspline para and model
  // compute the RA disp curve with Vkernel
  char foutnmR[100],foutnmL[100];

  for (int i=0;i<RAmodel.ngroup;i++){
	RAmodel.groups[i].flagcpttype=2;//compute RAdisp with Vkernel, AZdisp with Lkernel
  }
  compute_RAdisp(RAmodel,RApara,refmodel,refpara,Vkernel,Lkernel,1,1);

  sprintf(foutnmR,"RAdisp_R_Vkernel.txt");
  sprintf(foutnmL,"RAdisp_L_Vkernel.txt");
  write_RAdisp(foutnmR,foutnmL,RAmodel);  

  return 1;
}
//-----------------------------------


