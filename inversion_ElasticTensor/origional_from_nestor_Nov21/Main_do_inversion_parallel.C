// pay attention to the included CALMineos, the "computeDisp4IsoVs" version is used to compute disp correlated with Isotropic Vs, not anisotropic model in which Vsv&Vsh are different
//


#include<iostream>
#include<algorithm>
#include<vector>
#include<cmath>
#include<fstream>
#include <Eigen/Core>
#include<string>

using namespace std;
using namespace Eigen;

#include"./string_split.C"
#include"./generate_Bs.C"
#include"./gen_random.C"
#include"./INITstructure.h"
#include"CALpara_isolay.C"
#include"./CALgroup_smooth.C"
#include"CALmodel_LVZ_ET.C"
//#include"CALforward_Mineos_readK.C_changeING"
//#include"CALforward_Mineos_readK.C_bak"
#include"CALforward_Mineos_readK_parallel.C_bak"
#include "./ASC_rw.C"
#include "./BIN_rw.C"
#include "CALinv_isolay_rf_parallel.C"

#define _USE_MATH_DEFINES

/*
#include "./BIN_rw.C"
#include "./CALinv_isolay.C"
#define PI 3.14159265
*/



int main(int argc, char *argv[])
{
int npoint,Rsurflag,Lsurflag,Rmonoc,Lmonoc,PosAnic,iitercri1,iitercri2,ijumpcri1,ijumpcri2;
int i,j,k,isoflag,Nprem,k1,k2;
int AziampRsurflag,AziphiRsurflag,AziampLsurflag,AziphiLsurflag;
double bestmisfit, misfit, L,misfitcri;
float depcri1,depcri2,qpcri,qscri,lon,lat;
char inponm[100],Rgpindir[100],Rphindir[100],Lphindir[100],Lgpindir[100],kernelnmR[100],kernelnmL[100];
vector<string> AziampRdispnm,AziphiRdispnm,AziampLdispnm,AziphiLdispnm;
vector<string> Rdispnm,Ldispnm;
char nodeid[5],str[150],modnm[100],Lparanm[100],Rparanm[100],PREMnm[100],fparanm[100];
time_t start;
modeldef model0,model1,ttmodel,modelref,modelnew;
modeldef modeltemp,modelavg1;
paradef para0,para1,para2,pararef,paranew,tpara,paraavg1;
FILE *inpo,*fkernel,*fmisfit;
vector<vector<double> >  PREM;
vector<vector<vector<double> > > Vkernel,Lkernel;
vector<int> Lvmono,Lvgrad,Rvmono,Rvgrad,Vposani,idlst;
vector<paradef> paralst;
vector<double> parastd,LoveRAparastd;
vector<vector<double> > LoveAZparastd;
char modnm1[100],fRdispnm[100],fLdispnm[100],modnm2[100],fAZRdispnm[100],fAZLdispnm[100];
char outinitnm[100],lay[20],dirlay[50],tmpstr[500],fbinnm1[200],fbinnm2[200];
float inpamp,inpphi;

if(argc!=9){
printf("Usage: xx 1]input_point_file 2]output_dir_name 3]input_vsv_dir_name 4]Rphindir 5]Rgpindir 6]Lphindir 7]Lgindir\n");
printf("1] node_id node_lon node_lat\n3]dir/initmod/vsv_node_lon_lat.mod\n");
exit(0);
}
sprintf(dirlay,argv[2]);

  //----------------PARAMETERS-----------------------------------------
  isoflag=1; //isoflag==1: Vsv=Vsh, isoflag==0: Vsv!=Vsh
  Rsurflag=1; //surflag==1: open phase only. surfalg ==3 open phase and group, surflag==2: open group only
  Lsurflag=1;
  AziampRsurflag=1;
  AziphiRsurflag=1;
  AziampLsurflag=0;
  AziphiLsurflag=0;
  inpamp=0.25; //the weight of the azi_aniso disp curve, amp part (0~1)
  inpphi=0.25; //the weight of the azi_aniso disp curve, ang part (0-1)
  //the weight of iso dispersion curve is 1-inpamp-inpphi  
  iitercri1=100000;//12000 (mod1, 1cstlay)
  iitercri2=15000;
  ijumpcri1=8;
  ijumpcri2=5;
  depcri1=20.0;
  depcri2=80.0;
  qpcri=900.;//900.;
  qscri=250.;
  Rmonoc=1;
  Lmonoc=1;
  PosAnic=1;
  //Rvmono.push_back(0);
  Rvmono.push_back(1);
  //Rvmono.push_back(2);
  Lvmono.push_back(1);
  //Lvmono.push_back(2);
  //Rvgrad.push_back(0);
  Rvgrad.push_back(1);//Rvgrad.push_back(2);
  //Lvgrad.push_back(0);
  Lvgrad.push_back(1);
  Vposani.push_back(1);
  Vposani.push_back(2);
  //Vposani.push_back(1);Vposani.push_back(2);
  k1=0;k2=1;
  //----------------------------------------------------------------------

  //sprintf(PREMnm,"/home/jiayi/progs/jy/Mineos/Mineos-Linux64-1_0_2/DEMO/models/prem_noocean.txt");
  sprintf(PREMnm,"/home/jiayi/progs/jy/Mineos/Mineos-Linux64-1_0_2/DEMO/models/ak135_iso_nowater.txt");
  sprintf(inponm,argv[1]);
  sprintf(Rphindir,argv[4]);
  sprintf(Rgpindir,argv[5]);
  sprintf(Lphindir,argv[6]);
  sprintf(Lgpindir,argv[7]);
  sprintf(fparanm,argv[8]);

  readPREM(PREMnm,PREM,Nprem);

  sprintf(tmpstr,"if [ ! -d %s ]; then mkdir %s; fi",dirlay,dirlay);
  system(tmpstr);
  sprintf(tmpstr,"if [ ! -d %s/initmod ]; then mkdir %s/initmod; fi",dirlay,dirlay);
  system(tmpstr);
sprintf(tmpstr,"if [ ! -d %s/binmod ]; then mkdir %s/binmod; fi",dirlay,dirlay);
  system(tmpstr);

  //---------------------------------------------------------
  if((inpo=fopen(inponm,"r"))==NULL){
	printf("Cannot open points file %s!\n",inponm);exit(0);
  }
  npoint=0;

  while(1){
    if(fscanf(inpo,"%s %f %f",&nodeid[0],&lon,&lat)==EOF)
	break;
    npoint++;
    printf("Begin to work on point %d: id=%s lon=%f lat=%f\n",npoint,nodeid,lon,lat);

    start=time(0);
    bestmisfit=1e10;

    Rdispnm.clear();
    sprintf(str,"%s/disp_%.1f_%.1f.txt",Rphindir,lon,lat);
    Rdispnm.push_back(str);

    Ldispnm.clear();
    sprintf(str,"%s/disp_%.1f_%.1f.txt",Lphindir,lon,lat);
    Ldispnm.push_back(str);

    AziampRdispnm.clear();
    sprintf(str,"%s/aziamp_%.1f_%.1f.txt",Rphindir,lon,lat);
    AziampRdispnm.push_back(str);

    AziphiRdispnm.clear();
    sprintf(str,"%s/aziphi_%.1f_%.1f.txt",Rphindir,lon,lat);
    AziphiRdispnm.push_back(str);

    AziampLdispnm.clear();
    sprintf(str,"%s/aziamp_%.1f_%.1f.txt",Lphindir,lon,lat);
    AziampLdispnm.push_back(str);

    AziphiLdispnm.clear();
    sprintf(str,"%s/aziphi_%.1f_%.1f.txt",Lphindir,lon,lat);
    AziphiLdispnm.push_back(str);



    //-----the final output model, that will serve as input model for the next step 
    sprintf(outinitnm,"%s/initmod/ani_%s_%.1f_%.1f.mod",dirlay,nodeid,lon,lat);
    //-----the outpur binary files, all the accepted models during inversion
    sprintf(fbinnm1,"%s/binmod/ani_%s_%.1f_%.1f.bin",dirlay,nodeid,lon,lat);
    sprintf(fbinnm2,fbinnm1);//overwrite the first binary output


    //-----starting model-----
    sprintf(modnm,argv[3]);
    //---------------------------------------------------------
    initmodel(model0);
    initmodel(modelref);

    initpara(para0);
    initpara(para1);
    initpara(pararef);
  
    readdisp(model0,Rdispnm,Ldispnm,AziampRdispnm,AziphiRdispnm,AziampLdispnm,AziphiLdispnm,Rsurflag,Lsurflag,AziampRsurflag,AziphiRsurflag,AziampLsurflag,AziphiLsurflag); 
    //==check===
    //printf("test-- the number of AZ data: AZRamp.npper=%d AZRamp.pvel.size()=%d AZRamp.pvelo.size=%d\n",model0.data.AziampRdisp.npper,model0.data.AziampRdisp.pvel.size(),model0.data.AziampRdisp.pvelo.size());
    printf("the Rayleigh wave data:\n");
    for(i=0;i<model0.data.Rdisp.npper;i++){
    	printf("T=%g vel=%g unc=%g\n",model0.data.Rdisp.pper[i],model0.data.Rdisp.pvelo[i],model0.data.Rdisp.unpvelo[i]);
    }
    printf("the Rayleigh Azi wave data:\n");
    for(i=0;i<model0.data.AziampRdisp.npper;i++){
    	printf("T=%g vel=%g unc=%g\n",model0.data.AziampRdisp.pper[i],model0.data.AziampRdisp.pvelo[i],model0.data.AziampRdisp.unpvelo[i]);
    }
    for(i=0;i<model0.data.AziphiRdisp.npper;i++){
    	printf("T=%g vel=%g unc=%g\n",model0.data.AziphiRdisp.pper[i],model0.data.AziphiRdisp.pvelo[i],model0.data.AziphiRdisp.unpvelo[i]);
    }
    //
    
    readmodAniso(model0,modnm);// both m.g.LV/Rv are filled regardless of flags. (readin iso model)
    printf("#########group[0].vsvvalue[0]=%g vshvalue[0]=%g\n",model0.groups[0].vsvvalue[0],model0.groups[0].vshvalue[0]);
    //===check===
    for(i=0;i<model0.ngroup;i++){
    	printf("group %d\n",i);
	for(j=0;j<model0.groups[i].np;j++){
		printf("  np %d\n vsv=%g vsh=%g vpv=%g vph=%g eta=%g theta=%g phi=%g vpvs=%g\n",j,model0.groups[i].vsvvalue[j],model0.groups[i].vshvalue[j],model0.groups[i].vpvvalue[j],model0.groups[i].vphvalue[j],model0.groups[i].etavalue[j],model0.groups[i].thetavalue[j],model0.groups[i].phivalue[j],model0.groups[i].vpvs);
	}
    
    }
    //
    readpara(para0,fparanm);


    mod2para(model0,para0,para1);////fill both para.R/Lpara0 (they could be inequal if Rf*Lf>0, they are equal if Rf*Lf=0)

    checkParaModel(para1,model0);
    //===check===
    printf("check readpara & mod2para----\n");
    printf("inpara.flag=%d\n",para0.flag);
    for(i=0;i<para1.npara;i++){
	printf("parameter %d\n",i);
	printf(" value=%g  dv=%g ng=%g nv=%g pflag=%g LVflag=%g RWflag=%g LWflag=%g AZflag=%g\n",para1.parameter[i],para1.para0[i][2],para1.para0[i][4],para1.para0[i][5],para1.para0[i][6],para1.para0[i][7],para1.para0[i][8],para1.para0[i][9],para1.para0[i][10]);
    }
    //
    int flagupdaterho=0;
    Vpara2Lovepara(para1,model0,flagupdaterho);//there is para2mod inside the code

    int flagreadVkernel=1,flagreadLkernel=0;

    if(model0.flag==0){updatemodel(model0,flagupdaterho);}
    compute_dispMineos(model0,PREM,Nprem,Rsurflag,Lsurflag,0);

    //exit(0);
    printf("@@@ check, main, from initial model\n");
    for(j=0;j<model0.data.Rdisp.npper;j++)printf("  @@@ check, T=%g,vin=%g vold=%g\n",model0.data.Rdisp.pper[j],model0.data.Rdisp.pvelo[j],model0.data.Rdisp.pvel[j]);

    printf("@@@ check, do V kernel ====\n");
    //---obtain V kernel ---
    sprintf(kernelnmR,"%s/VkernelRp1ani_%s_%.1f_%.1f.txt",dirlay,nodeid,lon,lat);
    sprintf(kernelnmL,"%s/VkernelLp1ani_%s_%.1f_%.1f.txt",dirlay,nodeid,lon,lat);

    if(flagreadVkernel==1){
    	if((read_kernel(para1,model0,Vkernel,kernelnmR,kernelnmL,Rsurflag,Lsurflag,PREM,Nprem))==0){
		 printf ("#####!! read_kernel failed\n");
    	 	 sprintf(str,"echo point %d: id=%s lon=%f lat=%f >> point_rdKernel_failed_Ani.txt",npoint,nodeid,lon,lat);
         	 system(str);
		 continue;
     	} // if readkernel()==0    
    }//if flagreadkernel==1
    else{
    	compute_Vkernel(para1,model0,Vkernel,PREM,Nprem,Rsurflag,Lsurflag,flagupdaterho);
	cout<<"check finish compute_Vkernel\n";
	write_kernel(Vkernel,model0,para1,kernelnmR,kernelnmL,Rsurflag,Lsurflag);
	cout<<"check finish write_Vkernel\n";
    }//else flagreadkernel==1
    //====


    printf("@@@ check, do L kernel ====\n");
    //---obtain Love kernel ---
    sprintf(kernelnmR,"%s/LkernelRp1ani_%s_%.1f_%.1f.txt",dirlay,nodeid,lon,lat);
    sprintf(kernelnmL,"%s/LkernelLp1ani_%s_%.1f_%.1f.txt",dirlay,nodeid,lon,lat);
    if(flagreadLkernel==1){
    	if((read_kernel(para1,model0,Lkernel,kernelnmR,kernelnmL,Rsurflag,Lsurflag,PREM,Nprem))==0){
		 printf ("#####!! read_kernel failed\n");
    	 	 sprintf(str,"echo point %d: id=%s lon=%f lat=%f >> point_rdKernel_failed_Ani.txt",npoint,nodeid,lon,lat);
         	 system(str);
		 continue;
     	} // if readkernel()==0    
    	
    }
    else{
    	Vkernel2Lkernel(para1,model0,Vkernel,Lkernel,flagupdaterho);
    	write_kernel(Lkernel,model0,para1,kernelnmR,kernelnmL,Rsurflag,Lsurflag);
    }

/*    
 int d = Math.abs(a - b) % 360;
     int r = d > 180 ? 360 - d : d;
*/
    printf("test-- cpt misfit\n");
    //somethihng to do, 1st, in generating the ang disp curve, make it smooth
    //something to do, in computing the misfit for the angle phi, need to take the period of the angle into account. modify the compute_misfitDISP in CALmodel_LVZ_ET.C
    compute_misfitDISP(model0,Rsurflag,Lsurflag,AziampRsurflag,AziampLsurflag,AziphiRsurflag,AziphiLsurflag,inpamp,inpphi);
    //is the AZ disp filled? so are they taken into account in the misfit cpt? A: yes, both amp and angle are 0

    modelref=model0;
    pararef=para1;

    float theta;
    printf("test-- do inv\n");
    if((do_inv(2,-2,paralst,pararef,modelref,Rvmono,Lvmono,Rvgrad,Lvgrad,PREM,Vkernel,Lkernel,k1,k2,start,isoflag,Rsurflag,Lsurflag,AziampRsurflag,AziampLsurflag,AziphiRsurflag,AziphiLsurflag,Nprem,Rmonoc,Lmonoc,PosAnic,Vposani,iitercri1,ijumpcri1,fbinnm1,inpamp,inpphi,flagupdaterho))==0){
    	sprintf(str,"echo %s %.1f %.1f >> point_do_inv_failed.txt",nodeid,lon,lat);
	system(str);
    	continue;
    }
    //---then get the avg model for the second kernel computation---
    idlst.clear();
    if((para_avg(paralst,tpara,parastd,LoveRAparastd,LoveAZparastd,idlst))==0){cout<<"#### in para_avg,incorrect paralst.size()\n";exit(0);}
    para2mod(tpara,model0,modelavg1);//will the averaged vpv~eta be different from that computed from avearge_vsv,vsh??
    mod2para(modelavg1,tpara,paraavg1);
    updatemodel(modelavg1,flagupdaterho);


    //------------------write out model_average--------------------------------------------------



    sprintf(fRdispnm,"%s/Rdisp_%.0f_%s_%.1f_%.1f.txt",theta,dirlay,nodeid,lon,lat);
    sprintf(fLdispnm,"%s/Ldisp_%.0f_%s_%.1f_%.1f.txt",theta,dirlay,nodeid,lon,lat);
    sprintf(fAZRdispnm,"%s/AZRdisp_%.0f_%s_%.1f_%.1f.txt",theta,dirlay,nodeid,lon,lat);
    sprintf(fAZLdispnm,"%s/AZLdisp_%.0f_%s_%.1f_%.1f.txt",theta,dirlay,nodeid,lon,lat);

    get_misfitKernel(modelavg1,paraavg1,modelref,pararef,Vkernel,Lkernel,Rsurflag,Lsurflag,AziampRsurflag,AziampLsurflag,AziphiRsurflag,AziphiLsurflag,inpamp,inpphi,flagupdaterho);
    write_ASC(modelavg1,paraavg1,fRdispnm,fLdispnm,fAZRdispnm,fAZLdispnm,Rsurflag,Lsurflag,AziampRsurflag, AziampLsurflag,AziphiRsurflag, AziphiLsurflag);

    write_initmodAniso(outinitnm,modelavg1);

    sprintf(fRdispnm,"%s/RdispM_%.0f_%s_%.1f_%.1f.txt",theta,dirlay,nodeid,lon,lat);
    sprintf(fLdispnm,"%s/LdispM_%.0f_%s_%.1f_%.1f.txt",theta,dirlay,nodeid,lon,lat);
    printf("misfit: %8g\n Rmisfit: iso=%8g AZamp=%8g AZphi=%8g\nLmisfit: iso=%8g AZamp=%8g AZphi=%8g\n",modelavg1.data.misfit,modelavg1.data.Rdisp.misfit,modelavg1.data.AziampRdisp.misfit,modelavg1.data.AziphiRdisp.misfit,modelavg1.data.Ldisp.misfit,modelavg1.data.AziampLdisp.misfit,modelavg1.data.AziphiLdisp.misfit);
 
    //after applying the get_misfitMineosRA, since there is Lpara2Vpara inside this function, the Vpra becomes effective TI parameters (theta=0 phi=0, vsv~eta all changed)
    get_misfitMineosRA(modelavg1,paraavg1,PREM,Nprem,Rsurflag,Lsurflag,flagupdaterho);
    write_ASC(modelavg1,paraavg1,fRdispnm,fLdispnm,fAZRdispnm,fAZLdispnm,Rsurflag,Lsurflag,0,0,0,0);
    printf("\n---\nmisfit(only accounts RA from Mineos): %8g\n Rmisfit: iso=%8g AZamp=%8g AZphi=%8g\nLmisfit: iso=%8g AZamp=%8g AZphi=%8g\n",modelavg1.data.misfit,modelavg1.data.Rdisp.misfit,modelavg1.data.AziampRdisp.misfit,modelavg1.data.AziphiRdisp.misfit,modelavg1.data.Ldisp.misfit,modelavg1.data.AziampLdisp.misfit,modelavg1.data.AziphiLdisp.misfit);


/*
    //---then compute kernel2
    Vkernel.clear();Lkernel.clear();
    compute_Vkernel(paraavg1,modelavg1,Vkernel,PREM,Nprem,Rsurflag,Lsurflag,flagupdaterho);
    sprintf(kernelnmR,"%s/VkernelRp2ani_%s_%.1f_%.1f.txt",dirlay,nodeid,lon,lat);
    sprintf(kernelnmL,"%s/VkernelLp2ani_%s_%.1f_%.1f.txt",dirlay,nodeid,lon,lat);
    write_kernel(Vkernel,modelavg1,paraavg1,kernelnmR,kernelnmL,Rsurflag,Lsurflag);

    Vkernel2Lkernel(paraavg1,modelavg1,Vkernel,Lkernel,flagupdaterho);
    sprintf(kernelnmR,"%s/LkernelRp2ani_%s_%.1f_%.1f.txt",dirlay,nodeid,lon,lat);
    sprintf(kernelnmL,"%s/LkernelLp2ani_%s_%.1f_%.1f.txt",dirlay,nodeid,lon,lat);
    write_kernel(Lkernel,modelavg1,paraavg1,kernelnmR,kernelnmL,Rsurflag,Lsurflag);

    write_initmodAniso(outinitnm,modelavg1);

@@    //--do second inversion---
    if((do_inv(2,-2,))==0){
    	sprintf(str,"echo %s %.1f %.1f >> point_do_inv_failed.txt",nodeid,lon,lat);
	system(str);
    	continue;  
    }
    vector<paradef>(paralst).swap(paralst);//shrink size

    idlst.clear();
    if((para_avg())==0){}
    para2mod();
    mod2para();
    write_initmodAniso();
    
    get_misfitRAMineos();
    write_ASC();

    get_misfitKernel();
    write_ASC();
*/

    /*
    perturb_para();
    //para2mod();
    Vpara2Lovepara();

    for(i=0;i<para1.Rnpara;i++)printf("i=%d paraR=%g paraL=%g\n",i,para1.Rparameter[i],para1.Lparameter[i]);//--test--
    //updatemodelTibet(model0,depcri1,depcri2,qpcri,qscri);
    //compute_dispMineos(model0,PREM,Nprem,1,1);
    model1=model0;
    //para2mod(para1,model0,model1);
    get_misfitMineos(model1,para1,1,1,depcri1,depcri2,qpcri,qscri,PREM,Nprem);//Lp2Vp is inside it
    //------------------write out model_average--------------------------------------------------

    sprintf(modnm1,"%s/AnimodM_%s_%.1f_%.1f.txt",dirlay,nodeid,lon,lat);
    sprintf(modnm2,"%s/AnimodK_%s_%.1f_%.1f.txt",dirlay,nodeid,lon,lat);
    sprintf(Rdispnm1,"%s/AniRdispM_%s_%.1f_%.1f.txt",dirlay,nodeid,lon,lat);
    sprintf(Rdispnm2,"%s/AniRdispK_%s_%.1f_%.1f.txt",dirlay,nodeid,lon,lat);
    sprintf(Ldispnm1,"%s/AniLdispM_%s_%.1f_%.1f.txt",dirlay,nodeid,lon,lat);
    sprintf(Ldispnm2,"%s/AniLdispK_%s_%.1f_%.1f.txt",dirlay,nodeid,lon,lat);

    write_ASC(model1,para1,modnm1,Rdispnm1,Ldispnm1,1,1);


    vector<vector<vector<double> > >().swap(kernel1);//free
    vector<vector<vector<double> > >().swap(kernel2);//free
    vector<paradef>().swap(paralst);//free
    vector<double>().swap(Rparastd);//free
    vector<double>().swap(Lparastd);//free
  }//while 1
  vector<vector<double> >().swap(PREM);//free
  fclose(fmisfit); 
*/
  }//while1
  return 1;
}
//main
