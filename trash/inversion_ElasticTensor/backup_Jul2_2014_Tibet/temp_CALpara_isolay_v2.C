/*#include<iostream>
#include<vector>
#include<fstream>
#include<algorithm>
#include"INITstructure.h"
#include"/home/jiayi/progs/jy/HEAD/head_c++/gen_random.C"
*/
using namespace std;
/*========CONTENT============
int initpara(paradef &para)
//int readpara_single(vector<vector<string>> &para0, int &npara,const char* fname,int surflag)
int readpara(paradef &para,const char* Rfname, const char* Lfname, int Rsurflag, int Lsurflag)
//int gen_newpara_single ( vector<vector<double> > space1, vector<double > &parameter, int npara,int pflag)
//int gen_newpara(paradef inpara, paradef &outpara, int pflag,int isoflag ,vector<float> &Viso)
int mod2para(modeldef &model, paradef &inpara, paradef &outpara)
int para2mod(paradef para, modeldef inmodel, modeldef &outmodel)

//============================
*/
//class paracal{
//public:
	int initpara(paradef &para)
	{
	  //para.Rnpara=0;para.Lnpara=0;
	  para.npara=0;
     	  para.L=0.;
	  para.misfit=0.;
	  para.flag=0;
	  //para.Rparameter.clear();para.Rpara0.clear();para.Rspace1.clear();
	  //para.Lparameter.clear();para.Lpara0.clear();para.Lspace1.clear();
	  para.parameter.clear();para.LoveRAparameter.clear();
	  para.para0.clear();para.LoveAZparameter.clear();para.space1.clear();
 	  /* vector<double>().swap(para.Rparameter);vector<double>().swap(para.Lparameter);
	  vector<vector<double> >().swap(para.Rspace1);vector<vector<double> >().swap(para.Lspace1);
	  vector<vector<string> >().swap(para.Rpara0);vector<vector<string> >().swap(para.Lpara0);
	*/
	  return 1;
	}
//-----------------------------------------------------	 
	int checkParaModel(paradef para, modeldef model){
	  //check if the para satisfies certain criteria: If the model.group.flagcpttype=2 or 4 (part or all of the forward cpt will be through LoveParameters); 
	  //then the type of model must be layer (i.e., flag=1); 
	  //all np*5 (vsv~eta) must be in the parameter list, in order to produce partial derivative for all the np*5 parameters, which is required for the computation of Love parameter partial derivatives
	  int i,j,N,count;
	  int ppflagid,nvid,ngid;
	  int p0,p5,p6;

	  for(i=0;i<model.ngroup;i++){
	    //if(model.groups[i].flagcpttype!=2 and model.groups[i].flagcpttype!=4)continue;
	    N=model.groups[i].np*7;
	    count=0;
	    ppflagid=-1;nvid=-1,ngid=-1;

	    if(model.groups[i].flag==3 and model.groups[i].flagcpttype!=1){
	      printf("### checkParaModel. group%d, if groups is described by gradient, then the forward cpt must be through Vkernel and Vpara!\n",i);
	      exit(0);
	    }

	    for(j=0;j<para.npara;j++){
	      if((int)para.para0[j][4]!=i)continue;
	      p0=(int)para.para0[j][0];    
	      if(p0==0){	 

		//---- check the order of the layer     
		p5=(int)para.para0[j][5];
		if(p5>nvid){ppflagid=-1;}//reach the next layer in this group
	        else if (p5<nvid){
			printf("### checkParaModel, in group %d, we require the layer number nv (teh 6th colunm) is increading!\n",i);
			exit(0);
		}//if	
		nvid=p5;

		//---- check the order within each layer
		p6=(int)para.para0[j][6];     
	        if (p6<=ppflagid){
	      	  printf("### checkParaModel, in group %d, layer %d,  we require the para's ppflag (the 7th column) is ordered in an increasing order!\n",i,nvid);
		  exit(0);
	        }//if
	        ppflagid=p6;
	      }//if p0==0
	  
	      //---- count the total number of parameters
	      if((int)para.para0[j][6]>7)continue;
	      //parameter belongs to group i, and is one of the five parameters (vsv~eta)
	      count++;
	    }//for j <npara

	    if((model.groups[i].flagcpttype-2)*(model.groups[i].flagcpttype-4)==0 and count!=N){
	      printf("#### checkParaModel, the number of parameters(vsv~eta,theta,phi) in group %d is %d, but should be %d to enable the Love_para forward computation!\n",i,count,N);
	      exit(0);
	    }//if count

	  }//for i
	  return 1;
	}  //checkParaModel;


//-----------------------------------------------------	 
	//int readpara(vector<vector<string> > &para0, int &npara,const char* fname)
	int readpara(vector<vector<double> > &para0, int &npara,const char* fname)
	{  
	  fstream mff;
	  int k,i;
	  string line;
	  vector<string> v;
	  vector<double> vd;

	  //if (surflag>0){
	  mff.open(fname);
	  if ( not mff.is_open()){cout<<"########para file"<<fname<<"does not exist!\n exit!!\n";exit (0);}
	  k=0;
	  para0.clear();
	  while(getline(mff,line))
	  {
		v.clear();
		Split(line,v," ");
		//v.push_back("0"); //add a column indicating if this parameter belongs to vsv~eta and belongs to a groups whose flagcpttype==2 or ==4 (part or all of foward cpt is through the Love parameters); initial value is '0', modified in the checkParaModel function

		vd.clear();
		for(i=0;i<v.size();i++){
		  vd.push_back(atof(v[i].c_str()));
		}//for i
		if((int)vd[0]==0 and v.size()<7){
		  printf("#### readpara, wrong column #(should be >=7,1 flag,2 dv_type,3 dv,4 sigma,5 ng,6 nv,7 para_type,(8 LVflag)) in the line %d\n",k+1);
		  exit(0);
		}
		//--make the width of para0 to be 11; 11 flags for each para, see INITstructure.h for detail
		for(i=v.size();i<11;i++){
		  vd.push_back(0.);
		}
		if((int)vd[0]==1)vd[6]=9.;//the ppflag, indicating the type of parameter
		else if((int)vd[0]==-1)vd[6]=12.;
		printf("@@@ check, readpara, width of of para0 is %d, should be 8!\n",vd.size());
		k=k+1;
		para0.push_back(vd);	
	  }
	  mff.close();
	  npara=k;
	  //}//if surflag		
	  return 1;
	}//readpara
//----------------------------------------------------
/*	int readpara(paradef &para,const char* vsfname, const char* vpfname, const char* etafname, const char* thetafname,  const char* phifname,int Rsurflag, int Lsurflag)	
	{
	  readpara_single(para.vsvpara0,para.vsvnpara,vsfname,Rsurflag);
	  readpara_single(para.vshpara0,para.vshnpara,vsfname,Lsurflag);
	  readpara_single(para.vpvpara0,para.vpvnpara,vpfname,Rsurflag);
	  readpara_single(para.vphpara0,para.vphnpara,vpfname,Rsurflag);
	  readpara_single(para.etapara0,para.etanpara,etafname,Rsurflag);
	  readpara_single(para.thetapara0,para.thetanpara,thetafname,Rsurflag);
	  readpara_single(para.phipara0,para.phinpara,phifname,Rsurflag);

	  printf("-----check readpara: vsvnpara=%d vpvnpara=%d\n",para.vsvnpara,para.vpvnpara);
	  return 1;
	}//readpara
*/
//-----------------------------------------------------	 
// randomly generating new parameters, uniform or normal distribution
/*	int gen_newpara_single (  vector<vector<double> > space1, vector<double> &parameter, int npara,int pflag)
	//vector<double> gen_newpara_single ( const vector<vector<double> > &space1, vector<double> parameter, int npara,int pflag)
	{
	  int i,flag;
	  double newv,sigma,mean,perc;

	  //cout<<"gen_newpara! pflag="<<pflag<<endl;
	  
	  if(pflag==0)//uniform
	  {
		for(i=0;i<npara;i++)
			{
			  if(space1[i][2]<0.001)continue;
			  newv=gen_random_unif01()*(space1[i][1]-space1[i][0])+space1[i][0];//gen_random_unif01 => [0,1)
			//  printf("$$$$ gen_newpara unif:%d %g %g %g\n",i,inpara.space1[i][0],inpara.space1[i][1],newv);
			  parameter[i]=newv;
			}
	  }//if		  	  
	  else if(pflag==1)//normal
	  {
		for(i=0;i<npara;i++)
		{
		  flag=0;
		  mean=parameter[i];
		  sigma=space1[i][2];
		  if(sigma<0.001)continue;
		  while(flag<1)
		  {
			newv=gen_random_normal(mean,sigma); // normal distribution
			if(newv>space1[i][1] or newv<space1[i][0])continue;
			else flag=2;
		  }//while flag
		  parameter[i]=newv;
		}//fori

	  }//else if 1
	  else if(pflag==2)//set the perturbed parameter to have the same purterbation. e.g. uniform anisotropy in the mantle
          {
                perc=gen_random_unif01();
                for(i=0;i<npara;i++)                
		{
		//cout<<"test--- i="<<i<<endl;
		  if (space1[i][2]<0.001){continue;} // skip the para don't require any purterbation
		  newv=perc*(space1[i][1]-space1[i][0])+space1[i][0];
		  //printf("f2 i=%d old %g new %g perc %g\n",i,outpara.parameter[i],newv,perc);
		  parameter[i]=newv; 
               }

          }//else if 2
	  else if(pflag==3)//set the perturbed parameter to have the same purterbation. e.g. uniform anisotropy in the mantle
	  {
		//printf("pglag=%d,npara=%d\n",pflag,inpara.npara);
		flag=0;
		for (i=0;i<npara;i++)
		{
		 if (space1[i][2]<0.001)continue;
		 mean=parameter[i];
		 sigma=space1[i][2];
		 break;
		}
		if (i==npara){
		 // printf("gen_newPara, pflag==3, no para need to be purterbed! skip!\n");
		  return 1;
		  //return parameter;
		}
                while(flag<1)
                  {
                        newv=gen_random_normal(mean,sigma); // normal distribution
                        if(newv>space1[i][1] or newv<space1[i][0])continue;
                        else flag=2;
                  }//while flag
		perc=newv/mean;
		for(i=0;i<npara;i++)
		{
		  if (space1[i][2]<0.0001)continue; // skip the para don't require any purterbation
		  newv=parameter[i]*perc;
	//	  printf("f3 i=%d old %g new %g perc %g\n",i,outpara.parameter[i],newv,perc);
		  parameter[i]=newv;
		}
	  }//else if 3
	  else{printf("### wrong pflag for gen_newpara!!\n");exit(0);}
	  return 1;
	  //return parameter;
	}//gen new para single
//-----------------------------------------------------
	int gen_newpara_single_L (  vector<vector<double> > space1,vector<double> Rparameter, vector<double> &Lparameter, int npara,int pflag,vector<float> Viso)
	//vector<double> gen_newpara_single ( const vector<vector<double> > &space1, vector<double> parameter, int npara,int pflag)
	//========pay attention to the unit of Viso value. (e.g., if Viso[i] is anisotropy, does it mean ani=(vsh-vsv)*2/(vsv+vsh) or ani*100?)
	{
	  int i,flag;
	  double newv,sigma,mean,perc,f;

	  //cout<<"gen_newpara! pflag="<<pflag<<endl;
	  
	  if(pflag==0)//uniform
	  {
		for(i=0;i<npara;i++)
			{
			  if(Viso[i]>999.){Lparameter[i]=Rparameter[i];}
			  else if (Viso[i]<-999.){
			  if(space1[i][2]<0.001)continue;
			  newv=gen_random_unif01()*(space1[i][1]-space1[i][0])+space1[i][0];//gen_random_unif01 => [0,1)
			//  printf("$$$$ gen_newpara unif:%d %g %g %g\n",i,inpara.space1[i][0],inpara.space1[i][1],newv);
			  Lparameter[i]=newv;}
			  else{
			      //constant velocity difference-----
			      //Lparameter[i]=Rparameter[i]+Viso[i];
			      //constant anisotropy instead of constant Velocity difference---------
			      Lparameter[i]=Rparameter[i]*(1+Viso[i]/200.)/(1-Viso[i]/200.);
			  }
			}
	  }//if		  	  
	  else if(pflag==1)//normal
	  {
		for(i=0;i<npara;i++)
		{
		  if(Viso[i]>999.){Lparameter[i]=Rparameter[i];}
	  	  else if(Viso[i]<-999.)
		  { flag=0;
		    mean=Lparameter[i];
		    sigma=space1[i][2];
		    if(sigma<0.001)continue;
		    while(flag<1)
		    {
			newv=gen_random_normal(mean,sigma); // normal distribution
			if(newv>space1[i][1] or newv<space1[i][0])continue;
			else flag=2;
		    }//while flag
		    Lparameter[i]=newv;}//else if
		  else	  
		  {	  //constant velocity differnce -------
			  //Lparameter[i]=Rparameter[i]+Viso[i];
			  //constant anisotropy --------
			  Lparameter[i]=Rparameter[i]*(1+Viso[i]/200.)/(1-Viso[i]/200.);
		  }//else
		}//fori

	  }//else if 1	
	 else if(pflag==2)//set the perturbed parameter to have the same purterbation. e.g. uniform anisotropy in the mantle
          {
                perc=gen_random_unif01();
                for(i=0;i<npara;i++)                
		{
		  //if(Viso[i]>999.){Lparameter[i]=Rparameter[i];}
		  //else if(Viso[i]<-999.){Lparameter[i]=Rparameter[i]+Viso[i];}
		  //else{
		    if (space1[i][2]<0.001){continue;} // skip the para don't require any purterbation
		    newv=perc*(space1[i][1]-space1[i][0])+space1[i][0];
		    //printf("f2 i=%d old %g new %g perc %g\n",i,outpara.parameter[i],newv,perc);
		    Lparameter[i]=newv; 
		  //}//else
                }//for i

          }//else if 2
	  else if(pflag==3)//set the perturbed parameter to have the same purterbation. e.g. uniform anisotropy in the mantle
	  {
		//printf("pglag=%d,npara=%d\n",pflag,inpara.npara);
		flag=0;
		for (i=0;i<npara;i++)
		{
		 
		 if (space1[i][2]<0.001)continue;
		 mean=Lparameter[i];
		 sigma=space1[i][2];
		 break;
		}
		if (i==npara){
		 // printf("gen_newPara, pflag==3, no para need to be purterbed! skip!\n");
		  return 1;
		  //return parameter;
		}
                while(flag<1)
                  {
                        newv=gen_random_normal(mean,sigma); // normal distribution
                        if(newv>space1[i][1] or newv<space1[i][0])continue;
                        else flag=2;
                  }//while flag
		perc=newv/mean;
		for(i=0;i<npara;i++)
		{
		  if (space1[i][2]<0.0001)continue; // skip the para don't require any purterbation
		  newv=Lparameter[i]*perc;
	//	  printf("f3 i=%d old %g new %g perc %g\n",i,outpara.parameter[i],newv,perc);
		  Lparameter[i]=newv;
		}
	  }//else if 3
	  else{printf("### wrong pflag for gen_newpara!!\n");exit(0);}
	return 1;
	}//int gen_newpara_single_L 
//-----------------------------------------------------	 

	int gen_newpara(paradef inpara, paradef &outpara, int pflag, int isoflag , const vector<float> &Viso)
	{//change the p.parameter[] 
	  vector<double> parameter;
	  vector<vector<double> > space1;
	  double cc;
	  double tmppara[20];
	  int i;
	  outpara=inpara;
	    if(isoflag>0){//isotropic case, Rparameter=Lparameter
		parameter=inpara.Rparameter;
		space1=inpara.Rspace1;
		gen_newpara_single(space1,parameter,inpara.Rnpara,pflag);
		outpara.Rparameter=parameter;
		outpara.Lparameter=parameter;
	    }//isoflag
	    else{//anisotropy case
		//cout<<"test hi anisotropy para\n";
		parameter.clear();
		parameter=outpara.Lparameter;
		gen_newpara_single(outpara.Rspace1,outpara.Rparameter,outpara.Rnpara,pflag);
		// something interesting here, seems the outpara.Lparameter cannot be written twice~, which means I canNOT do this: first gen_newpara, then change the parameter according to Viso
		// don't understand why, I need to lean more about C++
		gen_newpara_single_L(outpara.Lspace1,outpara.Rparameter,outpara.Lparameter,outpara.Lnpara,pflag,Viso);	
 	    }//else isoflag
	    vector<vector<double> >().swap(space1);
	    vector<double>().swap(parameter);
  	    return 1;
	}//gen_newpara

*/

//-----------------------------------------------------	 
/*
the inpara tells us what para need to be purterbed.
inpara.para0[i][j], i=npara.
outpara.para0=inpara.para0

j:	0		1		2		3		4		5
    0-value/Bcoeff	1-[i][2]=dv     dv or dv*100	sigma		group id	value id
    1-gp thickness	else...=dv*100  based on [i][1]
   -1- vpvs											  model
	==								==		==      ====>>>> outpara.parameter[] (value/thik/vpvs), length=npara=#of lines in file para.in(=>inpara.para0)
													 outpara.parameter.push_back(model.groups[ng].velue[nv]/thick/vpvs)
			-------------------------------------					====>>>> outpara.space1[min,max,sigma], length=npara
*/

	int mod2para(modeldef &model, paradef &inpara, paradef &outpara)
	{
	  int i,ng,nv,p0,p1,p3,p4,p5,p6,p7;
	  double tv,tmin,tmax,sigma,p2;
	  vector<double> vt(3,0.);

	  outpara=inpara;
	  outpara.parameter.clear();
	  outpara.L=model.data.L;
	  outpara.misfit=model.data.misfit;


	  for(i=0;i<inpara.npara;i++){
	  	p0=(int)inpara.para0[i][0]; // flag; 0--value/Bcoeff 1--gp thickness -1--vpvs
	  	p1=(int)inpara.para0[i][1]; // type of perturbation; 1--abs dv; else--percentation
	  	p2=(int)inpara.para0[i][2]; // perturbatin; dv
		ng=(int)inpara.para0[i][4]; // grounp number; 
		p6=(int)inpara.para0[i][6]; // ppflag; 1--vsv 2--vsh 3--vpv 4--vph 5--eta 6--theta 7--phi  8-rho; 9-h; 10-dVsvcos; 11-dVsvsin; 12-vpvs;

		if(p0==0){
		  p5=(int)inpara.para0[i][5];
		  nv=p5;
		  if(p6==1)
		  	tv=model.groups[ng].vsvvalue[nv];
		  else if (p6==2)
		  	tv=model.groups[ng].vshvalue[nv];
		  else if (p6==3)
		  	tv=model.groups[ng].vpvvalue[nv];
		  else if (p6==4)
		  	tv=model.groups[ng].vphvalue[nv];
		  else if (p6==5)
		  	tv=model.groups[ng].etavalue[nv];
		  else if (p6==6)
		  	tv=model.groups[ng].thetavalue[nv];
		  else if (p6==7)
		  	tv=model.groups[ng].phivalue[nv];
		  else if (p6==9)
			tv=model.groups[ng].rhovalue[nv];
		  else {	
			printf("### mod2para, wrong number for the 7th column in para.in!");
			exit(0);}
		}// p0==0
		else if (p0==1){
			tv=model.groups[ng].thick;
		}
		else if (p0==-1){
			tv=model.groups[ng].vpvs;
		}
		else{
			printf("#### mod2para, wrong number for the 1st column in para.in!;");
			exit(0);
		}
		outpara.parameter.push_back(tv);
	  }	  
	  if(inpara.flag<1){ //initial flag
		for(i=0;i<outpara.npara;i++){
			//---fill the flags in para0
			if(model.groups[ng].flagcpttype==2 or model.groups[ng].flagcpttype==4){
		  	  outpara.para0[i][7]=-1;//the LVflag, 1--use Vpara&Vkernel to do forward computation; -1--use Lovepara&Lovekernel to do forward computation
			}
			else{
			  outpara.para0[i][7]=1;
			}
			if(p6==1){
			  outpara.para0[i][8]=outpara.para0[i][9]=0;
			  outpara.para0[i][10]=2;
			}
			else if(p6==2){
			  outpara.para0[i][9]=1;
			  outpara.para0[i][10]=4;
			}
			else if(p6==3){
			  outpara.para0[i][8]=1;
			}
			else if(p6==4 or p6==5){
			  outpara.para0[i][8]=1;
			  outpara.para0[i][10]=2;
			}
			else if(p6==9){
			  outpara.para0[i][8]=1;
			  outpara.para0[i][9]=1;
			}
			else if(p6==10 or p6==11){
			  outpara.para0[i][8]=2;
			  outpara.para0[i][10]=2;
			}

			//---fill the para.space1
			sigma=outpara.para0[i][3]; // if don't want to perturb this parameter, then set sigma=0.
			if(p1==1){
				tmin=tv-p2;tmax=tv+p2;
			}
			else{
				tmin=tv*(1-p2/100.);tmax=tv*(1+p2/100.);
				printf("npara%d, v=%g [%g,%g]\n",i,tv,tmin,tmax);
			}
			if(p6<6){ // velocity or eta, should >0
				tmin=max(0.,tmin);
				tmax=max(0.,tmax);
				tmax=max(tmin+0.001,tmax);
				// there was a constraint on the sedimental velocity, don't know why need it, so it's not added here;
			}
			vt[0]=tmin;vt[1]=tmax;vt[2]=sigma;
			outpara.space1.push_back(vt);
	  	}// if inpara.flag<1
	  
	  }// for i
	  outpara.flag=1;
	
	  return 1;
	}// mod2para

//-----------------------------------------------------	 
/*
combine information in para.para0[] and para.parameter[], we can fill model.groups.value[]/thick/vpvs
for i<para.npara
        para.para0[i][j], var=para.parameter[i]
        j:      0                       4                       5
                flag                    group id
                0-value/Bcoeff  ----------------------------->  value id,nv     ===> model.groups[ng].value[nv]=para.parameter[i]
                1-gp thickness                                                  ===> model.groups[ng].thick=para.parameter[i]
                -1- vpvs                                                        ===> model.groups[ng].vpvs=para.parameter[i]


*/
	int para2mod(paradef para, modeldef inmodel, modeldef &outmodel)
	{
	  int i, p0,p1,ng,p6,nv,pflag;
 	  double newv;
	  int np,flagvpv,flagvph,flagvsv,flagvsh;
	  flagvpv=flagvph=flagvsv=flagvsh=-1;
	   
	  outmodel=inmodel;
	  for(i=0;i<para.npara;i++){
	  	p0=(int)para.para0[i][0];//flag, explanations are in mod2para;
	  	p1=(int)para.para0[i][1];
	  	ng=(int)para.para0[i][4];
	  	p6=(int)para.para0[i][6];//ppflag, explanations are in mod2para;
		
		newv=para.parameter[i];

		if(p0==0){ //value
		  nv=(int)para.para0[i][5];
		  if(p6==1)
		  	  {outmodel.groups[ng].vsvvalue[nv]=newv;flagvsv=1;}
		  else if (p6==2)
		  	  {outmodel.groups[ng].vshvalue[nv]=newv;flagvsh=1;}
		  else if (p6==3)
		          {outmodel.groups[ng].vpvvalue[nv]=newv;flagvpv=1;}
		  else if (p6==4)
		  	  {outmodel.groups[ng].vphvalue[nv]=newv;flagvph=1;}
		  else if (p6==5)
			  outmodel.groups[ng].etavalue[nv]=newv;
		  else if (p6==6)
			  outmodel.groups[ng].thetavalue[nv]=newv;
		  else if (p6==7)
			  outmodel.groups[ng].phivalue[nv]=newv;
		}//if p0==0

		else if(p0==1){// groups thickness
			//outmodel.tthick=outmodel.tthick-outmodel.gorups[ng].thick+newv; // the tthick is changing
			outmodel.groups[outmodel.ngroup-1].thick=outmodel.groups[outmodel.ngroup-1].thick+(outmodel.groups[ng].thick-newv);//change the thickness of the last group while keeping the tthcik the same
			outmodel.groups[ng].thick=newv;
		
		}// else if p0==1
		
		else if (p0==-1){//vpvs
			outmodel.groups[ng].vpvs=newv;
		}
		else{cout<<"########wrong!!!!! in para2mode"<<endl;exit(0);}
	  }// for i


	  //---depending on the pflag of each group, fill all the unfilled vectors among those vsv...phi 7 vectors with appropriate values derived from certain relationships between vp&vs;

	  for(ng=0;ng<outmodel.ngroup;ng++){
		np=outmodel.groups[ng].np;
		pflag=outmodel.groups[ng].pflag;
	  	if(pflag<7){
			for(i=0;i<np;i++)
				{outmodel.groups[ng].thetavalue[i]=outmodel.groups[ng].phivalue[i]=0.;}
			
			if(pflag<5){
			  for(i=0;i<np;i++)outmodel.groups[ng].etavalue[i]=1.0;
			  
			  if(pflag==3){
			    for(i=0;i<np;i++){
			      if(flagvpv>0)outmodel.groups[ng].vphvalue[i]=outmodel.groups[ng].vpvvalue[i];
			      else if(flagvph>0)outmodel.groups[ng].vpvvalue[i]=outmodel.groups[ng].vphvalue[i];
			    }
			  }//pflag==3
			  
			  else if(pflag==2){
			      for(i=0;i<np;i++){
				      outmodel.groups[ng].vphvalue[i]=outmodel.groups[ng].vpvvalue[i]=outmodel.groups[ng].vsvvalue[i]*outmodel.groups[ng].vpvs;

			      }  
			  }//pflag==2

			  else if(pflag==1){
			  	if(flagvsv>0)outmodel.groups[ng].vshvalue[i]=outmodel.groups[ng].vsvvalue[i];
				else if(flagvsh>0)outmodel.groups[ng].vsvvalue[i]=outmodel.groups[ng].vshvalue[i];
				outmodel.groups[ng].vphvalue[i]=outmodel.groups[ng].vpvvalue[i]=outmodel.groups[ng].vsvvalue[i]*outmodel.groups[ng].vpvs;
			  }//pflag==1 

			}//pflag<5		
		} // if pflag<7	  
	  } // for ng

	  return 1;
	}// para2mod
	 
	
//-----------------------------------------------------	 
/*        int get_misfitDISPTibet(paradef &para, modeldef &inmodel,modeldef &outmodel,int cflag,float depcri1,float depcri2,float qpcri,float qscri,double permin,double permax)
        {
          paradef temppara;
          temppara=para;
          para2mod(para,inmodel,outmodel);
          updatemodelTibet(outmodel,depcri1,depcri2,qpcri,qscri);
          compute_disp(outmodel,cflag);
          compute_misfitDISP(outmodel);
          compute_misfitDISP2(outmodel,permin,permax);
          mod2para(outmodel,temppara,para);
          return 1;
        }//get misfit
*/
//----------------------------------------------------- 
	int para_avg(vector<paradef> &paralst, paradef &paraavg,vector<double> &Rparastd,vector<double> &Lparastd, vector<int> &idlst){
	//get the average parameters from a list of acceptable paras
	    int i,j,k,size,Ngood;
	    double mismin=1e10;

	    initpara(paraavg);
	    size=paralst.size();
	    if(size<1) return 0;
		
	    paraavg=paralst[0];
	    Rparastd.clear();Lparastd.clear();
/*	    initpara(paraavg);
	    paraavg.Rpara0=paralst[0].Rpara0;
	    paraavg.Lpara0=paralst[0].Lpara0;
	    size=paralst.size();
	    if(size<1) return 0;

	    for(j=0;j<paralst[0].Rnpara){//Rnpara==Lnpara
		paraavg.Rparameter.push_back(0.);
		paraavg.Lparameter.push_back(0.);
		Rparastd.push_back(0.);
		Lparastd.push_back(0.);
	    }//for j
*/

	    for(j=0;j<paralst[0].Rnpara;j++){//Rnpara==Lnpara
                Rparastd.push_back(0.);
                Lparastd.push_back(0.);
		paraavg.Rparameter[j]=paraavg.Lparameter[j]=0.;
            }//for j

	    for(i=0;i<size;i++)
		if(paralst[i].misfit<mismin) mismin=paralst[i].misfit;
	    //min(2*Mis_min,Mis_min+0.5)
	    cout<<"mismin="<<mismin<<endl;
	    if(mismin>0.5) mismin=mismin*2;
	    else mismin=mismin+0.5;
	    //  mismin=mismin+0.5;
	    // mismin=1e10;
	    Ngood=0;
	    for(i=0;i<size;i++){
		if(paralst[i].misfit<mismin){
		  Ngood++;
		  idlst.push_back(i);
		  for(j=0;j<paralst[0].Rnpara;j++){
			paraavg.Rparameter[j]=paraavg.Rparameter[j]+paralst[i].Rparameter[j];
			paraavg.Lparameter[j]=paraavg.Lparameter[j]+paralst[i].Lparameter[j];}
		}//if
	    }//for i
	   for(i=0;i<paralst[0].Rnpara;i++){
		paraavg.Rparameter[i]=paraavg.Rparameter[i]/Ngood;
		paraavg.Lparameter[i]=paraavg.Lparameter[i]/Ngood;
	   }//for i
	   cout<<"Ngood="<<Ngood<<endl; 
	   printf("size =%d, Ngood=%d mismin=%g\n",size,Ngood,mismin);
	   //if(Ngood>size){printf("WRONG HERE! Ngood > size!!\n");exit(0);}
	   for(i=0;i<paralst[0].Rnpara;i++){
		for(j=0;j<Ngood;j++){
		    k=idlst[j];
		    Rparastd[i]=Rparastd[i]+pow(paralst[k].Rparameter[i]-paraavg.Rparameter[i],2);
		    Lparastd[i]=Lparastd[i]+pow(paralst[k].Lparameter[i]-paraavg.Lparameter[i],2);
		}
		Rparastd[i]=sqrt(Rparastd[i]/Ngood);
		Lparastd[i]=sqrt(Lparastd[i]/Ngood);
	    }//for i
	 
	    return 1;
	}//paramodel_avg
//----------------------------------------------------------------
	int para_best(vector<paradef> paralst,paradef &parabest){
	  int size,i,idbest;
	  double misfit,misfitbest;
  	 
  	  misfitbest=1e10;
	  size=paralst.size();
	  for(i=0;i<size;i++){
		misfit=paralst[i].misfit;
		if(misfit<misfitbest){misfitbest=misfit;idbest=i;}
	  }//for i	
	  parabest=paralst[idbest];
	}//para_best

//-----------------------------------------------------	
	int model_avg_sub(vector<vector<float> > &vlst,vector<vector<float> > &stdlst,vector<float> &hlst,vector<modeldef> &modlst,int ng,float h0,float h,float dth)
	{
	  //h0=0;h=Hsed;hstd=hsedstd;dth=0.4;
	  //vlst[ithick]=[vsv,vsh,vsvmin,vsvmax,vshmin,vshmax];
	  vector<float> tv,dep2lst,dep1lst;
	  vector<int> iddep1lst;
	  float th,tth,fm1,fm2,dep;
	  float vsv,vsh,dep1,dep2,vsvmin,vsvmax,vshmin,vshmax;
	  int i,j,Ncount,size,k;
	  vector<vector<float> > tvlst;

	  hlst.clear();vlst.clear();stdlst.clear();tvlst.clear();
		  
	  size=modlst.size();
	
	  for(i=0;i<size;i++){	  
	    for(k=0;k<=ng;k++){dep2=dep2+modlst[i].groups[k].thick;}
	    dep2lst.push_back(dep2);
	    dep1lst.push_back(0.);
	    iddep1lst.push_back(0);
	  }
	  for(th=h0;th<h;th=th+dth){//thickness
		printf("th=%g h=%g\n",th,h);//---test---
		vsv=0.;vsh=0.;Ncount=0;fm1=0.;fm2=0.;
		tvlst.clear();
		vsvmin=1e10;vsvmax=-1;
		vshmin=1e10;vshmax=-1;
		for(i=0;i<size;i++){//iterate within model list
		  //dep2=0.;
		  //for(k=0;k<=ng;k++){dep2=dep2+modlst[i].groups[k].thick;}
		  dep2=dep2lst[i];
		  //dep1=dep2-modlst[i].groups[ng].thick;//dep1--upper of the layer; dep2--lower of the layer
		  if(th>dep2)continue;//tth=dep2; //tth = min(dep2,th) ####if th>dep2, continue, skip present mod, jumpt to next mod??
		  else tth=th;
		  
		  dep1=0;
		  //for(j=0;j<modlst[i].laym0.nlayer;j++){
		  //	dep1=dep1+modlst[i].laym0.thick[j];
		  for(j=iddep1lst[i];j<modlst[i].laym0.nlayer;j++){
			dep1=dep1lst[i]+modlst[i].laym0.thick[j];
			if(dep1>=tth ){
				vsv=vsv+modlst[i].laym0.vsv[j];
				vsh=vsh+modlst[i].laym0.vsh[j];
				Ncount++;	
				if (modlst[i].laym0.vsv[j]>vsvmax)vsvmax=modlst[i].laym0.vsv[j];
				if(modlst[i].laym0.vsv[j]<vsvmin)vsvmin=modlst[i].laym0.vsv[j];
				if(modlst[i].laym0.vsh[j]>vshmax)vshmax=modlst[i].laym0.vsh[j];
				if(modlst[i].laym0.vsh[j]<vshmin)vshmin=modlst[i].laym0.vsh[j];
				tv.clear();tv.push_back(modlst[i].laym0.vsv[j]);tv.push_back(modlst[i].laym0.vsh[j]);
				tvlst.push_back(tv);
				dep1lst[i]=dep1-modlst[i].laym0.thick[j];//
				iddep1lst[i]=j;//
				break;
			}//if dep>tth
		  }//for j Nlay
		}//for i Nmodel,size
		vsv=vsv/Ncount;
		vsh=vsh/Ncount;
		tv.clear();
		tv.push_back(vsv);tv.push_back(vsh);tv.push_back(vsvmin);tv.push_back(vsvmax);tv.push_back(vshmin);tv.push_back(vshmax);
		vlst.push_back(tv);
		hlst.push_back(th);
		for(i=0;i<Ncount;i++){ 
			fm1=fm1+pow(tvlst[i][0]-vsv,2);
			fm2=fm2+pow(tvlst[i][1]-vsv,2);
		}
		fm1=sqrt(fm1/Ncount);
		fm2=sqrt(fm2/Ncount);
		tv.clear();
		tv.push_back(fm1);tv.push_back(fm2);
		stdlst.push_back(tv);
	  }//for th
        vector<vector<float> >().swap(tvlst);
	vector<float>().swap(tv);
	vector<float>().swap(dep2lst);
	vector<float>().swap(dep1lst);
	vector<int>().swap(iddep1lst);
	return 1;
	}//mod avg sub
//------------------
 	int integrate_mod(vector<vector<float> > &outVsv,vector<vector<float> > &outVsh, vector<float> &hlst,vector<vector<float> > &vlst,vector<vector<float> > &stdlst)
	{ int i;vector<float> tv;
          for(i=0;i<vlst.size();i++){
                tv.clear();
                tv.push_back(hlst[i]);tv.push_back(vlst[i][0]);tv.push_back(stdlst[i][0]);tv.push_back(vlst[i][2]);tv.push_back(vlst[i][3]);
                outVsv.push_back(tv);
                tv.clear();
                tv.push_back(hlst[i]);tv.push_back(vlst[i][1]);tv.push_back(stdlst[i][1]);tv.push_back(vlst[i][4]);tv.push_back(vlst[i][5]);
                outVsh.push_back(tv);}
	  vector<float>().swap(tv);
	  return 1;
	}//integrate_mod
//-----------------------------------------------------------------
	int connect_modmin(vector<vector<float> > &outVsv,vector<vector<float> > &outVsh,vector<float> &hlst,vector<vector<float> > &vlst,vector<vector<float> > &stdlst,float h1,float h2)
	{
	  int i;vector<float> tv;
	  float h,v;
 	  for(i=0;i<vlst.size();i++){
	  	h=hlst[i];
		if (h<h1 or h>=h2)continue;
		tv.clear();
		v=vlst[i][0]-stdlst[i][0];
		tv.push_back(h);tv.push_back(v);
		outVsv.push_back(tv);
		tv.clear();
		v=vlst[i][1]-stdlst[i][1];
		tv.push_back(h);tv.push_back(v);
		outVsh.push_back(tv);
	  }
	  vector<float>().swap(tv);
	  return 1;
	}//connect_modmin
//-----------------------------------------------------------------
	int connect_modmid(vector<vector<float> > &outVsv,vector<vector<float> > &outVsh,vector<float> &hlst,vector<vector<float> > &vlst,vector<vector<float> > &stdlst,float h1,float h2)
	{
	  int i;vector<float> tv;
	  float h,v;
 	  for(i=0;i<vlst.size();i++){
	  	h=hlst[i];
		if (h<h1 or h>=h2)continue;
		tv.clear();
		v=vlst[i][0];
		tv.push_back(h);tv.push_back(v);
		outVsv.push_back(tv);
		tv.clear();
		v=vlst[i][1];
		tv.push_back(h);tv.push_back(v);
		outVsh.push_back(tv);
	  }
	  vector<float>().swap(tv);
	  return 1;
	}//connect_modmin
//-----------------------------------------------------------------
//-----------------------------------------------------------------
	int connect_modmax(vector<vector<float> > &outVsv,vector<vector<float> > &outVsh,vector<float> &hlst,vector<vector<float> > &vlst,vector<vector<float> > &stdlst,float h1,float h2)
	{
	  int i;vector<float> tv;
	  float h,v;
 	  for(i=0;i<vlst.size();i++){
	  	h=hlst[i];
		if (h<h1 or h>=h2)continue;
		tv.clear();
		v=vlst[i][0]+stdlst[i][0];
		tv.push_back(h);tv.push_back(v);
		outVsv.push_back(tv);
		tv.clear();
		v=vlst[i][1]+stdlst[i][1];
		tv.push_back(h);tv.push_back(v);
		outVsh.push_back(tv);
	  }
	  vector<float>().swap(tv);
	  return 1;
	}//connect_modmin
//-----------------------------------------------------------------
//-----------------------------------------------------------------
	int write_modavg( const char *fvsvnm, const char *fvshnm, vector<vector<float> > &Vsv1,vector<vector<float> > &Vsv2,vector<vector<float> > &Vsv3,vector<vector<float> > &Vsh1,vector<vector<float> > &Vsh2,vector<vector<float> > &Vsh3, float Hsed,float Hsedstd,float Hmoho,float Hmohostd )
	{
	    FILE *ftmp;
	    int i,j;
	    //cout<<"BEGIN!\n";//--test--
	    //printf("fvsvnm=%s fvshnm=%s, size, vsv1=%d vsv2=%d vsv3=%d vsh1=%d vsh2=%d vsh3=%d, Hsee=%g,Hsedstd=%g,Hmoho=%g,Hmohostd=%g\n",fvsvnm,fvshnm,Vsv1.size(),Vsv2.size(),Vsv3.size(),Vsh1.size(),Vsh2.size(),Vsh3.size());//---test---
	    if((ftmp=fopen(fvsvnm,"w"))==NULL){printf("####in write-modavg, canot open file %s to write!\n",fvsvnm);exit(0);};
	    //cout<<"hi1\n";//--test--
	    fprintf(ftmp,"sedi %8g %8g\n",Hsed,Hsedstd);
	    //cout<<"hi2\n";//--test--
	    fprintf(ftmp,"moho %8g %8g\n",Hmoho,Hmohostd);
	    //cout<<"hi3\n";//--test--
	    for(i=0;i<Vsv1.size();i++){
	       //fprintf(ftmp,"%8g %8g %8g %8g\n",Vsv1[i][0],Vsv1[i][1],Vsv2[i][1],Vsv3[i][1]);// if everything works fine, the hlst in Vsv123 should be the same
	       fprintf(ftmp,"%8g %8g %8g %8g %8g %8g\n",Vsv1[i][0],Vsv1[i][1],Vsv2[i][0],Vsv2[i][1],Vsv3[i][0],Vsv3[i][1]);// if everything works fine, the hlst in Vsv123 should be the same
	    }
	    fclose(ftmp);
	    //cout<<"hi4\n";//--test--	   
	    if((ftmp=fopen(fvshnm,"w"))==NULL){printf("####in write-modavg, canot open file %s to write!\n",fvshnm);exit(0);};
	    //cout<<"hi5\n";//--test--
	    fprintf(ftmp,"sedi %8g %8g\n",Hsed,Hsedstd);
	    //cout<<"hi6\n";//--test--
	    fprintf(ftmp,"moho %8g %8g\n",Hmoho,Hmohostd);
	    //cout<<"hi7\n";//--test--
	    for(i=0;i<Vsh1.size();i++){
	       fprintf(ftmp,"%8g %8g %8g %8g %8g %8g\n",Vsh1[i][0],Vsh1[i][1],Vsh2[i][0],Vsh2[i][1],Vsh3[i][0],Vsh3[i][1]);// if everything works fine, the hlst in Vsv123 should be the same
	    }
	    fclose(ftmp);
	    return 1;
	}//write_modavg
//--------------------------------------------------
	int model_avg( const char *fvsvnm, const char *fvshnm, vector<paradef> &paralst, modeldef modelref,float depcri1,float depcri2,float qpcri, float qscri, vector<int> &idlst)
	{
	  int size,i,n,k;
 	  paradef para;
	  modeldef mod;
	  char vsvnm[200],vshnm[200];
	  float Hsed=0.,Hmoho=0.,Hsedstd,Hmohostd,fm1=0.,fm2=0.;
	  vector<vector<float> > v1lst,v2lst,v3lst,std1lst,std2lst,std3lst;
	  vector<vector<float> > Vsvmin,Vsvmid,Vsvmax,Vshmin,Vshmid,Vshmax;
	  vector<float> h1lst,h2lst,h3lst,tv;
	  time_t start;
	  vector<modeldef> modlst;
	  cout<<"---begin! "<<time(0)<<endl;start=time(0);
	  size = idlst.size();  
	  if(size<1){printf("##### in mod_avg, idlst is empty!\n");exit(0);}
	  for (i=0;i<size;i++){
		k=idlst[i];
		para=paralst[k];
		para2mod(para,modelref,mod);
		updatemodelTibet(mod,depcri1,depcri2,qpcri,qscri);
		modlst.push_back(mod);
		Hsed=Hsed+mod.groups[0].thick;
		Hmoho=Hmoho+mod.groups[1].thick;
	  }//for i
	  Hsed=Hsed/size;
	  Hmoho=Hmoho/size;
	  for(i=0;i<size;i++){
	 	fm1=fm1+pow(modlst[i].groups[0].thick-Hsed,2);
	 	fm2=fm2+pow(modlst[i].groups[1].thick-Hmoho,2);
	  }
	  Hsedstd=sqrt(fm1/size);
	  Hmohostd=sqrt(fm2/size);
	  cout<<"----get_avg_sub---"<<time(0)-start<<"\n";//---test----
	  model_avg_sub(v1lst,std1lst,h1lst,modlst,0,0,Hsed+Hsedstd,0.5);
	  cout<<"finishsed\n";//--test
	  model_avg_sub(v2lst,std2lst,h2lst,modlst,1,Hsed-Hsedstd,Hmoho+Hmohostd,1.);
	  cout<<"finish cst\n";//---test
	  model_avg_sub(v3lst,std3lst,h3lst,modlst,2,Hmoho-Hmohostd,150,5.);
	  cout<<"finish mant\n";//---test
          vector<modeldef>().swap(modlst);
	  
	  cout<<"---connect model min/mid/max "<<time(0)-start<<"\n";//---test---
	  connect_modmin(Vsvmin,Vshmin,h1lst,v1lst,std1lst,0,Hsed+Hsedstd);
	  connect_modmin(Vsvmin,Vshmin,h2lst,v2lst,std2lst,Hsed+Hsedstd,Hmoho+Hmohostd);
	  connect_modmin(Vsvmin,Vshmin,h3lst,v3lst,std3lst,Hmoho+Hmohostd,200);
	  
	  connect_modmid(Vsvmid,Vshmid,h1lst,v1lst,std1lst,0,Hsed);
	  connect_modmid(Vsvmid,Vshmid,h2lst,v2lst,std2lst,Hsed,Hmoho);
	  connect_modmid(Vsvmid,Vshmid,h3lst,v3lst,std3lst,Hmoho,200);
	  
	  connect_modmax(Vsvmax,Vshmax,h1lst,v1lst,std1lst,0,Hsed-Hsedstd);
	  connect_modmax(Vsvmax,Vshmax,h2lst,v2lst,std2lst,Hsed-Hsedstd,Hmoho-Hmohostd);
	  connect_modmax(Vsvmax,Vshmax,h3lst,v3lst,std3lst,Hmoho-Hmohostd,200);
	  
          vector<float>().swap(h1lst);
          vector<float>().swap(h2lst);
          vector<float>().swap(h3lst);
          vector<float>().swap(tv);
          vector<vector<float> >().swap(v1lst);
          vector<vector<float> >().swap(v2lst);
          vector<vector<float> >().swap(v3lst);
          vector<vector<float> >().swap(std1lst);
          vector<vector<float> >().swap(std2lst);
          vector<vector<float> >().swap(std3lst);

	  cout<<"----write model avg----"<<time(0)-start<<"\n";//----test----
	  write_modavg(fvsvnm,fvshnm,Vsvmin,Vsvmid,Vsvmax,Vshmin,Vshmid,Vshmax,Hsed,Hsedstd,Hmoho,Hmohostd);

          vector<vector<float> >().swap(Vsvmin);
          vector<vector<float> >().swap(Vsvmid);
          vector<vector<float> >().swap(Vsvmax);
          vector<vector<float> >().swap(Vshmin);
          vector<vector<float> >().swap(Vshmid);
          vector<vector<float> >().swap(Vshmax);


	  return 1;
	}// model_avg
//-----------------------------------
