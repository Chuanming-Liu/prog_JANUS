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
int readpara_single(vector<vector<string>> &para0, int &npara,const char* fname,int surflag)
int readpara(paradef &para,const char* Rfname, const char* Lfname, int Rsurflag, int Lsurflag)
int gen_newpara_single ( vector<vector<double> > space1, vector<double > &parameter, int npara,int pflag)
int gen_newpara(paradef inpara, paradef &outpara, int pflag,int isoflag ,vector<float> &Viso)
int mod2para(modeldef &model, paradef &inpara, paradef &outpara,int Rflag, int Lflag)
int para2mod(paradef para, modeldef inmodel, modeldef &outmodel)

//============================
*/
//class paracal{
//public:
	int initpara(paradef &para)
	{
	  para.Rnpara=0;para.Lnpara=0;
	  para.L=0.;
	  para.misfit=0.;
	  para.flag=0;
	  //para.Rparameter.clear();para.Rpara0.clear();para.Rspace1.clear();
	  //para.Lparameter.clear();para.Lpara0.clear();para.Lspace1.clear();
	  vector<double>().swap(para.Rparameter);vector<double>().swap(para.Lparameter);
	  vector<vector<double> >().swap(para.Rspace1);vector<vector<double> >().swap(para.Lspace1);
	  vector<vector<string> >().swap(para.Rpara0);vector<vector<string> >().swap(para.Lpara0);
	  return 1;
	}
//-----------------------------------------------------	 
	int readpara_single(vector<vector<string> > &para0, int &npara,const char* fname,int surflag)
	{  
	  fstream mff;
	  int k,i;
	  string line;
	  vector<string> v;

	  if (surflag>0){
	  mff.open(fname);
	  if ( not mff.is_open()){cout<<"########para file"<<fname<<"does not exist!\n exit!!\n";exit (0);}
	  k=0;
	  para0.clear();
	  while(getline(mff,line))
	  {
		v.clear();
		Split(line,v," ");
		k=k+1;
		para0.push_back(v);	
	  }
	  mff.close();
	  npara=k;
	  }//if surflag		
	  return 1;
	}//readpara_single
//----------------------------------------------------
	int readpara(paradef &para,const char* Rfname, const char* Lfname, int Rsurflag, int Lsurflag)	
	{
	  readpara_single(para.Rpara0,para.Rnpara,Rfname,Rsurflag);
	  readpara_single(para.Lpara0,para.Lnpara,Lfname,Lsurflag);
	  printf("-----check readpara: Rnpara=%d Lnpara=%d\n",para.Rnpara,para.Lnpara);
	  return 1;
	}//readpara
//-----------------------------------------------------	 
// randomly generating new parameters, uniform or normal distribution
	int gen_newpara_single (  vector<vector<double> > space1, vector<double> &parameter, int npara,int pflag)
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
		cout<<"test--- i="<<i<<endl;
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
		//outpara.Lnpara=outpara.Rnpara;
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
	}



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
	int mod2para(modeldef &model, paradef &inpara, paradef &outpara,int Rflag, int Lflag)
	{//both Lpara and Rpara are filled in. ****** is this an appropriate way ?? *********
	 //L and R have the same thick/vpvs value
	  int i,ngr,ngl,nvr,nvl;
	  double tvr,tvl,tminr,tminl,tmaxr,tmaxl,sigmar,sigmal;
	  int p0l,p1l,p3l,p4l,p5l,p0r,p1r,p3r,p4r,p5r;
	  double p2l,p2r;
	  vector<double> vtrr(3,0.),vtrl(3,0.);
	  //*******shall I consider the isoflag ?? **********maybe not?? it is considered in para2mod function
	  if(Rflag<=0){inpara.Rpara0=inpara.Lpara0;}
 	  else if (Lflag<=0){inpara.Lpara0=inpara.Rpara0;}

	  outpara=inpara;
	  outpara.Lparameter.clear();outpara.Rparameter.clear();
	//  outpara.Lspace1.clear();outpar.Rspace1.clear();
	  outpara.L=model.data.L;
	  outpara.misfit=model.data.misfit;


	  for(i=0;i<inpara.Rnpara;i++)
	  {
	  	p0r=atoi(inpara.Rpara0[i][0].c_str());p0l=atoi(inpara.Lpara0[i][0].c_str());//same, sign
		p1r=atoi(inpara.Rpara0[i][1].c_str());p1l=atoi(inpara.Lpara0[i][1].c_str());//same, sign
		p2r=atof(inpara.Rpara0[i][2].c_str());p2l=atof(inpara.Lpara0[i][2].c_str());//diff, perturbation
		ngr=atoi(inpara.Rpara0[i][4].c_str());ngl=atoi(inpara.Lpara0[i][4].c_str());// group id
		if (p0r==0)//value. vs or Bspline coeff
			{
			 p5r=atoi(inpara.Rpara0[i][5].c_str());p5l=atoi(inpara.Lpara0[i][5].c_str());//same
			 nvr=p5r;nvl=p5l;
			 tvr=model.groups[ngr].Rvalue[nvr];
			 tvl=model.groups[ngl].Lvalue[nvl];}
		else if(p0r==1)// thickness of group
			{tvr=model.groups[ngr].thick;tvl=tvr;}
		else if(p0r==-1)// vpvs
			{tvr=model.groups[ngr].vpvs;tvl=tvr;}
		else
			{cout<<"#####Wrong!!! mode2para\n";
			exit (0);}
		outpara.Rparameter.push_back(tvr);
		outpara.Lparameter.push_back(tvl);
		if(inpara.flag<1) //initial flag;
		{
		  sigmar=atof(inpara.Rpara0[i][3].c_str());
		  sigmal=atof(inpara.Lpara0[i][3].c_str());
		  if(p1r==1)
			{
			  tminr=tvr-p2r;tminl=tvl-p2l;
			  tmaxr=tvr+p2r;tmaxl=tvl+p2l;
			}
		  else
			{
			  tminr=tvr-tvr*p2r/100.;tminl=tvl-tvl*p2l/100;
			  tmaxr=tvr+tvr*p2r/100.;tmaxl=tvl+tvl*p2l/100;
		          printf("-L:%8g %8g %8g %8g -R:%8g %8g %8g %8g\n",tvl,tminl,tmaxl,p2l,tvr,tminr,tmaxr,p2r);
			}
		  tminr=max(0.,tminr);tmaxr=max(0.,tmaxr);
		  tminl=max(0.,tminl);tmaxl=max(0.,tmaxl);
		  tmaxr=max(tminr+0.0001,tmaxr);tmaxl=max(tminl+0.0001,tmaxl);//************how about zero perturbation?? A: if i remembered right, if the sigma is set to be 0., then the perturbation would be zero.
		  if(p0r+i==0 and p5r==0) //p0=i=p5=0 => upper sediment
			{tminr=max(0.2,tminr); tmaxr=max(0.5,tmaxr);
			 tminl=max(0.2,tminl); tmaxl=max(0.5,tmaxl);} 
		  vtrr[0]=tminr; vtrr[1]=tmaxr; vtrr[2]=sigmar;
		  vtrl[0]=tminl; vtrl[1]=tmaxl; vtrl[2]=sigmal;
		  outpara.Rspace1.push_back(vtrr);
		  outpara.Lspace1.push_back(vtrl);
		  outpara.flag=1;
		}//if flag
	  }//for
	 return 1;
	}//mod2para
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
	//int para2mod(paradef &para, modeldef &inmodel, modeldef &outmodel)
	int para2mod(paradef para, modeldef inmodel, modeldef &outmodel)
	{
	  double newv,newvl;
	  int ng,nv,p0,i,ngl,nvl;
	  outmodel=inmodel;

	  for(i=0;i<para.Rnpara;i++){//Rnpara=Lnpara!!
                  p0=atoi(para.Rpara0[i][0].c_str());
                  newv=para.Rparameter[i];
		  newvl=para.Lparameter[i];
                  ng=atoi(para.Rpara0[i][4].c_str());
                  ngl=atoi(para.Lpara0[i][4].c_str());
                  if(p0==0){//value
                         nv=atoi(para.Rpara0[i][5].c_str());
                         nvl=atoi(para.Lpara0[i][5].c_str());
                         outmodel.groups[ng].Rvalue[nv]=newv;
                         outmodel.groups[ngl].Lvalue[nvl]=newvl; }//if
                  else if(p0==1){//group thickness
			//****** if in some cases you need to purterbing the thickness while inverting Love wave only.(i.e. the thickness in Rpara is constant, Rflag=0), need to add one more sentence considering the Rflag==0 case.***********OR maybe just make the Rpara.in have non-zero thickness space, then the perturbed thickness in R would be transfered to L (Mar. 27,2012)//
			outmodel.tthick=outmodel.tthick-outmodel.groups[ng].thick+newv; //there was a bug here, fixed on Feb 20.
                        outmodel.groups[ng].thick=newv;}
                  else if(p0==-1){//vpvs
                        outmodel.groups[ng].vpvs=newv;}
                  else {cout<<"########wrong!!!!! in para2mode"<<endl;exit(0);}
                }//for i        

	  return 1;
	}//para2mod

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
	   printf("size =%d, Ngood=%d misfit_cri=%g\n",size,Ngood,mismin);
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
	    cout<<"BEGIN!\n";//--test--
	    printf("fvsvnm=%s fvshnm=%s, size, vsv1=%d vsv2=%d vsv3=%d vsh1=%d vsh2=%d vsh3=%d, Hsee=%g,Hsedstd=%g,Hmoho=%g,Hmohostd=%g\n",fvsvnm,fvshnm,Vsv1.size(),Vsv2.size(),Vsv3.size(),Vsh1.size(),Vsh2.size(),Vsh3.size());//---test---
	    if((ftmp=fopen(fvsvnm,"w"))==NULL){printf("####in write-modavg, canot open file %s to write!\n",fvsvnm);exit(0);};
	    cout<<"hi1\n";//--test--
	    fprintf(ftmp,"sedi %8g %8g\n",Hsed,Hsedstd);
	    cout<<"hi2\n";//--test--
	    fprintf(ftmp,"moho %8g %8g\n",Hmoho,Hmohostd);
	    cout<<"hi3\n";//--test--
	    for(i=0;i<Vsv1.size();i++){
	       //fprintf(ftmp,"%8g %8g %8g %8g\n",Vsv1[i][0],Vsv1[i][1],Vsv2[i][1],Vsv3[i][1]);// if everything works fine, the hlst in Vsv123 should be the same
	       fprintf(ftmp,"%8g %8g %8g %8g %8g %8g\n",Vsv1[i][0],Vsv1[i][1],Vsv2[i][0],Vsv2[i][1],Vsv3[i][0],Vsv3[i][1]);// if everything works fine, the hlst in Vsv123 should be the same
	    }
	    fclose(ftmp);
	    cout<<"hi4\n";//--test--	   
	    if((ftmp=fopen(fvshnm,"w"))==NULL){printf("####in write-modavg, canot open file %s to write!\n",fvshnm);exit(0);};
	    cout<<"hi5\n";//--test--
	    fprintf(ftmp,"sedi %8g %8g\n",Hsed,Hsedstd);
	    cout<<"hi6\n";//--test--
	    fprintf(ftmp,"moho %8g %8g\n",Hmoho,Hmohostd);
	    cout<<"hi7\n";//--test--
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
