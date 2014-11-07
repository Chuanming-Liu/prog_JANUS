/*#include<fstream>
#include<iostream>
#include<vector>
#include<math.h>
#include<algorithm>
*/
//#include"INITstructure.h"
//#include"CALgroup.C"

//#include"/home/jiayi/progs/jy/HEAD/head_c++/string_split.C"
// this version,  in 'goodmodel', verified the gradient of the model. 
// this rf version, the updatemodelRF and compute_misfit (which takes the misfit of both DISP and RF into account) are added
using namespace std;
/*=====CONTENT=======
int initdisp (dispdef &disp)
int initmodel(modeldef &model)
int readdisp_single(dispdef &disp,vector<string> name, int surflag)
int readdisp(modeldef &model,vector<string> Rname,vector<string> Lname,int Rsurflag,int Lsurflag)
***surflag 0no/1phase/2group/3ph&gp
int readrf(modeldef &model,char *name)
int readmodIso(modeldef &model,const char *name)
int readmodAniso(modeldef &model,const char *name)
int updatemodelTibet(modeldef &model,float depcri1,float depcri2,float qpcri,float qscri)
int updatemodelRF(modeldef &model,float depcri1,float depcri2,float qpcri,float qscri)
int compute_misfitDISP_single(dispdef &disp,double &tempv1,double &tempv2,double &tS)
int compute_misfit(modeldef &model,int Rflag, int Lflag,float inp)
int goodmodel( modeldef &model, vector<int> vmono, vector<int> vgrad,int Rflag, int Lflag, int isoflag)
int model_avg
//===================
*/
//class modelcal{
//public:
	int initdisp (dispdef &disp)
	{
	  disp.npper=0;
	  disp.pL=1.;  //?????
	  disp.pmisfit=0.;
	  disp.fphase=0;
	  disp.ngper=0;
	  disp.gL=1.; //?????
	  disp.gmisfit=0.;
	  disp.fgroup=0;
	  disp.L=0.;
	  disp.misfit=0.;
	  return 1;	
	}
	int initmodel(modeldef &model)
	{
	  model.laym0.nlayer=0;
	  model.data.rf.nrfo=0;
	  model.data.rf.rt=0;
	  model.data.rf.L=0.;
	  model.data.rf.misfit=0.;
	  initdisp(model.data.Rdisp);
	  initdisp(model.data.Ldisp);
	  model.data.p=0.;
	  model.data.L=0.;
	  model.data.misfit=0.;
	  model.ngroup=0;
	  model.flag=0;
	  model.cc=2;
	  model.tthick=0.;//total thickness or all groups; computed when readmodel;
	  model.groups.clear();
	  model.laym0.vsv.clear();model.laym0.vsh.clear();model.laym0.vp.clear();model.laym0.vpvs.clear();model.laym0.rho.clear();model.laym0.qs.clear();model.laym0.qp.clear();model.laym0.thick.clear();
	  model.data.Rdisp.pper.clear();model.data.Rdisp.pvelo.clear();model.data.Rdisp.pvel.clear();model.data.Rdisp.unpvelo.clear();model.data.Rdisp.gper.clear();model.data.Rdisp.gvelo.clear();model.data.Rdisp.gvel.clear();model.data.Rdisp.ungvelo.clear();model.data.Rdisp.period1.clear();
	  model.data.Ldisp.pper.clear();model.data.Ldisp.pvelo.clear();model.data.Ldisp.pvel.clear();model.data.Ldisp.unpvelo.clear();model.data.Ldisp.gper.clear();model.data.Ldisp.gvelo.clear();model.data.Ldisp.gvel.clear();model.data.Ldisp.ungvelo.clear();model.data.Ldisp.period1.clear();
	  
	  return 1;
	}
	
//-----------------------------------------------------	
	int readdisp_single(dispdef &disp,vector<string> name, int surflag)//use vector for name, since the name could be 1 name or two names
	{
	  //FILE *ff;
	  ifstream mff;
	  double tmpper,tmpvel,tmpunc;
	  int i,size,flag=0;
	  string line,tmpnm;
 	  char ll[150];
          vector<string> v;
	  vector<double> vpper,vgper,vp(100);
	  vector<double>::iterator it;
	  if(disp.npper>0 or disp.ngper>0)
		{
		  cout<< "#########disp data existed already!\n exit!"<<endl;
		  exit(0);
		}
	  if(disp.fgroup>0 or disp.fphase>0)
		{ cout<< "##########disp data existed already!\n exit!"<<endl;
                  exit(0);}

	  disp.pper.clear(); disp.pvelo.clear();disp.unpvelo.clear();
	  if(surflag==1 or surflag==3) //surflag==1: open phase only. surfalg ==3 open phase and group, surflag==2: open group only
	  {
	    //cout<<"read phase disp!\n";
	    i=0;
	    mff.open(name[0].c_str());
            if ( not mff.is_open()){cout<<"###########disp phase file"<<name[0]<<"does not exist!\n exit!!\n";exit(0);}
            while(getline(mff,line))
                {
		  v.clear();
                  Split(line,v," ");
		  size=v.size();
		  
	         if(size>2)
		 {
		 disp.pper.push_back(atof(v[0].c_str()));
		 disp.pvelo.push_back(atof(v[1].c_str()));
		 disp.unpvelo.push_back(atof(v[2].c_str()));
		 }
		 else if(size >1)
		 {
		 disp.pper.push_back(atof(v[0].c_str()));
                 disp.pvelo.push_back(atof(v[1].c_str()));
		 disp.unpvelo.push_back(0.02);
		 }
		 else if(size>0)
		 {
		 disp.pper.push_back(atof(v[0].c_str()));
		 disp.pvelo.push_back(0.);
		 disp.unpvelo.push_back(0.02);
		 }
		 else
		 {
		 disp.pper.push_back(0.);
		 disp.pvelo.push_back(0.);
                 disp.unpvelo.push_back(0.02);
		 }
		  i=i+1;
		}//while
	    disp.npper=i;
	    disp.fphase=1;
	    mff.close();
	  }// 1 and 3
	  //============= read group ================
	  i=0;
	  disp.gper.clear(); disp.gvelo.clear();disp.ungvelo.clear();
	  if(surflag==2){tmpnm=name[0];flag=1;}
	  else if(surflag==3){tmpnm=name[1];flag=1;}
	  if(flag)
	  {
 	   mff.open(tmpnm.c_str());
	   //cout<<"read group diisp file! \n";
	   if ( not mff.is_open()){cout<<"####disp group file"<<tmpnm<<"does not exist!\n exit!!\n";exit(0);}
	   while(getline(mff,line))
		{
		// cout<<line<<endl;
		 v.clear();
		 Split(line,v," ");
		 size=v.size();
	         if(size>2)
		 {
		 disp.gper.push_back(atof(v[0].c_str()));
		 disp.gvelo.push_back(atof(v[1].c_str()));
		 disp.ungvelo.push_back(atof(v[2].c_str()));
		 }
		 else if(size >1)
		 {
		 disp.gper.push_back(atof(v[0].c_str()));
                 disp.gvelo.push_back(atof(v[1].c_str()));
		 disp.ungvelo.push_back(0.04);
		 }
		 else if(size>0)
		 {
		 disp.gper.push_back(atof(v[0].c_str()));
		 disp.gvelo.push_back(0.);
		 disp.ungvelo.push_back(0.04);
		 }
		 else
		 {
		 disp.gper.push_back(0.);
		 disp.gvelo.push_back(0.);
                 disp.ungvelo.push_back(0.04);
		 }
                i=i+1;
		}//while
	    disp.ngper=i;
	    disp.fgroup=1;
	    mff.close();
	  }//if flag
	  vpper=disp.pper;
	  vgper=disp.gper;
	  sort(vpper.begin(),vpper.end());
	  sort(vgper.begin(),vgper.end());
	  it=set_union(vpper.begin(),vpper.end(),vgper.begin(),vgper.end(),vp.begin());
	  size = int(it-vp.begin());
	  vp.resize(size);
	  disp.period1=vp;
	  return 1;
	}//readdisp_single
//----------------------------------------------------
	int readdisp(modeldef &model,vector<string> Rname,vector<string> Lname,int Rsurflag,int Lsurflag)
	{
	  //surflag 0--no  1--phase only 2--group only  3--both group&phase
	  if (Rsurflag>0)
		readdisp_single(model.data.Rdisp,Rname,Rsurflag);
	  if (Lsurflag>0)
		readdisp_single(model.data.Ldisp,Lname,Lsurflag);
	  return 1;
	}//readdisp
//-----------------------------------------------------	
	int readrf(modeldef &model,char *name)
	{
	  //FILE *ff;
	  fstream mff;
	  double tmpto,tmprfo,tmpunrfo;
	  int i=0,size;
	  string line;
	  vector<string> v;
	  if(model.data.rf.nrfo>0)
		{
		  cout<<"#########data existed in rf already!! exit!\n";
		  exit(0);
		}
	  mff.open(name);
	  if ( not mff.is_open()){cout<<"#########rf file"<<name<<"does not exist!\n exit!!\n";exit(0);}
	  model.data.rf.to.clear();model.data.rf.rfo.clear(); model.data.rf.unrfo.clear();
	while(getline(mff,line))
	{
	    v.clear();
	    Split(line,v," ");
	    size=v.size();
	    if(size>2)
	    {
		model.data.rf.to.push_back(atof(v[0].c_str()));
		model.data.rf.rfo.push_back(atof(v[1].c_str()));
		model.data.rf.unrfo.push_back(atof(v[2].c_str()));
	    }
	    else if(size>1)
	    {
		model.data.rf.rfo.push_back(atof(v[1].c_str()));
		model.data.rf.to.push_back(atof(v[0].c_str()));
		model.data.rf.unrfo.push_back(0.1);
	    }
	    else if(size>0)
	    {
		model.data.rf.rfo.push_back(atof(v[1].c_str()));
                model.data.rf.to.push_back(0.);
                model.data.rf.unrfo.push_back(0.1);
	    }
	   else
	   {
		model.data.rf.rfo.push_back(0.);
                model.data.rf.to.push_back(0.);
                model.data.rf.unrfo.push_back(0.1);
 	   }
	    i=i+1;	    
	  }//while
	  model.data.rf.rt=int(1./(model.data.rf.to[1]-model.data.rf.to[0]));
	  model.data.rf.nrfo=i;
	  mff.close();
	  return 1;
	}//readrf

//-----------------------------------------------------	
/*
fill model.groups[]

Input:
column:	0		1			2			3					....			N-1
	iid		groups[iid].flag	groups[iid].thick	np=groups[iid].value.size()		groups[iid].value[]	groups[iid].vpvs
			1-layered=>N=4+2np+1									(and groups[iid].ratio[]
			2-Bspline=>N=4+np+1									 if flag==1)
			4-gradient-------------------------------------> np=2					Value could be either V
			5-water layer ---------------------------------> np=1					or Bspline coefficient
						==
						sum of col2=> model.tthick
*/
	int readmodIso(modeldef &model,const char *name)
	{//Love and Rayleigh have the same input model;
	  string line;
	  int iid,size,i=0,tnp;
	  ifstream mff;
	  vector<string> v;
	  mff.open(name);
	  while(getline(mff,line))
		{
		 groupdef gp;
		 model.groups.push_back(gp);
	   	 model.groups[i].np=0;
		 model.groups[i].flag=-1;
		 model.groups[i].thick=0.;
		 model.groups[i].flagBs=-1;
	         model.groups[i].nlay=20;
		 model.groups[i].vpvs=1.75;
		 i++;
		}
	  model.ngroup=i;
	  mff.close(); //ANY OTHER WAY FOR GOING TO THE BEGINNING OF FILE?
	  mff.open(name);
	  if ( not mff.is_open()){cout<<"########model file"<<name<<"does not exist!\n exit!!\n";exit(0);}
	  while(getline(mff,line))
		{
		  v.clear();
		  Split(line,v," ");
		  size=v.size();
		  iid=atoi(v[0].c_str());
		  model.groups[iid].flag=atoi(v[1].c_str());//type of model: 1--layered 2--Bspline 4--gradient
		  model.groups[iid].thick=atof(v[2].c_str());//h of each group
		  tnp=atoi(v[3].c_str());
		  model.groups[iid].np=tnp;//#of parameters
		  //=========check the column # of input file====
		  if(model.groups[iid].flag == 4 and tnp!=2) //grad
			{cout<<"##########for gradient model, ONLY 2 values!\n";exit(0);}
		  else if ( (model.groups[iid].flag==1 and size != 4+2*tnp+1 ) or (model.groups[iid].flag==2 and size!=4+tnp+1))//lay or Bs
			{
			  cout<<"#########wrong input for model: "<<line<< "\t\t size: "<<size<<endl;
			  printf("model.group[%d].flag=%d\n",iid,model.groups[iid].flag);
			  //for(int k=0;k<size;k++)cout<<"//  "<<v[k]<<" "<<k;
			  exit(0);
			}
		  else if(model.groups[iid].flag==5 and tnp !=1)//water layer , a bug here, fixed on Dec 30,2011
			{cout<<"########for water model, ONLY 1 value!\n";exit(0);}
		  model.groups[iid].Rvalue.clear();model.groups[iid].Lvalue.clear();model.groups[iid].ratio.clear();

		  for(i=0;i< tnp;i++)
			{
			  model.groups[iid].Rvalue.push_back(atof(v[4+i].c_str()));
			  if( model.groups[iid].flag==1)//layered
				  {model.groups[iid].ratio.push_back(atof(v[4+i+tnp].c_str()));}
			}//for
		  model.groups[iid].Lvalue=model.groups[iid].Rvalue;//Love and Rayleigh have the same input model;
		  model.groups[iid].vpvs=(atof(v[size-1].c_str()));
		  model.tthick = model.tthick + model.groups[iid].thick;
		 // cout<<line<<endl;
		} //while getline
	  mff.close();
	  //if(model.ngroup>1)
	  //	{cout<<"flag0: "<<model.groups[0].flag<<"   flag1:"<<model.groups[1].flag<<endl;}
	  return 1;
	}//readmod


//----------------------------------------------------- 
	int readmodAniso(modeldef &model,const char *name){
	  string line;
	  int iid,size,i=0,tnp;
	  ifstream mff;
	  vector<string> v;
	  mff.open(name);
	  model.tthick=0.;// a bug here, fixed on Mar 29, 2012. (without this sentence, if the input model isn't initiated, then the tthick would be wrong.)
	  while(getline(mff,line))
		{
		 groupdef gp;
		 model.groups.push_back(gp);
	   	 model.groups[i].np=0;
		 model.groups[i].flag=-1;
		 model.groups[i].thick=0.;
		 model.groups[i].flagBs=-1;
	         model.groups[i].nlay=20;
		 model.groups[i].vpvs=1.75;
		 i++;
		}
	  model.ngroup=i;
	  mff.close(); //ANY OTHER WAY FOR GOING TO THE BEGINNING OF FILE?
	  mff.open(name);
	  if ( not mff.is_open()){cout<<"########model file"<<name<<"does not exist!\n exit!!\n";exit(0);}
	  while(getline(mff,line))
		{
		  v.clear();
		  Split(line,v," ");
		  size=v.size();
		  iid=atoi(v[0].c_str());
		  model.groups[iid].flag=atoi(v[1].c_str());//type of model: 1--layered 2--Bspline 4--gradient
		  model.groups[iid].thick=atof(v[2].c_str());//h of each group
		  tnp=atoi(v[3].c_str());
		  model.groups[iid].np=tnp;//#of parameters
		  //=========check the column # of input file====
		  if(model.groups[iid].flag == 4 and tnp!=2) //grad
			{cout<<"##########for gradient model, ONLY 2 values!\n";exit(0);}
		  else if ( (model.groups[iid].flag==1 and size != 4+3*tnp+1 ) or (model.groups[iid].flag==2 and size!=4+tnp*2+1))//lay or Bs
			{
			  cout<<"#########wrong input for model: "<<line<< "\t\t size: "<<size<<endl;
			  printf("model.group[%d].flag=%d\n",iid,model.groups[iid].flag);
			  //for(int k=0;k<size;k++)cout<<"//  "<<v[k]<<" "<<k;
			  exit(0);
			}
		  else if(model.groups[iid].flag==5 and tnp !=1)//water layer , a bug here, fixed on Dec 30,2011
			{cout<<"########for water model, ONLY 1 value!\n";exit(0);}
		  model.groups[iid].Rvalue.clear();model.groups[iid].Lvalue.clear();model.groups[iid].ratio.clear();

		  for(i=0;i< tnp;i++)
			{
			  model.groups[iid].Rvalue.push_back(atof(v[4+i*2].c_str()));
			  model.groups[iid].Lvalue.push_back(atof(v[4+i*2+1].c_str()));
			  if( model.groups[iid].flag==1)//layered
				  {model.groups[iid].ratio.push_back(atof(v[4+i+tnp*2].c_str()));}
			}//for
		  model.groups[iid].vpvs=(atof(v[size-1].c_str()));
		  model.tthick = model.tthick + model.groups[iid].thick;
		 // cout<<line<<endl;
		} //while getline
	  mff.close();
	  //if(model.ngroup>1)
	  //	{cout<<"flag0: "<<model.groups[0].flag<<"   flag1:"<<model.groups[1].flag<<endl;}
	  return 1;
	}//readmodAniso

//----------------------------------------------------- 
        int updatemodelTibet(modeldef &model,float depcri1,float depcri2,float qpcri,float qscri)
        {
          int i,j,tnlay=0;
	  double tRvalue,tLvalue,tAvalue,tthick,tvpvs,tdep=0.,tvsv,tvsh,tvp,tvs,trho,tqs,tqp;
	  model.laym0.vsv.clear();model.laym0.vsh.clear();model.laym0.vpvs.clear();model.laym0.vp.clear();model.laym0.rho.clear();model.laym0.qs.clear();model.laym0.qp.clear();model.laym0.thick.clear();
          for(i=0;i<model.ngroup;i++)
                {
		 updategroup(model.groups[i]);//based on g.Rv/Lv/t fill g.Rv1/Lv1/t1/Av1. Since Rv/Lv always exist regardless of flagLR, both Rv1/Lv1 are filled based on Rv/Lv.
		 
		 for(j=0;j<model.groups[i].nlay;j++)
			{
			  tRvalue=model.groups[i].Rvalue1[j];
			  tLvalue=model.groups[i].Lvalue1[j];
			  tAvalue=model.groups[i].Avalue1[j];
			  tthick=model.groups[i].thick1[j];
			  tdep=tdep+tthick;
			  if(model.groups[i].flag==5)//water layer
			  {
				tvpvs=-1.;
				tvsh=0.;
			        tvsv=0.;
				tvp=tRvalue;
				trho=1.02;
				tqs=10000.;
				tqp=57822.;
			  }
			  else if((i==0 and model.groups[0].flag!=5) or ( i==1 and model.groups[0].flag==5)) //layer1 sediment OR layer1 water, layer2 sediment
			  {
				tvpvs=model.groups[i].vpvs;
				tvsv=tRvalue;
				tvsh=tLvalue;
				tvs=tAvalue;
				tvp=tvpvs*tvs;
				tqp=160.; tqs=80.; trho=0.541+0.3601*tvp;
			  }
			  else//crust and mantle
			  {
				tvpvs=model.groups[i].vpvs;
				tvsv=tRvalue;tvsh=tLvalue;tvs=tAvalue;
				tvp=tvpvs*tvs;
				if(tdep<18.) // changed from 10 to 18 on Aug 12, 2012
					{tqp=1400.; tqs=600.;}//AK135
				else if(tdep<depcri1)
					{tqp=900.;tqs=400.;} // there was a bug here, fixed on Dec. 27,2011

				else if (tdep<depcri2 ) //depcri1<dep<depcri2 ==> the radial anisotropy zone
					{tqp=qpcri;tqs=qscri;}
				else if (tdep<80)// there was a bug here, added on Aug 9, 2012
                                        {tqp=900.;tqs=400.;}

				else 
					{tqp=200.;tqs=80.;} // there was a bug here, fixed on Dec. 27,2011
				if(tvp<7.5){trho=0.541+0.3601*tvp;} 
				else{trho=3.35;} //# Kaban, M. K et al. (2003), Density of the continental roots: Compositional and thermal contributions
			  }//else

			  model.laym0.vpvs.push_back(tvpvs);
			  model.laym0.vsv.push_back(tvsv);
			  model.laym0.vsh.push_back(tvsh);
			  model.laym0.vp.push_back(tvp);
			  model.laym0.rho.push_back(trho);
			  model.laym0.qs.push_back(tqs);
			  model.laym0.qp.push_back(tqp);
			  model.laym0.thick.push_back(tthick);
			}//for j
		  tnlay=tnlay+model.groups[i].nlay;
                }//for i

	  model.laym0.nlayer=tnlay;
	  model.flag=1;//updated, from layered vs,h get other layered para
	  return 1;
        }//updatemodel Tibet


//----------------------------------------------------- 
        int updatemodelRF(modeldef &model,float depcri1,float depcri2,float qpcri,float qscri)
        { // this is very similar to the updatemodelTibet; Only that here, the updategrouprf layerize the model in a different way (layer other than smooth model, so the last layer doesn't have a 0km thickness)
          int i,j,tnlay=0;
	  double tRvalue,tLvalue,tAvalue,tthick,tvpvs,tdep=0.,tvsv,tvsh,tvp,tvs,trho,tqs,tqp;
	  model.laym0.vsv.clear();model.laym0.vsh.clear();model.laym0.vpvs.clear();model.laym0.vp.clear();model.laym0.rho.clear();model.laym0.qs.clear();model.laym0.qp.clear();model.laym0.thick.clear();
          for(i=0;i<model.ngroup;i++)
                {
		 updategrouprf(model.groups[i]);//based on g.Rv/Lv/t fill g.Rv1/Lv1/t1/Av1. Since Rv/Lv always exist regardless of flagLR, both Rv1/Lv1 are filled based on Rv/Lv.
		 
		 for(j=0;j<model.groups[i].nlay;j++)
			{
			  tRvalue=model.groups[i].Rvalue1[j];
			  tLvalue=model.groups[i].Lvalue1[j];
			  tAvalue=model.groups[i].Avalue1[j];
			  tthick=model.groups[i].thick1[j];
			  tdep=tdep+tthick;
			  if(model.groups[i].flag==5)//water layer
			  {
				tvpvs=-1.;
				tvsh=0.;
			        tvsv=0.;
				tvp=tRvalue;
				trho=1.02;
				tqs=10000.;
				tqp=57822.;
			  }
			  else if((i==0 and model.groups[0].flag!=5) or ( i==1 and model.groups[0].flag==5)) //layer1 sediment OR layer1 water, layer2 sediment
			  {
				tvpvs=model.groups[i].vpvs;
				tvsv=tRvalue;
				tvsh=tLvalue;
				tvs=tAvalue;
				tvp=tvpvs*tvs;
				tqp=160.; tqs=80.; trho=0.541+0.3601*tvp;
			  }
			  else//crust and mantle
			  {
				tvpvs=model.groups[i].vpvs;
				tvsv=tRvalue;tvsh=tLvalue;tvs=tAvalue;
				tvp=tvpvs*tvs;
				if(tdep<18.) // changed from 10 to 18 on Aug 12, 2012
					{tqp=1400.; tqs=600.;}//AK135
				else if(tdep<depcri1)
					{tqp=900.;tqs=400.;} // there was a bug here, fixed on Dec. 27,2011

				else if (tdep<depcri2 ) //depcri1<dep<depcri2 ==> the radial anisotropy zone
					{tqp=qpcri;tqs=qscri;}
				else if (tdep<80)// there was a bug here, added on Aug 9, 2012
                                        {tqp=900.;tqs=400.;}

				else 
					{tqp=200.;tqs=80.;} // there was a bug here, fixed on Dec. 27,2011
				if(tvp<7.5){trho=0.541+0.3601*tvp;} 
				else{trho=3.35;} //# Kaban, M. K et al. (2003), Density of the continental roots: Compositional and thermal contributions
			  }//else

			  model.laym0.vpvs.push_back(tvpvs);
			  model.laym0.vsv.push_back(tvsv);
			  model.laym0.vsh.push_back(tvsh);
			  model.laym0.vp.push_back(tvp);
			  model.laym0.rho.push_back(trho);
			  model.laym0.qs.push_back(tqs);
			  model.laym0.qp.push_back(tqp);
			  model.laym0.thick.push_back(tthick);
			}//for j
		  tnlay=tnlay+model.groups[i].nlay;
                }//for i

	  model.laym0.nlayer=tnlay;
	  model.flag=1;//updated, from layered vs,h get other layered para
	  return 1;
        }//updatemodel RF

//-----------------------------------------------------	
	int compute_misfitDISP_single(dispdef &disp,double &tempv1,double &tempv2,double &tS)
	{
          int i,j,k;
          double tp,td,p;

          for(i=0;i<disp.npper;i++)
                {tempv1=tempv1+pow((disp.pvelo[i]-disp.pvel[i])/disp.unpvelo[i],2);}
          if(tempv1>0.)
                {
                  disp.pmisfit=sqrt(tempv1/disp.npper);
                  tS=tempv1;
                  if(tS>50.)tS=sqrt(tS*50);
                  if(tS>50.)tS=sqrt(tS*50);
                  disp.pL=exp(-0.5*tS);
                }//if

          for(i=0;i<disp.ngper;i++)
                {tempv2=tempv2+pow((disp.gvelo[i]-disp.gvel[i])/disp.ungvelo[i],2);}
          if(tempv2>0.)
                {
                  disp.gmisfit=sqrt(tempv2/disp.ngper);
                  tS=tempv2;
                  if(tS>50)tS=sqrt(tS*50);
                  if(tS>50)tS=sqrt(tS*50);
                  disp.gL=exp(-0.5*tS);
                }//if	
	}//compute_misfitDISP_single
//-----------------------------------------------------	
//	int compute_misfitDISP(modeldef &model,int Rflag, int Lflag)
	int compute_misfit(modeldef &model,int Rflag, int Lflag,float inp)
	{
	  double tmisfit1,tL1,tS1,tS2,tS;
	  double tempv1=0.,tempv2=0.,tempv3=0.,tempv4=0.,tempv5=0.;
	  int Rnp,Rng,Lnp,Lng,RFn;
	  Rnp=Rng=Lnp=Lng=0;
	  //printf("nT: nRp=%d nRg=%d nLp=%g nLg=%d\n",model.data.Rdisp.npper,model.data.Rdisp.ngper,model.data.Ldisp.npper,model.data.Ldisp.ngper);//---test
	  //--------------------
	  if (Rflag>0){compute_misfitDISP_single(model.data.Rdisp,tempv1,tempv2,tS1);Rnp=model.data.Rdisp.npper;Rng=model.data.Rdisp.ngper;}
	  tmisfit1=sqrt((tempv1+tempv2)/(model.data.Rdisp.npper+model.data.Rdisp.ngper));
	  //tL1=model.data.Rdisp.gL*model.data.Rdisp.pL;
	  tS=(tempv1+tempv2);
	  if(tS>50)tS=sqrt(tS*50);
	  if(tS>50)tS=sqrt(tS*50);
	  tL1=exp(-0.5*tS);

	  model.data.Rdisp.L=tL1;
	  model.data.Rdisp.misfit=tmisfit1;

	  //--------------------
	  if (Lflag>0){compute_misfitDISP_single(model.data.Ldisp,tempv3,tempv4,tS2);Lnp=model.data.Ldisp.npper;Lng=model.data.Ldisp.ngper;}
	  tmisfit1=sqrt((tempv3+tempv4)/(model.data.Ldisp.npper+model.data.Ldisp.ngper));
	  //tL1=model.data.Ldisp.gL*model.data.Ldisp.pL;
	  tS=(tempv3+tempv4);
	  if(tS>50)tS=sqrt(tS*50);
	  if(tS>50)tS=sqrt(tS*50);
	  tL1=exp(-0.5*tS);

	  model.data.Ldisp.L=tL1;
	  model.data.Ldisp.misfit=tmisfit1;

	  //-------------------------
	  int k=0;
	  if(model.data.rf.nrfo>0){
	  //RFn=model.data.rf.nn;
 	  for (int i=0;i<model.data.rf.nrfo;i++){
	     for (int j =0;j<model.data.rf.tn.size();j++){
		if(pow((model.data.rf.to[i]-model.data.rf.tn[j]),2)<0.0001 and model.data.rf.to[i]<=10. and model.data.rf.to[i]>=0.) //In case the obs and pred RF have different time interval. Use first 10 sec only
		{
			tempv5=tempv5+pow((model.data.rf.rfo[i]-model.data.rf.rfn[j])/model.data.rf.unrfo[i],2);
			k=k+1;
		}
	     }// for j
	  }//for i
	  //tS=tempv5/RFn;
	  tS=tempv5;
	  if (tS>50)tS=sqrt(tS*50);
	  if (tS>50)tS=sqrt(tS*50);
	  tmisfit1=sqrt(tempv5/k);
	  tL1=exp(-0.5*tS);

	  model.data.rf.L=tL1;
  	  model.data.rf.misfit=tmisfit1;
	  }// if nrfo>0
	  //-------------------------------------------------------

	  //*************how to compute data.L data.misfit when R and L are both presented?
	  tS=0.5*(inp*(tempv1+tempv2+tempv3+tempv4)+(1-inp)*tempv5);
	  if (tS>50.)tS=sqrt(tS*50);
	  if(tS>50.)tS=sqrt(tS*50);
	  model.data.L=exp(-0.5*tS);
	  //model.data.L=model.data.Rdisp.L*model.data.Ldisp.L;
	  //model.data.misfit=sqrt((tempv1+tempv2+tempv3+tempv4)/(model.data.Rdisp.npper+model.data.Rdisp.ngper+model.data.Ldisp.npper+model.data.Ldisp.ngper));	
	  model.data.misfit=inp*sqrt((tempv1+tempv2+tempv3+tempv4)/(Rnp+Rng+Lnp+Lng))+(1-inp)*tmisfit1;	
 	  //*****************
	  return 1;
	}//compute misfit
//-----------------------------------------------------	

//-----------------------------------------------------	
	int goodmodel( modeldef &model, vector<int> vmono, vector<int> vgrad,int Rflag, int Lflag, int isoflag)
	 {

	  // for anisotropic case, cannot check goodmodel for both R and L. need to run this twice, check L R seperately
	  int i,j,var=0;
	  vector<int>::iterator id;
	  double gradient;
	  if(isoflag>0 or Rflag>0){ // iso case, R and L have the same model.g.value1
	      for(i=0;i<model.ngroup-1;i++)//between layers
		{if(model.groups[i+1].Rvalue1[0]<model.groups[i].Rvalue1.back()) 
			{//cout<<"case 1\n"; //---test----
       			return 0;}
		}  	 
	      for(id=vmono.begin();id<vmono.end();id++)// monotonic change in group vmono[?]
              {
	    	j=*id;
	    	for(i=0;i<model.groups[j].nlay-1;i++)
		    { gradient=(model.groups[j].thick1[i])/(model.groups[j].Rvalue1[i]-model.groups[j].Rvalue1[i+1]);
		      /* this is for LVZ, require LVZ not too slow
		       * if(gradient>0. and gradient <70.)
		      {//printf("V1=%g h1=%g V2=%g h2=%g\n",model.groups[j].Rvalue1[i],model.groups[j].Rvalue1[i+1],model.groups[j].thick1[i],model.groups[j].thick1[i+1]);
		       return 0;}
		       */
			  
		      if(gradient>0){return 0;}	    
		    }//for i
	     }//for id
	      //--------------newly added, Jun 2, 2012--------
	      //require the slope in mantle part larger than 70.
	     for(i=0;i<model.groups[2].nlay-1;i++){
	         gradient=(model.groups[2].thick1[i])/(model.groups[2].Rvalue1[i]-model.groups[2].Rvalue1[i+1]);
		 if(gradient>0. and gradient <50.){
			 //printf("mantle grad !\n");
			 return 0;}
	     }
	     // require the vel at all depth (the maximum = g.tthcik, which comes from the input) to be smaller than 4.9 ---- revised on Aug 23, 2012------
	     for (i=0;i<model.groups[2].nlay;i++){
	     if(model.groups[2].Rvalue1[i]>4.9 or model.groups[2].Lvalue1[i]>4.9){
		     //printf("too large Vmoho\n");
		     return 0;}
	     }
	     //------------------------------------------------
	     for(id=vgrad.begin();id<vgrad.end();id++)//gradient check for the 1st two velue in group vgrad[?]
	      {
	        j=*id;
	        if(model.groups[j].Rvalue1[1]<model.groups[j].Rvalue1[0])
		{//cout<<"case 3\n"; //---test----
		    return 0;}
	      } 
	
	  }//if isoflag or Rflag>0

	  
	  else if (Lflag>0){
              for(i=0;i<model.ngroup-1;i++)//between layers
                {if(model.groups[i+1].Lvalue1[0]<model.groups[i].Lvalue1.back())
                 return 0;} 
              for(id=vmono.begin();id<vmono.end();id++)// monotonic change in group vmono[?]
              {
                j=*id;
                for(i=0;i<model.groups[j].nlay-1;i++)
		    { gradient=(model.groups[j].thick1[i])/(model.groups[j].Lvalue1[i]-model.groups[j].Lvalue1[i+1]);
		      if(gradient>0. and gradient <70.)return 0; }
	     }//for id
              for(id=vgrad.begin();id<vgrad.end();id++)//gradient check for the 1st two velue in group vgrad[?]
              {
                j=*id;
                if(model.groups[j].Lvalue1[1]<model.groups[j].Lvalue1[0])
                    return 0;
              }			
	  
	  }//if Lflag
	
	  return 1;
	}//goodmodel

//-----------------------------------------------------	
	int positiveAni( modeldef model, vector<int> vmono)
	{// require the L to be always faster than R in group vmono[?]
	  vector<int>::iterator id;
	  int i,j;
	  for(id=vmono.begin();id<vmono.end();id++){
		j=*id;
		for(i=0;i<model.groups[j].nlay;i++){
		  if(model.groups[j].Lvalue1[i]<model.groups[j].Rvalue1[i])
			{//printf("--@@@ng=%d nv=%d Lv1=%g Rv1=%g\n",j,i,model.groups[j].Lvalue1[i],model.groups[j].Rvalue1[i]);
			return 0;}
		}//for i
	  }//for id
	  return 1;
	}//positiveAni
	
