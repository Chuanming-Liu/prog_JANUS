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
// this version, on Apr 22, 2014, added constraint on vp/vs 

using namespace std;
/*=====CONTENT=======
int initdisp (dispdef &disp)
int initmodel(modeldef &model)
int readdisp_single(dispdef &disp,vector<string> name, int surflag)
int readdisp(modeldef &model,vector<string> Rname,vector<string> Lname,int Rsurflag,int Lsurflag)
***surflag 0no/1phase/2group/3ph&gp
int readrf(modeldef &model,char *name)
*int readmodIso(modeldef &model,const char *name)
int readmodAniso(modeldef &model,const char *name)
int updatemodelTibet(modeldef &model,float depcri1,float depcri2,float qpcri,float qscri)
int compute_misfitDISP_single(dispdef &disp,double &tempv1,double &tempv2,double &tS)
int compute_misfitDISP(modeldef &model,int Rflag, int Lflag, int AziampRflag, int AziampLflag, int AziphiRflag, int AziphiLflag, float inpamp, float inpphi )
int goodmodel( modeldef &model, vector<int> vmono, vector<int> vgrad,int Rflag, int Lflag, int isoflag)
int positiveAni( modeldef model, vector<int> vmono)
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
	  model.laym0.vsv.clear();
	  model.laym0.vsh.clear();
	  model.laym0.vpv.clear();
	  model.laym0.vph.clear();
	  model.laym0.eta.clear();
	  model.laym0.theta.clear();
	  model.laym0.phi.clear();
	  model.laym0.vpvs.clear();
	  model.laym0.rho.clear();
	  model.laym0.qs.clear();
	  model.laym0.qp.clear();
	  model.laym0.thick.clear();

	  model.data.Rdisp.pper.clear();
	  model.data.Rdisp.pvelo.clear();
	  model.data.Rdisp.pvel.clear();
	  model.data.Rdisp.unpvelo.clear();
	  model.data.Rdisp.gper.clear();
	  model.data.Rdisp.gvelo.clear();
	  model.data.Rdisp.gvel.clear();
	  model.data.Rdisp.ungvelo.clear();
	  model.data.Rdisp.period1.clear();

	  model.data.Ldisp.pper.clear();
	  model.data.Ldisp.pvelo.clear();
	  model.data.Ldisp.pvel.clear();
	  model.data.Ldisp.unpvelo.clear();
	  model.data.Ldisp.gper.clear();
	  model.data.Ldisp.gvelo.clear();
	  model.data.Ldisp.gvel.clear();
	  model.data.Ldisp.ungvelo.clear();
	  model.data.Ldisp.period1.clear();
	  
	  model.data.AziampRdisp.pper.clear();
	  model.data.AziampRdisp.pvelo.clear();
	  model.data.AziampRdisp.pvel.clear();
	  model.data.AziampRdisp.unpvelo.clear();
	  model.data.AziampRdisp.gper.clear();
	  model.data.AziampRdisp.gvelo.clear();
	  model.data.AziampRdisp.gvel.clear();
	  model.data.AziampRdisp.ungvelo.clear();
	  model.data.AziampRdisp.period1.clear();
	  
	  model.data.AziampLdisp.pper.clear();
	  model.data.AziampLdisp.pvelo.clear();
	  model.data.AziampLdisp.pvel.clear();
	  model.data.AziampLdisp.unpvelo.clear();
	  model.data.AziampLdisp.gper.clear();
	  model.data.AziampLdisp.gvelo.clear();
	  model.data.AziampLdisp.gvel.clear();
	  model.data.AziampLdisp.ungvelo.clear();
	  model.data.AziampLdisp.period1.clear();
	  
	  model.data.AziphiRdisp.pper.clear();
	  model.data.AziphiRdisp.pvelo.clear();
	  model.data.AziphiRdisp.pvel.clear();
	  model.data.AziphiRdisp.unpvelo.clear();
	  model.data.AziphiRdisp.gper.clear();
	  model.data.AziphiRdisp.gvelo.clear();
	  model.data.AziphiRdisp.gvel.clear();
	  model.data.AziphiRdisp.ungvelo.clear();
	  model.data.AziphiRdisp.period1.clear();
	  
	  model.data.AziphiLdisp.pper.clear();
	  model.data.AziphiLdisp.pvelo.clear();
	  model.data.AziphiLdisp.pvel.clear();
	  model.data.AziphiLdisp.unpvelo.clear();
	  model.data.AziphiLdisp.gper.clear();
	  model.data.AziphiLdisp.gvelo.clear();
	  model.data.AziphiLdisp.gvel.clear();
	  model.data.AziphiLdisp.ungvelo.clear();
	  model.data.AziphiLdisp.period1.clear();
	  
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
	int readdisp(modeldef &model,vector<string> Rname,vector<string> Lname,vector<string> AziampRname,vector<string> AziphiRname,vector<string> AziampLname,vector<string> AziphiLname,int Rsurflag,int Lsurflag, int AziampRsurflag, int AziphiRsurflag,int AziampLsurflag,int AziphiLsurflag)
	{
	  //surflag 0--no  1--phase only 2--group only  3--both group&phase
	  if (Rsurflag>0)
		readdisp_single(model.data.Rdisp,Rname,Rsurflag);
	  if (Lsurflag>0)
		readdisp_single(model.data.Ldisp,Lname,Lsurflag);
	  if(AziampRsurflag>0)
		readdisp_single(model.data.AziampRdisp,AziampRname,AziampRsurflag); 
	  if(AziphiRsurflag>0)
		readdisp_single(model.data.AziphiRdisp,AziphiRname,AziphiRsurflag); 
	  if(AziampLsurflag>0)
		readdisp_single(model.data.AziampLdisp,AziampLname,AziampLsurflag); 
	  if(AziphiLsurflag>0)
		readdisp_single(model.data.AziphiLdisp,AziphiLname,AziphiLsurflag); 
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
column:	0		1			2				3			4				5				....			N-1
	iid		groups[iid].flag	groups[iid].pflag	groups[iid].flagcpttype		groups[iid].thick	np=groups[iid].value.size()		groups[iid].value[]	groups[iid].vpvs
			1-layered=>N=nn+tnp*(pflag+1)+1	check the readmodAniso for detail								(and groups[iid].ratio[]
			2-Bspline=>N=nn+tnp*pflag+1													 if flag==1)
			4-gradient-------------------------------------> np=2											 Value could be either V
			5-water layer ---------------------------------> np=1											 or Bspline coefficient
						==
						sum of col2=> model.tthick
*/
// ******** a colum flagcpttype is added between pflag and thick; it indicates the type of forward computation that will be used for this group; 1--use Vkernel to do all cpt (model is TI or iso) 2--use Vkernel(for RA) and Lovekernel(for AZ) 3--use Vkernel(for RA) and Azikernel(for AZ) 4--use Lovekernel (for both RA and AZ)
// // Vkernel-- dC/dV, dC/dh;
// // Lovekernel-- dC/dX, (dC/dh); ==> right now, this can be used for layerized model only (flag==1)
// // Azikernel-- dC/dVcos, dC/dVsin;
//

//----------------------------------------------------- 
	int readmodAniso(modeldef &model,const char *name){
	  string line;
	  int iid,size,i=0,tnp,pflag,nn;
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
	  nn=6;//# of para before the value-para shows up
	  while(getline(mff,line))
		{
		  v.clear();
		  Split(line,v," ");
		  size=v.size();
		  iid=atoi(v[0].c_str());
		  for(int a=0;a<size;a++){//---test---
			printf("v[%d]=%s  ",a,v[a].c_str());
			}
		  printf("\n");//---test--
		  model.groups[iid].flag=atoi(v[1].c_str());//type of model: 1--layered 2--Bspline 4--gradient 5--water layer
		  pflag=model.groups[iid].pflag=atoi(v[2].c_str());// type of input model: if not specified, vsv=vsh,vpv=vph,eta=1, and vp=vs*vpvs; 1--vs; 2--vsv,vsh; 3--vsv,vsh,vp; 4--vsv,vsh,vpv,vph; 5--vsv,vsh,vpv,vph,eta;
		  model.groups[iid].flagcpttype=atoi(v[3].c_str());
		  model.groups[iid].thick=atof(v[4].c_str());//h of each group
		  if((model.groups[iid].flagcpttype==2 or model.groups[iid].flagcpttype==4) and (model.groups[iid].flag-1)*(model.groups[iid].flag-3)){// Lovepara/kernel, only for layered model or Bspline model that will be changed to point model(flag=3 NOT 2) ; BS
			printf("##### readmod, wrong flagcpttype vs. flag! Lovepara/kernel, only for layered or Bspline (flag=3(Bsp2Point) not flag=2(Bsp)) model\ngroup%d, flagcpttype=%d, flag=%d\n",iid,model.groups[iid].flagcpttype,model.groups[iid].flag);

			exit(0);
		  }

		  tnp=atoi(v[5].c_str());
		  model.groups[iid].np=tnp;//#of parameters
		  
		  //=========check the column # of input file====
		  if(model.groups[iid].flag == 4 and tnp!=2) //grad
			{cout<<"##########for gradient model, ONLY 2 values!\n";exit(0);}
		  else if ( (model.groups[iid].flag==1 and size != nn+(pflag+1)*tnp+1 ) or (model.groups[iid].flag==2 and size!=nn+tnp*pflag+1))//lay or Bs
			{
			  cout<<"#########wrong input for model: "<<line<< "\t\t size: "<<size<<endl;
			  printf("should be %d if flag==1, %d if flag==2, here flag=%d\n",nn+(pflag+1)*tnp+1,nn+tnp*pflag+1,model.groups[iid].flag);
			  printf("model.group[%d].flag=%d pflag=%d nn=%d tnp=%d\n",iid,model.groups[iid].flag,pflag,nn,tnp);
			  //for(int k=0;k<size;k++)cout<<"//  "<<v[k]<<" "<<k;
			  exit(0);
			}
		  else if(model.groups[iid].flag==5 and tnp !=1)//water layer , a bug here, fixed on Dec 30,2011
			{cout<<"########for water model, ONLY 1 value!\n";exit(0);}
		  model.groups[iid].vsvvalue.clear();
		  model.groups[iid].vshvalue.clear();
		  model.groups[iid].vpvvalue.clear();
		  model.groups[iid].vphvalue.clear();
		  model.groups[iid].etavalue.clear();
		  model.groups[iid].thetavalue.clear();
		  model.groups[iid].phivalue.clear();
		  model.groups[iid].rhovalue.clear();
		  model.groups[iid].ratio.clear();
		  model.groups[iid].vpvs=(atof(v[size-1].c_str()));
		  model.tthick = model.tthick + model.groups[iid].thick;
		  for(i=0;i<tnp;i++)
		 	{
			  printf("para%d\n",i);//--test--
				
			  if(pflag>0){ 	  	
			    model.groups[iid].vsvvalue.push_back(atof(v[nn+i*pflag].c_str()));

			    if(model.groups[iid].flag==1){//layered
				    model.groups[iid].ratio.push_back(atof(v[nn+i+tnp*pflag].c_str()));
			    }// if flag==1

			    if(pflag>1){
			      model.groups[iid].vshvalue.push_back(atof(v[nn+i*pflag+1].c_str()));

			      if(pflag>2){
			        model.groups[iid].vpvvalue.push_back(atof(v[nn+i*pflag+2].c_str()));
				
				if(pflag>3){
				  model.groups[iid].vphvalue.push_back(atof(v[nn+i*pflag+3].c_str()));
				
				  if(pflag>4){
				    model.groups[iid].etavalue.push_back(atof(v[nn+i*pflag+4].c_str()));
				      if(pflag>6){
					model.groups[iid].thetavalue.push_back(atof(v[nn+i*pflag+5].c_str()));
					model.groups[iid].phivalue.push_back(atof(v[nn+i*pflag+6].c_str()));				
				      }
				      else{
			    		model.groups[iid].thetavalue.push_back(0.);
			    		model.groups[iid].phivalue.push_back(0.);
				      }

				  }// if pflag>4
				  else{//pflag==4
				    model.groups[iid].etavalue.push_back(1.0);
				    model.groups[iid].thetavalue.push_back(0.);
			    	    model.groups[iid].phivalue.push_back(0.);
				  }// else pflag>4

				}// if pflag>3
				else{//pflag==3
				  model.groups[iid].vphvalue.push_back(model.groups[iid].vpvvalue[i]);
				  model.groups[iid].etavalue.push_back(1.0);
				  model.groups[iid].thetavalue.push_back(0.);
			    	  model.groups[iid].phivalue.push_back(0.);
				}//else pflag>3

			      }// if pflag>2
			      else{//pflag==2
			        model.groups[iid].vpvvalue.push_back(model.groups[iid].vsvvalue[i]*model.groups[iid].vpvs);
			        model.groups[iid].vphvalue.push_back(model.groups[iid].vsvvalue[i]*model.groups[iid].vpvs);
			        model.groups[iid].etavalue.push_back(1.0);
				model.groups[iid].thetavalue.push_back(0.);
			    	model.groups[iid].phivalue.push_back(0.);
			      }// else pflag>2

			    }//if pflag>1	  
			    else{//pflag==1
			      model.groups[iid].vshvalue.push_back(model.groups[iid].vsvvalue[i]);
			      model.groups[iid].vpvvalue.push_back(model.groups[iid].vsvvalue[i]*model.groups[iid].vpvs);
			      model.groups[iid].vphvalue.push_back(model.groups[iid].vsvvalue[i]*model.groups[iid].vpvs);
			      model.groups[iid].etavalue.push_back(1.0);
			      model.groups[iid].thetavalue.push_back(0.);
			      model.groups[iid].phivalue.push_back(0.);
			    }//else pflag>1
			  }//pflag>0

			  if(model.groups[iid].flag==5){//water layer
			  	model.groups[iid].rhovalue.push_back(1.02);}
			  else if((iid==0 and model.groups[0].flag!=5) or (iid==1 and model.groups[0].flag==5)){// layer1 sediment OR layer2 seidment (the layer1 is water)
				 model.groups[iid].rhovalue.push_back(0.541+0.3601*.5*(model.groups[iid].vpvvalue[i]+model.groups[iid].vphvalue[i])); }
			  else{// crust or mantle ########### this part may need modification. in the updatemodel, there is a vp<7.5km/s criteria, but here, since vpv,vph value can be B-spline, criteria is not clear;
			  	model.groups[iid].rhovalue.push_back(0.541+0.3601*.5*(model.groups[iid].vpvvalue[i]+model.groups[iid].vphvalue[i]));}
			  
			 }//for i<np

		 // cout<<line<<endl;
		} //while getline
	  mff.close();
	  //if(model.ngroup>1)
	  //	{cout<<"flag0: "<<model.groups[0].flag<<"   flag1:"<<model.groups[1].flag<<endl;}
	  return 1;
	}//readmodAniso

//----------------------------------------------------- 
        int updatemodel(modeldef &model, int flagupdaterho)
        {
	  //the updatemodel, the vel values in the layer is the value when the ET is flat, i.e., the effect of theta is not taken into account, it's not the vel of the effective TI medium, but the vel of the un-rotated medium
          int i,j,tnlay=0; //flagupdaterho;
	  double tvsvvalue,tvshvalue,tvpvvalue,tvphvalue,tetavalue,tthetavalue,tphivalue,trhovalue;
	  double tthick,tvpvs,tdep=0.,tvsv,tvsh,tvpv,tvph,teta,ttheta,tphi,trho,tqs,tqp;
	  model.laym0.vsv.clear();model.laym0.vsh.clear(); model.laym0.vpv.clear();model.laym0.vph.clear();model.laym0.eta.clear();model.laym0.theta.clear();model.laym0.phi.clear();
	  model.laym0.vpvs.clear();model.laym0.rho.clear();model.laym0.qs.clear();model.laym0.qp.clear();model.laym0.thick.clear();

	  //------------- chose if update rho based on vp or not; if not, the rho would be the rho computed from the initial model (from the readmodAniso step);
	  //flagupdaterho=1;
	  //-------------

	  for(i=0;i<model.ngroup;i++)
                {
		 updategroup(model.groups[i]);//based on g.Rv/Lv/t fill g.Rv1/Lv1/t1/Av1. Since Rv/Lv always exist regardless of flagLR, both Rv1/Lv1 are filled based on Rv/Lv.	
		 for(j=0;j<model.groups[i].nlay;j++)
			{
			  tvsvvalue=model.groups[i].vsvvalue1[j];
			  tvshvalue=model.groups[i].vshvalue1[j];
			  tvpvvalue=model.groups[i].vpvvalue1[j];
			  tvphvalue=model.groups[i].vphvalue1[j];
			  tetavalue=model.groups[i].etavalue1[j];
			  tthetavalue=model.groups[i].thetavalue1[j];
			  tphivalue=model.groups[i].phivalue1[j];
			  trhovalue=model.groups[i].rhovalue1[j];

			  tthick=model.groups[i].thick1[j];
			  tdep=tdep+tthick;
			  if(model.groups[i].flag==5)//water layer
			  {
				tvpvs=-1.;
				tvsh=0.;
			        tvsv=0.;
				tvpv=tvpvvalue;
				tvph=tvphvalue;
				teta=1.0;
				ttheta=0.0;
				tphi=0.0;
				trho=1.02;
				tqs=10000.;
				tqp=57822.;
			  }
			  else if((i==0 and model.groups[0].flag!=5) or ( i==1 and model.groups[0].flag==5)) //layer1 sediment OR layer1 water, layer2 sediment
			  {
				tvpvs=model.groups[i].vpvs;
				tvsv=tvsvvalue;
				tvsh=tvshvalue;
				tvpv=tvpvvalue;
				tvph=tvphvalue;
				teta=tetavalue;
				ttheta=tthetavalue;
				tphi=tphivalue;
				//tvp=tvpvs*tvs;
				tqp=160.; tqs=80.; trho=0.541+0.3601*(tvpv+tvph)/2.;
			  }
			  else//crust and mantle
			  {
				tvpvs=model.groups[i].vpvs;
				tvsv=tvsvvalue;
				tvsh=tvshvalue;
				tvpv=tvpvvalue;
				tvph=tvphvalue;
				teta=tetavalue;
				ttheta=tthetavalue;
				tphi=tphivalue;
				//trho=trhovalue;
				//tvp=tvpvs*tvsv; // modified/changed on Mar 27, 2012; I think the Vp/Vs from RF indicate Vpv/Vsv.
				if(tdep<18.) // changed from 10 to 18 on Aug 12, 2012
					{tqp=1400.; tqs=600.;}//AK135
				else if (tdep<80)// there was a bug here, added on Aug 9, 2012
                                        {tqp=900.;tqs=400.;}
				else 
					{tqp=200.;tqs=80.;} // there was a bug here, fixed on Dec. 27,2011
				if((tvpv+tvph)/2.0<7.5){trho=0.541+0.3601*(tvpv+tvph)/2.0;} 
				else{trho=3.35;} //# Kaban, M. K et al. (2003), Density of the continental roots: Compositional and thermal contributions
			  }//else

			  model.laym0.vpvs.push_back(tvpvs);
			  model.laym0.vsv.push_back(tvsv);
			  model.laym0.vsh.push_back(tvsh);
			  model.laym0.vpv.push_back(tvpv);
			  model.laym0.vph.push_back(tvph);
			  model.laym0.eta.push_back(teta);
			  model.laym0.theta.push_back(ttheta);
			  model.laym0.phi.push_back(tphi);
			  model.laym0.qs.push_back(tqs);
			  model.laym0.qp.push_back(tqp);
			  model.laym0.thick.push_back(tthick);
			  if(flagupdaterho)
			  	model.laym0.rho.push_back(trho);//########### RHO ###########
			  else
				model.laym0.rho.push_back(trhovalue);
			}//for j
		  tnlay=tnlay+model.groups[i].nlay;
                }//for i

	  model.laym0.nlayer=tnlay;
	  model.flag=1;//updated, from layered vs,h get other layered para
	  return 1;
        }//updatemodel 

//-----------------------------------------------------	
	int compute_misfitDISP_single(dispdef &disp,double &tempv,double &tS)
	{
          int i,j,k;
          double tp,td,p,ttS,tempv1=0.,tempv2=0.;

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

	  tempv=tempv1+tempv2;
	  disp.misfit=sqrt((tempv1+tempv2)/(disp.npper+disp.ngper));
	  ttS=(tempv1+tempv2);
	  if(ttS>50)ttS=sqrt(ttS*50);
	  if(ttS>50)ttS=sqrt(ttS*50);
	  disp.L=exp(-0.5*ttS);

	}//compute_misfitDISP_single
//-----------------------------------------------------	
	int compute_misfitDISP_single_phi(dispdef &disp,double &tempv,double &tS, int flag)
	{
          int i,j,k;
          double tp,td,p,ttS,tempv1=0.,tempv2=0.;
	  float T,phidiff1,phidiff2; 


	  //----IMPORTANT PARAMETER, the period of phi----set this according to the flag in forward computation (compute_AZdisp)
	  if(flag==1)T=180.;//for 2psi Rayleigh wave
	  else if (flag==2)T=90.;//for 4psi Love wave
	  else{printf("### compute_misfitDISP_single_phi, wrong flag, use 1 or 2!\n");exit(0);}

          for(i=0;i<disp.npper;i++)
                {
      		  phidiff1=fabs(disp.pvelo[i]-disp.pvel[i]);
		  while(phidiff1>T){phidiff1-=T;}
		  //--check @@
		  if(phidiff1<0 or phidiff1>T){printf("@@@ hey, wrong phidiff1(%g)\n",phidiff1);exit(0);}
		  //
		  //the max absolute value of ang difference is T/2
		  phidiff2=phidiff1>T/2?T-phidiff1:phidiff1;
		  tempv1=tempv1+pow(phidiff2/disp.unpvelo[i],2);}
          if(tempv1>0.)
                {
                  disp.pmisfit=sqrt(tempv1/disp.npper);
                  tS=tempv1;
                  if(tS>50.)tS=sqrt(tS*50);
                  if(tS>50.)tS=sqrt(tS*50);
                  disp.pL=exp(-0.5*tS);
                }//if

          for(i=0;i<disp.ngper;i++)
                {
		  phidiff1=fabs(disp.gvelo[i]-disp.gvel[i]);
		  while(phidiff1>T){phidiff1-=T;}
		  phidiff2=phidiff1>T/2?T-phidiff1:phidiff1;
		  tempv2=tempv2+pow(phidiff2/disp.ungvelo[i],2);}
          if(tempv2>0.)
                {
                  disp.gmisfit=sqrt(tempv2/disp.ngper);
                  tS=tempv2;
                  if(tS>50)tS=sqrt(tS*50);
                  if(tS>50)tS=sqrt(tS*50);
                  disp.gL=exp(-0.5*tS);
                }//if	

	  tempv=tempv1+tempv2;
	  disp.misfit=sqrt((tempv1+tempv2)/(disp.npper+disp.ngper));
	  ttS=(tempv1+tempv2);
	  if(ttS>50)ttS=sqrt(ttS*50);
	  if(ttS>50)ttS=sqrt(ttS*50);
	  disp.L=exp(-0.5*ttS);

	}//compute_misfitDISP_single_phi

	
//-----------------------------------------------------	
	int compute_misfitDISP(modeldef &model,int Rflag, int Lflag, int AziampRflag, int AziampLflag, int AziphiRflag, int AziphiLflag, float inpamp, float inpphi )
	{
	  double tmisfit1,tL1,tS1,tS2,tS,tS3,tS4,tS5,tS6;
	  double tempvR=0.,tempvRamp=0.,tempvL=0.,tempvLamp=0.,tempvRphi=0.,tempvLphi=0.;
	  int Rn,Rampn,Rphin,Ln,Lampn,Lphin;
	  Rn=Rampn=Rphin=Ln=Lampn=Lphin=0;
	  //printf("nT: nRp=%d nRg=%d nLp=%g nLg=%d\n",model.data.Rdisp.npper,model.data.Rdisp.ngper,model.data.Ldisp.npper,model.data.Ldisp.ngper);//---test
	  if (Rflag>0){compute_misfitDISP_single(model.data.Rdisp,tempvR,tS1);Rn=model.data.Rdisp.npper+model.data.Rdisp.ngper;}
	  /*tmisfit1=sqrt((tempv1+tempv2)/(model.data.Rdisp.npper+model.data.Rdisp.ngper));
	  //tL1=model.data.Rdisp.gL*model.data.Rdisp.pL;
	  tS=(tempv1+tempv2);
	  if(tS>50)tS=sqrt(tS*50);
	  if(tS>50)tS=sqrt(tS*50);
	  tL1=exp(-0.5*tS);
	  model.data.Rdisp.L=tL1;
	  model.data.Rdisp.misfit=tmisfit1;*/

	  if (Lflag>0){compute_misfitDISP_single(model.data.Ldisp,tempvL,tS2);Ln=model.data.Ldisp.npper+model.data.Ldisp.ngper;}
	  
	  if(AziampRflag>0){compute_misfitDISP_single(model.data.AziampRdisp,tempvRamp,tS3);Rampn=model.data.AziampRdisp.npper+model.data.AziampRdisp.ngper;}
	  if(AziampLflag>0){compute_misfitDISP_single(model.data.AziampLdisp,tempvLamp,tS4);Lampn=model.data.AziampLdisp.npper+model.data.AziampLdisp.ngper;}
	  if(AziphiRflag>0){compute_misfitDISP_single_phi(model.data.AziphiRdisp,tempvRphi,tS5,1);Rphin=model.data.AziphiRdisp.npper+model.data.AziphiRdisp.ngper;}
	  if(AziphiLflag>0){compute_misfitDISP_single_phi(model.data.AziphiLdisp,tempvLphi,tS6,2);Lphin=model.data.AziphiLdisp.npper+model.data.AziphiLdisp.ngper;}

	  //*************how to compute data.L data.misfit when R and L are both presented?
	  tS=0.5*((1.-inpamp-inpphi)*(tempvR+tempvL)+inpamp*(tempvRamp+tempvLamp)+inpphi*(tempvRphi+tempvLphi));
	  //tS=0.5*((1.-inpamp-inpphi)*(tempvR+tempvL*4)+inpamp*(tempvRamp+tempvLamp)+inpphi*(tempvRphi+tempvLphi)); //---test--- weitLoveMore
	  if (tS>50.)tS=sqrt(tS*50);
	  if(tS>50.)tS=sqrt(tS*50);
	  model.data.L=exp(-0.5*tS);
	  model.data.misfit=(1-inpamp-inpphi)*sqrt((tempvR+tempvL)/max(1,Rn+Ln))+inpamp*sqrt((tempvRamp+tempvLamp)/max(1,Rampn+Lampn))+inpphi*sqrt((tempvRphi+tempvLphi)/max(1,Rphin+Lphin));
 	  //*****************
	  return 1;
	}//compute misfitDISP
//-----------------------------------------------------	

//-----------------------------------------------------	


 int goodmodel( modeldef &model, vector<int> vmono, vector<int> vgrad,int Rflag, int Lflag, int isoflag)
	 {
	  // for anisotropic case, cannot check goodmodel for both R and L. need to run this twice, check L R seperately
	  int i,j;
	  double var=0;
	  vector<int>::iterator id;
	  double gradient;
	  if(isoflag>0 or Rflag>0){ // iso case, R and L have the same model.g.value1
	      for(i=0;i<model.ngroup-1;i++)//between groups, require positive vel jump
		{if(model.groups[i+1].vsvvalue1[0]<model.groups[i].vsvvalue1.back()) 
			{//cout<<"case 1\n"; //---test----
       			return 0;}
		}  	 

	      for(id=vmono.begin();id<vmono.end();id++)// monotonic change in group vmono[?]
              {
	    	j=*id;
	    	for(i=0;i<model.groups[j].nlay-1;i++)
		    { /*gradient=(model.groups[j].thick1[i])/(model.groups[j].vsvvalue1[i]-model.groups[j].vsvvalue1[i+1]);
		      if(gradient>0. and gradient <70.)
		      {//printf("V1=%g h1=%g V2=%g h2=%g\n",model.groups[j].vsvvalue1[i],model.groups[j].vsvvalue1[i+1],model.groups[j].thick1[i],model.groups[j].thick1[i+1]);
		       return 0;}*/
		
			if(model.groups[j].vsvvalue1[i]>model.groups[j].vsvvalue1[i+1]){return 0;}
		    }//for i
	     }//for id
		///*
	      //--------------revised; newly added, Jun 2, 2012--------
	      //require the slope in mantle part larger than 70.
	     for(i=0;i<model.groups[2].nlay-1;i++){
	         gradient=(model.groups[2].thick1[i])/(model.groups[2].vsvvalue1[i]-model.groups[2].vsvvalue1[i+1]);
		 if(gradient>0. and gradient <50.){
			 //printf("mantle grad !\n");
			 return 0;}
	     }
	     // require the vel at all depth (the maximum = g.tthcik, which comes from the input) to be smaller than 4.9 ---- revised on Aug 23, 2012------
		//*/
	     for (i=0;i<model.groups[2].nlay;i++){
	     if(model.groups[2].vsvvalue1[i]>4.9 or model.groups[2].vshvalue1[i]>4.9){
		     //printf("too large Vmoho\n");
		     return 0;}
	     }
		//*/
	     //------------------------------------------------
	     for(id=vgrad.begin();id<vgrad.end();id++)//gradient check for the 1st two velue in group vgrad[?]
	      {
	        j=*id;
	        if(model.groups[j].vsvvalue1[1]<model.groups[j].vsvvalue1[0])
		{ //cout<<"1st two value gradient\n"; //---test----
		    return 0;}
	      } 

	      ///*
	      //-------- revised; newly added, Apr 22, 2014. constrains on the Vp/Vs=(vph+vpv)/(vsh+vsv)
	      for(i=1;i<model.ngroup;i++){//skip the sedimentary layer; by default, the sediment is group0
		for(j=0;j<model.groups[i].nlay-1;j++){
			var=(model.groups[i].vphvalue1[j]+model.groups[i].vpvvalue1[j])/(model.groups[i].vshvalue1[j]+model.groups[i].vsvvalue1[j]);
			if(var<1.65 or var>1.85){
				//printf("reject ig=%d ilay=%d vp/vs=(%g+%g)/(%g+%g)=%g\n",i,j,model.groups[i].vphvalue1[j],model.groups[i].vpvvalue1[j],model.groups[i].vshvalue1[j],model.groups[i].vsvvalue1[j],var);
				return 0;}
	     	/*//---put upper limit on the amp of Vp and Vs anisotropy. //---test---
	     	if(fabs(2*(model.groups[i].vshvalue1[j]-model.groups[i].vsvvalue1[j])/(model.groups[i].vshvalue1[j]+model.groups[i].vsvvalue1[j]))>0.15 or fabs(2*(model.groups[i].vphvalue1[j]-model.groups[i].vpvvalue1[j])/(model.groups[i].vphvalue1[j]+model.groups[i].vpvvalue1[j]))>0.1){
			return 0;}*/
		}
	      }//vp/vs constraint
	     //*/
	  }//if isoflag or Rflag>0

	  
	  else if (Lflag>0){
              for(i=0;i<model.ngroup-1;i++)//between layers
                {if(model.groups[i+1].vshvalue1[0]<model.groups[i].vshvalue1.back())
                 return 0;} 
              for(id=vmono.begin();id<vmono.end();id++)// monotonic change in group vmono[?]
              {
                j=*id;
                for(i=0;i<model.groups[j].nlay-1;i++)
		    { gradient=(model.groups[j].thick1[i])/(model.groups[j].vshvalue1[i]-model.groups[j].vshvalue1[i+1]);
		      if(gradient>0. and gradient <70.)return 0; }
	     }//for id
              for(id=vgrad.begin();id<vgrad.end();id++)//gradient check for the 1st two velue in group vgrad[?]
              {
                j=*id;
                if(model.groups[j].vshvalue1[1]<model.groups[j].vshvalue1[0])
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
			if(model.groups[j].vshvalue1[i]<model.groups[j].vsvvalue1[i] or (model.groups[j].vphvalue1[i]<model.groups[j].vpvvalue1[i])){
			//if(model.groups[j].vshvalue1[i]<model.groups[j].vsvvalue1[i] ){ //rm the VpRA>0 constraint, ---test---
				//printf("--@@@ ng=%d, nv=%d, vsh=%g, vsv=%g\n",j,i,model.groups[j].vshvalue1[i],model.groups[j].vsvvalue1[i]);
				return 0;
			}
		}//for i
	  }//for id
	  return 1;
	}//positiveAni
	
