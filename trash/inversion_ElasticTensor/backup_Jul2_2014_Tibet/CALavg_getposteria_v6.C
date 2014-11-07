// this code reads the binary file, and compute the average model.
// the binary file contains all the accepted models, it was written during the inversion process.
// this version write out the vsv vsh result for both step1 and step2 do_inv. Also, write out prior and posteria distribution 
// attentin, this version is for the smooth model, the thickness of the last layer is 0 for smooth model, so in the subroutine model_avg_sub, (j=0;j<modlst[i].laym0.nlayer-1 [NOT nlayer!!!];j++)
//
//the connecting has problem dealing with layerized model, who has multiple values at one single depth. need to be modified!
// this version, output both extrinsic values, and apparent values
// this version ,modified on Jan 20, 2014, output the info for vsv,vsh,vpv,vph,eta,theta,phi; while previous version output vsv,vsh,and vsRA
// this version, modified on jan23, 2014, could deal with multi-phi situation. seperate the final average model into at most 2 groups based on the phi value distribution.
// this version, use para_avg_multiple_gp_v2.C instead of para_avg_multiple_gp.C, could also deal with both crust and mantle phi; and the mantle phi seperation can be diabled by setting idphiM<0

// this version, use the para_avg_multiple_gp_v3.C, which enables the computation of parabest for each phi group, so there are multiple parabest output
//
#include<iostream>
#include<algorithm>
#include<vector>
#include<cmath>
#include<fstream>
#include"./string_split.C"
#include"./generate_Bs.C"
#include"./gen_random.C"
#include"./INITstructure.h"
#include "init.C"
#include "./BIN_rw.C"
//#include "para_avg_multiple_gp.C"
#include "para_avg_multiple_gp_v3.C"
using namespace std;
//----------------------------------------------------- 
//----------------------------------------------------- 
	int para_avg_V(vector<paradef> &paralst, paradef &parabest,paradef &paraavg,vector<double> &parastd,vector<int> &idlst){
	//copied from CALpara_isolay.C, but here, only deal with the Vparameter, not the Love parameters
	//get the average parameters from a list of acceptable paras
	    int i,j,k,size,Ngood,ibadp,idmin;
	    double mismin=1e10;
 	    vector<int> badv;

	    initpara(paraavg);
	    initpara(parabest);
	    size=paralst.size();
	    if(size<1) return -1;
		
	    paraavg=paralst[0];
 	    parabest=paralst[0];
	    parastd.clear();

//what's the size of love para
	    for(j=0;j<paralst[0].npara;j++){
                parastd.push_back(0.);
		paraavg.parameter[j]=0.;
            }//for j

	    for(i=0;i<size;i++)
		{if(paralst[i].misfit<mismin) {idmin=i;mismin=paralst[i].misfit;}}
	    parabest=paralst[idmin];
	    //min(2*Mis_min,Mis_min+0.5)
	    printf("@@@ para_avg mismin=%g\n",mismin);
	    //if(mismin>1.5) mismin=mismin*2;
	    //else mismin=mismin+0.5;
	    mismin=mismin+0.5;
	    // mismin=1e10;
	    Ngood=0;
	    for(i=0;i<size;i++){
		if(paralst[i].misfit<mismin){
		  Ngood++;
		  idlst.push_back(i);
		  for(j=0;j<paralst[0].npara;j++)
			paraavg.parameter[j]+=paralst[i].parameter[j];
		}//if
	    }//for i
	   for(i=0;i<paralst[0].npara;i++){
		paraavg.parameter[i]/=Ngood;
	   }//for i
	   //cout<<"Ngood="<<Ngood<<endl; 
	   printf("@@@ para_avg size =%d, Ngood=%d mismin=%g\n",size,Ngood,mismin);
	   //if(Ngood>size){printf("WRONG HERE! Ngood > size!!\n");exit(0);}
	   for(i=0;i<paralst[0].npara;i++){
		for(j=0;j<Ngood;j++){
		    k=idlst[j];
		    parastd[i]+=pow(paralst[k].parameter[i]-paraavg.parameter[i],2);
		}
		parastd[i]=sqrt(parastd[i]/Ngood);
		if(parastd[i]>50){//indicating it's phi and has period problem that hasn't been taken into account
			badv.push_back(i);
		}
	    }//for i
	 
	    //========added Nov23, 2013; check the para_std, if std too large (>50), indicating it's phi and has period problem that hasn't been taken into account
	    for(i=0;i<badv.size();i++){
		ibadp=badv[i];
		paraavg.parameter[ibadp]=0.;
		for(j=0;j<Ngood;j++){
			k=idlst[j];
			if(paralst[k].parameter[ibadp]>90){paralst[k].parameter[ibadp]-=180;}
			paraavg.parameter[ibadp]+=paralst[k].parameter[ibadp];
		}
		paraavg.parameter[ibadp]/=Ngood;
	    }//for i <badv.size
	    for(i=0;i<badv.size();i++){
		ibadp=badv[i];
		parastd[ibadp]=0.;
		for(j=0;j<Ngood;j++){
			k=idlst[j];
	    		parastd[ibadp]+=pow(paralst[k].parameter[ibadp]-paraavg.parameter[ibadp],2);
		}
		parastd[ibadp]=sqrt(parastd[ibadp]/Ngood);
	    }

	    //return 1;
	    return idmin;
	}//para_avg_V
//-----------------------------------------------------

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
	  //vlst[ithick]=[vsv,vsh,vpv,vph,eta,theta,phi,  vsvmin~phimin,  vsvmax~phimax];
	  //stdlst[ithick]=[vsvstd~phistd]

	  vector<float> tv,dep2lst,dep1lst;
	  vector<int> iddep1lst;
	  float th,tth,fm1,fm2,dep;
	  float vsv,vsh,dep1,dep2,vsvmin,vsvmax,vshmin,vshmax;
	  float tvsh,tvsv,tani,ani,animin,animax;
	  float fm3;
	  char str[500];
	  int i,j,Ncount,size,k,flag;
	  vector<vector<float> > tvlst;
	  float vpvmin,vpvmax,vphmin,vphmax,etamin,etamax,thetamin,thetamax,phimin,phimax;
	  float tvpv,tvph,teta,ttheta,tphi,vpv,vph,eta,theta,phi;

	  hlst.clear();vlst.clear();stdlst.clear();tvlst.clear();
		  
	  size=modlst.size();
	
	  for(i=0;i<size;i++){	// obtain the end depth for this mod_avg
	    dep2=0.;  
	    for(k=0;k<=ng;k++){dep2=dep2+modlst[i].groups[k].thick;}
	    dep2lst.push_back(dep2);
	    dep1lst.push_back(0.);
	    iddep1lst.push_back(0);
	  }
	  flag=0;
	  for(th=h0;th<=h+dth;th=th+dth){//thickness
		//printf("th=%g h=%g dep1[0]=%g dep2[0]=%g\n",th,h,dep1lst[0],dep2lst[0]);//---test---
		if(th>h){//modified on Nov 24, 2013, allow the th to go very close to h 
			if(flag==0){flag=1;th=h-0.001;}
			else{break;}
		}
		if(fabs(th-h+dth)<1e-4)th=h;
		vsv=0.;vsh=0.;Ncount=0;fm1=0.;fm2=0.;ani=0.;fm3=0.;tth=0.;
		vpv=0.;vph=0.;eta=0.;theta=0.;phi=0.;
		tvlst.clear();
		vsvmin=1e10;vsvmax=-1;
		vshmin=1e10;vshmax=-1;
		animin=1e10;animax=-1;
		vpvmin=1e10;vpvmax=-1.;
		vphmin=1e10;vphmax=-1.;
		etamin=1e10;etamax=-1.;
		thetamin=1e10;thetamin=-1.;
		phimin=1e10;phimax=-1.;
		for(i=0;i<size;i++){//iterate within model list
		  //dep2=0.;
		  //for(k=0;k<=ng;k++){dep2=dep2+modlst[i].groups[k].thick;}
		  dep2=dep2lst[i];// end depth
		  //dep1=dep2-modlst[i].groups[ng].thick;//dep1--upper of the layer; dep2--lower of the layer
		  if(th>dep2)continue;//tth=dep2; //tth = min(dep2,th) ####if th>dep2, continue, skip present mod, jumpt to next mod??, this means th is outside the present layer of present model, so skip.
		  else tth=th;
		  dep1=0;
		  for(j=0;j<modlst[i].laym0.nlayer-1;j++){ // modeified on Aug 27, 2012. For the smooth model, the thickness of the last layer is 0 
		    dep1=dep1+modlst[i].laym0.thick[j];
		  //for(j=iddep1lst[i];j<modlst[i].laym0.nlayer;j++){
		//	dep1=dep1+dep1lst[i]+modlst[i].laym0.thick[j];
			if(dep1>=tth ){//  the interpolation, there was a bug here (for smooth model, though it works for layered model), changed on Aug27, 2012. another bug modified on Oct 8, 2013
				tvsv=(modlst[i].laym0.vsv[j+1]-modlst[i].laym0.vsv[j])/(modlst[i].laym0.thick[j])*(tth-dep1+modlst[i].laym0.thick[j])+modlst[i].laym0.vsv[j];
				tvsh=(modlst[i].laym0.vsh[j+1]-modlst[i].laym0.vsh[j])/(modlst[i].laym0.thick[j])*(tth-dep1+modlst[i].laym0.thick[j])+modlst[i].laym0.vsh[j];
				tvpv=(modlst[i].laym0.vpv[j+1]-modlst[i].laym0.vpv[j])/(modlst[i].laym0.thick[j])*(tth-dep1+modlst[i].laym0.thick[j])+modlst[i].laym0.vpv[j];
				tvph=(modlst[i].laym0.vph[j+1]-modlst[i].laym0.vph[j])/(modlst[i].laym0.thick[j])*(tth-dep1+modlst[i].laym0.thick[j])+modlst[i].laym0.vph[j];
				teta=(modlst[i].laym0.eta[j+1]-modlst[i].laym0.eta[j])/(modlst[i].laym0.thick[j])*(tth-dep1+modlst[i].laym0.thick[j])+modlst[i].laym0.eta[j];
				ttheta=(modlst[i].laym0.theta[j+1]-modlst[i].laym0.theta[j])/(modlst[i].laym0.thick[j])*(tth-dep1+modlst[i].laym0.thick[j])+modlst[i].laym0.theta[j];
				tphi=(modlst[i].laym0.phi[j+1]-modlst[i].laym0.phi[j])/(modlst[i].laym0.thick[j])*(tth-dep1+modlst[i].laym0.thick[j])+modlst[i].laym0.phi[j];	
				
				if(tvsv<0 or tvsv>5){
					printf("test--- th=%g, imod=%d, ilay=%d, tvsv=%g=[(%g-%g)/(%g)+%g]\n",th,i,j,tvsv,modlst[i].laym0.vsv[j+1],modlst[i].laym0.vsv[j],modlst[i].laym0.thick[j],modlst[i].laym0.vsv[j]);}
				tani=100*(tvsh-tvsv)/(sqrt((2*tvsv*tvsv+tvsh*tvsh)/3.0));
				vsv=vsv+tvsv;// modified on Aug27, 2012
				vsh=vsh+tvsh;//
				ani=ani+tani;
				vpv=vpv+tvpv;
				vph=vph+tvph;
				eta=eta+teta;
				theta=theta+ttheta;
				phi=phi+tphi;
				Ncount++;	
				if (tvsv>vsvmax)vsvmax=tvsv;
				if(tvsv<vsvmin)vsvmin=tvsv;
				if(tvsh>vshmax)vshmax=tvsh;
				if(tvsh<vshmin)vshmin=tvsh;
				if(tani>animax)animax=tani;if(tani<animin)animin=tani;
				if(tvpv>vpvmax)vpvmax=tvpv;if(tvpv<vpvmin)vpvmin=tvpv;
				if(tvph>vphmax)vphmax=tvph;if(tvph<vphmin)vphmin=tvph;
				if(teta>etamax)etamax=teta;if(teta<etamin)etamin=teta;
				if(ttheta>thetamax)thetamax=ttheta;if(ttheta<thetamin)thetamin=ttheta;
				if(tphi>phimax)phimax=tphi;if(tphi<phimin)phimin=tphi;
				//---tvlst: vsv,vsh ..., phi
				tv.clear();tv.push_back(tvsv);tv.push_back(tvsh);//tv.push_back(tani);
				tv.push_back(tvpv);tv.push_back(tvph);tv.push_back(teta);tv.push_back(ttheta);tv.push_back(tphi);
				tvlst.push_back(tv);
			//	dep1lst[i]=dep1-modlst[i].laym0.thick[j];//
			//	iddep1lst[i]=j;//
				break;
			}//if dep>tth
		  }//for j Nlay
		}//for i Nmodel,size
		vsv=vsv/Ncount;
		vsh=vsh/Ncount;
		ani=ani/Ncount;
		vpv=vpv/Ncount;
		vph/=Ncount;
		eta/=Ncount;
		theta/=Ncount;
		phi/=Ncount;
		if (isnan(vsv) or isnan(vsh) or isnan(ani) or isnan(vpv) or isnan(vph) or isnan(eta) or isnan(theta) or isnan(phi)){
			printf("Hey, NaN happen!!! something wrong!!\n Ncount=%d th=%g tth=%g h0=%g h=%g\n",Ncount,th,tth,h0,h);
			sprintf(str,"echo WRONG Ncount = %d >> point_finished_Ani.txt ",Ncount);
			system(str);
			//exit(0);
		}
		//---vlst: vsv~phi, vsvmin~phimin, vsvmax~phimax---
		tv.clear();
		//tv.push_back(vsv);tv.push_back(vsh);tv.push_back(vsvmin);tv.push_back(vsvmax);tv.push_back(vshmin);tv.push_back(vshmax);tv.push_back(ani);tv.push_back(animin);tv.push_back(animax);
		tv.push_back(vsv);tv.push_back(vsh);tv.push_back(vpv);tv.push_back(vph);tv.push_back(eta);tv.push_back(theta);tv.push_back(phi);
		tv.push_back(vsvmin);tv.push_back(vshmin);tv.push_back(vpvmin);tv.push_back(vphmin);tv.push_back(etamin);tv.push_back(thetamin);tv.push_back(phimin);
		tv.push_back(vsvmax);tv.push_back(vshmax);tv.push_back(vpvmax);tv.push_back(vphmax);tv.push_back(etamax);tv.push_back(thetamax);tv.push_back(phimax);

		vlst.push_back(tv);
		hlst.push_back(th);
	
		//---stdlst: vsvstd~phistd
		tv.clear();
		fm1=0.;fm2=0.;	
		for(i=0;i<Ncount;i++){ 
			fm1=fm1+pow(tvlst[i][0]-vsv,2);
			fm2=fm2+pow(tvlst[i][1]-vsh,2);//----there was a bug here, fixed on May 18, 2012
		}
		fm1=sqrt(fm1/Ncount);
		fm2=sqrt(fm2/Ncount);
		tv.push_back(fm1);tv.push_back(fm2);
		fm1=0.;fm2=0.;		
		for(i=0;i<Ncount;i++){ 
			fm1=fm1+pow(tvlst[i][2]-vpv,2);
			fm2=fm2+pow(tvlst[i][3]-vph,2);
		}
		fm1=sqrt(fm1/Ncount);
		fm2=sqrt(fm2/Ncount);
		tv.push_back(fm1);tv.push_back(fm2);
		fm1=0.;fm2=0.;fm3=0.;
		for(i=0;i<Ncount;i++){ 
			fm1=fm1+pow(tvlst[i][4]-eta,2);
			fm2=fm2+pow(tvlst[i][5]-theta,2);
			fm3=fm3+pow(tvlst[i][6]-phi,2);
		}
		fm1=sqrt(fm1/Ncount);
		fm2=sqrt(fm2/Ncount);
		fm3=sqrt(fm3/Ncount);
		tv.push_back(fm1);tv.push_back(fm2);tv.push_back(fm3);
		
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
/* 	int integrate_mod(vector<vector<float> > &outVsv,vector<vector<float> > &outVsh,vector<vector<float> > &outAni, vector<float> hlst,vector<vector<float> > vlst,vector<vector<float> > stdlst)
	{ int i;vector<float> tv;
          for(i=0;i<vlst.size();i++){
                tv.clear();
                tv.push_back(hlst[i]);tv.push_back(vlst[i][0]);tv.push_back(stdlst[i][0]);tv.push_back(vlst[i][2]);tv.push_back(vlst[i][3]);
                outVsv.push_back(tv);
                tv.clear();
                tv.push_back(hlst[i]);tv.push_back(vlst[i][1]);tv.push_back(stdlst[i][1]);tv.push_back(vlst[i][4]);tv.push_back(vlst[i][5]);
                outVsh.push_back(tv);
	  	tv.clear();
		tv.push_back(hlst[i]);tv.push_back(vlst[i][6]);tv.push_back(stdlst[i][2]);tv.push_back(vlst[i][7]);tv.push_back(vlst[i][8]);
		outAni.push_back(tv);	  	
	  }
	  vector<float>().swap(tv);
	  return 1;
	}//integrate_mod
*/
//-----------------------------------------------------------------
	int connect_modmin(vector<vector<float> > &out,vector<float> hlst,vector<vector<float> > vlst,vector<vector<float> > stdlst,float h1,float h2, int id)
	{// can choose output avg-std or the min; 
	 //the vlst has info of min (column 7~13), info of avg is column 0~6; info of max is 14~20
	  int i;vector<float> tv;
	  float h,v;
 	  for(i=0;i<vlst.size();i++){
	  	h=hlst[i];
		if (h<h1)continue;
		else if (h>=h2)break;
		tv.clear();
		v=vlst[i][id]-stdlst[i][id];//avg-std
		//v=vlst[i][7+id];//min
		tv.push_back(h);tv.push_back(v);
		out.push_back(tv);
	  }
	  vector<float>().swap(tv);
	  return 1;
	}//connect_modmin
//-----------------------------------------------------------------
	int connect_modmid(vector<vector<float> > &out,vector<float> hlst,vector<vector<float> > vlst,vector<vector<float> > stdlst,float h1,float h2, int id)
	{// modified on Aug 28, 2012. add the std into vector, and output it. previous version has some problems (the hlst for Vmin, Vmid, Vmax could be different, so when read the output file, obtaining unc using Vid-Vmin is incorrect!!)
	  int i;vector<float> tv;
	  float h,v;
 	  for(i=0;i<vlst.size();i++){
	  	h=hlst[i];
		if (h<h1)continue;
		else if (h>=h2)break;
		tv.clear();
		v=vlst[i][id];
		tv.push_back(h);tv.push_back(v);tv.push_back(stdlst[i][id]);
		out.push_back(tv);
	  }
	  vector<float>().swap(tv);
	  return 1;
	}//connect_modmin
//-----------------------------------------------------------------
//-----------------------------------------------------------------
	int connect_modmax(vector<vector<float> > &out,vector<float> hlst,vector<vector<float> > vlst,vector<vector<float> > stdlst,float h1,float h2, int id)
	{
	  int i;vector<float> tv;
	  float h,v;
 	  for(i=0;i<vlst.size();i++){
	  	h=hlst[i];
		if (h<h1)continue;
		else if (h>=h2) break;
		tv.clear();
		v=vlst[i][id]+stdlst[i][id];//avg+std
		//v=vlst[i][14+id];//max
		tv.push_back(h);tv.push_back(v);
		out.push_back(tv);
	  }
	  vector<float>().swap(tv);
	  return 1;
	}//connect_modmin
//-----------------------------------------------------------------
//-----------------------------------------------------------------
	int write_modavg( const char *foutnm,vector<vector<float> > &Vsv1,vector<vector<float> > &Vsv2,vector<vector<float> > &Vsv3,float Hsed,float Hsedstd,float Hmoho,float Hmohostd )
	{// attention, the hlst for Vmin, Vmid, Vmax could be different in many cases. So, use the associated hlst while plotting/reading ...
// this writting method has some problems, since the size of Vmin,Vmid,Vmax is different, I only write to the min(len(V)). In this case, the V with the largest size cannot be fully written out (cannot reach the maximum depth).
 	
	    //Vsv1[h,vmin];Vsv2[h,vmid,vstd];Vsv3[h,vmax]
	    FILE *ftmp;
	    int i,j,sizemin;
	    if((ftmp=fopen(foutnm,"w"))==NULL){printf("####in write-modavg, canot open file %s to write!\n",foutnm);exit(0);};
	    fprintf(ftmp,"sedi %8g %8g\n",Hsed,Hsedstd);
	    fprintf(ftmp,"moho %8g %8g\n",Hmoho,Hmohostd);
	    sizemin=min(Vsv1.size(),Vsv2.size());
	    if(sizemin>Vsv3.size())sizemin=Vsv3.size();
	    for(i=0;i<sizemin;i++){
	       //fprintf(ftmp,"%8g %8g %8g %8g\n",Vsv1[i][0],Vsv1[i][1],Vsv2[i][1],Vsv3[i][1]);
	       fprintf(ftmp,"%8g %8g %8g %8g %8g %8g %8g\n",Vsv1[i][0],Vsv1[i][1],Vsv2[i][0],Vsv2[i][1],Vsv3[i][0],Vsv3[i][1],Vsv2[i][2]);
	    }
	    fclose(ftmp);

	    return 1;
	}//write_modavg
//--------------------------------------------------
       //int model_avg2( const char *fvsvnm, const char *fvshnm,char *faninm, vector<modeldef> &modlstall, vector<int> &idlst)
	int model_avg2(vector<string> fnmlst, vector<modeldef> &modlstall, vector<int> &idlst)
	{
	  int size,i,n,k;
 	  paradef para;
	  modeldef mod;
	  char str[200];
	  float depmax; 
	  char vsvnm[200],vshnm[200];
	  float Hsed=0.,Hmoho=0.,Hmat=0.,Hsedstd,Hmohostd,fm1=0.,fm2=0.;
	  vector<vector<float> > v1lst,v2lst,v3lst,std1lst,std2lst,std3lst;
	  vector<vector<float> > Vsvmin,Vsvmid,Vsvmax,Vshmin,Vshmid,Vshmax;
	  vector<vector<float> > Animin,Animid,Animax;
	  vector<float> h1lst,h2lst,h3lst,tv;
	  vector<modeldef> modlst;
	  time_t start;
	  //cout<<"---begin! "<<time(0)<<endl;
	  start=time(0);
	  size = idlst.size();  
	  if(size<1){printf("##### in mod_avg, idlst is empty!\n");exit(0);}
	  //printf("ok1\n");//--test--
	  for (i=0;i<size;i++){
		k=idlst[i];
	 	modlst.push_back(modlstall[k]);
		Hsed=Hsed+modlst[i].groups[0].thick; // THIS IS thickness not depth
		Hmoho=Hmoho+modlst[i].groups[1].thick;// Hmoho changes from thickness to depth later
		Hmat=Hmat+modlst[i].groups[2].thick;
	  	//sprintf(str,"echo %g %g >> temp_h.txt\n",modlst[i].groups[0].thick+modlst[i].groups[1].thick+modlst[i].groups[2].thick,modlst[i].tthick);
		//system(str);
	  }//for i
	  Hsed=Hsed/size;
	  Hmoho=Hmoho/size;
	  Hmat=Hmat/size;
	  //printf("ok2\n");//--test--
	  for(i=0;i<size;i++){
	 	fm1=fm1+pow(modlst[i].groups[0].thick-Hsed,2);
	 	fm2=fm2+pow(modlst[i].groups[1].thick-Hmoho,2);
	  }
	  //printf("ok3\n");//--test--
	  Hsedstd=sqrt(fm1/size);
	  Hmohostd=sqrt(fm2/size);
	  Hmoho=Hmoho+Hsed; // *********  the previous Hmoho is the thickness of Moho not the depth of Moho; modified on Oct 10, 2012; 
	  depmax=Hmoho+Hmat;
	  printf("test-- Hsed=%g Hsedstd=%g\n",Hsed,Hsedstd);
	  //printf("ok4\n");//--test--
	  model_avg_sub(v1lst,std1lst,h1lst,modlst,0,0,Hsed+Hsedstd,0.05);
	  cout<<"finish sed\n";//--test
	  model_avg_sub(v2lst,std2lst,h2lst,modlst,1,Hsed-Hsedstd,Hmoho+Hmohostd,1.);
	  cout<<"finish cst\n";//---test
	  model_avg_sub(v3lst,std3lst,h3lst,modlst,2,Hmoho-Hmohostd,depmax,1.);
	  cout<<"finish mant\n";//---test
	  
	  cout<<"---connect model min/mid/max time_used="<<time(0)-start<<"\n";//---test---
	  for(i=0;i<7;i++){
	  printf("i=%d, fnm=%s\n",i,fnmlst[i].c_str());
	  Vsvmin.clear();Vsvmid.clear();Vsvmax.clear();
	  connect_modmin(Vsvmin,h1lst,v1lst,std1lst,0,Hsed+Hsedstd,i);
	  connect_modmin(Vsvmin,h2lst,v2lst,std2lst,Hsed+Hsedstd,Hmoho+Hmohostd,i);
	  connect_modmin(Vsvmin,h3lst,v3lst,std3lst,Hmoho+Hmohostd,depmax,i);
	  
	  connect_modmid(Vsvmid,h1lst,v1lst,std1lst,0,Hsed,i);
	  connect_modmid(Vsvmid,h2lst,v2lst,std2lst,Hsed,Hmoho,i);
	  connect_modmid(Vsvmid,h3lst,v3lst,std3lst,Hmoho,depmax,i);
	  
	  connect_modmax(Vsvmax,h1lst,v1lst,std1lst,0,Hsed-Hsedstd,i);
	  connect_modmax(Vsvmax,h2lst,v2lst,std2lst,Hsed-Hsedstd,Hmoho-Hmohostd,i);
	  connect_modmax(Vsvmax,h3lst,v3lst,std3lst,Hmoho-Hmohostd,depmax,i);
	  write_modavg(fnmlst[i].c_str(),Vsvmin,Vsvmid,Vsvmax,Hsed,Hsedstd,Hmoho,Hmohostd);
	  }//for i 		  

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

	  //cout<<"----write model avg----"<<time(0)-start<<"\n";//----test----
	  

          vector<vector<float> >().swap(Vsvmin);
          vector<vector<float> >().swap(Vsvmid);
          vector<vector<float> >().swap(Vsvmax);
          vector<vector<float> >().swap(Vshmin);
          vector<vector<float> >().swap(Vshmid);
          vector<vector<float> >().swap(Vshmax);
          vector<vector<float> >().swap(Animin);
          vector<vector<float> >().swap(Animid);
          vector<vector<float> >().swap(Animax);
	  vector<modeldef>().swap(modlst);

	  return 1;
	}// model_avg2
//-----------------------------------
int get_posteriaDist(double depmin, double depmax, double depstep, vector<double> vsv, vector<double> vsh, vector<double> vpv, vector<double> vph, vector<double> eta, vector<double> theta, vector<double> phi, vector<double> thick, FILE *out)
{//for each input model, write out the vsv~phi value from depmin(1st para) to depmax(2ed para) at depdx(3rd) interval.
  int i,size;
  double dep,h;
  vector<double> hlst;
  double tvsv,tvsh,tvpv,tvph,teta,ttheta,tphi;

  size=vsv.size();
  dep=0.;
  hlst.push_back(dep);// there was a bug here, fixed on Aug 27, 2012; The first node is at the surface, depth=0.;
  /*
    *|
     |__
        *|
	 |__
	    *
   */
  for (i=0;i<size;i++){
  	dep=dep+thick[i];
	hlst.push_back(dep);
  }
  if (depmin<hlst[0] or depmax>hlst[size-1]){
    printf("#### the depmin (%g) depmax(%g) is outside the reference H range (%g ~ %g)\n",depmin,depmax,hlst[0],hlst[size-1]);
    return 0;  
  }

  for (h=depmin;h<=depmax;h=h+depstep){
    for(i=1;i<size;i++){
	if (hlst[i]>h){
		//printf("h=%g hlst[i]=%g\n",h,hlst[i]);//--test--
		tvsv=(vsv[i]-vsv[i-1])/(hlst[i]-hlst[i-1])*(h-hlst[i-1])+vsv[i-1];
		tvsh=(vsh[i]-vsh[i-1])/(hlst[i]-hlst[i-1])*(h-hlst[i-1])+vsh[i-1];
		tvpv=(vpv[i]-vpv[i-1])/(hlst[i]-hlst[i-1])*(h-hlst[i-1])+vpv[i-1];
		tvph=(vph[i]-vph[i-1])/(hlst[i]-hlst[i-1])*(h-hlst[i-1])+vph[i-1];
		teta=(eta[i]-eta[i-1])/(hlst[i]-hlst[i-1])*(h-hlst[i-1])+eta[i-1];
		// for theta&phi, the average isn't appropriate, so, the angels near the discountinuity doesn't make sence---need modification
		ttheta=(theta[i]-theta[i-1])/(hlst[i]-hlst[i-1])*(h-hlst[i-1])+theta[i-1];
		tphi=(phi[i]-phi[i-1])/(hlst[i]-hlst[i-1])*(h-hlst[i-1])+phi[i-1];
		fprintf(out,"%8g %8g %8g %8g %8g %8g %8g ",tvsv,tvsh,tvpv,tvph,teta,ttheta,tphi);break;}
    }
  }

  return 1;
}//get_psteriaDist

//-----------------------------------
//-----------------------------------
int main(int argc, char *argv[])
{
// this version, we will calculate the uncertainty of anisotropy------	
  vector<modeldef> modelall;
  vector<paradef> paraall;
  //vector<double> parastd;
  vector<int> signall,iiterall,iaccpall,idlst;
  char inponm[300],dirbin[300],fbname[200],foutnm[300],fvsvnm[300],fvshnm[300],nodeid[10],outdir[300],str[200];
  char fvpvnm[300],fvphnm[300],fetanm[300],fthetanm[300],fphinm[300],fparanm[300],fmodBnm[300];
  char faninm[300];
  //char fnmlst[7][300];
  vector<string> fnmlst;
  FILE *inpo, *out;
  int npoint,i,RAflag,idmin;
  float tvsv,tvsh,tani;
  vector<float> anilst;	  
  paradef paraavg,parabest,parastd;
  vector<paradef> parabestlst;
  vector<int> idminlst;
  float lon,lat;  

  if(argc!=5){
    printf("Usage: XX 1]point_file 2]bin_dir 3]out_dir\n1]point_file: nodeid lon lat\n2]bin_dir: bindir/ani(vsv)_nodeid_lon_lat.bin\n3]out_dir: outdir/vsv(vsh)_lon_lat.txt\n4]compute_avg for flat vel(1) or titled vel(2; the vel from the effective TI medium)\n");
    exit(0);
  }
  sprintf(inponm,argv[1]);
  sprintf(dirbin,argv[2]);
  sprintf(outdir,argv[3]);
  RAflag=atoi(argv[4]);

  printf("inponm=%s\ndirbin=%s\noutdir=%s\n",inponm,dirbin,outdir);

  if((inpo=fopen(inponm,"r"))==NULL){
        printf("Cannot open points file %s!\n",inponm);exit(0);
  }
  sprintf(str,"if [ ! -d %s ]; then mkdir %s;fi",outdir,outdir);
  system(str);

  npoint=0;
  //int ktopo;
  while(1){
    if(fscanf(inpo,"%s %f %f",&nodeid[0],&lon,&lat)==EOF)
    //if(fscanf(inpo,"%s %f %f %d",&nodeid[0],&lon,&lat,&ktopo)==EOF)
        break;
    npoint++;
    printf("=====begin working on point: %d %s %f %f===\n",npoint,nodeid,lon,lat);
    initpara(paraavg);
    initpara(parabest);
    if(RAflag==1){
    sprintf(fbname,"%s/ani_%s_%.1f_%.1f.bin",dirbin,nodeid,lon,lat);
    sprintf(foutnm,"%s/post_%.1f_%.1f.txt",outdir,lon,lat);
    sprintf(fvsvnm,"%s/vsv_%.1f_%.1f.txt",outdir,lon,lat);
    sprintf(fvshnm,"%s/vsh_%.1f_%.1f.txt",outdir,lon,lat);
    sprintf(faninm,"%s/ani_%.1f_%.1f.txt",outdir,lon,lat);
    sprintf(fvpvnm,"%s/vpv_%.1f_%.1f.txt",outdir,lon,lat);
    sprintf(fvphnm,"%s/vph_%.1f_%.1f.txt",outdir,lon,lat);
    sprintf(fetanm,"%s/eta_%.1f_%.1f.txt",outdir,lon,lat);
    sprintf(fthetanm,"%s/theta_%.1f_%.1f.txt",outdir,lon,lat);
    sprintf(fphinm,"%s/phi_%.1f_%.1f.txt",outdir,lon,lat);
    sprintf(fparanm,"%s/para_%.1f_%.1f.txt",outdir,lon,lat);
    sprintf(fmodBnm,"%s/modB_%.1f_%.1f.txt",outdir,lon,lat);
    }
    else if (RAflag==2){
    sprintf(fbname,"%s/ani_%s_%.1f_%.1f.bin_effTI",dirbin,nodeid,lon,lat);
    sprintf(foutnm,"%s/post_%.1f_%.1f.txt_effTI",outdir,lon,lat);
    sprintf(fvsvnm,"%s/vsv_%.1f_%.1f.txt_effTI",outdir,lon,lat);
    sprintf(fvshnm,"%s/vsh_%.1f_%.1f.txt_effTI",outdir,lon,lat);
    sprintf(faninm,"%s/ani_%.1f_%.1f.txt_effTI",outdir,lon,lat);
    sprintf(fvpvnm,"%s/vpv_%.1f_%.1f.txt_effTI",outdir,lon,lat);
    sprintf(fvphnm,"%s/vph_%.1f_%.1f.txt_effTI",outdir,lon,lat);
    sprintf(fetanm,"%s/eta_%.1f_%.1f.txt_effTI",outdir,lon,lat);
    sprintf(fthetanm,"%s/theta_%.1f_%.1f.txt_effTI",outdir,lon,lat);
    sprintf(fphinm,"%s/phi_%.1f_%.1f.txt_effTI",outdir,lon,lat);
    sprintf(fparanm,"%s/para_%.1f_%.1f.txt_effTI",outdir,lon,lat);
    sprintf(fmodBnm,"%s/modB_%.1f_%.1f.txt_effTI",outdir,lon,lat);
    }
    else{
      printf("### CALavg inproper value for the RAflag (should be 1 or 2)\n");
      return 0;
    }
    fnmlst.push_back(fvsvnm);fnmlst.push_back(fvshnm);fnmlst.push_back(fvpvnm);fnmlst.push_back(fvphnm);fnmlst.push_back(fetanm);fnmlst.push_back(fthetanm);fnmlst.push_back(fphinm);
    ifstream mff(fbname);
    if (! mff.good()){
	printf("Cannot access file %s\n",fbname);
    	sprintf(str,"echo %s %.1f %.1f >> points_CannotAccess.txt",nodeid,lon,lat);
	continue;
    }
    mff.close();
    printf("test-- begin read_bin\n");
    read_bin(modelall,paraall,fbname,signall,iiterall,iaccpall);
    printf("test-- finish read_bin\n");
    printf("test-- begin para_avg\n");

    int idphi,idphiM;
    idphi=6;
    //idphiM=paraall[0].npara-1;
    idphiM=-1;
    vector<paradef> paraavglst,parastdlst;
    vector<vector<int> > idlstlst;
    vector<string> fnmlsttmp;
    idminlst=para_avg_multiple_gp(idphi,idphiM,paraall,parabestlst,paraavglst,parastdlst,idlstlst,1);
    if(idminlst.size()<1){cout<<"### in para_avg, incorrect paralst.size()\n";exit(0);}
    for(int ig=0;ig<idlstlst.size();ig++){
	parabest=parabestlst[ig];
	idlst=idlstlst[ig];
    	printf("size=%d\n",idlst.size());
	paraavg=paraavglst[ig];
	parastd=parastdlst[ig];
	sprintf(str,"%s_phigp%d",fparanm,ig);
	if((out=fopen(str,"w"))==NULL){printf("### Cannot open %s to write!!\n",fparanm);exit(0);}
	for(i=0;i<paraavg.npara;i++){ //npara, para_avg, para_avg_std, parabest;
		fprintf(out,"%5d %8.4f %8.4f %8.4f\n",i,paraavg.parameter[i],parastd.parameter[i],parabest.parameter[i]);
	}
	fclose(out);
	fnmlsttmp.clear();
	for(i=0;i<fnmlst.size();i++){
		sprintf(str,"%s_phigp%d",fnmlst[i].c_str(),ig);
		fnmlsttmp.push_back(str);
	}
	model_avg2(fnmlsttmp,modelall,idlst);
	printf("Read bin in all phi_group%d, model_size=%d\n",ig,modelall.size());
	sprintf(str,"%s_phigp%d",foutnm,ig);
	if((out=fopen(str,"w"))==NULL){printf("### Cannot open %s to write!!\n",foutnm);exit(0);};
	for(int j=0;j<idlst.size();j++){
		i=idlst[j];
		if(signall[i]==0){fprintf(out,"prior ");}
	        else if(signall[i]==1){
		fprintf(out,"posteria "); }
		else{printf("#### problem with the readin sign, it's not 1 or 0, but=%d",signall[i]);exit(0);}
		if((get_posteriaDist(10.,110.,5.,modelall[i].laym0.vsv,modelall[i].laym0.vsh,modelall[i].laym0.vpv,modelall[i].laym0.vph,modelall[i].laym0.eta,modelall[i].laym0.theta,modelall[i].laym0.phi,modelall[i].laym0.thick,out))==0){continue;}
		fprintf(out," g.thick %8g %8g %8g ALLmisfit %8g  Rmisfit %8g %8g %8g Lmisfit %8g %8g %8g L %8g\n",modelall[i].groups[0].thick,modelall[i].groups[1].thick,modelall[i].groups[2].thick,modelall[i].data.misfit,modelall[i].data.Rdisp.misfit,modelall[i].data.AziampRdisp.misfit,modelall[i].data.AziphiRdisp.misfit,modelall[i].data.Ldisp.misfit,modelall[i].data.AziampLdisp.misfit,modelall[i].data.AziphiLdisp.misfit,modelall[i].data.L);
	}
	fclose(out);

   	//---write out the best fitting model v(h) of this phi group---- (not para)
	idmin=idminlst[ig];
    	sprintf(str,"%s_phigp%d",fmodBnm,ig);
    	if((out=fopen(str,"w"))==NULL){printf("### Cannot open %s to write!!\n",fmodBnm);exit(0);};
    	float h=0.;
    	printf("idmin=%d,misfit=%g\n",idmin,paraall[idmin].misfit);//---check---
    	for(i=0;i<modelall[idmin].laym0.nlayer;i++){
		fprintf(out,"%8g  %8g %8g %8g %8g %8g %8g %8g\n",h,modelall[idmin].laym0.vsv[i],modelall[idmin].laym0.vsh[i],modelall[idmin].laym0.vpv[i],modelall[idmin].laym0.vph[i],modelall[idmin].laym0.eta[i],modelall[idmin].laym0.theta[i],modelall[idmin].laym0.phi[i]);
		h+=modelall[idmin].laym0.thick[i];
    	}
    	fclose(out);

   }//for ig


/*
    idmin=para_avg_V(paraall,parabest,paraavg,parastd,idlst);
    if(idmin<0){cout<<"### in para_avg, incorrect paralst.size()\n";exit(0);}
    printf("size=%d\n",idlst.size());
    if((out=fopen(fparanm,"w"))==NULL){printf("### Cannot open %s to write!!\n",fparanm);exit(0);}
    for(i=0;i<paraavg.npara;i++){ //npara, para_avg, para_avg_std, parabest;
        fprintf(out,"%5d %8.4f %8.4f %8.4f\n",i,paraavg.parameter[i],parastd[i],parabest.parameter[i]);
    }
    fclose(out);

    model_avg2(fnmlst,modelall,idlst); 
    printf("Read bin in all model_size=%d\n",modelall.size());
    if((out=fopen(foutnm,"w"))==NULL){printf("### Cannot open %s to write!!\n",foutnm);exit(0);};
    for(int j=0;j<idlst.size();j++){
	i=idlst[j];    
	if(signall[i]==0){fprintf(out,"prior ");} 
	else if(signall[i]==1){	
	fprintf(out,"posteria "); }   
	else{printf("#### problem with the readin sign, it's not 1 or 0, but=%d",signall[i]);exit(0);}

	//--- for each (accepted[if signall=1] or searched[signall=0]) model, write out the vsv,vsh value from depmin(1st para) to depmax(2ed para) at depdx(3rd) interval. 
    	if((get_posteriaDist(10.,110.,5.,modelall[i].laym0.vsv,modelall[i].laym0.vsh,modelall[i].laym0.vpv,modelall[i].laym0.vph,modelall[i].laym0.eta,modelall[i].laym0.theta,modelall[i].laym0.phi,modelall[i].laym0.thick,out))==0){continue;}
	fprintf(out," g.thick %8g %8g %8g ALLmisfit %8g  Rmisfit %8g %8g %8g Lmisfit %8g %8g %8g L %8g\n",modelall[i].groups[0].thick,modelall[i].groups[1].thick,modelall[i].groups[2].thick,modelall[i].data.misfit,modelall[i].data.Rdisp.misfit,modelall[i].data.AziampRdisp.misfit,modelall[i].data.AziphiRdisp.misfit,modelall[i].data.Ldisp.misfit,modelall[i].data.AziampLdisp.misfit,modelall[i].data.AziphiLdisp.misfit,modelall[i].data.L);
    }
    fclose(out);
	*/
    vector<modeldef>().swap(modelall);
    vector<paradef>().swap(paraall);
    vector<int>().swap(signall);
    vector<int>().swap(iiterall);
    vector<int>().swap(iaccpall);
    vector<int>().swap(idlst);
    
  }//while 1
 
  return 1;
}




