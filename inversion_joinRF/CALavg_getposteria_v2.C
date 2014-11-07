// this code reads the binary file, and compute the average model.
// the binary file contains all the accepted models, it was written during the inversion process.
// this version write out the vsv vsh result for both step1 and step2 do_inv. Also, write out prior and posteria distribution 
// attentin, this version is for the smooth model, the thickness of the last layer is 0 for smooth model, so in the subroutine model_avg_sub, (j=0;j<modlst[i].laym0.nlayer-1 [NOT nlayer!!!];j++)
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
/*
#include"./CALgroup_smooth.C"
#include"./CALgroup_rf.C"
#include"./CALmodel_LVZ_rf.C"
#include"./CALpara_isolay_v2.C""
#include "./CALrf.C"
#include "./CALMineos_readK_smooth_rf_dpv2_parallel.C"
#include "./ASC_rw_rf.C"
#include "./CALinv_isolay_rf_parallel.C"
*/

using namespace std;

	int para_avg(vector<paradef> &paralst, paradef &paraavg,vector<double> &Rparastd,vector<double> &Lparastd, vector<int> &idlst){
	//get the average parameters from a list of acceptable paras
	    int i,j,k,size,Ngood;
	    double mismin=1e10;
	  
	    size=paralst.size();
	    if(size<1) return 0;
	    paraavg=paralst[0];	
	    Rparastd.clear();Lparastd.clear();
	    initpara(paraavg);
	    paraavg.Rpara0=paralst[0].Rpara0;
	    paraavg.Lpara0=paralst[0].Lpara0;
	    size=paralst.size();
	    if(size<1) return 0;

	    for(j=0;j<paralst[0].Rnpara;j++){//Rnpara==Lnpara
		paraavg.Rparameter.push_back(0.);
		paraavg.Lparameter.push_back(0.);
		Rparastd.push_back(0.);
		Lparastd.push_back(0.);
	    }//for j


	    for(j=0;j<paralst[0].Rnpara;j++){//Rnpara==Lnpara
                Rparastd.push_back(0.);
                Lparastd.push_back(0.);
		paraavg.Rparameter[j]=paraavg.Lparameter[j]=0.;
            }//for j

	    for(i=0;i<size;i++)
		if(paralst[i].misfit<mismin) mismin=paralst[i].misfit;
	    //min(2*Mis_min,Mis_min+0.5)
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
	  float tvsh,tvsv,tani,ani,animin,animax;
	  float fm3;
	  char str[500];
	  int i,j,Ncount,size,k;
	  vector<vector<float> > tvlst;

	  hlst.clear();vlst.clear();stdlst.clear();tvlst.clear();
		  
	  size=modlst.size();
	
	  for(i=0;i<size;i++){	// obtain the end depth for this mod_avg
	    dep2=0.;  
	    for(k=0;k<=ng;k++){dep2=dep2+modlst[i].groups[k].thick;}
	    dep2lst.push_back(dep2);
	    dep1lst.push_back(0.);
	    iddep1lst.push_back(0);
	  }
	  for(th=h0;th<h;th=th+dth){//thickness
		//printf("th=%g h=%g dep1[0]=%g dep2[0]=%g\n",th,h,dep1lst[0],dep2lst[0]);//---test---
		vsv=0.;vsh=0.;Ncount=0;fm1=0.;fm2=0.;ani=0.;fm3=0.;tth=0.;
		tvlst.clear();
		vsvmin=1e10;vsvmax=-1;
		vshmin=1e10;vshmax=-1;
		animin=1e10;animax=-1;
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
				if(tvsv<0 or tvsv>5){
					printf("test--- th=%g, imod=%d, ilay=%d, tvsv=%g=[(%g-%g)/(%g)+%g]\n",th,i,j,tvsv,modlst[i].laym0.vsv[j+1],modlst[i].laym0.vsv[j],modlst[i].laym0.thick[j],modlst[i].laym0.vsv[j]);}
				tani=100*(tvsh-tvsv)/(sqrt((2*tvsv*tvsv+tvsh*tvsh)/3.0));
				vsv=vsv+tvsv;// modified on Aug27, 2012
				vsh=vsh+tvsh;//
				ani=ani+tani;
				Ncount++;	
				if (tvsv>vsvmax)vsvmax=tvsv;
				if(tvsv<vsvmin)vsvmin=tvsv;
				if(tvsh>vshmax)vshmax=tvsh;
				if(tvsh<vshmin)vshmin=tvsh;
				if(tani>animax)animax=tani;if(tani<animin)animin=tani;
				tv.clear();tv.push_back(tvsv);tv.push_back(tvsh);tv.push_back(tani);
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
		if (isnan(vsv) or isnan(vsh) or isnan(ani)){
			printf("Hey, NaN happen!!! something wrong!!\n Ncount=%d th=%g tth=%g h0=%g h=%g\n",Ncount,th,tth,h0,h);
			sprintf(str,"echo WRONG Ncount = %d >> point_finished_Ani.txt ",Ncount);
			system(str);
			//exit(0);
		}
		tv.clear();
		tv.push_back(vsv);tv.push_back(vsh);tv.push_back(vsvmin);tv.push_back(vsvmax);tv.push_back(vshmin);tv.push_back(vshmax);tv.push_back(ani);tv.push_back(animin);tv.push_back(animax);
		vlst.push_back(tv);
		hlst.push_back(th);
		for(i=0;i<Ncount;i++){ 
			fm1=fm1+pow(tvlst[i][0]-vsv,2);
			fm2=fm2+pow(tvlst[i][1]-vsh,2);//----there was a bug here, fixed on May 18, 2012
			fm3=fm3+pow(tvlst[i][2]-ani,2);
		}
		fm1=sqrt(fm1/Ncount);
		fm2=sqrt(fm2/Ncount);
		fm3=sqrt(fm3/Ncount);
		tv.clear();
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
 	int integrate_mod(vector<vector<float> > &outVsv,vector<vector<float> > &outVsh,vector<vector<float> > &outAni, vector<float> hlst,vector<vector<float> > vlst,vector<vector<float> > stdlst)
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
//-----------------------------------------------------------------
	int connect_modmin(vector<vector<float> > &outVsv,vector<vector<float> > &outVsh,vector<vector<float> > &outAni,vector<float> hlst,vector<vector<float> > vlst,vector<vector<float> > stdlst,float h1,float h2)
	{
	  int i;vector<float> tv;
	  float h,v;
 	  for(i=0;i<vlst.size();i++){
	  	h=hlst[i];
		if (h<h1)continue;
		else if (h>=h2)break;
		tv.clear();
		v=vlst[i][0]-stdlst[i][0];
		tv.push_back(h);tv.push_back(v);
		outVsv.push_back(tv);
		tv.clear();
		v=vlst[i][1]-stdlst[i][1];
		tv.push_back(h);tv.push_back(v);
		outVsh.push_back(tv);
		tv.clear();
		v=vlst[i][6]-stdlst[i][2];
		tv.push_back(h);tv.push_back(v);
		outAni.push_back(tv);
	  }
	  vector<float>().swap(tv);
	  return 1;
	}//connect_modmin
//-----------------------------------------------------------------
	int connect_modmid(vector<vector<float> > &outVsv,vector<vector<float> > &outVsh,vector<vector<float> > &outAni,vector<float> &hlst,vector<vector<float> > &vlst,vector<vector<float> > &stdlst,float h1,float h2)
	{// modified on Aug 28, 2012. add the std into vector, and output it. previous version has some problems (the hlst for Vmin, Vmid, Vmax could be different, so when read the output file, obtaining unc using Vid-Vmin is incorrect!!)
	  int i;vector<float> tv;
	  float h,v;
 	  for(i=0;i<vlst.size();i++){
	  	h=hlst[i];
		if (h<h1)continue;
		else if (h>=h2)break;
		tv.clear();
		v=vlst[i][0];
		tv.push_back(h);tv.push_back(v);tv.push_back(stdlst[i][0]);
		outVsv.push_back(tv);
		tv.clear();
		v=vlst[i][1];
		tv.push_back(h);tv.push_back(v);tv.push_back(stdlst[i][1]);
		outVsh.push_back(tv);
		tv.clear();
		v=vlst[i][6];
		tv.push_back(h);tv.push_back(v);tv.push_back(stdlst[i][2]);
		outAni.push_back(tv);
	  }
	  vector<float>().swap(tv);
	  return 1;
	}//connect_modmin
//-----------------------------------------------------------------
//-----------------------------------------------------------------
	int connect_modmax(vector<vector<float> > &outVsv,vector<vector<float> > &outVsh,vector<vector<float> > &outAni,vector<float> &hlst,vector<vector<float> > &vlst,vector<vector<float> > &stdlst,float h1,float h2)
	{
	  int i;vector<float> tv;
	  float h,v;
 	  for(i=0;i<vlst.size();i++){
	  	h=hlst[i];
		if (h<h1)continue;
		else if (h>=h2) break;
		tv.clear();
		v=vlst[i][0]+stdlst[i][0];
		tv.push_back(h);tv.push_back(v);
		outVsv.push_back(tv);
		tv.clear();
		v=vlst[i][1]+stdlst[i][1];
		tv.push_back(h);tv.push_back(v);
		outVsh.push_back(tv);
		tv.clear();
		v=vlst[i][6]+stdlst[i][2];
		tv.push_back(h);tv.push_back(v);
		outAni.push_back(tv);
	  }
	  vector<float>().swap(tv);
	  return 1;
	}//connect_modmin
//-----------------------------------------------------------------
//-----------------------------------------------------------------
	int write_modavg( const char *fvsvnm, const char *fvshnm, char *faninm,vector<vector<float> > &Vsv1,vector<vector<float> > &Vsv2,vector<vector<float> > &Vsv3,vector<vector<float> > &Vsh1,vector<vector<float> > &Vsh2,vector<vector<float> > &Vsh3, vector<vector<float> > &Ani1,vector<vector<float> > &Ani2,vector<vector<float> > &Ani3,float Hsed,float Hsedstd,float Hmoho,float Hmohostd )
	{// attention, the hlst for Vmin, Vmid, Vmax could be different in many cases. So, use the associated hlst while plotting/reading ...
// this writting method has some problems, since the size of Vmin,Vmid,Vmax is different, I only write to the min(len(V)). In this case, the V with the largest size cannot be fully written out (cannot reach the maximum depth).
 
	    FILE *ftmp;
	    int i,j,sizemin;
	    if((ftmp=fopen(fvsvnm,"w"))==NULL){printf("####in write-modavg, canot open file %s to write!\n",fvsvnm);exit(0);};
	    fprintf(ftmp,"sedi %8g %8g\n",Hsed,Hsedstd);
	    fprintf(ftmp,"moho %8g %8g\n",Hmoho,Hmohostd);
	    sizemin=min(Vsv1.size(),Vsv2.size());
	    if(sizemin>Vsv3.size())sizemin=Vsv3.size();
	    for(i=0;i<sizemin;i++){
	       //fprintf(ftmp,"%8g %8g %8g %8g\n",Vsv1[i][0],Vsv1[i][1],Vsv2[i][1],Vsv3[i][1]);
	       fprintf(ftmp,"%8g %8g %8g %8g %8g %8g %8g\n",Vsv1[i][0],Vsv1[i][1],Vsv2[i][0],Vsv2[i][1],Vsv3[i][0],Vsv3[i][1],Vsv2[i][2]);
	    }
	    fclose(ftmp);
	    if((ftmp=fopen(fvshnm,"w"))==NULL){printf("####in write-modavg, canot open file %s to write!\n",fvshnm);exit(0);};
	    fprintf(ftmp,"sedi %8g %8g\n",Hsed,Hsedstd);
	    fprintf(ftmp,"moho %8g %8g\n",Hmoho,Hmohostd);
	    for(i=0;i<sizemin;i++){
	       fprintf(ftmp,"%8g %8g %8g %8g %8g %8g %8g\n",Vsh1[i][0],Vsh1[i][1],Vsh2[i][0],Vsh2[i][1],Vsh3[i][0],Vsh3[i][1],Vsh2[i][2]);
	    }
	    fclose(ftmp);
	    if((ftmp=fopen(faninm,"w"))==NULL){printf("####in write-modavg, canot open file %s to write!\n",faninm);exit(0);};
	    fprintf(ftmp,"sedi %8g %8g\n",Hsed,Hsedstd);
	    fprintf(ftmp,"moho %8g %8g\n",Hmoho,Hmohostd);
	    for(i=0;i<sizemin;i++){
		    fprintf(ftmp,"%8g %8g %8g %8g %8g %8g %8g\n",Ani1[i][0],Ani1[i][1],Ani2[i][0],Ani2[i][1],Ani3[i][0],Ani3[i][1],Ani2[i][2]);
	    }		    
	    fclose(ftmp);
	    return 1;
	}//write_modavg
//--------------------------------------------------
       int model_avg2( const char *fvsvnm, const char *fvshnm,char *faninm, vector<modeldef> &modlstall, vector<int> &idlst)
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
	  cout<<"---begin! "<<time(0)<<endl;start=time(0);
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
	  model_avg_sub(v1lst,std1lst,h1lst,modlst,0,0,Hsed+Hsedstd,0.2);
	  cout<<"finish sed\n";//--test
	  model_avg_sub(v2lst,std2lst,h2lst,modlst,1,Hsed-Hsedstd,Hmoho+Hmohostd,1.);
	  cout<<"finish cst\n";//---test
	  model_avg_sub(v3lst,std3lst,h3lst,modlst,2,Hmoho-Hmohostd,depmax,5.);
	  cout<<"finish mant\n";//---test
	  
	  cout<<"---connect model min/mid/max time_used="<<time(0)-start<<"\n";//---test---
	  connect_modmin(Vsvmin,Vshmin,Animin,h1lst,v1lst,std1lst,0,Hsed+Hsedstd);
	  connect_modmin(Vsvmin,Vshmin,Animin,h2lst,v2lst,std2lst,Hsed+Hsedstd,Hmoho+Hmohostd);
	  connect_modmin(Vsvmin,Vshmin,Animin,h3lst,v3lst,std3lst,Hmoho+Hmohostd,depmax);
	  
	  connect_modmid(Vsvmid,Vshmid,Animid,h1lst,v1lst,std1lst,0,Hsed);
	  connect_modmid(Vsvmid,Vshmid,Animid,h2lst,v2lst,std2lst,Hsed,Hmoho);
	  connect_modmid(Vsvmid,Vshmid,Animid,h3lst,v3lst,std3lst,Hmoho,depmax);
	  
	  connect_modmax(Vsvmax,Vshmax,Animax,h1lst,v1lst,std1lst,0,Hsed-Hsedstd);
	  connect_modmax(Vsvmax,Vshmax,Animax,h2lst,v2lst,std2lst,Hsed-Hsedstd,Hmoho-Hmohostd);
	  connect_modmax(Vsvmax,Vshmax,Animax,h3lst,v3lst,std3lst,Hmoho-Hmohostd,depmax);
	  
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
	  write_modavg(fvsvnm,fvshnm,faninm,Vsvmin,Vsvmid,Vsvmax,Vshmin,Vshmid,Vshmax,Animin,Animid,Animax,Hsed,Hsedstd,Hmoho,Hmohostd);

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
int get_posteriaDist(double depmin, double depmax, double depstep, vector<double> &vsv, vector<double> &vsh, vector<double> &thick, FILE *out)
{
  int i,size;
  double dep,h;
  vector<double> hlst;
  double tvsv,tvsh;

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
		fprintf(out,"%8g %8g ",tvsv,tvsh);break;}
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
  vector<double> Rparastd,Lparastd;
  vector<int> signall,iiterall,iaccpall,idlst;
  char inponm[300],dirbin[300],fbname[200],foutnm[300],fvsvnm[300],fvshnm[300],nodeid[10],outdir[300],str[200];
  char fbname1[300],fvsvnm1[300],fvshnm1[300],foutnm1[300];
  char faninm[300];
  FILE *inpo, *out;
  int npoint,i;
  float tvsv,tvsh,tani;
  vector<float> anilst;	  
  paradef paraavg;
  float lon,lat;  

  if(argc!=4){
    printf("Usage: XX 1]point_file 2]bin_dir 3]out_dir\n1]point_file: nodeid lon lat\n2]bin_dir: bindir/ani(vsv)_nodeid_lon_lat.bin\n3]out_dir: outdir/vsv(vsh)_lon_lat.txt\n");
    exit(0);
  }
  sprintf(inponm,argv[1]);
  sprintf(dirbin,argv[2]);
  sprintf(outdir,argv[3]);

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
    sprintf(fbname,"%s/ani_%s_%.1f_%.1f.bin",dirbin,nodeid,lon,lat);
    sprintf(foutnm,"%s/post_%.1f_%.1f.txt",outdir,lon,lat);
   sprintf(fvsvnm,"%s/vsv_%.1f_%.1f.txt",outdir,lon,lat);
    sprintf(fvshnm,"%s/vsh_%.1f_%.1f.txt",outdir,lon,lat);
    sprintf(faninm,"%s/ani_%.1f_%.1f.txt",outdir,lon,lat);

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
    if((para_avg(paraall,paraavg,Rparastd,Lparastd,idlst))==0){cout<<"### in para_avg, incorrect paralst.size()\n";exit(0);}
    printf("size=%d\n",idlst.size());
    model_avg2(fvsvnm,fvshnm,faninm,modelall,idlst); 

    printf("Read bin in all model_size=%d\n",modelall.size());
    if((out=fopen(foutnm,"w"))==NULL){printf("### Cannot open %s to write!!\n",foutnm);exit(0);};
    for(int j=0;j<idlst.size();j++){
	i=idlst[j];    
	if(signall[i]==0){fprintf(out,"prior ");} 
	else if(signall[i]==1){	
	fprintf(out,"posteria "); }   
	else{printf("#### problem with the readin sign, it's not 1 or 0, but=%d",signall[i]);exit(0);}

    	if((get_posteriaDist(10.,110.,5.,modelall[i].laym0.vsv,modelall[i].laym0.vsh,modelall[i].laym0.thick,out))==0){continue;}
	fprintf(out,"%8g %8g %8g %8g %8g %8g\n",modelall[i].groups[0].thick,modelall[i].groups[1].thick,modelall[i].data.misfit,modelall[i].data.Rdisp.misfit,modelall[i].data.Ldisp.misfit,modelall[i].data.L);
    }
    fclose(out);
/* ************************************************************8
    //------------add the step1 distribution information -------------------
    printf("---begin the step1 info---\n");
    modelall.clear();paraall.clear();
    signall.clear();iiterall.clear();iaccpall.clear();Rparastd.clear();Lparastd.clear();idlst.clear();
    initpara(paraavg);
    
    ifstream mff1(fbname1);
    if(! mff1.good()){
    	printf("Cannot access file %s\n",fbname1);
	sprintf(str,"echo %s %.1f %.1f >> points_CannotAccess.txt",nodeid,lon,lat);
	continue;
    }
    mff1.close();
    
    read_bin(modelall,paraall,fbname1,signall,iiterall,iaccpall);
    printf("paraall.size()=%d fbname1=%s\n",paraall.size(),fbname1);
    if((para_avg(paraall,paraavg,Rparastd,Lparastd,idlst))==0){cout<<"### in para_avg_step1, incorrect paralst.size()\n";exit(0);}
    printf("size_step1=%d\n",idlst.size());
    model_avg2(fvsvnm1,fvshnm1,faninm1,modelall,idlst);
    printf("step1 Read bin in all model_size=%d\n",modelall.size());
    if((out=fopen(foutnm1,"w"))==NULL){printf("### Cannot open %s to write!!\n",foutnm1);exit(0);}
    for(i=0;i<modelall.size();i++){

        if(signall[i]==0){fprintf(out,"prior ");} 
	else if(signall[i]==1){	
	fprintf(out,"posteria "); }   
	else{printf("#### problem with the readin sign, it's not 1 or 0, but=%d",signall[i]);exit(0);}

	if((get_posteriaDist(10.,90.,10.,modelall[i].laym0.vsv,modelall[i].laym0.vsh,modelall[i].laym0.thick,out))==0){continue;}
	fprintf(out,"%8g %8g %8g %8g %8g %8g\n",modelall[i].groups[0].thick,modelall[i].groups[1].thick,modelall[i].data.misfit,modelall[i].data.Rdisp.misfit,modelall[i].data.Ldisp.misfit,modelall[i].data.L);
    }
    fclose(out);

*/  

    vector<modeldef>().swap(modelall);
    vector<paradef>().swap(paraall);
    vector<double>().swap(Rparastd);
    vector<double>().swap(Lparastd);
    vector<int>().swap(signall);
    vector<int>().swap(iiterall);
    vector<int>().swap(iaccpall);
    vector<int>().swap(idlst);
    
  }//while 1
 
  return 1;
}




