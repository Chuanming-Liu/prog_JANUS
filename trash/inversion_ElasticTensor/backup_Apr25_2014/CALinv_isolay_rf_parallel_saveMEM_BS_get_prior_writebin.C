// this rf version, the get_misfitKernel is slightly different
// right now, the isoflag is useless, since it's not used in the gen_newpara function
// the cpu time counted is strange, not sure what the problem is, but the wall time is working properly
//-------
#include <sys/time.h>
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}
//-------
int temp_writepara(FILE *f, paradef para,modeldef model, int iaccp, int ithread, int flag)
{
//
// (1)ithread (2)flag (3)iaccp  (4)'L' (5)m.d.L  (6)'RL' (7)m.d.R.L  (8)m.d.AZampR.L (9)m.d.AZphiR.L (10)'LL' (11)m.d.L.L  (12)m.d.AZampL.L  (13)m.d.AZphiL.L  (14)'misfit' (15)m.d.misfit  (16)'Rm' (17)m.d.R.misfit  (18)m.d.AZampR.misfit  (19)m.d.AZphiR.misfit  (20)'Lm' (21)m.d.L.m   (22)m.d.AZampL.m  (23)m.d.AZphiL.m	
	int i;
	fprintf(f,"%d %d %d L %8g  RL %8g %8g %8g LL %8g %8g %8g  misfit %8g Rm %8g %8g %8g Lm %8g %8g %8g\t\t",ithread,flag,iaccp,model.data.L,model.data.Rdisp.L,model.data.AziampRdisp.L,model.data.AziphiRdisp.L,model.data.Ldisp.L, model.data.AziampLdisp.L,model.data.AziphiLdisp.L,  model.data.misfit, model.data.Rdisp.misfit, model.data.AziampRdisp.misfit, model.data.AziphiRdisp.misfit,model.data.Ldisp.misfit, model.data.AziampLdisp.misfit,model.data.AziphiLdisp.misfit);
	for(i=0;i<para.npara;i++)
	{ fprintf(f,"%8g ",para.parameter[i]);
	}
	//---check---write out the Love wave 4psi signal amp, write peak and avg of the 1st 4 values
	double max=-1.;
	double sum=0.;
	int n=4;
	for(i=0;i<model.data.AziampLdisp.npper;i++){
		if(model.data.AziampLdisp.pvel[i]/model.data.Ldisp.pvel[i]*200>max)max=model.data.AziampLdisp.pvel[i]/model.data.Ldisp.pvel[i]*200;
		if(i<n)sum+=model.data.AziampLdisp.pvel[i]/model.data.Ldisp.pvel[i]*200;
	}
	fprintf(f,"\t\t %g %g",max,sum/n);
	//----------
	fprintf(f,"\n");
	return 1;
}//temp_writepara


//int do_inv_BS(int id,double misfitcri,vector<paradef> &aralst,paradef refparaBS,modeldef refmodelBS,paradef refpara,modeldef refmodel,vector<int> Rvmono,vector<int> Lvmono,vector<int> Rvgrad, vector<int> Lvgrad,vector<vector<double> > PREM, vector<vector<vector<double> > > Vkernel, vector<vector<vector<double> > > Lkernel,int k1,int k2,time_t start,int isoflag, int Rsurflag, int Lsurflag, int AziampRflag, int AziampLflag, int AziphiRflag, int AziphiLflag,int Nprem,int Rmonoc,int Lmonoc, int PosAni, vector<int> Vposani, int iitercri,int ijumpcri,char *fbinnm,float inpamp, float inpphi, int flagupdaterho){
int do_inv_BS(const int num_thread,const int id,const double misfitcri, vector<paradef> &paralst,const paradef refparaBS,const modeldef refmodelBS,const paradef refpara,const modeldef refmodel,const vector<int> Rvmono,const vector<int> Lvmono,const vector<int> Rvgrad,const vector<int> Lvgrad,const vector<vector<double> > PREM,const  vector<vector<vector<double> > > Vkernel, const vector<vector<vector<double> > > Lkernel,const int k1,const int k2,const time_t start,const int isoflag,const int Rsurflag, const int Lsurflag,const int AziampRflag,const int AziampLflag,const int AziphiRflag,const int AziphiLflag,const int Nprem,const int Rmonoc,const int Lmonoc,const int PosAni, const vector<int> Vposani, const int iitercri,const int ijumpcri,const char *fbinnm,const float inpamp, const float inpphi, const int flagupdaterho){
  int tflag,iiter,iaccp,ibad;
  double oldL,oldRL,oldLL,newL,newRL,newLL,oldmisfit,newmisfit;
  double prob,prandom;
  char str[500];
  time_t now,dtime;
  const time_t start2=time(0);
  const double start_walltime0=time(NULL);
  const double start_cputime0=clock();
  const double start_walltime=omp_get_wtime();//get_wall_time();
  const double start_cputime=get_cpu_time();
  modeldef ttmodel,tmodel,tempmod;
  paradef para1,para2,temppara,tpara;
  ofstream outbin,outbinRA;
  if ( id >1 ){
  outbin.open(fbinnm,ios_base::out|ios_base::binary);
  sprintf(str,"%s_effTI",fbinnm);
  outbinRA.open(str,ios_base::out|ios_base::binary);
  }
  //int key=1;

  //---get the phi id for crust, the group1's phi_para
  int idphiC;
  for(int i=0;i<refpara.npara;i++){
      if((int)refpara.para0[i][6]==7 and (int)refpara.para0[i][4]==1){idphiC=i;break;}
  }//for
  printf("inv, crust idphi=%d\n",idphiC);
  
  
  //const int accpcri=2000;//3000;
  const int accpcri=5000*2;//3000;
  const int Njump=ijumpcri; //set this to the times of the max_thread_num to max the efficiency of the code; Also, need to pay attetion if Njump is too small to be enough for the Monte-Carlo search. (usually, i require Njump>=5)
  int iloop=0;

  printf("num of jump=%d\n",Njump);
  iiter=iaccp=ibad=0;

  int countitt=0,countacc=0;
  //vector<paradef> tparalst;
  //-----
  if(num_thread>0){
  omp_set_dynamic(0);     // Explicitly disable dynamic teams
  omp_set_num_threads(num_thread); // Use 4 threads for all consecutive parallel regions
  }
  //while(iaccp<accpcri){
  while(iloop<10 or iaccp<accpcri){
  //while(iaccp<accpcri and iloop<4){
  //maybe, just use Njump jumps all the time. eacho jump itterate iitercri times before it terminates. No break out of the code is necessary then.
  //num_threads(Njump)
  //#pragma omp parallel for default(none) shared(start_walltime0,start_cputime0,start_walltime,start_cputime,paralst,flagupdaterho,iloop,start,k1,k2,Njump,iiter,iaccp,dtime,start2,accpcri,iitercri,misfitcri,id,key,PosAni,isoflag,inpamp,inpphi,refparaBS,refmodelBS,refpara,refmodel,Rmonoc,Lmonoc,Rvgrad,Lvgrad,Vposani,Rvmono,Lvmono,Lsurflag,Rsurflag,AziampRflag,AziampLflag,AziphiRflag,AziphiLflag,Vkernel,Lkernel,outbin,outbinRA) private(tpara,para2,para1,ttmodel,tmodel,oldLL,oldRL,oldL,countitt,countacc,ibad,tparalst,tflag,newL,newRL,newLL,newmisfit,now,prob,oldmisfit,prandom) 
  #pragma omp parallel default(none) shared(idphiC,outbin,outbinRA,paralst,iloop,iiter,iaccp) private(tpara,para2,para1,ttmodel,tmodel,oldLL,oldRL,oldL,countitt,countacc,ibad,tflag,newL,newRL,newLL,newmisfit,now,prob,oldmisfit,prandom) 
  {
  printf("threads=%d\n",omp_get_num_threads());//--check---
  #pragma omp for schedule(dynamic,1)
  for(int i=0;i<Njump;i++){
  	//--temp---
  	FILE *ftemp;
	char str[100];
	sprintf(str,"temp_LMisfit_id%d_loop%d_thread%d.txt",id,iloop,i);
	if((ftemp=fopen(str,"w"))==NULL){
		printf("Cannot open temp_LMsifit file to write! %d\n",i);
		exit(0);
	}
	paradef RApara,para1BS,para2BS;
	modeldef RAmodel,model1;
  	//-----

	int ijump=i;  
	//tparalst.clear();

	tpara=refpara;
	para1=refpara;
	para2=refpara;
	oldLL=oldRL=oldL=0.;
	countacc=0;

	para1BS=refparaBS;
	para2BS=refparaBS;

	//--for each jump, first get the starting para
	gen_newpara(refparaBS, para1BS,k1);
	/*--check--
	for(int m=0;m<para1BS.npara;m++){
		int tng=(int)para1BS.para0[m][4];
		if(tng==2){printf("para%d ppflag=%d, nv=%d %g->%g\n",m,(int)para1BS.para0[m][6],(int)para1BS.para0[m][5],refparaBS.parameter[m],para1BS.parameter[m]);}
	}
	*/
 	/*--check---
 	int j;
    	printf("\n\nbefore Bsp2P\n");
    	ttmodel=refmodelBS;
    	for(i=0;i<ttmodel.ngroup;i++){
       	 	printf("refBSmodel group%d\n",i);
 	       	for(j=0;j<ttmodel.groups[i].np;j++){
                	printf("\tvsv=%.2f vsh=%.2f vpv=%.2f vph=%.2f eta=%.2f vpvs=%.2f RAvs=%.2f RAvp=%.2f\n",ttmodel.groups[i].vsvvalue[j],ttmodel.groups[i].vshvalue[j],ttmodel.groups[i].vpvvalue[j],ttmodel.groups[i].vphvalue[j],ttmodel.groups[i].etavalue[j],ttmodel.groups[i].vpvvalue[j]/ttmodel.groups[i].vsvvalue[j],(ttmodel.groups[i].vshvalue[j]-ttmodel.groups[i].vsvvalue[j])/(ttmodel.groups[i].vshvalue[j]+ttmodel.groups[i].vsvvalue[j])*50,(ttmodel.groups[i].vphvalue[j]-ttmodel.groups[i].vpvvalue[j])/(ttmodel.groups[i].vphvalue[j]+ttmodel.groups[i].vpvvalue[j])*50);
        	}
 	}
    	
    	printf("\n\nbefore Bsp2P\n");
    	ttmodel=refmodel;
    	for(i=0;i<ttmodel.ngroup;i++){
       	 	printf("refPmodel group%d\n",i);
 	       	for(j=0;j<ttmodel.groups[i].np;j++){
                	printf("\tvsv=%.2f vsh=%.2f vpv=%.2f vph=%.2f eta=%.2f vpvs=%.2f RAvs=%.2f RAvp=%.2f\n",ttmodel.groups[i].vsvvalue[j],ttmodel.groups[i].vshvalue[j],ttmodel.groups[i].vpvvalue[j],ttmodel.groups[i].vphvalue[j],ttmodel.groups[i].etavalue[j],ttmodel.groups[i].vpvvalue[j]/ttmodel.groups[i].vsvvalue[j],(ttmodel.groups[i].vshvalue[j]-ttmodel.groups[i].vsvvalue[j])/(ttmodel.groups[i].vshvalue[j]+ttmodel.groups[i].vsvvalue[j])*50,(ttmodel.groups[i].vphvalue[j]-ttmodel.groups[i].vpvvalue[j])/(ttmodel.groups[i].vphvalue[j]+ttmodel.groups[i].vpvvalue[j])*50);
        	}
 	}
    	*/
	Bsp2Point(refmodelBS,para1BS,model1,para1,flagupdaterho); //para1BS->para1;
	para2mod(para1,model1,ttmodel);
 	/*--check---
    	printf("\n\nafter Bsp2P\n");
    	for(i=0;i<ttmodel.ngroup;i++){
       	 	printf("Pointmodel group%d\n",i);
 	       	for(j=0;j<ttmodel.groups[i].np;j++){
                	printf("\tvsv=%.2f vsh=%.2f vpv=%.2f vph=%.2f eta=%.2f vpvs=%.2f RAvs=%.2f RAvp=%.2f\n",ttmodel.groups[i].vsvvalue[j],ttmodel.groups[i].vshvalue[j],ttmodel.groups[i].vpvvalue[j],ttmodel.groups[i].vphvalue[j],ttmodel.groups[i].etavalue[j],ttmodel.groups[i].vpvvalue[j]/ttmodel.groups[i].vsvvalue[j],(ttmodel.groups[i].vshvalue[j]-ttmodel.groups[i].vsvvalue[j])/(ttmodel.groups[i].vshvalue[j]+ttmodel.groups[i].vsvvalue[j])*50,(ttmodel.groups[i].vphvalue[j]-ttmodel.groups[i].vpvvalue[j])/(ttmodel.groups[i].vphvalue[j]+ttmodel.groups[i].vpvvalue[j])*50);
        	}
 	}
    	exit(0);
    	*/
	/*--check--
	printf("\n\nafter P2M npara, BS=%d, Point=%d\n",para1BS.npara,para1.npara);
	for(int m=0;m<para1.npara;m++){
		int tng=(int)para1.para0[m][4];
		if(tng==2){printf("para%d ppflag=%d, nv=%d %g->%g\n",m,(int)para1.para0[m][6],(int)para1.para0[m][5],refpara.parameter[m],para1.parameter[m]);}
	}
	exit(0);
	*/

	//gen_newpara(refpara,para1,k1); //refpara-->para1
	//para2mod(para1,refmodel,ttmodel);// para--> model value, also, some values (e.g., vpv~eta) maybe updated according to other values (vsv,vsh)
	//mod2para(ttmodel,tpara,para1);// passing the updated values to para (those were updated in model, not in para; e.g. the vph-vpv and eta are scaled by vsh-vsv,), the process of passing back (mod-->para) for the changed para values has been added to para2mod 

	// a question: after getting a new group of para. while checking if the model satisfy some criteria (e.g. monotoniclly increasing vel), shall I use the unrotated value(gen_newpara-->para2mod-->mod2para-->goodmodel) or rotated effective TI value (gen_newpara-->Vpara2Lovepara-->Lovepara2Vpara-->para2mod-->mod2para-->goodmodel)?
	
	/*--check--@@@
	for(int m=0;m<tpara.npara;m++){
		printf("tpara%d=%g para1_=%g\n",m,tpara.parameter[m],para1.parameter[m]);
	}
	for(int m=0;m<ttmodel.ngroup;m++){
		printf("group%d, theta=%g phi=%g eta=%g\n",m,ttmodel.groups[m].thetavalue[0],ttmodel.groups[m].phivalue[0],ttmodel.groups[m].etavalue[0]);
	}
	exit(0);
	*/
	
	if(PosAni>0){//require the L always faster than R
	    //cout<<"check posAni\n";//test---
            updatemodel(ttmodel,flagupdaterho);
	    //the updatemodel, the values in the layer is the value when the ET is flat, i.e., the effect of theta is not taken into account, it's not the vel of the effective TI medium, but the vel of the un-rotated medium
            ibad=0;
	    while(positiveAni(ttmodel,Vposani)==0){
		ibad++;
		gen_newpara(para2BS, para1BS,k1);
		//para2mod(para1BS,refmodelBS,ttmodel);
		Bsp2Point(refmodelBS,para1BS,model1,para1,flagupdaterho);// para1BS->para1;
		para2mod(para1,model1,ttmodel);
		//gen_newpara(para2,para1,k1);
		//para2mod(para1,refmodel,ttmodel);
		//mod2para(ttmodel,tpara,para1);	
		updatemodel(ttmodel,flagupdaterho);
		if(ibad>10000){printf("######POSITIVE ANISO MODEL CANNOT BE SATISFIED !!!! ijump=%d\n",i);break;}
	    }	
	}//if PosAni>0
	/*--check--
	printf("posAni jump%d; ibad=%d time used walltime=%g sec; cputime=%g sec\n",i,ibad,time(NULL)-start_walltime0,(clock()-start_cputime0)/CLOCKS_PER_SEC);
	printf("posAni jump%d time used walltime=%g sec; cputime=%g sec\n",i,omp_get_wtime()-start_walltime,get_cpu_time()-start_cputime);
	*/
	
	if(Rmonoc+Lmonoc>0){
	    //cout<<"check monoc\n";//test---
	    if(PosAni<=0){updatemodel(ttmodel,flagupdaterho);}
	    ibad=0;
	    while(goodmodel(ttmodel,Rvmono,Rvgrad,Rmonoc,0,isoflag) ==0 or goodmodel(ttmodel,Lvmono,Lvgrad,0,Lmonoc,0)==0) 
		{ibad++;
		 gen_newpara(para2BS, para1BS,k1);
		 //para2mod(para1BS,refmodelBS,ttmodel);
		 Bsp2Point(refmodelBS,para1BS,model1,para1,flagupdaterho);// para1BS->para1;
		 para2mod(para1,model1,ttmodel);
		 //gen_newpara(para2,para1,k1);
		 //para2mod(para1,refmodel,ttmodel);
		 //mod2para(ttmodel,tpara,para1);
		 updatemodel(ttmodel,flagupdaterho);
		 if(ibad>10000){printf("##### ibad=%d GOOD MODEL CANNOT BE SATISFIED UNDER MONOC==1!! ijump=%d\n",ibad,ijump);break;}
		}//while	  
	}//if monoc  
	/*--check--
	printf("monoc jump%d; ibad=%d time used walltime=%g sec; cputime=%g sec\n",i,ibad,time(NULL)-start_walltime0,(clock()-start_cputime0)/CLOCKS_PER_SEC);
	printf("monoc jump%d time used walltime=%g sec; cputime=%g sec\n",i,omp_get_wtime()-start_walltime,get_cpu_time()-start_cputime);
	*/
	//printf("\n\n\ntest--begin get misfit Kernel\n");
	get_misfitKernel(ttmodel,para1,refmodel,refpara,Vkernel,Lkernel,Rsurflag,Lsurflag,AziampRflag,AziampLflag,AziphiRflag,AziphiLflag,inpamp,inpphi,flagupdaterho);
	para2=para1;
	para2BS=para1BS;
	//printf("test--finish get misfit Kernel\n\n\n");
	//printf("test-- Rnpara=%d Lnpara=%d\n",para1.Rnpara,para1.Lnpara);
	//----then begin to sample the model space----------------INVERSION--------------------------------
	for(countitt=0;countitt<iitercri;countitt++){
	  	iiter++;
		if(iiter%10000==0){
			printf("iiter=%d iaccp=%d dtime=%.1f ijump=%d iloop=%d;  ",iiter,iaccp,(float)(time(0)-start),i,iloop);
			printf("inv time used wtime=%.2fs cputime=%.2fs\n",omp_get_wtime()-start_walltime,get_cpu_time()-start_cputime);
		}
  		//---------
		gen_newpara(para2BS, para1BS,k2);
		Bsp2Point(refmodelBS,para1BS,model1,para1,flagupdaterho);// para1BS->para1;
		para2mod(para1,model1,ttmodel);
		//gen_newpara(para2,para1,k2); //in:para2 out:para1
		//para2mod(para1,refmodel,ttmodel);
		//mod2para(ttmodel,tpara,para1);
		//para2=para1;
  		//-----------------Criteria 1, goodmodel------------------------
         	 if(PosAni>0){//require the L always faster than R
	   		 //cout<<"check 2posani | ";//test---
            		updatemodel(ttmodel,flagupdaterho);
            		if(positiveAni(ttmodel,Vposani)==0){
		    		tflag=1;
				//printf("not posAni model %d\n",iiter);
		    		continue;
            		}

          	}//if PosAni>0
	
	 	if(Rmonoc+Lmonoc>0){
            		if(PosAni<=0){updatemodel(ttmodel,flagupdaterho);}
	    		if(goodmodel(ttmodel,Rvmono,Rvgrad,Rmonoc,0,isoflag)==0 or goodmodel(ttmodel,Lvmono,Lvgrad,0,Lmonoc,0)==0){
	    		tflag=2;	
			//printf("not good model %d\n",iiter);
    	    		continue;}
        	}//if monoc 

		get_misfitKernel(ttmodel,para1,refmodel,refpara,Vkernel,Lkernel,Rsurflag,Lsurflag,AziampRflag,AziampLflag,AziphiRflag,AziphiLflag,inpamp,inpphi,flagupdaterho);
		newL=ttmodel.data.L;//para1.L;
		newRL=ttmodel.data.Rdisp.L;
		newLL=ttmodel.data.Ldisp.L;
		newmisfit=para1.misfit; 

		/*---check---@@@
		if(ttmodel.groups[1].thetavalue[0]>30){
			printf("@@@-- theta=%g AZampRdisp:\n",ttmodel.groups[1].thetavalue[0]);
			for(int m=0;m<ttmodel.data.AziampRdisp.npper;m++){
				printf("%g R: %8g vs (input) %8g L: %8g vs %8g\n",ttmodel.data.AziampRdisp.pper[m],ttmodel.data.AziampRdisp.pvel[m],ttmodel.data.AziampRdisp.pvelo[m],ttmodel.data.AziampLdisp.pvel[m],ttmodel.data.AziampLdisp.pvelo[m]);
			}
 	   		printf("\n---\nmisfit: %8g\n Rmisfit: iso=%8g AZamp=%8g AZphi=%8g\nLmisfit: iso=%8g AZamp=%8g AZphi=%8g\n",ttmodel.data.misfit,ttmodel.data.Rdisp.misfit,ttmodel.data.AziampRdisp.misfit,ttmodel.data.AziphiRdisp.misfit,ttmodel.data.Ldisp.misfit,ttmodel.data.AziampLdisp.misfit,ttmodel.data.AziphiLdisp.misfit);
			exit(0);
		}
		*/

		//----------Criteria 2, acceptable L-----------------------------------
		/*
		if(misfitcri<-1){
	  	    if(newL<oldL){
	    		prob=(oldL-newL)/oldL;
	   		prandom=gen_random_unif01();
	    		if(prandom<prob) { 
		    		tflag=3;
		    		continue; }
	  	    }	

		}
		else {
			if(newmisfit>misfitcri)continue;
		}
		*/

		//----------Pass critero, accpet para----------------------------------
		para2BS=para1BS;
		para2=para1;
		oldL=newL;
		oldRL=newRL;
		oldLL=newLL;
		oldmisfit=newmisfit;
		//SAVEMEM tparalst.push_back(para1);
		//printf("test-- accp Rnpara=%d,%d Lnpara=%d,%d\n",para1.Rnpara,para1.Rparameter.size(),para1.Lnpara,para1.Lparameter.size());
		if(PosAni+Lmonoc+Rmonoc<1){// no PosAni or monoc requirement, so model hasn't been updated
			updatemodel(ttmodel,flagupdaterho);
			
		}
		countacc++;		
		//printf("write out\n");//===test
		temp_writepara(ftemp,para1,ttmodel,countacc,i,1);
		#pragma omp critical (updateflag_write_model)
		{
		  //printf("Hey, write model iaccp %d\n",iaccp);
		  //paralstBS.push_back(para1BS);
		  iaccp++;
		  if(id>1){	 
		  ///paralst.push_back(para1);
		  write_bin(ttmodel,outbin,para1,1,i,iaccp);
		  RApara=para1;
		  tmodel=ttmodel;
		  Lovepara2Vpara(RApara,tmodel);
		  para2mod_static(RApara,tmodel,RAmodel);
		  updatemodel(RAmodel,flagupdaterho);
		  
		  RApara.parameter[idphiC]=para1.parameter[idphiC];//################# this is a test ######
		  write_bin(RAmodel,outbinRA,RApara,1,i,iaccp); 
		  RApara.parameter[idphiC]=0.; //################# this is a test ######
		  }
		  else{//--this part is just a test,---check---
			///paralst.push_back(para1BS);
		  }
		  //now=time(0);
		  //dtime=dtime+now-start2;
		  //if((iaccp>1000 and (ijump>ijumpcri) and (iiter>1000000 or iaccp>accpcri)) or dtime>1800.){
			//key=-1;
			//#pragma omp flush(key,iaccp,dtime)
		  //}//if
		}//critical (updateflag)
		
	}//for count

	/* SAVE MEM
	if(id>1){
		#pragma omp critical (write_model)
		{
			for(int j=0;j<countacc;j++){
				write_bin(tmodellst[j],outbin,tparalst[j],1,i,j);
				RApara=tparalst[j];
				tmodel=tmodellst[j];
				Lovepara2Vpara(RApara,tmodel);			
				para2mod_static(RApara,tmodel,RAmodel);
				updatemodel(RAmodel,flagupdaterho);
				write_bin(RAmodel,outbinRA,RApara,1,i,j);
				//what's the layered vel i'm writting? the vel when the E.T. is flat
				//i,j--> ijump,iaccp
			}
		}//critical(write_model)
	}*/
  //----temp---
  fclose(ftemp);
  }//for i	
  }//pragma
  iloop++;
  }//while
  printf("finish do_inv iiter=%d, iaccp=%d iloop=%d\n",iiter,iaccp,iloop);
  //if(iaccp<accpcri)return 0;
  if (id>1){
  outbin.close();
  outbinRA.close();
  }
  return 1;
}
