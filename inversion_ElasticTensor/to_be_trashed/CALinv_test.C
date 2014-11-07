// this rf version, the get_misfitKernel is slightly different
// right now, the isoflag is useless, since it's not used in the gen_newpara function
//
int temp_writepara(FILE *f, paradef para,modeldef model, int iaccp, int ithread, int flag)
{
//
// (1)ithread (2)flag (3)iaccp  (4)'L' (5)m.d.L  (6)'RL' (7)m.d.R.L  (8)m.d.AZampR.L (9)m.d.AZphiR.L (10)'LL' (11)m.d.L.L  (12)m.d.AZampL.L  (13)m.d.AZphiL.L  (14)'misfit' (15)m.d.misfit  (16)'Rm' (17)m.d.R.misfit  (18)m.d.AZampR.misfit  (19)m.d.AZphiR.misfit  (20)'Lm' (21)m.d.L.m   (22)m.d.AZampL.m  (23)m.d.AZphiL.m	
	int i;
	fprintf(f,"%d %d %d L %8g  RL %8g %8g %8g LL %8g %8g %8g  misfit %8g Rm %8g %8g %8g Lm %8g %8g %8g\t\t",ithread,flag,iaccp,model.data.L,model.data.Rdisp.L,model.data.AziampRdisp.L,model.data.AziphiRdisp.L,model.data.Ldisp.L, model.data.AziampLdisp.L,model.data.AziphiLdisp.L,  model.data.misfit, model.data.Rdisp.misfit, model.data.AziampRdisp.misfit, model.data.AziphiRdisp.misfit,model.data.Ldisp.misfit, model.data.AziampLdisp.misfit,model.data.AziphiLdisp.misfit);
	for(i=0;i<para.npara;i++)
	{ fprintf(f,"%8g ",para.parameter[i]);
	}
	fprintf(f,"\n");
	return 1;
}//temp_writepara


int do_inv(int id,double misfitcri,vector<paradef> &paralst,paradef refpara,modeldef refmodel,vector<int> Rvmono,vector<int> Lvmono,vector<int> Rvgrad, vector<int> Lvgrad,vector<vector<double> > PREM, vector<vector<vector<double> > > Vkernel, vector<vector<vector<double> > > Lkernel,int k1,int k2,time_t start,int isoflag, int Rsurflag, int Lsurflag, int AziampRflag, int AziampLflag, int AziphiRflag, int AziphiLflag,int Nprem,int Rmonoc,int Lmonoc, int PosAni, vector<int> Vposani, int iitercri,int ijumpcri,char *fbinnm,float inpamp, float inpphi, int flagupdaterho){
  int tflag,key,iiter,iaccp,ibad;
  double oldL,oldRL,oldLL,newL,newRL,newLL,oldmisfit,newmisfit;
  double prob,prandom;
  char str[500];
  time_t now,dtime,start2;
  modeldef ttmodel,tmodel,tempmod;
  paradef para1,para2,temppara,tpara;
  ofstream outbin(fbinnm,ios_base::out|ios_base::binary);
  sprintf(str,"%s_effTI",fbinnm);
  ofstream outbinRA(str,ios_base::out|ios_base::binary);

  int accpcri=2000;//3000;
  int Njump=ijumpcri; //set this to the times of the max_thread_num to max the efficiency of the code; Also, need to pay attetion if Njump is too small to be enough for the Monte-Carlo search. (usually, i require Njump>=5)
  int iloop=0;

  printf("num of threads=%d\n",Njump);
  key=1;iiter=iaccp=ibad=0;
//  oldL=refmodel.data.L;
//  oldmisfit=refmodel.data.misfit;



  start2=time(0);
  int countitt=0,countacc=0;
  vector<paradef> tparalst;
  vector<modeldef> tmodellst;
  //-----
  while(iaccp<accpcri and iloop<2){
  //maybe, just use Njump jumps all the time. eacho jump itterate iitercri times before it terminates. No break out of the code is necessary then.
  //num_threads(Njump)
  #pragma omp parallel for default(none) shared(flagupdaterho,iloop,start,k1,k2,Njump,iiter,iaccp,dtime,start2,accpcri,iitercri,misfitcri,id,key,PosAni,isoflag,inpamp,inpphi,refpara,refmodel,Rmonoc,Lmonoc,Rvgrad,Lvgrad,Vposani,Rvmono,Lvmono,Lsurflag,Rsurflag,AziampRflag,AziampLflag,AziphiRflag,AziphiLflag,Vkernel,Lkernel,outbin,outbinRA) private(tpara,para2,para1,ttmodel,tmodel,oldLL,oldRL,oldL,countitt,countacc,ibad,tparalst,tmodellst,tflag,newL,newRL,newLL,newmisfit,now,prob,oldmisfit,prandom) 
  for(int i=0;i<Njump;i++){
  	//--temp---
  	FILE *ftemp;
	char str[100];
	sprintf(str,"temp_LMisfit_id%d_loop%d_thread%d.txt",id,iloop,i);
	if((ftemp=fopen(str,"w"))==NULL){
		printf("Cannot open temp_LMsifit file to write! %d\n",i);
		exit(0);
	}
	paradef RApara;
	modeldef RAmodel;
  	//-----

	int ijump=i;  
	tparalst.clear();
	tmodellst.clear();

	tpara=refpara;
	para1=refpara;
	para2=refpara;
	oldLL=oldRL=oldL=0.;
	countacc=0;

	//--for each jump, first get the starting para
	gen_newpara(refpara,para1,k1); //refpara-->para1
	para2mod(para1,refmodel,ttmodel);// para--> model value, also, some values (e.g., vpv~eta) maybe updated according to other values (vsv,vsh)
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
		gen_newpara(para2,para1,k1);
		para2mod(para1,refmodel,ttmodel);
		//mod2para(ttmodel,tpara,para1);	
		updatemodel(ttmodel,flagupdaterho);
		if(ibad>10000){printf("######POSITIVE ANISO MODEL CANNOT BE SATISFIED !!!! ijump=%d\n",i);break;}
	    }	
	}//if PosAni>0
	
	if(Rmonoc+Lmonoc>0){
	    //cout<<"check monoc\n";//test---
	    if(PosAni<=0){updatemodel(ttmodel,flagupdaterho);}
	    ibad=0;
	    while(goodmodel(ttmodel,Rvmono,Rvgrad,Rmonoc,0,isoflag) ==0 or goodmodel(ttmodel,Lvmono,Lvgrad,0,Lmonoc,0)==0) 
		{ibad++;
		 gen_newpara(para2,para1,k1);
		 para2mod(para1,refmodel,ttmodel);
		 //mod2para(ttmodel,tpara,para1);
		 updatemodel(ttmodel,flagupdaterho);
		 if(ibad>10000){printf("##### ibad=%d GOOD MODEL CANNOT BE SATISFIED UNDER MONOC==1!! ijump=%d\n",ibad,ijump);break;}
		}//while	  
	}//if monoc  
	//printf("\n\n\ntest--begin get misfit Kernel\n");
	get_misfitKernel(ttmodel,para1,refmodel,refpara,Vkernel,Lkernel,Rsurflag,Lsurflag,AziampRflag,AziampLflag,AziphiRflag,AziphiLflag,inpamp,inpphi,flagupdaterho);
	para2=para1;
	//printf("test--finish get misfit Kernel\n\n\n");
	//printf("test-- Rnpara=%d Lnpara=%d\n",para1.Rnpara,para1.Lnpara);
	//----then begin to sample the model space----------------INVERSION--------------------------------
	for(countitt=0;countitt<iitercri;countitt++){
	  	iiter++;
		if(iiter%2000==0){
			printf("iiter=%d iaccp=%d dtime=%f ijump=%d iloop=%d\n",iiter,iaccp,(float)(time(0)-start),i,iloop);
		}
  		//---------
		gen_newpara(para2,para1,k2); //in:para2 out:para1
		para2mod(para1,refmodel,ttmodel);
		//mod2para(ttmodel,tpara,para1);
		//para2=para1;
  		//-----------------Criteria 1, goodmodel------------------------
         	 if(PosAni>0){//require the L always faster than R
	   		 //cout<<"check 2posani | ";//test---
            		updatemodel(ttmodel,flagupdaterho);
            		if(positiveAni(ttmodel,Vposani)==0){
		    		tflag=1;
		    		continue;
            		}

          	}//if PosAni>0
	
	 	if(Rmonoc+Lmonoc>0){
            		if(PosAni<=0){updatemodel(ttmodel,flagupdaterho);}
	    		if(goodmodel(ttmodel,Rvmono,Rvgrad,Rmonoc,0,isoflag)==0 or goodmodel(ttmodel,Lvmono,Lvgrad,0,Lmonoc,0)==0){
	    		tflag=2;	
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

		//----------Pass critero, accpet para----------------------------------

		para2=para1;
		oldL=newL;
		oldRL=newRL;
		oldLL=newLL;
		oldmisfit=newmisfit;
		tparalst.push_back(para1);
		//printf("test-- accp Rnpara=%d,%d Lnpara=%d,%d\n",para1.Rnpara,para1.Rparameter.size(),para1.Lnpara,para1.Lparameter.size());
		if(PosAni+Lmonoc+Rmonoc<1){// no PosAni or monoc requirement, so model hasn't been updated
			updatemodel(ttmodel,flagupdaterho);
			
		}
		tmodellst.push_back(ttmodel);
		countacc++;		

		temp_writepara(ftemp,para1,ttmodel,countacc,i,1);
		#pragma omp critical (updateflag)
		{

		  paralst.push_back(para1);
		  iaccp++;
		  //now=time(0);
		  //dtime=dtime+now-start2;
		  //if((iaccp>1000 and (ijump>ijumpcri) and (iiter>1000000 or iaccp>accpcri)) or dtime>1800.){
			//key=-1;
			//#pragma omp flush(key,iaccp,dtime)
		  //}//if
		}//critical (updateflag)
		
	}//for count

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
	}
  //----temp---
  fclose(ftemp);
  }//for i	
  iloop++;
  }
  printf("finish do_inv iiter=%d, iaccp=%d iloop=%d\n",iiter,iaccp,iloop);
  //if(iaccp<accpcri)return 0;
  return 1;
  outbin.close();
  outbinRA.close();
}
