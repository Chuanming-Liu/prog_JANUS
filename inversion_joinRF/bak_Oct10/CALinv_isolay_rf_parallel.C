// this rf version, the get_misfitKernel is slightly different
//
int temp_writepara(FILE *f, paradef para,modeldef model, int iaccp, int ithread, int flag)
{
	int i;
	fprintf(f,"%d %d %d L %8g %8g %8g  misfit %8g %8g %8g\t",ithread,flag,iaccp,model.data.L,model.data.Rdisp.L,model.data.rf.L,model.data.misfit,model.data.Rdisp.misfit,model.data.rf.misfit);
	for(i=0;i<para.Rnpara;i++)
	{ fprintf(f,"%8g ",para.Rparameter[i]);
	}
	fprintf(f,"\n");
	return 1;
}//temp_writepara


int do_inv(int id,double misfitcri,vector<paradef> &paralst,paradef refpara,modeldef refmodel,vector<int> Rvmono,vector<int> Lvmono,vector<int> Rvgrad, vector<int> Lvgrad,vector<vector<double> > PREM, vector<vector<vector<double> > > kernel,int k1,int k2,time_t start,int isoflag, int Rsurflag, int Lsurflag,int Nprem,float depcri1,float depcri2,float qscri,float qpcri,int Rmonoc,int Lmonoc, int PosAni, vector<int> Vposani,vector<float> &Viso, int iitercri,int ijumpcri,char *fbinnm,float inp){
  int tflag,key,iiter,iaccp,ibad;
  double oldL,oldRL,oldLL,newL,newRL,newLL,oldmisfit,newmisfit;
  double prob,prandom;
  time_t now,dtime,start2;
  modeldef ttmodel,tempmod;
  paradef para1,para2,temppara;
  ofstream outbin(fbinnm,ios_base::out|ios_base::binary);
  int accpcri=3000;
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
  while(iaccp<accpcri and iloop<5){
  //maybe, just use Njump jumps all the time. eacho jump itterate iitercri times before it terminates. No break out of the code is necessary then.
  //num_threads(Njump)
  #pragma omp parallel for default(none) shared(paralst,Viso,iloop,start,k1,k2,Njump,iiter,iaccp,dtime,start2,accpcri,iitercri,misfitcri,id,key,PosAni,isoflag,depcri1,depcri2,qpcri,qscri,inp,refpara,refmodel,Rmonoc,Lmonoc,Rvgrad,Lvgrad,Vposani,Rvmono,Lvmono,Lsurflag,Rsurflag,kernel,outbin) private(para2,para1,ttmodel,oldLL,oldRL,oldL,countitt,countacc,ibad,tparalst,tmodellst,tflag,newL,newRL,newLL,newmisfit,now,prob,oldmisfit,prandom) 
  for(int i=0;i<Njump;i++){
  	//--temp---
  	FILE *ftemp;
	char str[100];
	sprintf(str,"temp_LMisfit_id%d_loop%d_thread%d.txt",id,iloop,i);
	if((ftemp=fopen(str,"w"))==NULL){
		printf("Cannot open temp_LMsifit file to write! %d\n",i);
		exit(0);
	}
  	//-----

	int ijump=i;  
	tparalst.clear();
	tmodellst.clear();

	para1=refpara;
	para2=refpara;
	oldLL=oldRL=oldL=0.;
	countacc=0;

	//--for each jump, first get the starting para
	gen_newpara(refpara,para1,k1,isoflag,Viso); //refpara-->para1
	para2mod(para1,refmodel,ttmodel);
	if(PosAni>0){//require the L always faster than R
	    //cout<<"check posAni\n";//test---
            updatemodelTibet(ttmodel,depcri1,depcri2,qpcri,qscri);
            ibad=0;
	    while(positiveAni(ttmodel,Vposani)==0){
		ibad++;
		gen_newpara(para2,para1,k1,isoflag,Viso);
		para2mod(para1,refmodel,ttmodel);
		updatemodelTibet(ttmodel,depcri1,depcri2,qpcri,qscri);
		if(ibad>10000){printf("######POSITIVE ANISO MODEL CANNOT BE SATISFIED !!!! ijump=%d\n",i);break;}
	    }	
	}//if PosAni>0
	
	if(Rmonoc+Lmonoc>0){
	    //cout<<"check monoc\n";//test---
	    if(PosAni<=0){updatemodelTibet(ttmodel,depcri1,depcri2,qpcri,qscri);}
	    ibad=0;
	    while(goodmodel(ttmodel,Rvmono,Rvgrad,Rmonoc,0,isoflag)*goodmodel(ttmodel,Lvmono,Lvgrad,0,Lmonoc,0)==0)
		{ibad++;
		 gen_newpara(para2,para1,k1,isoflag,Viso);
		 para2mod(para1,refmodel,ttmodel);
		 updatemodelTibet(ttmodel,depcri1,depcri2,qpcri,qscri);
		 if(ibad>10000){printf("##### ibad=%d GOOD MODEL CANNOT BE SATISFIED UNDER MONOC==1!! ijump=%d\n",ibad,ijump);break;}
		}//while	  
	}//if monoc  
	// compute_dispKernel(model1,para1,kernel);
	//printf("test-- get misfit Kernel\n");
	get_misfitKernel(ttmodel,para1,refmodel,refpara,kernel,Rsurflag,Lsurflag,depcri1,depcri2,qpcri,qscri,inp);	  
	para2=para1;
	//printf("test--finish get misfit Kernel\n");
	//printf("test-- Rnpara=%d Lnpara=%d\n",para1.Rnpara,para1.Lnpara);
	//----then begin to sample the model space----------------INVERSION--------------------------------
	for(countitt=0;countitt<iitercri;countitt++){
	  	iiter++;
		if(iiter%1000==0){
			printf("iiter=%d iaccp=%d dtime=%f ijump=%d iloop=%d\n",iiter,iaccp,(float)(time(0)-start),i,iloop);
		}
  		//---------
		gen_newpara(para2,para1,k2,isoflag,Viso); //in:para2 out:para1
		para2mod(para1,refmodel,ttmodel);

		//para2=para1;
  		//-----------------Criteria 1, goodmodel------------------------
         	 if(PosAni>0){//require the L always faster than R
	   		 //cout<<"check 2posani | ";//test---
            		updatemodelTibet(ttmodel,depcri1,depcri2,qpcri,qscri);
            		if(positiveAni(ttmodel,Vposani)==0){
		    		tflag=1;
		    		continue;
            		}

          	}//if PosAni>0
	
	 	if(Rmonoc+Lmonoc>0){
            		if(PosAni<=0){updatemodelTibet(ttmodel,depcri1,depcri2,qpcri,qscri);}
	    		if(goodmodel(ttmodel,Rvmono,Rvgrad,Rmonoc,0,isoflag)*goodmodel(ttmodel,Lvmono,Lvgrad,0,Lmonoc,0)==0){
	    		tflag=2;	
    	    		continue;}
        	}//if monoc 

		get_misfitKernel(ttmodel,para1,refmodel,refpara,kernel,Rsurflag,Lsurflag,depcri1,depcri2,qpcri,qscri,inp);
		newL=ttmodel.data.L;//para1.L;
		newRL=ttmodel.data.Rdisp.L;
		newLL=ttmodel.data.Ldisp.L;
		newmisfit=para1.misfit; 

		//----------Criteria 2, acceptable L-----------------------------------
		if(misfitcri<-1){
	/*	  if(newRL<oldRL)
		  {//continue;
	  	   prob=(oldRL-newRL)/oldRL;
	 	   prandom=gen_random_unif01();
	  	   if(prandom<prob)continue;
	 	  }
		*/
	      /*  if(newLL<oldLL){
	  	//printf("### newll<oldLL, %g<%g\n",newLL,oldLL);//---test---
	 	prob=(oldLL-newLL)/oldLL;
	  	prandom=gen_random_unif01();
	  	if(prandom<prob)continue;
	  	}
		*/	
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
			updatemodelTibet(ttmodel,depcri1,depcri2,qpcri,qscri);
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
  if(iaccp<accpcri)return 0;
  return 1;

}
