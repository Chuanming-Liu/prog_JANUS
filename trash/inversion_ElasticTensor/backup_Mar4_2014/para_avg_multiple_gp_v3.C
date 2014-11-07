// this is a function that's similar to the para_avg, but this function could handle the case of multiple peak.
// decide if there is multiple peaks(by defulat <=2 peaks) in the phi(strike) parameter. if there is, then group the parameters into 2 groups. one group with phi~phi1, ther other group with phi~phi1+90, and Npk=2; if there isn't, Npk=1
// do normal para_avg for each para group (#of para gp is Npk), and return avg, std

// this version could handle multiple(<=2) phi group in both crust and mantle. The final output has Nphi_group_crust*Nphi_group_mantle  paraavgs
// the seperation of mantle group can be disabled by setting idphiM<0

// this version, get parabest for each phi group
//
float convert(float vin,float vref,float T){
  float v;
  v=vin;
  //while(vref>T)vref-=T;
  //while(vref<0)vref+=T;

  while(v-vref>0.5*T)v-=T;
  while(v-vref<-0.5*T)v+=T;
  return v;
}//convert
//--------------------------

vector<int> seperate_gp(vector<double> &vlst,int &Ngp){
  //this is used to group the parameters according to the phi(strke) value in the crust (HOW ABOUT MANTLE?)
  // input the para list (a list of phi) that will be judged during the grouping process. 
  //return two index lists, each list belongs to one phi group
  //
  // 
  //#############PARAMETER
  float distcri=5;//if |peak-avg|>distcri, then we think there are multiple groups of phi value
  float dv,dist,tv,dist1,dist2;
  float T=180.; //period of phi
  double avg1,avg2,pk;
  int i,j,k,Nv,nmax,Nbin;  
  vector<int> binlst,indexflaglst;  

  Nv=vlst.size();
  //--compute two kind of average values; avg1=sum(v)/N; avg2=sum(v if v<T/2; v-T if v>T/2)/N
  avg1=0.;avg2=0.;
  for(i=0;i<Nv;i++){
	if(vlst[i]<0 or vlst[i]>T){
		printf("#### seperate_gp the value (%g) in input list is outside the given range[0,%g]",vlst[i],T);
		exit(0);
	}
	avg1+=vlst[i];
	if(vlst[i]<T/2)avg2+=vlst[i];
	else{avg2+=(vlst[i]-T);}
  }//for i<Nv
  avg1/=Nv;
  avg2/=Nv;  

  //-- compute the peak value; bin the list and get the value correlated with most data
  dv=5.;//bin width
  Nbin=int(T/dv);
  if(int(T*10)%int(dv*10)!=0){
	printf("### seperate_gp, the dv(%g) is not good! mod(T,dv)=mod(%g,%g)!=0\n",dv,T,dv);
	exit(0);
  }
  
  for(i=0;i<Nbin;i++){
	binlst.push_back(0);
  }
  for(i=0;i<Nv;i++){
	k=int(floor(vlst[i]/dv));
	if(k==Nbin){printf("#### seperate_gp, strange case, vlst[%d]=%g,k=%d==%d\n",i,vlst[i],k,Nbin);exit(0);}
	binlst[k]++;
  }  
  //get the value with max binlst value
  k=0;
  nmax=-100;
  for(i=0;i<Nbin;i++){
	if(binlst[i]>nmax){nmax=binlst[i];k=i;}
  }
  pk=(k+0.5)*dv;//the mid value of that bin; normally, I think pk is correlated with the group with c<=0, but I might be wrong
  printf("hey, peak value = %g avg=%g & %g===\n",pk,avg1,avg2);//---check---

  //-- choose between avg1, avg2 based on their dist to pk, and get the |avg-pk|
  avg1=convert(avg1,pk,T);
  avg2=convert(avg2,pk,T);
  dist1=fabs(avg1-pk);
  dist2=fabs(avg2-pk);
  dist=dist1<dist2?dist1:dist2;
  
  //--decide if there is multiple peaks in the vlst based on the value of dist
  if(dist>distcri)Ngp=2;
  else Ngp=1;

  printf("dist=%g distcri=%g, Ngp=%d\n",dist,distcri,Ngp);//--check--
  //--group the index of vlst into Ngp group(s)
  indexflaglst.clear();
  if(Ngp==1){//only one group
	for(i=0;i<Nv;i++){
		indexflaglst.push_back(1);//this lst tells if this para belongs to gp1 or 2
		vlst[i]=convert(vlst[i],pk,T);		
	}//for i<Nv
  }//if Ngp==1
  else{//two groups, one with center value around pk, the other around pk+T*0.5 (pk+90)
 	for(i=0;i<Nv;i++){
		tv=convert(vlst[i],pk,T);
		if(fabs(tv-pk)<0.23*T){//### 0.23 is an arbitrary number
			//this para belongs to the groups with center value ~pk
			vlst[i]=tv;
			indexflaglst.push_back(1);	
			continue;
		}
		tv=convert(vlst[i],pk+0.5*T,T);
		if(fabs(tv-pk-0.5*T)<0.23*T){//belongs to the gp2
			//printf("%g->%g\n",vlst[i],tv);//---check--
			vlst[i]=tv;
			indexflaglst.push_back(2);
			continue;
		}
		indexflaglst.push_back(0);
		printf("---value %g, no gp (peak~%g and %g) wants it\n",vlst[i],pk,pk+0.5*T);//---check---
	}//for i<Nv
  }//else two groups

  if(indexflaglst.size()!=Nv){printf("something wrong, size doesn't match! %d!=%d\n",indexflaglst.size(),Nv);exit(0);}//---check---
  

  return indexflaglst;
}//seperate_gp
//--------------------------
vector<double> compute_paraavg(vector<paradef> paralst,vector<int> idlst, vector<double> &stdlst, int id ){
        int i,j,k;
        int Npara,Ngood;
        vector<double> vavglst;

        Npara=paralst[0].npara;
        Ngood=idlst.size();

	//printf("begin compute_paraavg %d %d\n",paralst[0].parameter.size(),paralst[0].LoveRAparameter.size());//--check--
        if((id-1)*(id-2)!=0){printf("### compute_paraavg, wrong value for id! should be 1 or 2\n");exit(0);}

        stdlst.clear();
        for(i=0;i<Npara;i++){
                vavglst.push_back(0.);
                stdlst.push_back(0.);
        }

	
        for(j=0;j<Npara;j++){
                for(i=0;i<Ngood;i++){
                        k=idlst[i];
                        if(id==1)vavglst[j]+=paralst[k].parameter[j];
                        else if(id==2)vavglst[j]+=paralst[k].LoveRAparameter[j];
                }
                vavglst[j]/=Ngood;
        }//for j<Npara

        for(j=0;j<Npara;j++){
                for(i=0;i<Ngood;i++){
                        k=idlst[i];
                        if(id==1)stdlst[j]+=pow(paralst[k].parameter[j]-vavglst[j],2);     
                        else {stdlst[j]+=pow(paralst[k].LoveRAparameter[j]-vavglst[j],2);}
                }
                stdlst[j]=sqrt(stdlst[j]/Ngood);
		//---check--
		if(stdlst[j]>50){
                        printf("###compute_paraavg, something wrong(?),para%d, avg=%g,large std %g\n",j,vavglst[j],stdlst[j]);
			//exit(0);
		}
		//
	}//for j<Npara                                                                              
  //printf("finish compute_paraavg\n");//--check--
  return vavglst;
}//compute_paraavg
//--------------------------
vector<vector<double> > compute_paraavgAZ(vector<paradef> paralst,vector<int> idlst, vector<vector<double> > &stdlst){
        int i,j,k;
        int Npara,Ngood;
        vector<double> LAZp0(2,0.);
	vector<vector<double> > vavglst;

        Npara=paralst[0].npara;
        Ngood=idlst.size();

	//printf("begin compute_paraavg %d %d\n",paralst[0].parameter.size(),paralst[0].LoveRAparameter.size());//--check--

        stdlst.clear();
        for(i=0;i<Npara;i++){
                stdlst.push_back(LAZp0);
		vavglst.push_back(LAZp0);
        }

	
        for(j=0;j<Npara;j++){
                for(i=0;i<Ngood;i++){
                        k=idlst[i];
                        vavglst[j][0]+=paralst[k].LoveAZparameter[j][0];
			vavglst[j][1]+=paralst[k].LoveAZparameter[j][1];                   
                }
                vavglst[j][0]/=Ngood;
		vavglst[j][0]/=Ngood;
        }//for j<Npara

        for(j=0;j<Npara;j++){
                for(i=0;i<Ngood;i++){
                        k=idlst[i];
			stdlst[j][0]+=pow(paralst[k].LoveAZparameter[j][0]-vavglst[j][0],2);     
			stdlst[j][1]+=pow(paralst[k].LoveAZparameter[j][1]-vavglst[j][1],2);
			
                }
                stdlst[j][0]=sqrt(stdlst[j][0]/Ngood);
		stdlst[j][1]=sqrt(stdlst[j][1]/Ngood);
		
	}//for j<Npara       
                                                                   
  //printf("finish compute_paraavg AZ\n");//--check--
  return vavglst;
}//compute_paraavgAZ
//--------------------------
vector<int> para_avg_multiple_gp(int idphi,int idphiM, vector<paradef> &paralst, vector<paradef> &parabestlst, vector<paradef> &paraavglst, vector<paradef> &parastdlst, vector<vector<int> > &idlstlst, int flag){
  // in this function, if idphiM<0 then, won't do mantle group seperation based on mantle phi
  // flag indicate if average is for para.parameter(flag=1) or for para.parameter/LoveRAparameter/LoveAZparameter (flag=3)
  int i,j,k,size,idmin,Ngp,NgpM,igp,pflag,Ngood;
  vector<int> indexflaglst,indexflaglstM,idlst,idlstnew,idminlst;
  vector<double> philst,philstM;  
  double mismin;  
  //vector<double> parastd;
  paradef paraavg,parastd;

  if((flag-1)*(flag-3)!=0){
	printf("### para_avg_multiple_gp, flag should be 1(compute avg for para.parameter), or 3 (compute avg for para.parameter, para.LoveRAparameter, and para.LoveAZparameter), but input is %d\n",flag);
	exit(0);
  }
  //printf("hey begin!\n");//--check--
  idlstlst.clear();
  size=paralst.size();
  if(size<1) return idminlst;
 
  paraavg=paralst[0];
  //parabest=paralst[0];
  parastd=paralst[0];
  //parastd.clear();
 
  //--initialize the paralst and parastd --> will be handled by the function compute_paraavg

  //--get the ipara for the phi
  /*for(i=0;i<paralst[0].npara;i++){
	printf("para0.size=%d\n",paralst[0].para0.size());//--check--
	igp=(int)paralst[0].para0[i][4];
	pflag=(int)paralst[0].para0[i][6];
	if(igp==1 and pflag==7){//since the phi is constant in crust, just take the 1st phi
		idphi=i;
		break;
	}
  }//for i
  */
  printf("idphi=%d\n",idphi);//--check--

  //--get the smallest misfit and get the index lst for model within misfit criteria-- this is just a groups of model with acceptable misfit; then after group seperation, we will re-select the model with higher criteria;(two criteria, one in seperate_gp[based on its dist to the avg phi value], one in the later part of this subroutine[based on the mismin of that phi gorup])
  mismin=1e10;
  for(i=0;i<size;i++){
	if(paralst[i].misfit<mismin){idmin=i;mismin=paralst[i].misfit;}
  }
  //parabest=paralst[idmin];
  printf("@@@ para_avg mismin=%g\n",mismin);
  //mismin=mismin+0.5;
  mismin=(mismin*2>(mismin+1.0))?mismin*2:(mismin+1.0); //arbitrary selection criteria
  Ngood=0;
  idlst.clear();
  for(i=0;i<size;i++){
	if(paralst[i].misfit<mismin){
		Ngood++;
		idlst.push_back(i);
	}
  }

  //---get philst(only from mod with small misfit), and call function seperate_gp to group the philst
  philst.clear();
  for(i=0;i<Ngood;i++){
	k=idlst[i];
	philst.push_back(paralst[k].parameter[idphi]);
  }  

  Ngp=0;
  indexflaglst=seperate_gp(philst,Ngp);
  
  //--- put the modified (+/-T) phi back into the paralst; and compute the avg for all parameters for each gp(seperated based on philst)
  for(i=0;i<Ngood;i++){
	k=idlst[i];
	/*--check---
	if(fabs(paralst[k].parameter[idphi]-philst[i])>1){
		printf("  %g->%g\n",paralst[k].parameter[idphi],philst[i]);
	}*/
	paralst[k].parameter[idphi]=philst[i];
  }
  //


  //## if want to seperate the model also based on mantle phi,then
  if(idphiM>=0){ // do do group seperation for mantle
  philstM.clear();
  for(i=0;i<Ngood;i++){
	k=idlst[i];
	philstM.push_back(paralst[k].parameter[idphiM]);
  }
  NgpM=0;
  indexflaglstM=seperate_gp(philstM,NgpM);
  for(i=0;i<Ngood;i++){
	k=idlst[i];
	paralst[k].parameter[idphiM]=philstM[i];
  }
  }//if idphiM>0 
  else{// do not do group seperation for mantle
	indexflaglstM.clear();
	NgpM=1;
	for(i=0;i<Ngood;i++)indexflaglstM.push_back(1);
  }

  for(i=0;i<Ngp;i++){
	for(j=0;j<NgpM;j++){
  		idlstnew.clear();
		//---find the parabest for this phi group
		mismin=1e10;
		idmin=-1;
		for(k=0;k<Ngood;k++){
			if(paralst[idlst[k]].misfit<mismin and indexflaglst[k]==i+1 and paralst[idlst[k]].misfit<mismin){
				mismin=paralst[idlst[k]].misfit;
				idmin=idlst[k];
			}
		}
		idminlst.push_back(idmin);
		parabestlst.push_back(paralst[idmin]);
		mismin=mismin+0.5;

  		//--prepare the new idlst
  		for(k=0;k<Ngood;k++){
			if(paralst[idlst[k]].misfit<mismin and indexflaglst[k]==i+1 and indexflaglstM[k]==j+1)idlstnew.push_back(idlst[k]);
		}//for k<Ngood
		printf("====the %dth group has %d models,mismin(after+0.5)=%g=====\n",i*Ngp+j,idlstnew.size(),mismin);//--check--
		if(flag==1)paraavg.parameter=compute_paraavg(paralst,idlstnew,parastd.parameter,1);
		else if(flag==3){
			paraavg.parameter=compute_paraavg(paralst,idlstnew,parastd.parameter,1);
			paraavg.LoveRAparameter=compute_paraavg(paralst,idlstnew,parastd.LoveRAparameter,2);
			paraavg.LoveAZparameter=compute_paraavgAZ(paralst,idlstnew,parastd.LoveAZparameter);
		}
		if(idlstnew.size()<2)continue;	
		idlstlst.push_back(idlstnew);
    		paraavglst.push_back(paraavg);
		parastdlst.push_back(parastd);
  	}//for j<NgpM
  }//for i<Ngp

  return idminlst;
  //return idmin;
}//para_avg_multiple_gp
