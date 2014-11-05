/*========CONTENT============
int initpara(paradef &para)
int checkParaModel(paradef para, modeldef model)
int readpara(vector<vector<double> > &para0, const char* fname)
int mod2para(modeldef &model, paradef &inpara, paradef &outpara)
int para2mod(paradef para, modeldef inmodel, modeldef &outmodel)
int gen_newpara_single (  vector<vector<double> > space1, vector<double> &parameter, int npara,int pflag)
int gen_newpara(paradef inpara, paradef &outpara, int pflag )
int para_avg(vector<paradef> &paralst, paradef &paraavg,vector<double> &Rparastd,vector<double> &Lparastd, vector<int> &idlst)
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
	  int i,j,N,N2,count,countsigma;
	  int ppflagid,nvid,ngid;
	  int p0,p5,p6,ng;
	  // (int)para.para0[j][4]:ng
	  // (int)para.para0[j][5]:nv
	  // (int)para.para0[j][6]:type of para
 	  // use after mod2para
	  if(para.flag<1){
	  	printf("Use checkParaModel after the para.space1 is filled, i.e., after mod2para\n");
	  	exit(0);
	  }

	  for(i=0;i<model.ngroup;i++){
	    //if(model.groups[i].flagcpttype!=2 and model.groups[i].flagcpttype!=4)continue;
	    N=model.groups[i].np*7; //number of para required
	    N2=model.groups[i].np*5; //number of Vkerne need  to compute
	    count=0;
	    countsigma=0;
	    ppflagid=-1;nvid=-1,ngid=-1;

	    if(model.groups[i].flag==4 and model.groups[i].flagcpttype!=1){
	      printf("### checkParaModel. group%d, if groups is described by gradient, then the forward cpt must be through Vkernel and Vpara!\n",i);
	      exit(0);
	    }


	    for(j=0;j<para.npara;j++){
	      if((int)para.para0[j][4]!=i)continue;

	      p0=(int)para.para0[j][0];    
     	      p5=(int)para.para0[j][5];//nv
	      p6=(int)para.para0[j][6];//flag of para     
	      if(p0==0){//para is value/Bcoeff	 

		// this criteria must be used when computing the partial derivatives, but during the computation, if partial derivatives are read from outside, this isn't required; probably need to be modified ...
	        if(p6==1 and model.groups[i].flagcpttype==3 and para.space1[j][2]*(para.space1[j][2]-0.001)<0. ){
			printf("### checkParaModel, in group %d, whose flagcpttype=3, the partial derivative of the vsv must be computed, and the sigma should be >0(freely perturb) or <-1(vsv stay constant)!\n",i);
			exit(0);
		}
		

		if(para.space1[j][2]<-1. and (p6-2)*(p6-6)*(p6-7)*(p6-10)*(p6-11)==0 ){//this theta or phi or AZcos or AZsin will be scaled according to other values(the 1st value in that group) during the para2mod process
			//also,  vsh can be sacled according to the 1st value in that group; to set the same amount of VsRA; if vsv has space<-1, it mean this vsv won't be changed
			if(p5<1){
				printf("### checkParaModel, in group %d, only the theta,phi,AZcos, AZsin for the layer(nv)>=1 (start from 0) can be scaled! here nv=%d\n",i,p5);
				exit(0);
			}
		
		}
		      
		//---- check the order of the layer     
		if(p5>nvid){ppflagid=-1;}//reach the next layer in this group
	        else if (p5<nvid){
			printf("### checkParaModel, in group %d, we require the layer number nv (teh 6th colunm) is increading!\n",i);
			exit(0);
		}//if	
		nvid=p5;

		//---- check the order within each layer
		if(p6<1 or p6>11){
			printf("inproper ppflag(%d,should be within 1-11) for the %dth parameter in para.in input\n",p6,i);
			exit(0);
		}
		if((p6==10 or p6==11) and model.groups[i].flagcpttype!=3){
			printf("inproper ppflag, ppflag=%d, but m.g.flagcpttype!=3 but =%d\n",p6, model.groups[i].flagcpttype);
			exit(0);
		}
		if (p6<=ppflagid){
	      	  printf("### checkParaModel, in group %d, layer %d,  we require the para's ppflag (the 7th column) is ordered in an increasing order!\n",i,nvid);
		  exit(0);
	        }//if
	        ppflagid=p6;
	      }//if p0==0
	  
	      //---- count the total number of parameters
	      if(p6>7)continue;
	      //parameter belongs to group i, and is one of the five parameters (vsv~eta)
	      count++;
	      if(p6<6 and para.space1[j][2]*(para.space1[j][2]+1)>0)//if sigma >0 or sigma<-1, its Vkernel will be computed
	      {countsigma++;}
	    }//for j <npara

	    if((model.groups[i].flagcpttype-2)*(model.groups[i].flagcpttype-4)==0){
		if(count!=N){
      		    printf("#### checkParaModel, the number of parameters(vsv~eta,theta,phi) in group %d is %d, but should be %d to enable the Love_para forward computation!\n",i,count,N);
	      	    exit(0);}
	        if(countsigma!=N2){
		    printf("#### checkParaModel, the number of parameters with sigma>0 or <-1 (i.e. will compute Vkernel for this para) in group %d is %d, but should be %d to enable right Love_para forward computation!\n",i,countsigma,N2);
		    exit(0);
		}
	    }//if count

	  }//for i
	  return 1;
	}  //checkParaModel;


//-----------------------------------------------------	 
//	0		1	2	3	4   5	6	7		8			9		10	
//	flag		dv_type	dv	sigma	ng  nv	ppflag	LVflag		RW_flag			LW_flag		AZ_flag
//	0--value	1--dv				1-vsv	1--use Vkernel	1--for RAdisp cpt	0--not for	0/2/4 	
//	1--gp_thick	else dv*100			2-vsh 	to do forward	2--for Azidisp cpt	   L wave	indicating cos(X*phi)
//	-1--vpvs					3-vpv	computation	0--not for Rayleigh	1/-1--is for
//							4-vph			   wave
//							5-eta	-1--use Lkernel
//							6-theta
//							7-phi
//							8-rho
//							9-h
//							10-dVsv_cos
//							11-dVsv_sin
//							12-vpvs
//
	//int readpara(vector<vector<string> > &para0, int &npara,const char* fname)
	int readpara(paradef &para, const char* fname)
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
	  para.para0.clear();
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
		if((int)vd[0]==1)vd[6]=9.;//thickness; vd[6] is the ppflag, indicating the type of parameter
		else if((int)vd[0]==-1)vd[6]=12.;//vpvs
		//printf("@@@ check, readpara, width of of para0 is %d, should be 11!\n",vd.size());
		k=k+1;
		para.para0.push_back(vd);	
	  }
	  mff.close();
	  para.npara=k;
	  //}//if surflag		
	  return 1;
	}//readpara
//----------------------------------------------------

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
	  	p2=(float)inpara.para0[i][2]; // perturbatin; dv
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
		  else if (p6==6 or p6==10)
		  	tv=model.groups[ng].thetavalue[nv];
		  else if (p6==7 or p6==11)
		  	tv=model.groups[ng].phivalue[nv];
		  else if (p6==9)
			tv=model.groups[ng].rhovalue[nv];
		  else {	
			printf("### mod2para, wrong number for the 7th column in para.in! ppflag=%d\n",p6);
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

		if(inpara.flag<1){
			//---fill the para.space1
			sigma=outpara.para0[i][3]; // if don't want to perturb this parameter, then set sigma=0.
			if(p1==1){
				tmin=tv-p2;tmax=tv+p2;
				if(fabs(tmin-tmax)<1E-5 and sigma>1E-5){
					printf("### inproper para.in &mod set up for the %dth parameter(starts from 0), min,max=[%g,%g] RESET IT\n",i,tmin,tmax);
					exit(0);
				}
			}
			else{
				tmin=tv*(1-p2/100.);tmax=tv*(1+p2/100.);
				if(fabs(tmin-tmax)<1E-5 and sigma>1E-5){
					printf("### inproper para.in &mod set up for the %dth parameter(starts from 0), min,max=[%g,%g] RESET IT\n",i,tmin,tmax);
					exit(0);
				}
			}
			if(p6<8){ // velocity or eta or theta or phi, should >0
				tmin=max(0.,tmin);
				tmax=max(0.,tmax);
				tmax=max(tmin+0.001,tmax);
				// there was a constraint on the sedimental velocity, don't know why need it, so it's not added here;
			}
			vt[0]=tmin;vt[1]=tmax;vt[2]=sigma;
			//---test---
			//if(p6==7){vt[0]=90.;vt[1]=180.;}
			//if(p6==6){vt[0]=60.;vt[1]=90.;}
			outpara.space1.push_back(vt);
			printf("@@@ npara%d, v=%g [%g,%g] p6=%d\n",i,tv,vt[0],vt[1],p6);
		}//inpara.flag<1
	  }//for i<npara	  


	  if(inpara.flag<1){ //initial flag
		for(i=0;i<outpara.npara;i++){
			ng=(int)inpara.para0[i][4];
			p6=(int)inpara.para0[i][6];
			//---fill the flags in para0
			if(model.groups[ng].flagcpttype==2 or model.groups[ng].flagcpttype==4){
		  	  outpara.para0[i][7]=-1;//the LVflag, 1--use Vpara&Vkernel to do forward computation; -1--use Lovepara&Lovekernel to do forward computation
			}
			else{
			  outpara.para0[i][7]=1;
			}
			//printf("@@@ check, mod2para, the %dth LVflag=%g ng=%d,flagcpttype=%d \n",i,outpara.para0[i][7],ng,model.groups[ng].flagcpttype);
			if(p6==1){
			  outpara.para0[i][8]=1;
			  outpara.para0[i][9]=-1;
			  outpara.para0[i][10]=2;
			}
			else if(p6==2){
			  outpara.para0[i][8]=1;
			  outpara.para0[i][9]=-1;
			  outpara.para0[i][10]=4;
			}
			else if(p6==3){
			  outpara.para0[i][8]=1;
			}
			else if(p6==4  or p6==5){
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
	int para2mod_static(paradef para, modeldef inmodel, modeldef &outmodel)
	{
	  //this is different from the para2mod defined below this function; it won't modify the para anymore, it just use the given para, and transfer it to model, regardless if space1[i][2]<-1 or not;
	  // this is used in computing Vkernel!	
	  int i, p0,p1,ng,p6,nv,pflag;
 	  double newv;
	  int np,flagvpv,flagvph,flagvsv,flagvsh;
	  flagvpv=flagvph=flagvsv=flagvsh=-1;
	  float factor;

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
		  else if (p6==3){//vpv
			  outmodel.groups[ng].vpvvalue[nv]=newv;flagvpv=1;}
		  else if (p6==4){//vph
			  outmodel.groups[ng].vphvalue[nv]=newv;flagvph=1;}
		  else if (p6==5){//eta
			   outmodel.groups[ng].etavalue[nv]=newv;}
		  else if (p6==6 or p6==10){//theta
			  outmodel.groups[ng].thetavalue[nv]=newv;
		  	  }//p6=6 or 10	
		  else if (p6==7 or p6==11){//phi
			  outmodel.groups[ng].phivalue[nv]=newv;
		  	  }//if p6==7 or 11
		}//if p0==0

		else if(p0==1){// groups thickness
			outmodel.tthick=outmodel.tthick-outmodel.groups[ng].thick+newv; // the tthick is changing
			//outmodel.groups[outmodel.ngroup-1].thick=outmodel.groups[outmodel.ngroup-1].thick+(outmodel.groups[ng].thick-newv);//change the thickness of the last group while keeping the tthcik the same // find a bug, if ngroup<3, this choice would cause problem!
			outmodel.groups[ng].thick=newv;
		
		}// else if p0==1
		
		else if (p0==-1){//vpvs
			outmodel.groups[ng].vpvs=newv;
		}
		else{cout<<"########wrong!!!!! in para2mode"<<endl;exit(0);}
	  }// for i


	  return 1;
	}// para2mod_static
	 


//-----------------------------------------------------	 

	int para2mod(paradef &para, modeldef inmodel, modeldef &outmodel)
	{
	  //for a parameteri, if its space[i][2] (i.e., the sigma) <-1,then, use the default RAvp, RAvs, eta relationship to scale vpv,vph,eta according to vsv,vsh;
	  //ALSO, since I add & in front of para, the para will be changed accordingly (change to the scaled value, so don't need mod2para after para2mod in order to make sure para and mod have the same value)
	  //RAvs=(vsh-vsv)/((vsv+vsh)/2)
	  //vpv=vsv*vpvs
	  //A=0.25*RAvs
	  //vph=(1+A)/(1-A)*vpv
	  //eta=1-4.2*RAvs 	
	  //
	  //another way to scale eta, make the media to have purely ellipsoidal anisotropy
	  //eta=0.5*(1+(vpv^2-2*vsv^2)/(vph^2-2*vsv^2))
	  //
	  //for vsh, if the sigma<-1; then scale it according to the vsv in the current layer & the RAvs value in the 1st layer of that group
	  //set them to have the same VsRA (now I use the simplified computation in which RAvs=(vsh-vsv)/(vsv+vsh)/2)
	  //
	  //
	  int i, p0,p1,ng,p6,nv,pflag;
 	  double newv;
	  int np,flagvpv,flagvph,flagvsv,flagvsh;
	  flagvpv=flagvph=flagvsv=flagvsh=-1;
	  float factor;

	  outmodel=inmodel;
	  for(i=0;i<para.npara;i++){
	  	p0=(int)para.para0[i][0];//flag, explanations are in mod2para;
	  	p1=(int)para.para0[i][1];
	  	ng=(int)para.para0[i][4];
	  	p6=(int)para.para0[i][6];//ppflag, explanations are in mod2para;
		
		newv=para.parameter[i];

		if(p0==0){ //value
		  nv=(int)para.para0[i][5];
		  if(p6==1){
				outmodel.groups[ng].vsvvalue[nv]=newv;flagvsv=1;}
		  else if (p6==2){//vsh
			  if(para.space1[i][2]<-1){//scale it according to the VsRA in the layer 1 of this gorup
			  	factor=(outmodel.groups[ng].vshvalue[0]-outmodel.groups[ng].vsvvalue[0])/(outmodel.groups[ng].vshvalue[0]+outmodel.groups[ng].vsvvalue[0]); //0.5*RAvs
				outmodel.groups[ng].vshvalue[nv]=(1+factor)/(1-factor)*outmodel.groups[ng].vsvvalue[nv];
				para.parameter[i]=outmodel.groups[ng].vshvalue[nv];	
			  }
	       		  else		  
			  	outmodel.groups[ng].vshvalue[nv]=newv;
			  flagvsh=1;}
		  else if (p6==3){//vpv
			  if(para.space1[i][2]<-1.)//use the RAvp, eta, RAvs relation
				{outmodel.groups[ng].vpvvalue[nv]=outmodel.groups[ng].vpvs*outmodel.groups[ng].vsvvalue[nv];
				 para.parameter[i]=outmodel.groups[ng].vpvvalue[nv];
				}
			  else	  
			   	{outmodel.groups[ng].vpvvalue[nv]=newv;}
			  flagvpv=1;
			  }//p6==3

		  else if (p6==4){//vph
			  if(para.space1[i][2]<-1.)//use the RAvp, eta, RAvs relation
			  {
			  	if(ng<2){//only do scaling in the crust (sed+crystalline crust) part
				factor=(outmodel.groups[ng].vshvalue[nv]-outmodel.groups[ng].vsvvalue[nv])/(outmodel.groups[ng].vshvalue[nv]+outmodel.groups[ng].vsvvalue[nv])*0.5; //0.25*RAvs (VpRA=0.5VsRA)
				//factor=(outmodel.groups[ng].vshvalue[nv]-outmodel.groups[ng].vsvvalue[nv])/(outmodel.groups[ng].vshvalue[nv]+outmodel.groups[ng].vsvvalue[nv]); //0.5*RAvs (VpRA=VsRA)
				factor=(1+factor)/(1-factor);	  	
				outmodel.groups[ng].vphvalue[nv]=factor*outmodel.groups[ng].vpvvalue[nv];}
				else{outmodel.groups[ng].vphvalue[nv]=outmodel.groups[ng].vpvvalue[nv];}
				para.parameter[i]=outmodel.groups[ng].vphvalue[nv];
			   }
			   else	
			  	{outmodel.groups[ng].vphvalue[nv]=newv;} 
			   flagvph=1;	
		  	   }//p6==4

		  else if (p6==5){//eta
			   if(para.space1[i][2]<-1.)//use the RAvp, eta, RAvs relation
			   {
				  if(ng<2){
			   	  //factor=(outmodel.groups[ng].vshvalue[nv]-outmodel.groups[ng].vsvvalue[nv])/(outmodel.groups[ng].vshvalue[nv]+outmodel.groups[ng].vsvvalue[nv])*2;//RAvs
				  //outmodel.groups[ng].etavalue[nv]=1.0-4.2*factor;
				  //c==0, elliptical case
 				  factor=(pow(outmodel.groups[ng].vpvvalue[nv],2)-2*pow(outmodel.groups[ng].vsvvalue[nv],2))/(pow(outmodel.groups[ng].vphvalue[nv],2)-2*pow(outmodel.groups[ng].vsvvalue[nv],2));
				  outmodel.groups[ng].etavalue[nv]=0.5*(1+factor);
				  }
				  else{outmodel.groups[ng].etavalue[nv]=1.0;}
   				  para.parameter[i]=outmodel.groups[ng].etavalue[nv];	
			   }
		  	   else
			  	{outmodel.groups[ng].etavalue[nv]=newv;}
			   }//p6==5
		  else if (p6==6 or p6==10){//theta
			  if(para.space1[i][2]<-1)//set the theta to the value of the 1st layer in that group
			  {
			  	outmodel.groups[ng].thetavalue[nv]=outmodel.groups[ng].thetavalue[0];
				para.parameter[i]=outmodel.groups[ng].thetavalue[nv];
			  }
			  else
			  	{outmodel.groups[ng].thetavalue[nv]=newv;}
		  	  }//p6=6 or 10	
		  else if (p6==7 or p6==11){//phi
			  if(para.space1[i][2]<-1)
			  {
			  	outmodel.groups[ng].phivalue[nv]=outmodel.groups[ng].phivalue[0];
				para.parameter[i]=outmodel.groups[ng].phivalue[nv];
			  }	  
			  else
			  	{outmodel.groups[ng].phivalue[nv]=newv;}
		  	  }//if p6==7 or 11
		}//if p0==0

		else if(p0==1){// groups thickness
			outmodel.tthick=outmodel.tthick-outmodel.groups[ng].thick+newv; // the tthick is changing
			//outmodel.groups[outmodel.ngroup-1].thick=outmodel.groups[outmodel.ngroup-1].thick+(outmodel.groups[ng].thick-newv);//change the thickness of the last group while keeping the tthcik the same // find a bug, if ngroup<3, this choice would cause problem!
			outmodel.groups[ng].thick=newv;
		
		}// else if p0==1
		
		else if (p0==-1){//vpvs
			outmodel.groups[ng].vpvs=newv;
		}
		else{cout<<"########wrong!!!!! in para2mode"<<endl;exit(0);}
	  }// for i


	  //--the following value filling process have been performed in the readmodAniso function.
	  //---depending on the pflag of each group, fill all the unfilled vectors among those vsv...phi 7 vectors with appropriate values derived from certain relationships between vp&vs;
	/*
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
	*/
	  return 1;
	}// para2mod
	 


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
		  //printf("test-- gen_newpara i=%d (out of %d-1) mean=%g sigma=%g space=[%g~%g]\n",i,npara,mean,sigma,space1[i][0],space1[i][1]);
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

	int gen_newpara(paradef inpara, paradef &outpara, int pflag )
	{//change the p.parameter[]

	 outpara=inpara;
	 gen_newpara_single(inpara.space1,outpara.parameter,inpara.npara,pflag);

	  return 1;
	}

//----------------------------------------------------- 
	int para_avg(vector<paradef> &paralst, paradef &parabest, paradef &paraavg,vector<double> &parastd,vector<double> &LoveRAparastd, vector<vector<double> > &LoveAZparastd,vector<int> &idlst){
	//get the average parameters from a list of acceptable paras
	//also, get the best parameter (min misfit)
	    int i,j,k,size,Ngood,ibadp,idmin;
	    double mismin=1e10;
	    vector<double> LAZp0(2,0.);
	    vector<int> badv;

	    initpara(paraavg);
  	    initpara(parabest);
	    size=paralst.size();
	    if(size<1) return 0;
		
	    paraavg=paralst[0];
	    parabest=paralst[0];
	    parastd.clear();LoveRAparastd.clear();LoveAZparastd.clear();

//what's the size of love para
	    for(j=0;j<paralst[0].npara;j++){// the size of Love para is the same as V para
                parastd.push_back(0.);
                LoveRAparastd.push_back(0.);
		LoveAZparastd.push_back(LAZp0);
		paraavg.parameter[j]=paraavg.LoveRAparameter[j]=paraavg.LoveAZparameter[j][0]=paraavg.LoveAZparameter[j][1]=0.;
            }//for j

	    for(i=0;i<size;i++)
		{if(paralst[i].misfit<mismin) {idmin=i;mismin=paralst[i].misfit;}}
	    //min(2*Mis_min,Mis_min+0.5)
	    parabest=paralst[idmin];
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
		  for(j=0;j<paralst[0].npara;j++){
			paraavg.parameter[j]+=paralst[i].parameter[j];
			paraavg.LoveRAparameter[j]+=paralst[i].LoveRAparameter[j];
			paraavg.LoveAZparameter[j][0]+=paralst[i].LoveAZparameter[j][0];
			paraavg.LoveAZparameter[j][1]+=paralst[i].LoveAZparameter[j][1];}
		}//if
	    }//for i
	   for(i=0;i<paralst[0].npara;i++){
		paraavg.parameter[i]/=Ngood;
		paraavg.LoveRAparameter[i]/=Ngood;
		paraavg.LoveAZparameter[i][0]/=Ngood;
		paraavg.LoveAZparameter[i][1]/=Ngood;
	   }//for i
	   //cout<<"Ngood="<<Ngood<<endl; 
	   printf("@@@ para_avg size =%d, Ngood=%d mismin=%g\n",size,Ngood,mismin);
	   //if(Ngood>size){printf("WRONG HERE! Ngood > size!!\n");exit(0);}
	   for(i=0;i<paralst[0].npara;i++){
		for(j=0;j<Ngood;j++){
		    k=idlst[j];
		    parastd[i]+=pow(paralst[k].parameter[i]-paraavg.parameter[i],2);
		    LoveRAparastd[i]+=pow(paralst[k].LoveRAparameter[i]-paraavg.LoveRAparameter[i],2);
		    LoveAZparastd[i][0]+=pow(paralst[k].LoveAZparameter[i][0]-paraavg.LoveAZparameter[i][0],2);
		    LoveAZparastd[i][1]+=pow(paralst[k].LoveAZparameter[i][1]-paraavg.LoveAZparameter[i][1],2);
		}
		parastd[i]=sqrt(parastd[i]/Ngood);
		LoveRAparastd[i]=sqrt(LoveRAparastd[i]/Ngood);
		LoveAZparastd[i][0]=sqrt(LoveAZparastd[i][0]/Ngood);
		LoveAZparastd[i][1]=sqrt(LoveAZparastd[i][1]/Ngood);
		if (parastd[i]>50){//indicating it's phi and has period problem that hasn't been taken into account
			badv.push_back(i);
		}
	    }//for i
	    //========added Nov23, 2013; check the para_std, if std too large (>50), indicating it's phi and has period problem that hasn't been taken into account
	    for(i=0;i<badv.size();i++){
		ibadp=badv[i];
		paraavg.parameter[ibadp]=0.;//paraavg.LoveRAparameter[ibadp]=paraavg.LoveAZparameter[ibadp][0]=paraavg.LoveAZparameter[ibadp][1]=0.;
		for(j=0;j<Ngood;j++){
			k=idlst[j];
			if(paralst[k].parameter[ibadp]>90){//the phi should be the same among parameter, LoveRApara.. LoveAZpara..
				paralst[k].parameter[ibadp]-=180;
				//paralst[k].LoveRAparameter[ibadp]-=180;
				//paralst[k].LoveAZparameter[ibadp][0]-=180;
				//paralst[k].LoveAZparameter[ibadp][1]-=180;
			}
			paraavg.parameter[ibadp]+=paralst[k].parameter[ibadp];
		}//for j
		paraavg.parameter[ibadp]/=Ngood;
		paraavg.LoveRAparameter[ibadp]=paraavg.parameter[ibadp];
		paraavg.LoveAZparameter[ibadp][0]=paraavg.parameter[ibadp];
		paraavg.LoveAZparameter[ibadp][1]=paraavg.parameter[ibadp];
	    }//for i<badv.size
	    for(i=0;i<badv.size();i++){
		ibadp=badv[i];
		parastd[ibadp]=0.;
		for(j=0;j<Ngood;j++){
			k=idlst[j];
			parastd[ibadp]+=pow(paralst[k].parameter[ibadp]-paraavg.parameter[ibadp],2);
		}//for j
		parastd[ibadp]=sqrt(parastd[ibadp]/Ngood);
		LoveRAparastd[ibadp]=parastd[ibadp];
		LoveAZparastd[ibadp][0]=parastd[ibadp];
		LoveAZparastd[ibadp][1]=parastd[ibadp];
	    }//for i
	    return 1;
	}//para_avg
//-----------------------------------------------------
