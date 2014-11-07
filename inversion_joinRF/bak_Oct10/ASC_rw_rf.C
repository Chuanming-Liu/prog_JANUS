using namespace std;
int write_ASCdisp_single(dispdef disp, FILE *outdisp){
  int i;
  if(disp.npper>0)
  {
   if(disp.ngper>0)//write both
      {for(i=0;i<disp.npper;i++) // ATTENTION!! Here, by default, pper and gper are equal!!!! Otherwise, the writing part need to be changed!!!!
        fprintf(outdisp,"%5g ph %10g %10g %10g gp %10g %10g %10g\n",disp.pper[i],disp.pvelo[i],disp.pvel[i],disp.unpvelo[i],disp.gvelo[i],disp.gvel[i],disp.ungvelo[i]);}//if ngper>0
   else
      {for(i=0;i<disp.npper;i++) fprintf(outdisp,"%5g ph %10g %10g %10g gp 0 0 0\n",disp.pper[i],disp.pvelo[i],disp.pvel[i],disp.unpvelo[i]); }
 }//if npper>0
 else if (disp.ngper>0)
      {for(i=0;i<disp.ngper;i++)fprintf(outdisp,"%5g ph 0 0 0 gp %10g %10g %10g\n",disp.gper[i],disp.gvelo[i],disp.gvel[i],disp.ungvelo[i]);}

return 1;
}
//-------------------------------------------------------------------
int write_ASC_rf(modeldef &model,paradef &para,char *namemod,char *Rnamedisp,char *Lnamedisp,char *rfname,int Rsurflag, int Lsurflag)
{
  FILE *outmod,*Routdisp,*Loutdisp,*rfout;
  int i;
  if((outmod=fopen(namemod,"w"))==NULL)
  {
   printf("Cannot open file to write %s!!!\n",namemod);
   exit(0);
  }
  if((Routdisp=fopen(Rnamedisp,"w"))==NULL)
  {
   printf("Cannot open file to write %s!!!\n",Rnamedisp);
   exit(0);
  }
  if((Loutdisp=fopen(Lnamedisp,"w"))==NULL)
  {
   printf("Cannot open file to write %s!!!\n",Lnamedisp);
   exit(0);
  }
  if((rfout=fopen(rfname,"w"))==NULL)
  {
   printf("Cannot open file to write %s!!!\n",rfname);
   exit(0);
  }	  
  fprintf(outmod,"misfit:%g disp- %g %g %g %g rf- %g  L: %g disp- %g %g %g %g rf- %g\n",model.data.misfit,model.data.Rdisp.pmisfit,model.data.Rdisp.gmisfit,model.data.Ldisp.pmisfit,model.data.Ldisp.gmisfit,model.data.rf.misfit,model.data.L,model.data.Rdisp.pL,model.data.Rdisp.gL,model.data.Ldisp.pL,model.data.Ldisp.gL,model.data.rf.L);
  fprintf(outmod,"%d %d %d %d\n",para.Rnpara,para.Lnpara,model.ngroup,model.laym0.nlayer);
  for(i=0;i<para.Rnpara;i++)
        fprintf(outmod,"%10g %10g",para.Rparameter[i],para.Lparameter[i]);
  fprintf(outmod,"\n");
  for(i=0;i<model.laym0.nlayer;i++)
        fprintf(outmod,"%10g %10g %10g %10g %10g %10g %10g %10g\n",model.laym0.thick[i],model.laym0.vsv[i],model.laym0.vsh[i],model.laym0.vp[i],model.laym0.rho[i],model.laym0.qs[i],model.laym0.qp[i],model.laym0.vpvs[i]);

  if(Rsurflag>0){printf("Rsurflag>0, npper=%d ngper=%d\n",model.data.Rdisp.npper,model.data.Rdisp.ngper);write_ASCdisp_single(model.data.Rdisp, Routdisp);}
  if(Lsurflag>0){write_ASCdisp_single(model.data.Ldisp, Loutdisp);}

  for(i=0;i<model.data.rf.nrfo;i++){
  	for(int j=0;j<model.data.rf.tn.size();j++){
	  //In case the obs and pred RF have different time interval. Use first 10 sec only
	  if(pow((model.data.rf.to[i]-model.data.rf.tn[j]),2)<0.0001 and model.data.rf.to[i]<=10. and model.data.rf.to[i]>=0.){
		fprintf(rfout,"%5g %10g %10g %10g\n",model.data.rf.to[i],model.data.rf.rfo[i],model.data.rf.rfn[j],model.data.rf.unrfo[i]);		
	  }
	}
  
  }
  fclose(outmod);
  fclose(Routdisp);
  fclose(Loutdisp);
  fclose(rfout);
  return 1;
}


//-------------------------------------------------------------------
int write_initmodIso(char *foutnm,modeldef model){
//write initial model that would be used as input model for MC inversion
  int i,j;
  FILE *fout;
  if((fout=fopen(foutnm,"w"))==NULL){cout<<"### cannot open file to write "<<foutnm<<endl;exit(0);}
  for(i=0;i<model.ngroup;i++){
        fprintf(fout,"%5d %5d %8g %5d",i,model.groups[i].flag,model.groups[i].thick,model.groups[i].np);
        for(j=0;j<model.groups[i].np;j++)
          fprintf(fout," %8g",model.groups[i].Rvalue[j]);
	// there was a bug here, fixed on May 22,2012
	if(model.groups[i].flag==1){//layered model
		for(j=0;j<model.groups[i].np;j++)
			fprintf(fout," %8g",model.groups[i].ratio[j]);}
	fprintf(fout," %8g\n",model.groups[i].vpvs);
  }//for i
  fclose(fout);
  return 1;
}

//-------------------------------------------------------------------
int write_initmodAniso(char *foutnm,modeldef model){
//write initial model that would be used as input model for MC inversion
  int i,j;
  FILE *fout;
  if((fout=fopen(foutnm,"w"))==NULL){cout<<"### cannot open file to write "<<foutnm<<endl;exit(0);}
  for(i=0;i<model.ngroup;i++){
        fprintf(fout,"%5d %5d %8g %5d",i,model.groups[i].flag,model.groups[i].thick,model.groups[i].np);
        for(j=0;j<model.groups[i].np;j++)
          fprintf(fout," %8g %8g",model.groups[i].Rvalue[j],model.groups[i].Lvalue[j]);
	// there was a bug here, fixed on May 22,2012
	if(model.groups[i].flag==1){//layered model
		for(j=0;j<model.groups[i].np;j++)
			fprintf(fout," %8g",model.groups[i].ratio[j]);}
        fprintf(fout," %8g\n",model.groups[i].vpvs);
  }//for i
  fclose(fout);
  return 1;
}

