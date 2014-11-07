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
	
