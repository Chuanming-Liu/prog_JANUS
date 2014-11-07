def readfile (file):
	deplst=[];vlst=[];unclst=[]
	i=0
	for line in open(file):
		i+=1
		if('nan' in line):
                        print line
                        continue
                if(i==2):
			Hmoho=float(l[1]),Hmohounc=float(l[2])
		elif(i>2):	
			deplst.append(float(l[2]))
			vlst.append(float(l[3]))
			unclst.append(float(l[6]))
	return Hmoho,Hmohounc,deplst,vlst,unclst
#======================	
def get_aniC_aniM_fromV(fvsv,fvsh):
        #computes the crustal average ani (average over depmin~moho), and mantle ani; (this is only for the mantle with constant ani); this computes ani from vsv and vsh, not directly read in ani
        i=0
        sum=0;sumunc=0;dep1=0;ani1=0;unc1=0;H=0
	
	Hmohovsv,Hmohouncvsv,depvsvlst,vsvlst,uncvsvlst=readfile(fvsv)
	Hmohovsh,Hmohouncvsh,depvshlst,vshlst,uncvshlst=readfile(fvsh)

	if(fabs(Hmohovsv-Hmohovsh)>0.5*(Hmohouncvsv+Hmohouncvsh)):
		print "Hey, moho is different in fvsv and fvsh:\n\t %s~%g,\n\t %s~%g\n"%(fvsv,Hmohovsv,fvsh,Hmohovsh)
		sys.exit()

	depmin=Hmohovsv*0.25; ################### the start depth for computing anisotropy average
	depmax=Hmohovsv; ################### the end depth for computing anisotropy average
	for i in range(len(depvsvlst)):
		if(fabs(depvsvlst[i]-depvshlst[i])>1e-5):
			print "inconsistent depth! %g vs. %g\n"%(depvsvlst[i],depvshlst[i])
			sys.exit()
		vsv=vsvlst[i];vsh=vshlst[i];uncvsv=uncvsvlst[i];uncvsh=uncvshlst[i]
		dep=depvshlst[i];
		if(dep<depmin):
			ani= (vsh-vsv)*200/(vsv+vsh)
			A=(vsh-vsv);B=(vsv+vsh);
			dA=sqrt(uncvsh**2+uncvsv**2);dB=dA
			unc=fabs(200/B**2*(B*dA-A*dB))
			print dep,ani,unc
			ani1=ani;unc1=unc;dep1=dep;
		elif(dep<Hmoho):
			ani=(vsh-vsv)*200/(vsv+vsh)
			A=(vsh-vsv);B=(vsv+vsh);
                        dA=sqrt(uncvsh**2+uncvsv**2);dB=dA
                        unc=fabs(200/B**2*(B*dA-A*dB))
			sum=sum+(ani+ani1)*(dep-dep1)*0.5
			sumunc=sumunc+(unc+unc1)*(dep-dep1)*0.5
			H=H+(dep-dep1)
			ani1=ani;unc1=unc;dep1=dep;
		else:# since mantle has constant ani, so only take the 1st value is enough.
			aniavg=sum/H;
			uncavg=sumunc/H;
			animantle=(vsh-vsv)*200/(vsv+vsh);
			A=(vsh-vsv);B=(vsv+vsh);
                        dA=sqrt(uncvsh**2+uncvsv**2);dB=dA
                        uncmantle=fabs(200/B**2*(B*dA-A*dB))
			break
	return aniavg,uncavg,animantle,uncmantle	
#===================
