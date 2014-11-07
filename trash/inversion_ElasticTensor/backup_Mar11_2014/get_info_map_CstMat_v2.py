# this is used get the map view of some information 
# for every point, get the info of that para
# theta, phi in the crust; phi and theta in the mantle
# this version, use vsv&vsh file to compute ani value. Also, it reads in the multiple group (_phigp0/1) files

import os
import os.path
import sys
from math import *

############################################
fout=sys.argv[1]
invid = sys.argv[2]
flagBest = int(sys.argv[3]) # write out the best-fitting para (1) or average para (0); the para only indicate the theta&phi or cos&sin; not the amplitude of radial anisotropy
surffix=sys.argv[4]
indir = "/projects/jixi7887/work/US/inv_ET_BS"
indirdata = "/lustre/janus_scratch/jixi7887/US/inv_ET_BS/inv_%s"%(invid)
#####
fpoint = "%s/data_disp_Jan7/point_info_v2.txt"%(indir)
#fpoint = "%s/point_temp1.txt"%(indir)
#####
paraidlst=[5,6,75,76] #theta,phi in the crust, and theta, phi in the mantle
idAZ=[2,3] # tells which paraidlst[id] is the thetaM and phiM
#fout = "./para_map.txt"
#invid=int(sys.argv[2]) #fpara="%s/%s_inv_%d/bin_avg/para_%.1f_%.1f.txt"%(indir,name,invid,lon,lat)
paranmlst=['thetaC','phiC','thetaM','phiM'] #name of each para
flagCS2AP=0 # (1)indicate the read in is cos&sin, need to change them into amp&phi; (0)otherwise, indicate the read in is the value we want, no further computation is needed
############################################

#---------------------------------------------------
def readfile (file):
        deplst=[];vlst=[];unclst=[]
        i=0
        for line in open(file):
		l=line.rstrip().split()
                i+=1
                if('nan' in line):
                        print line
                        continue
                if(i==2):
                        Hmoho=float(l[1]);Hmohounc=float(l[2])
                elif(i>2):
                        deplst.append(float(l[2]))
                        vlst.append(float(l[3]))
                        unclst.append(float(l[6]))
        return Hmoho,Hmohounc,deplst,vlst,unclst

#---------------------------------------------------
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
		vs=sqrt((2*vsv**2+vsh**2)/3.)
                if(dep<depmin):
                        ani= (vsh-vsv)*200/(vsv+vsh)
			#ani=(vsh-vsv)*100/vs;
                        A=(vsh-vsv);B=(vsv+vsh);
                        dA=sqrt(uncvsh**2+uncvsv**2);dB=dA
                        unc=fabs(200/B**2*(B*dA-A*dB))
                        #print dep,ani,unc
                        ani1=ani;unc1=unc;dep1=dep;
                elif(dep<Hmohovsv):
                        ani=(vsh-vsv)*200/(vsv+vsh)
			#ani=(vsh-vsv)*100/vs;
                        A=(vsh-vsv);B=(vsv+vsh);
                        dA=sqrt(uncvsh**2+uncvsv**2);dB=dA
                        unc=fabs(200/B**2*(B*dA-A*dB))
                        sum=sum+(ani+ani1)*(dep-dep1)*0.5
                        sumunc=sumunc+(unc+unc1)*(dep-dep1)*0.5
                        H=H+(dep-dep1)
                        ani1=ani;unc1=unc;dep1=dep;
                else:# since mantle has constant ani, so only take the 1st value is enough.############ may need modification ##
                        aniavg=sum/H;
                        uncavg=sumunc/H;
                        animantle=(vsh-vsv)*200/(vsv+vsh);
			#animantle=(vsh-vsv)*100/vs;
                        A=(vsh-vsv);B=(vsv+vsh);
                        dA=sqrt(uncvsh**2+uncvsv**2);dB=dA
                        uncmantle=fabs(200/B**2*(B*dA-A*dB))
                        break
        return aniavg,uncavg,animantle,uncmantle        
#---------------------------------------------------
def get_aniC_aniM(fani):
        #computes the crustal average ani (average over depmin~moho), and mantle ani; (this is only for the mantle with constant ani)
        i=0
        sum=0;sumunc=0;dep1=0;ani1=0;unc1=0;H=0
        for line in open(fani):
                i=i+1
                l=line.rstrip().split()
                if 'nan' in line:
                        print line
                        continue
                if (i==2):
                        Hmoho=float(l[1]) #depth of moho
                        depmin=Hmoho/4
                elif(i>2):
                        dep=float(l[2])
                        if(dep<depmin):
                                ani=float(l[3])
                                unc=float(l[6])
                                ani1=ani;unc1=unc;dep1=dep;
                        elif(dep<=Hmoho):
                                ani=float(l[3])
                                unc=float(l[6])
                                sum=sum+(ani+ani1)*(dep-dep1)/2.
                                sumunc=sumunc+(unc+unc1)*(dep-dep1)/2.
                                H=H+(dep-dep1)
                                ani1=ani;unc1=unc;dep1=dep;
                        else:# since mantle has constant ani, so only take the 1st value is enough.
                                #print "H=%f moho=%f"%(H,Hmoho)
                                aniavg=sum/H;
                                uncavg=sumunc/H;
                                animantle=float(l[3])
                                uncmantle=float(l[6])
                                break
        return aniavg,uncavg,animantle,uncmantle

#---------------------------------------------------
def get_misfit(fani):
	i=0
	for line in open(fani):
		i=i+1
		if(i==2):
			l=line.rstrip().split()
			misfit = float(l[1])
			break
			
	return misfit

#----Main
print "writting %s"%fout
out=open(fout,"w")
stnmlst=[];lonlst=[];latlst=[]
Npoint=0
for line in open(fpoint):
	l=line.rstrip().split()
	stnmlst.append(l[0])
	lonlst.append(float(l[1]))
	latlst.append(float(l[2]))
	Npoint=Npoint+1

#I09A_inv_2/bin_avg/para_-118.0_44.0.txt:     0   1.4080   0.0000
Nparaid=len(paraidlst)
#paravaluelst2=[]
#parastdlst2=[]
for i in range(Npoint):
	paravaluelst1=[]
	parastdlst1=[]

	stnm=stnmlst[i]
	lon=lonlst[i]
	lat=latlst[i]
	name = "%s_%.1f_%.1f"%(stnm,lon,lat)
	fpara="%s/%s_inv_%s/bin_avg/para_%.1f_%.1f.txt%s"%(indirdata,name,invid,lon,lat,surffix)
	fvsv="%s/%s_inv_%s/bin_avg/vsv_%.1f_%.1f.txt%s"%(indirdata,name,invid,lon,lat,surffix)
	fvsh="%s/%s_inv_%s/bin_avg/vsh_%.1f_%.1f.txt%s"%(indirdata,name,invid,lon,lat,surffix)
	if (not( os.path.exists(fpara) and os.path.exists(fvsv) and os.path.exists(fvsh))):
		print "file %s doesn't exist! skip\n"%(fpara)
		continue
	count=0
	for line in open(fpara):
		l=line.rstrip().split()
		ipara=int(l[0])
		if ipara not in paraidlst:
			continue
		if ( flagBest == 0 ):
			paravaluelst1.append(float(l[1])) #average para
		else:
			paravaluelst1.append(float(l[3])) #best fitting para
		parastdlst1.append(float(l[2]))				
		count=count+1
		if(count==Nparaid):
			break;
  	if ( flagCS2AP==1): # indicate the read in is cos&sin, need to change them into amp&phi; otherwise, indicate the read in is the value we want, no further computation is needed
	  for i in range(Nparaid):
		if i == idAZ[0]:
			cos=paravaluelst1[i]
			unccos=parastdlst1[i]
		elif i==idAZ[1]:		
			sin=paravaluelst1[i]
			uncsin=parastdlst1[i]
	  os.system("/home/jixi7887/progs/jy/inversion_ElasticTensor/cs2ap %f %f %f %f > temp.txt"%(cos,sin,unccos,uncsin))	
	  for line in open("temp.txt"): #amp AZ uncamp uncAZ
		l=line.rstrip().split("_")
		amp=float(l[0]);AZ=float(l[1]);uncamp=float(l[2]);uncAZ=float(l[3])
	  paravaluelst1[idAZ[0]]=amp
	  parastdlst1[idAZ[0]]=uncamp
	  paravaluelst1[idAZ[1]]=AZ
	  parastdlst1[idAZ[1]]=uncAZ

	#paravaluelst2.append(paravaluelst1)
	#parastdlst2.append(parastdlst1)
	
	out.write("%s %7.1f %7.1f "%(stnm,lon,lat))
	for i in range(Nparaid):
		#---change the phiC and phiM to [0,180] range
		Tphi=180
		if(paranmlst[i]=="phiC" or paranmlst[i]=="phiM"):
			while(paravaluelst1[i]>Tphi):
				paravaluelst1[i]-=Tphi
		out.write(" %5s %8.4f %8.4f "%(paranmlst[i],paravaluelst1[i],parastdlst1[i]))

	#---compute average ani for crust, mantle
	# ../I11A_inv_2/bin_avg/ani_-116.0_44.0.txt
	#sedi 0.238023 2.27762e-05
	#moho  33.1773 0.00370676
       	#0        0        0        0        0        0        0
    	#0.05        0     0.05        0     0.05        0        0
	"""
	fani="%s/%s_inv_%s/bin_avg/ani_%.1f_%.1f.txt"%(indirdata,name,invid,lon,lat)
	aniavg,uncavg,animantle,uncmantle=get_aniC_aniM(fani)
	out.write(" aniC %8.4f %8.4f aniM %8.4f %8.4f"%(aniavg,uncavg,animantle,uncmantle))

	fani="%s/%s_inv_%s/bin_avg/ani_%.1f_%.1f.txt_effTI"%(indirdata,name,invid,lon,lat)
	aniavg,uncavg,animantle,uncmantle=get_aniC_aniM(fani)
	out.write(" aniCeffTI %8.4f %8.4f aniMeffTI %8.4f %8.4f"%(aniavg,uncavg,animantle,uncmantle))
	"""
	aniavg,uncavg,animantle,uncmantle=get_aniC_aniM_fromV(fvsv,fvsh)
	out.write(" aniC %8.4f %8.4f aniM %8.4f %8.4f"%(aniavg,uncavg,animantle,uncmantle))

        fvsv="%s/%s_inv_%s/bin_avg/vsv_%.1f_%.1f.txt_effTI%s"%(indirdata,name,invid,lon,lat,surffix)
        fvsh="%s/%s_inv_%s/bin_avg/vsh_%.1f_%.1f.txt_effTI%s"%(indirdata,name,invid,lon,lat,surffix)
        aniavg,uncavg,animantle,uncmantle=get_aniC_aniM_fromV(fvsv,fvsh)
        out.write(" aniCeffTI %8.4f %8.4f aniMeffTI %8.4f %8.4f"%(aniavg,uncavg,animantle,uncmantle))


	#out.write("\n")

	#../L32A_inv_8/Animod_0_L32A_-98.0_42.0.txt
	fani = "%s/%s_inv_%s/Animod_0_%s_%.1f_%.1f.txt%s"%(indirdata,name,invid,stnm,lon,lat,surffix)
	if (os.path.exists(fani)):
		misfit=get_misfit(fani)
		out.write(" misfit_avg %8.4f "%(misfit))
	else:
		out.write(" misfit_avg nan ")
		
	fani = "%s/%s_inv_%s/AnimodB_0_%s_%.1f_%.1f.txt"%(indirdata,name,invid,stnm,lon,lat)
	if (os.path.exists(fani)):
		misfitB=get_misfit(fani)
		out.write(" misfit_best %8.4f "%(misfitB))
	else:
		out.write(" misfit_best nan ")
	out.write("\n")	
	#out.write(" misfit_avg %8.4f misfit_best %8.4f\n"%(misfit,misfitB))
out.close()



