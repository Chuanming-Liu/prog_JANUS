# this is used get the map view of some information 
# for every point, get the info of that para
# theta, phi in the crust; phi and theta in the mantle
# this

import os
import os.path
import sys

############################################
fout=sys.argv[1]
invid = sys.argv[2]
flagBest = int(sys.argv[3]) # write out the best-fitting para (1) or average para (0); the para only indicate the theta&phi or cos&sin; not the amplitude of radial anisotropy
indir = "/projects/jixi7887/work/US/inv_ET_BS"
indirdata = "/lustre/janus_scratch/jixi7887/US/inv_ET_BS/inv_%s"%(invid)
plotflag=0
if (plotflag==1): #the single point testing situation
	fpoint = "%s/point_info_AGU.txt"%(indir)
else: # the multipoint inversion situation
	fpoint = "%s/data_disp_Jan7/point_info_v2.txt"%(indir)
	#fpoint = "%s/point_temp_even.txt"%(indir)
paraidlst=[5,6,75,76] #theta,phi in the crust, and theta, phi in the mantle
idAZ=[2,3] # tells which paraidlst[id] is the thetaM and phiM
#fout = "./para_map.txt"
#invid=int(sys.argv[2]) #fpara="%s/%s_inv_%d/bin_avg/para_%.1f_%.1f.txt"%(indir,name,invid,lon,lat)
paranmlst=['thetaC','phiC','thetaM','phiM'] #name of each para
flagCS2AP=0 # (1)indicate the read in is cos&sin, need to change them into amp&phi; (0)otherwise, indicate the read in is the value we want, no further computation is needed
############################################

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
	if (plotflag==1):
		name = stnm
	else:
		#name = stnm
		name = "%s_%.1f_%.1f"%(stnm,lon,lat)
	fpara="%s/%s_inv_%s/bin_avg/para_%.1f_%.1f.txt%s"%(indirdata,name,invid,lon,lat,surffix)
	if (not os.path.exists(fpara)):
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
		out.write(" %5s %8.4f %8.4f "%(paranmlst[i],paravaluelst1[i],parastdlst1[i]))

	#---compute average ani for crust, mantle
	"""
	 ../I11A_inv_2/bin_avg/ani_-116.0_44.0.txt
	sedi 0.238023 2.27762e-05
	moho  33.1773 0.00370676
       	0        0        0        0        0        0        0
    	0.05        0     0.05        0     0.05        0        0
	"""
	fani="%s/%s_inv_%s/bin_avg/ani_%.1f_%.1f.txt"%(indirdata,name,invid,lon,lat)
	aniavg,uncavg,animantle,uncmantle=get_aniC_aniM(fani)
	out.write(" aniC %8.4f %8.4f aniM %8.4f %8.4f"%(aniavg,uncavg,animantle,uncmantle))

	fani="%s/%s_inv_%s/bin_avg/ani_%.1f_%.1f.txt_effTI"%(indirdata,name,invid,lon,lat)
	aniavg,uncavg,animantle,uncmantle=get_aniC_aniM(fani)
	out.write(" aniCeffTI %8.4f %8.4f aniMeffTI %8.4f %8.4f"%(aniavg,uncavg,animantle,uncmantle))
	#out.write("\n")

	#../L32A_inv_8/Animod_0_L32A_-98.0_42.0.txt
	fani = "%s/%s_inv_%s/Animod_0_%s_%.1f_%.1f.txt"%(indirdata,name,invid,stnm,lon,lat)
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



