# this is used get the map view of some information 
# for every point, get the info of that para
# theta, phi in the crust; amp and phi in the mantle

import os
import os.path
import sys

############################################
indir = "/projects/jixi7887/work/US/inv_ET"
plotflag=1
if (plotflag==1): #the single point testing situation
	fpoint = "%s/point_info_AGU.txt"%(indir)
else: # the multipoint inversion situation
	fpoint = "%s/data_disp_v2/point_info_v2.txt"%(indir)
paraidlst=[7,8,35,36] #theta,phi in the crust, and AZcos, AZsin in the mantle
idAZ=[2,3] # tells which paraidlst[id] is the AZcos and AZsin
#fout = "./para_map.txt"
fout=sys.argv[1]
invid=int(sys.argv[2])
paranmlst=['thetaC','phiC','ampM','phiM'] #name of each para
############################################

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
		name = "%s_%.1f_%.1f"%(stnm,lon,lat)
	fpara="%s/%s_inv_%d/bin_avg/para_%.1f_%.1f.txt"%(indir,name,invid,lon,lat)
	#fpara="%s/%s_%.1f_%.1f_inv_%d/bin_avg/para_%.1f_%.1f.txt"%(indir,stnm,lon,lat,invid,lon,lat)
	if (not os.path.exists(fpara)):
		print "file %s doesn't exist! skip\n"%(fpara)
		continue
	count=0
	for line in open(fpara):
		l=line.rstrip().split()
		ipara=int(l[0])
		if ipara not in paraidlst:
			continue
		paravaluelst1.append(float(l[1])) #average para
		#paravaluelst1.append(float(l[3])) #best fitting para
		parastdlst1.append(float(l[2]))				
		count=count+1
		if(count==Nparaid):
			break;
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
	fani="%s/%s_inv_%d/bin_avg/ani_%.1f_%.1f.txt"%(indir,name,invid,lon,lat)
	aniavg,uncavg,animantle,uncmantle=get_aniC_aniM(fani)
	out.write(" aniC %8.4f %8.4f aniM %8.4f %8.4f"%(aniavg,uncavg,animantle,uncmantle))

	fani="%s/%s_inv_%d/bin_avg/ani_%.1f_%.1f.txt_effTI"%(indir,name,invid,lon,lat)
	aniavg,uncavg,animantle,uncmantle=get_aniC_aniM(fani)
	out.write(" aniCeffTI %8.4f %8.4f aniMeffTI %8.4f %8.4f"%(aniavg,uncavg,animantle,uncmantle))
	#out.write("\n")

	#../L32A_inv_8/Animod_0_L32A_-98.0_42.0.txt
	fani = "%s/%s_inv_%d/Animod_0_%s_%.1f_%.1f.txt"%(indir,name,invid,stnm,lon,lat)
	if (os.path.exists(fani)):
		misfit=get_misfit(fani)
		out.write(" misfit_avg %8.4f "%(misfit))
	else:
		out.write(" misfit_avg nan ")
		
	fani = "%s/%s_inv_%d/AnimodB_0_%s_%.1f_%.1f.txt"%(indir,name,invid,stnm,lon,lat)
	if (os.path.exists(fani)):
		misfitB=get_misfit(fani)
		out.write(" misfit_best %8.4f "%(misfitB))
	else:
		out.write(" misfit_best nan ")
	out.write("\n")	
	#out.write(" misfit_avg %8.4f misfit_best %8.4f\n"%(misfit,misfitB))
out.close()



