#!/usr/bin/python
# find the end perid according to the SNR distribution, the use this as the end T to do AFTAN. ie, each path has its own end T.
# The eT cannot be bigger than eTmax; track from the tail of snr curve, find the period T with snr > snr_cri, set eT=T ;and for those path with low SNR (snr < 5 even at eTmin), set the eT to be eTmin => This may still introduce some 2 pi ambiguity.

import sys
import string
import os
from math import *
#from numpy import *
##################### parameter ###########################
snr_cri = 5.0
eTmax = 150.0 # the maxT of model is 50sec  the minT of model is 8sec
eTmin = 50.0
bT = 14.0
fstalst = sys.argv[1] #only read the 1st column stnm
fevtlst = sys.argv[2]#only read the 1st column evtnm
datadirmain = sys.argv[3]#datadirmain/year/evnm/evnm.stnm.LHT.sac(_snr.txt) , in the code i'll set datadir = datadirmain/year
outdir =  sys.argv[4]#store the aftan input information file
piO4= sys.argv[5] # the phase shift
preddispdir=sys.argv[6] #preddispdir/evnm/evnm.stnm.dat is the predicted dispersion curve
###########################################################
os.system("mkdir -p %s"%(outdir))
stnm = [];
Nst=0
for line in open(fstalst):
        l=line.rstrip().split()
   #if l[0]=="AHCHZ" or l[0]=="YNMEL":
	if l[0] in stnm:
		print "repeat station %s"%l[0]
		continue
	stnm.append(l[0])
	Nst=Nst+1
evnm=[];
evyear=[]
Nev=0
for line in open(fevtlst):
	l=line.rstrip().split()
	if l[0] in evnm:
		print "repeat event %s"%(l[0])
		continue
	evnm.append(l[0])
	evyear.append(l[0][0:4])
	Nev=Nev+1

print "read in %d stations, %d events"%(Nst,Nev)
for ne in range(Nev):
	evt=evnm[ne]
	fout = open("%s/Infosaclst_eT_%s_snr%d_T%d_%d.txt"%(outdir,evt,int(snr_cri),int(eTmin),int(eTmax)),"w")
	print evt,ne
	datadir="%s/%s"%(datadirmain,evyear[ne])
	for ns in range(Nst):
		sta = stnm[ns];
		##### depending on the data path, need to reset this part ###########
		cornm = "%s.%s.LHT.sac"%(evt,sta)
		fin = "%s/%s/%s_snr.txt"%(datadir,evt,cornm)
		preddispnm="%s/%s/%s.%s.dat"%(preddispdir,evt,evt,sta)
		perlst=[];snrlst=[]
		Nsnr=0 #check the existance of both SAC_s and snr files, though extract info from the snr file only.
		#print "%s/%s/%s"%(datadir,evt,cornm),fin
		#sys.exit()
		cornm="%s/%s/%s"%(datadir,evt,cornm)
		#if ( not os.path.exists("%s/%s/%s"%(datadir,evt,cornm))) or ( not os.path.exists(fin)): 
		if ( not os.path.exists(cornm)) or ( not os.path.exists(fin)): 
			continue;
		if (not os.path.exists(preddispnm)):
			print "Cannot find preddisp for %s.%s.dat"%(evt,sta)
			continue
		#"""
		for line in open(fin):
			l=line.rstrip().split()
			perlst.append(float(l[0]))
			snrlst.append(float(l[1]))
			Nsnr=Nsnr+1
		flag = 0
		for k in range(Nsnr-1,-1,-1):
			if snrlst[k]<snr_cri:
				continue;
			if perlst[k]>eTmax:
				eT=eTmax;flag=1;
			elif perlst[k] > eTmin:
				eT=ceil(perlst[k]);flag=1;
			break;
		#%g 1.5 5.0 14.0 100.0 40.0 1.0 1.0 0.2 1.0 0 %s %s # cv's Rayleigh_eqk para
		if flag == 0:
                        #fout.write("3 1.5 5 %f %f  20 1 0.5 0.2 2 %s\n"%(bT,eTmin,cornm))
                        fout.write("%s 1.5 5 %f %f  40.0 1.0 1.0 0.2 1.0 0 %s %s\n"%(piO4,bT,eTmin,cornm,preddispnm))
			continue;
                fout.write("%s 1.5 5 %f %f  40.0 1.0 1.0 0.2 1.0 0 %s %s\n"%(piO4,bT,eT,cornm,preddispnm))
		#fout.write("3 1.5 5 %f %f 20 1 0.5 0.2 2 %s\n"%(bT,eT,cornm))
		#"""
	fout.close()				
			
