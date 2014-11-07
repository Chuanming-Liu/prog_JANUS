# this is used to fit the observed V(theta) with 1psi/2psi/4psi based on users choice.
# d=G*M --> M=inv(G.T*G)*G.T*d
# Cov(M)=inv(G.T,G), sqrt of the diagonal terms are the std of model parameters.

import sys
import string
import os
from numpy import * 

def readdata(indata):
 flagdata=0
 # read data -----
 lonlst=[]
 latlst=[]
 numlst=[]
 datalst=[] #2d array
 azilst=[]
 unclst=[]
 Npoint=0
 for line in open(indata):
	l=line.rstrip().split()
	if(flagdata==0): # read information
		Npoint=Npoint+1
		lon=float(l[0]);lat=float(l[1]);num=float(l[2])
		i=0
		flagdata=1
		#--change lon to -180~180
		if(lon>180):
			lon=lon-360
		lonlst.append(lon); latlst.append(lat); numlst.append(num)
		tdata=[]
		tazi=[]
		tunc=[]
		#print lon,lat,num
	else: # read data
		azi=float(l[0]);vel=float(l[1]);unc=float(l[2])
		i=i+1
		if(unc>2): #err
			print "### hey, are you reading info instead of data?!"
			print line
		if(i==num):
			flagdata=0 #finish reading data for this point
			datalst.append(tdata)
			azilst.append(tazi)
			unclst.append(tunc)
		tdata.append(vel)
		tazi.append(azi)
		tunc.append(unc)

 print "read in %d point"%Npoint
 return lonlst,latlst,numlst,datalst,azilst,unclst,Npoint

#===============================
def do_lstsq(data,azi,unc,id): # do the lstsq fitting, id is the choice of fitting
	N=len(unc)
	W=eye(N)*1./unc #weighting matrix
	G12=vstack([cos(azi),sin(azi),cos(2*azi),sin(2*azi),ones(N)]).T #the G matrix for 1psi+2psi fitting.
	G24=vstack([cos(2*azi),sin(2*azi),cos(4*azi),sin(4*azi),ones(N)]).T #the G matrix for 2psi+4psi fitting.
	G2=vstack([cos(2*azi),sin(2*azi),ones(N)]).T #the G matrix for 2psi fitting.
	G4=vstack([cos(4*azi),sin(4*azi),ones(N)]).T #the G matrix for 4psi fitting.
	G14=vstack([cos(azi),sin(azi),cos(4*azi),sin(4*azi),ones(N)]).T # the G matrix for 1psi+4psi fitting
	if id == '12psi':
		G=G12
	elif id == '24psi':
		G=G24
	elif id == '2psi':
		G=G2
	elif id == '4psi':
		G=G4
	elif id == '14psi':
		G=G14
	else:
		print "## Hey wrong input for id, type either 12psi, 24psi, 2psi, or 4psi!"
		sys.exit()
	
	data_w=dot(W,data)
	G_w=dot(W,G)
	M_w=linalg.lstsq(G_w,data_w)
	Cov_w=linalg.pinv(dot(G_w.T,G_w)) #pseudo inverse
	#print W,"\n",Cov_w
	#print linalg.pinv(dot(G.T,G))
	UncM_w=sqrt(Cov_w.diagonal())
	#print UncM_w
	return M_w[0],M_w[1],UncM_w	#return model, residual,and model_unc
#========= compute syn-data based on the inverted model; and also convert model to psi and amp  ===
def transfer(Ac,As,uncAc,uncAs):
	amp=sqrt(Ac**2+As**2)
	if(Ac**2+As**2<1E-30): #Ac=As=0
		phi=NaN
	elif (As**2<1E-30):#As=0
		if(Ac>0):
			phi=0.5*pi
		else:
			phi=1.5*pi
	elif(Ac**2<1E-30):#Ac=0
		if(As>0):
			phi=0.
		else:
			phi=pi
	else:	
		phi=arctan(Ac/As)
		#!!!!!!!! amp can be negative !!!!!!!!!!
        	if(As/cos(phi)<0.):
          		phi=phi+pi
	uncamp=sqrt((Ac*uncAc)**2+(As*uncAs)**2)/amp
	uncphi=1/(1+(Ac/As)**2)*sqrt((uncAc/As)**2+(Ac*uncAs/As/As)**2) #(rad)
	#uncamp=abs((Ac*uncAc+As*uncAs)/amp)
	#uncphi=abs((Ac*uncAs+As*uncAc)/amp/amp)
	return amp,phi,uncamp,uncphi

def do_syn(M,U,id,azi):
	# since fast axis = phi/X, so unc_of_FA=unc_phi/X 
	#t=array(range(0,361))*deg2rad
	t=azi
	var=[]
	varFA=[] #in this list, the directino is the fast axis direction FA=T/4-phi/X
	T1=2*pi;T2=pi;T4=pi/2;
	if (id=='12psi'):
		syn=M[0]*cos(t)+M[1]*sin(t)+M[2]*cos(t*2)+M[3]*sin(t*2)+M[4]
		psi1A,psi1t,unc1A,unc1t=transfer(M[0],M[1],U[0],U[1])
		psi2A,psi2t,unc2A,unc2t=transfer(M[2],M[3],U[2],U[3])
		psi2t=psi2t/2.
		unc2t=unc2t/2.
		
		var=[M[4],psi1A,psi1t,psi2A,psi2t]
		varFA=[M[4],psi1A,T1/4-psi1t,psi2A,T2/4-psi2t]
		varunc=[U[4],unc1A,unc1t,unc2A,unc2t]
		syn=psi1A*sin(t+psi1t)+psi2A*sin(2*(t+psi2t))+M[4]
	elif (id=='24psi'):
		syn=M[0]*cos(t*2)+M[1]*sin(t*2)+M[2]*cos(t*4)+M[3]*sin(t*4)+M[4]
		psi2A,psi2t,unc2A,unc2t=transfer(M[0],M[1],U[0],U[1])
		psi4A,psi4t,unc4A,unc4t=transfer(M[2],M[3],U[2],U[3])
		psi2t=psi2t/2.
		psi4t=psi4t/4.
		unc2t=unc2t/2.
		unc4t=unc4t/4.		

		var=[M[4],psi2A,psi2t,psi4A,psi4t]
		varFA=[M[4],psi2A,T2/4-psi2t,psi4A,T4/4-psi4t]
		varunc=[U[4],unc2A,unc2t,unc4A,unc4t]
		syn=psi2A*sin(2*(t+psi2t))+psi4A*sin(4*(t+psi4t))+M[4]
	elif (id =='2psi'):
		syn1=M[0]*cos(t*2)+M[1]*sin(t*2)+M[2]
		psi2A,psi2t,unc2A,unc2t=transfer(M[0],M[1],U[0],U[1])
		psi2t=psi2t/2.
		unc2t=unc2t/2.

		var=[M[2],psi2A,psi2t]
		varFA=[M[2],psi2A,T2/4-psi2t]
		varunc=[U[2],unc2A,unc2t]
		syn=psi2A*sin(2*(t+psi2t))+M[2]
		"""
		print syn1,"\n",syn
		print M[0],M[1],M[2]
		print psi2A,psi2t,psi2t*rad2deg,M[2]
		"""
	elif (id =='4psi'):
		syn=M[0]*cos(t*4)+M[1]*sin(t*4)+M[2]
		psi4A,psi4t,unc4A,unc4t=transfer(M[0],M[1],U[0],U[1])
		psi4t=psi4t/4.
		unc4t=unc4t/4.

		var=[M[2],psi4A,psi4t]
		varFA=[M[2],psi4A,T4/4-psi4t]
		varunc=[U[2],unc4A,unc4t]
		syn=psi4A*sin(4*(t+psi4t))+M[2]
	elif (id =='14psi'):
		syn=M[0]*cos(t)+M[1]*sin(t)+M[2]*cos(t*4)+M[3]*sin(t*4)+M[4]
		psi1A,psi1t,unc1A,unc1t=transfer(M[0],M[1],U[0],U[1])
		psi4A,psi4t,unc4A,unc4t=transfer(M[2],M[3],U[2],U[3])
		psi4t=psi4t/4.
		unc4t=unc4t/4

		var=[M[4],psi1A,psi1t,psi4A,psi4t]
		varFA=[M[4],psi1A,T1/4-psi1t,psi4A,T4/4-psi4t]
		varunc=[U[4],unc1A,unc1t,unc4A,unc4t]
		syn=psi1A*sin(t+psi1t)+psi4A*sin(4*(t+psi4t))+M[4]
	return syn,var,varunc,varFA
		


#================ Main ===========
if (len(sys.argv)!=5):
	print "Input:\n 1) input data file (e.g., 10s_iso_ani_v1.ani)\n 2)output map file name (writes the fit_model at every point on the map\n 3)output point file directory (it creates single file for every point on the map!! Type NO to turn off this option)\n 4) choice of fit, 2psi, 4psi, 12psi, or 24psi"
	sys.exit()
#------------attention that psi is different from fast axis! 
# Xpsi: fitting_curve=A*sin(X*(t+Xpsi))+C ==> fast_axis = pi/2/X-Xpsi=T/4-Xpsi
# this can be understood if you draw a picture; +Xpsi moves the sin curve left by Xpsi, and 1/4 period of the sin curve is pi/2/X; the positive peak was at pi/2/X,then is moved left by Xpsi --> the positive peak is at pi/2/X=Xpsi

indata=sys.argv[1] # input data file
outpsimap=sys.argv[2] #output map file: )lon )lat )num )chi2 )V0 )amp )amp_unc  )psi(in deg) )unc_psi ...
outpointdir=sys.argv[3] 
id=sys.argv[4]
#output information for each point:
#V0 amp1 psi1 amp2 psi2 ... { can be used to compute syn vel(t)=SUM_X[ampX*sin(X*(t+psiX))]+V0 }
#azi vel unc {from the read in data}

"""
input data format:

-123.600000 44.400000 8 	#lon lat num_of azi_bin
110.000000 3.267035 0.011641	#azi vel unc
130.000000 3.245972 0.006339
...
"""
if (outpointdir !="NO"):
	os.system("if [ ! -d %s ]; then mkdir %s; fi"%(outpointdir,outpointdir))

# read in data ---
if ( not os.path.exists(indata)):
	print "input file %s does not exist!!!"%(indata)
	sys.exit()
lonlst,latlst,numlst,datalst,azilst,unclst,Npoint=readdata(indata)

global deg2rad,rad2deg
deg2rad=pi/180.
rad2deg=180./pi

outmap=open(outpsimap,"w")
for i in range(Npoint):
	if(i%500==0):
		print "point %d"%(i)
	lon=lonlst[i]; lat=latlst[i]; num=numlst[i]
	if (num<4): # too few points
		continue
	data=array(datalst[i]) #
	azi=array(azilst[i])*deg2rad
	unc=array(unclst[i]) #
	#do inversion---------
	model,residual,Munc=do_lstsq(data,azi,unc,id)
	#if Munc[0]<0:
	#	print "point %g_%g gives a sigular matrix, no meanful unc is computed!\n"%(lon,lat)
	#transfer parameters to other format-----
	syn,var,varunc,varFA=do_syn(model,Munc,id,azi)

	#--- filter the input data, only use data that's not far from syn (dist<2*data_unc)
	"""
	dataflt=[];uncflt=[];aziflt=[];
	for j in range(int(num)):
		if(abs(data[j]-syn[j])<2*unc[j]):
			dataflt.append(data[j])
			uncflt.append(unc[j])
			aziflt.append(azi[j])
	num=len(dataflt)
	if len(dataflt)<4: #too few points
		continue
	dataflt=array(dataflt);uncflt=array(uncflt);aziflt=array(aziflt)
	model,residual,Munc=do_lstsq(dataflt,aziflt,uncflt,id)
	syn,var,varunc,varFA=do_syn(model,Munc,id,aziflt) # here I use azi instead of aziflt in purpose, keep the number of syn the same as the # of origional input data, otherwise, could use aziflt to get the reduced chi2
	"""
	#compute misfit ----------
	variance=sum((syn-data)**2)/num
	chi2=sum(((syn-data)/unc)**2)/num #chi2/num--> reduced chi2

	#write out result -----------
	# the map file:lon lat num chi2 v0 dV0 Apsi1 dApsi1 FA1 dFA1 (Apsi2 dApsi2...)  (deg)
	# the point file: Line1: psitype V0 Apsi1 phi/X1 (Apsi2 phi/X2) (rad, here it is psiX, NOT FA)
	#		  Line 2~end: azimuth vel unc (input data, deg or rad)

	str = " "
	str2= " "
	for j in range(len(var)):
		if (j==0 or j%2!=0): #1st, or odd number, V0 or amplitude
			str="%s %5g %5g"%(str,var[j],varunc[j])
			str2="%s %5g"%(str2,var[j])
		else: # even number other than 0, theta, rad --> deg
			# for 2psi, T=pi; if fit 1psi, then change T to 2pi
			while(var[j]>pi):
				var[j]=var[j]-pi
			while(var[j]<0):
				var[j]=var[j]+pi
			while(varFA[j]>pi):
				varFA[j]=varFA[j]-pi
			while(varFA[j]<0):
				varFA[j]=varFA[j]+pi
			str="%s %5g %5g"%(str,varFA[j]*rad2deg,varunc[j]*rad2deg)
			str2="%s %5g"%(str2,var[j])
	#---write map file ---------------------
	#---lon, lat, azi_num, chi2, v0, dv0, A, dA, psi, dpsi, (in degree)...
	#print lon,lat,num,variance,chi2
	outmap.write("%5g %5g %5g %5g %s\n"%(lon,lat,num,chi2,str))

	#---write file for each point ---------- the azi[X] are in rad	(no unc info)
	if (outpointdir !="NO"):
		outpoint=open("%s/point_%.1f_%.1f.txt"%(outpointdir,lon,lat),"w")
		outpoint.write("%s %s\n"%(id,str2))  
		for j in range(int(num)):
			# for 2psi, T=pi; if fit 1psi, then change T to 2pi
			while(azi[j]>2*pi):
				azi[j]=azi[j]-2*pi
			while(azi[j]<0):
				azi[j]=azi[j]+2*pi
			outpoint.write("%5g %5g %5g\n"%(azi[j]*rad2deg,data[j],unc[j]))
		outpoint.close()
outmap.close()	


