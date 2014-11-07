# from the output of get_info_map_CstMat_v2/3.py compare the output for gp0 and gp1, then get the point that has only 1 group of result 

import sys
import os

if(len(sys.argv)!=4):
	print "Input filegp0, filegp1, name of output file"
	sys.exit()

filegp0 = sys.argv[1]
filegp1 = sys.argv[2]
foutnm = sys.argv[3]

fY0N1="%s_Y0N1.txt"%foutnm
fY1N0="%s_Y1N0.txt"%foutnm
outY0N1=open(fY0N1,"w")
outY1N0=open(fY1N0,"w")
print "write diff point files: %s %s"%(fY0N1,fY1N0)
#N02D  -123.0    41.0  thetaC  19.0325   8.7723   phiC 104.4161  15.4731  thetaM  19.0273   4.3816   phiM  64.3324   9.0248  aniC   3.8231   2.0187 aniM   6.5007   2.5857 aniCeffTI   3.8613   2.0078 aniMeffTI   6.5244   2.5785 misfit_avg   2.1982  misfit_best nan 

plst0=[]
for line in open(filegp0):
	l=line.rstrip().split()
	str=" %s %s %s "%(l[0],l[1],l[2])
	plst0.append(str)

plst1=[]
for line in open(filegp1):
	l=line.rstrip().split()
	str=" %s %s %s "%(l[0],l[1],l[2])
	plst1.append(str)

iY0N1=0
for i in range(len(plst0)):
	str0=plst0[i]
	if str0 not in plst1:
		outY0N1.write("%s\n"%str0)
		iY0N1+=1
iY1N0=0
for i in range(len(plst1)):
	str1=plst1[i]
	if str1 not in plst0:
		outY1N0.write("%s\n"%str1)
		iY1N0+=1
print "%d points in file0 not file1, and %d in file1 not file0"%(iY0N1,iY1N0)

outY0N1.close()
outY1N0.close()
if(iY0N1==0):
	os.system("rm -f %s"%fY0N1)
if(iY1N0==0):
	os.system("rm -f %s"%fY1N0)

