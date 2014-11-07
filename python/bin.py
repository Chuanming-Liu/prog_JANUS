import math;
import sys;
import string;

#if (len(sys.argv)< 5):
#    print "input [in.f] [n.col] [width] [min] [max]";
if (len(sys.argv)< 3):
    print "input [in.f] [n.col] [width]";
    sys.exit();

ncol = int(sys.argv[2]); #starts from 0
width = float(sys.argv[3]);
#tmin = float(sys.argv[4]);
#tmax = float(sys.argv[5]);
tmin=1E10
tmax=-1E10
for line in open(sys.argv[1]):
	l=line.rstrip().split()
	tvalue=float(l[ncol])
	if tvalue<tmin:
		tmin=tvalue
	if tvalue>tmax:
		tmax=tvalue
tmin=int(tmin/width)*width
n_bin = int((tmax-tmin)/width-0.5)+1;

mins = [];
maxs = [];
mids = [];
nn = [];
np = [];
for i in range(n_bin):
    mins.append(tmin + i*width);
    maxs.append(tmin + (i+1)*width);
    mids.append(tmin + (i+0.5)*width);
    nn.append(0);
    np.append(0);

n_t = 0;
for l1 in open(sys.argv[1]):
    l2 = l1.rstrip().split();
    tvalue = float(l2[ncol]); 
    for i in range(n_bin):
        if tvalue >= mins[i] and tvalue < maxs[i]:
             nn[i] = nn[i]+1;
             n_t = n_t + 1;
        else:
             continue;

for i in range(n_bin):
    t_r = 100.*float(nn[i])/float(n_t);
    np[i] = t_r;
    print mins[i],maxs[i],mids[i],nn[i],np[i];
