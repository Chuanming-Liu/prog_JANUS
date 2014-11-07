#compute mean and std from a distribution
# do a reselction
import sys;
import math;
import string;

def get_mean(array):
	return sum(array)/len(array);
def get_std(array,mean):
	N=len(array)
	sum=0.
	for i in range(N):
		sum=(array[i]-mean)**2+sum
	return math.sqrt(sum/N)

if (len(sys.argv) != 3):
	print "input [file] [colume]";

cid = int(sys.argv[2]);

value = [];
for l1 in open(sys.argv[1],"r"):
	l1 = l1.rstrip();
	l2 = l1.split();
	if (len(l2) < cid):
		print "no %d th colum!!" % cid;
		continue;
        if float(l2[cid-1]) > -10000000000000. and float(l2[cid-1]) < 10000000000000:
		value.append(float(l2[cid-1]));

m1 = get_mean(value) #numpy.mean(value);
std1 = get_std(value,m1) #numpy.std(value);
value1 = [];
for vv in value:
	if (vv < 3*std1 + m1 and vv > m1 - 3*std1):
		value1.append(vv);

mean=get_mean(value1)
std = get_std(value1,mean)
#print value
if (len(value) != 0):
	#print numpy.mean(value1), numpy.std(value1), 0.5*(max(value1) + min(value1)), 0.5* (max(value1) - min(value1)),sum(value1)/len(value1);
	print mean, std, 0.5*(max(value1) + min(value1)), 0.5* (max(value1) - min(value1)),sum(value1)/len(value1);
else:
	print "NaN Nan Nan Nan";
