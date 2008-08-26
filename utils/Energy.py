#!/usr/bin/env python
from numpy import *
import sys

def corr(i,x,mean,var):
    N=len(x)
    if var==0:#if the variance is 0 return an effectively infinity corr
        return 1e100
    corr=1.0/var*1.0/(N-i)*sum((x[0:N-i]-mean)*(x[i:N]-mean))
    return corr

def Stats(x):
    N=len(x)
    mean=sum(x)/(N+0.0)
    xSquared=x*x
    var=sum(xSquared)/(N+0.0)-mean*mean
    i=0          
    tempC=0.5
    kappa=0.0
    while (tempC>0 and i<(N-1)):
        kappa=kappa+2.0*tempC
        i=i+1
        tempC=corr(i,x,mean,var)
    if kappa == 0.0:
        kappa = 1.0
    Neff=(N+0.0)/(kappa+0.0)
    error=sqrt(var/Neff)
    return (mean,var,error,kappa)

def Block (x, blockfactor):
    N = len(x)
    Nblocks = N / blockfactor
    xb = zeros(Nblocks)
    for i in range(0,Nblocks):
        start = i*blockfactor
        end   = (i+1)*blockfactor
        xb[i] = sum(x[start:end])/(blockfactor+0.0)
    return xb

def MeanErrorString (mean, error):
     if (mean!=0.0):
          meanDigits = math.floor(math.log(abs(mean))/math.log(10))
     else:
          meanDigits=2
     if (error!=0.0):
          rightDigits = -math.floor(math.log(error)/math.log(10))+1
     else:
          rightDigits=2
     if (rightDigits < 0):
          rightDigits = 0
     formatstr = '%1.' + '%d' % rightDigits + 'f'
     meanstr  = formatstr % mean
     errorstr = formatstr % error
     return meanstr + ' +/- ' + errorstr
#     return (meanstr, errorstr)


file = open (sys.argv[1], 'r')
names = file.readline().split()
file.close()

s = loadtxt(sys.argv[1])
c = s.shape[0]
if (len(sys.argv) > 3):
    s = s / float(sys.argv[3])
if len(sys.argv) > 2:
    first = int(sys.argv[2])
else:
    first = 20

for i in range(2,len(names)):
    n = names[i];
    data = s[first:c,i-1];
    (avg, var, error, kapp) = Stats(data)
    print '%-20s = %s' % (n, MeanErrorString(avg,error))
    
E = s[first:c,1];
Vloc  = s[first:c,2];
KE    = s[first:c,3];
E_E   = s[first:c,4];
locPP = s[first:c,5];
NLPP  = s[first:c,6];
IonIon = s[first:c,7];



#avg = mean (E);
#var = mean ((E-avg)*(E-avg));
#err = sqrt(var/c);

N = len(E)
#print "Blockfactor      Std Error"
#errlist = []
#for factor in range(10,N/2,10):
#    Eblock = Block(E,factor)
#    E2 = Eblock*Eblock
#    avg = sum (Eblock)/(len(Eblock)+0.0)
#    var = var=sum((Eblock-avg)*(Eblock-avg))/(len(Eblock)+0.0)
#    error = sqrt (var/(len(Eblock)+0.0))
#    errlist.append(error)
#    print "   %4d      %8.6f" % (factor, error)

(avg, var, error, kapp) = Stats(E)
#error = max(errlist)
print "E     =  " + MeanErrorString(avg,error)
(avg, var, error, kapp) = Stats(Vloc)
print "Vloc  =  " + MeanErrorString(avg,error)
(avg, var, error, kapp) = Stats(KE)
print "KE    =  " + MeanErrorString(avg,error)
(avg, var, error, kapp) = Stats(E_E)
print "E_E   =  " + MeanErrorString(avg,error)
(avg, var, error, kapp) = Stats(locPP)
print "locPP = " + MeanErrorString(avg,error)
(avg, var, error, kapp) = Stats(NLPP)
print "NLPP  = " + MeanErrorString(avg,error)
(avg, var, error, kapp) = Stats(IonIon)
print "IonIon  = " + repr(avg)



