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
     return (meanstr, errorstr)
#     return (meanstr, errorstr)



file = open (sys.argv[1], 'r')
names = file.readline().split()
file.close()

M = int(sys.argv[3])
print '# Forward-walking distance is ' + repr(M)

numat = 0
for n in names:
    if (n[0:4] == 'FW_F'):
        d = n.split('_')
        m = int(d[-1])
        if (m == M):
            iat = int(d[-3])
            dim = int(d[-2])
            numat = max(iat+1,numat)
collist = [ [ [], [], [] ] ]
for iat in range(0,numat-1):
    collist.append([ [],[],[] ])

col = 0
for n in names:
    if (n[0:4] == 'FW_F'):
        d = n.split('_')
        m = int(d[-1])
        if (m == M):
            iat = int(d[-3])
            dim = int(d[-2])
            collist[iat][dim].append(col-1)
    
    if(n[0:4] == 'F_AA'):
        d = n.split('_')
        iat = int(d[-2])
        dim = int(d[-1])
        if (iat < numat):
            collist[iat][dim].append(col-1)        
    col = col + 1

s = loadtxt(sys.argv[1])
c = s.shape[0]
if len(sys.argv) > 2:
    first = int(sys.argv[2])
else:
    first = 20

print '# Atom               Force              +/-             Error'
for iat in range(0,numat):
    means  = []
    errors = []
    for dim in range(0,3):
        Ftot = 0.0*array(c-first)
        for col in collist[iat][dim]:
#            print 'iat = %d  dim = %d  col=%d' % (iat, dim, col)
            Ftot = Ftot + s[first:c+1,col]
        (avg, var, error, kapp) = Stats(Ftot)
        means.append (avg)
        errors.append(error)
#        print '%-20s = %s' % (n, MeanErrorString(avg,error))
    print '%4d  %10.5f %10.5f %10.5f    %10.5f %10.5f %10.5f' \
        % (iat+1,means[0], means[1], means[2], errors[0], errors[1], errors[2])
