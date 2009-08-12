#!/usr/bin/env python
from pylab import *
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

s = loadtxt(sys.argv[1])
c = s.shape[0]
scale = 1.0
if (len(sys.argv) > 3):
    scale = 1.0 / float(sys.argv[3])
if len(sys.argv) > 2:
    first = int(sys.argv[2])
else:
    first = 20
    
E      = s[first:c,1]*scale;
Var    = s[first:c,2]*scale;
Weight = s[first:c,3];
Pop    = s[first:c,4];
Etrial = s[first:c,5]*scale;
step = range(first,c,1)

N = len(E)
subplots_adjust(hspace=0.0000001, wspace=0.0000001)
ax1 = subplot(211)
title(sys.argv[1])
plot (step, E,'b');
myformatter = matplotlib.ticker.ScalarFormatter(False)
gca().axes.yaxis.set_major_formatter(myformatter)
hold(True);
plot(step, Etrial, 'r');
ylabel ('Energy (Hartrees)');
ax2 = subplot(2, 1, 2, sharex=ax1)
plot(step, Pop, 'r');
ylabel ('Population (walkers)')
yticklabels = ax2.get_yticklabels();
ytx = yticks()[0].tolist()
ytx.pop()
ytx2 = array(ytx)
yticks(ytx2)
setp(ax1.get_xticklabels(), visible=False)
xlabel('DMC step')
show();



