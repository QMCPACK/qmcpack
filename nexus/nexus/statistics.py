
import numpy as np


def theil_sen(x,y):
    x = np.asarray(x)
    y = np.asarray(y)
    assert len(x)==x.size
    assert len(y)==y.size
    assert len(x)==len(y)
    x = x.ravel()
    y = y.ravel()
    ms = []
    for i,(xi,yi) in enumerate(zip(x,y)):
        for xj,yj in zip(x[i+1:],y[i+1:]):
            ms.append((yi-yj)/(xi-xj))
    m = np.median(ms)
    bs = []
    for xi,yi in zip(x,y):
        bs.append(yi-m*xi)
    b = np.median(bs)
    return m,b
#end def theil_sen


def acf_autocorr_time(data_samples):
    if data_samples.ndim>1:
        assert np.max(data_samples.shape)==data_samples.size
        data_samples = data_samples.ravel()
    x = data_samples
    assert x.ndim==1,'data array must be 1-dimensional'
    N     = len(x)
    mean  = x.mean()
    var   = x.var()
    i     = 0          
    tempC = 0.5
    kappa = 0.0
    mtmp  = mean
    if np.abs(var)<1e-15:
        kappa = 1.0
    else:
        ovar=1.0/var
        while (tempC>0 and i<(N-1)):
            kappa=kappa+2.0*tempC
            i=i+1
            tempC = ovar/(N-i)*np.sum((x[0:N-i]-mtmp)*(x[i:N]-mtmp))
        kappa = max(kappa,1.0)
    t_auto = kappa
    return t_auto
#end def acf_autocorr_time


def reblocked_autocorr_time(data_samples,min_blocks=10,plot=False,show=False):
    if data_samples.ndim>1:
        assert np.max(data_samples.shape)==data_samples.size
        data_samples = data_samples.ravel()
    assert data_samples.ndim==1,'data array must be 1-dimensional'
    assert len(data_samples)>0,'data array must not be empty'
    assert min_blocks>=1,'minimum number of blocks must be at least one'
    if len(data_samples)==1:
        t_auto = 1.
        return t_auto
    nblocks = len(data_samples)
    nreblock_max = int(np.floor(nblocks/min_blocks))
    block_lens  = []
    data_means = []
    data_errs  = []
    # length 1 "reblocking"
    data_errs1 = data_samples.std(axis=0).ravel()/np.sqrt(nblocks)
    if np.all(data_errs1==0):
        return 1.
    block_lens.append(1)
    data_errs.append(data_errs1)
    # block lengths > 1
    for nrb in range(2,nreblock_max+1):
        ds = []
        for n in range(0,nblocks,nrb):
            if n+nrb>nblocks:
                break
            ds.append(data_samples[n:n+nrb].mean(axis=0))
        ds = np.array(ds)
        block_lens.append(nrb)
        data_errs.append(ds.std(axis=0).ravel()/np.sqrt(len(ds)))
    block_lens = np.array(block_lens)
    data_errs = np.array(data_errs)
    for bl in range(len(block_lens)):
        data_errs[bl] /= data_errs1
    dem = data_errs.mean(axis=1)
    des = data_errs.std(axis=1)
    assert len(dem)==len(block_lens)
    if len(block_lens)>1:
        p = theil_sen(block_lens,dem)
        m,b = p
        if m>0:
            err_max = np.polyval(p,[block_lens[-1]])[0]
        else:
            err_max = np.median(dem)
    else:
        err_max = dem[0]
    t_auto = float((err_max)**2)
    if plot:
        import matplotlib.pyplot as plt
        plt.figure(tight_layout=True)
        plt.errorbar(block_lens,dem,des,fmt='b.-')
        if len(block_lens)>1:
            plt.plot(block_lens,np.polyval(p,block_lens),'r--')
        plt.axhline(err_max,color='k')
        plt.xlabel('reblocking factor')
        plt.ylabel('errorbar')
        plt.title('t_auto = {}'.format(t_auto))
        if show:
            plt.show()
    return t_auto
#end def reblocked_autocorr_time


def series_stats(x,t_auto=None):
    if t_auto is None:
        #t_auto = autocorr_time(x)
        t_auto = 1.0
    N        = len(x)
    N_eff    = N/t_auto
    x_mean   = np.mean(x)
    x_stderr = np.std(x)/np.sqrt(N_eff)
    return x_mean,x_stderr,t_auto
#end def series_stats
