
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


def acf_autocorr_time(x):
    if x.ndim>1:
        assert np.max(x.shape)==x.size
        x = x.ravel()
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


def reblocked_autocorr_time(x,min_blocks=10,plot=False,show=False):
    if x.ndim>1:
        assert np.max(x.shape)==x.size
        x = x.ravel()
    assert x.ndim==1,'data array must be 1-dimensional'
    assert len(x)>0,'data array must not be empty'
    assert min_blocks>=1,'minimum number of blocks must be at least one'
    if len(x)==1:
        t_auto = 1.
        return t_auto
    nblocks = len(x)
    nreblock_max = int(np.floor(nblocks/min_blocks))
    block_lens  = []
    data_means = []
    data_errs  = []
    # length 1 "reblocking"
    data_errs1 = x.std(axis=0).ravel()/np.sqrt(nblocks)
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
            ds.append(x[n:n+nrb].mean(axis=0))
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


def reblocked_autocorr_time2(x,confidence=0.99,plot=False,show=False):
    """Estimate the integrated autocorrelation time by data blocking.

    This implementation follows the blocking transformation and automatic
    stopping test of Flyvbjerg and Petersen.  The returned value is the
    variance-inflation factor ``N*Var(mean)/Var(data)``, consistent with an
    effective sample size of ``N/t_auto``.
    """
    from scipy.stats import chi2

    x = np.asarray(x,dtype=float)
    if x.ndim>1 and np.max(x.shape)==x.size:
        x = x.ravel()
    if x.ndim!=1:
        raise ValueError('data array must be 1-dimensional')
    if len(x)==0:
        raise ValueError('data array must not be empty')
    if not np.all(np.isfinite(x)):
        raise ValueError('data array must contain only finite values')
    if not 0.<confidence<1.:
        raise ValueError('confidence must be between zero and one')
    if len(x)==1 or x.var()==0.:
        return 1.

    # The pairwise blocking transformation requires a power-of-two length.
    nlevel = int(np.floor(np.log2(len(x))))
    nused  = 2**nlevel
    xb = x[:nused].copy()

    block_sizes = 2**np.arange(nlevel,dtype=int)
    block_counts = nused//block_sizes
    variances = np.empty(nlevel,dtype=float)
    lag1_covariances = np.empty(nlevel,dtype=float)

    for level in range(nlevel):
        nblock = len(xb)
        delta = xb-xb.mean()
        variances[level] = np.mean(delta**2)
        lag1_covariances[level] = np.sum(delta[:-1]*delta[1:])/nblock
        xb = .5*(xb[0::2]+xb[1::2])

    # Under the null hypothesis that a blocking level is uncorrelated, this
    # statistic is asymptotically chi-square distributed.  The first level
    # that passes the test supplies the variance of the sample mean.
    correlations = np.divide(
        lag1_covariances,
        variances,
        out=np.zeros_like(lag1_covariances),
        where=variances>0.)
    weights = 2.**np.arange(nlevel,0,-1)
    test_statistics = np.cumsum((correlations**2*weights)[::-1])[::-1]
    critical_values = chi2.ppf(confidence,np.arange(nlevel,0,-1))
    passing_levels = np.flatnonzero(test_statistics<critical_values)
    if len(passing_levels)>0:
        selected_level = int(passing_levels[0])
    else:
        import warnings
        selected_level = nlevel-1
        warnings.warn(
            'blocking analysis did not reach an uncorrelated level; '
            'the autocorrelation-time estimate may be unreliable',
            RuntimeWarning,
            stacklevel=2)

    tau_levels = block_sizes*variances/variances[0]
    t_auto = max(float(tau_levels[selected_level]),np.finfo(float).eps)

    if plot:
        import matplotlib.pyplot as plt
        tau_errors = tau_levels*np.sqrt(2./(block_counts-1))
        plt.figure(tight_layout=True)
        plt.errorbar(block_sizes,tau_levels,yerr=tau_errors,fmt='b.-')
        plt.axvline(block_sizes[selected_level],color='k',linestyle='--')
        plt.axhline(t_auto,color='r')
        plt.xscale('log',base=2)
        plt.xlabel('block length')
        plt.ylabel('integrated autocorrelation time')
        plt.title('t_auto = {}'.format(t_auto))
        if show:
            plt.show()
    return t_auto
#end def reblocked_autocorr_time2


def integrated_autocorr_time(x, c=5.0):
    """Estimate integrated autocorrelation time with Geyer's IMS method.

    Autocorrelations are computed with an FFT.  Geyer's initial positive
    sequence of adjacent autocorrelation pairs is then made non-increasing
    with a linear-time pool-adjacent-violators algorithm.  This estimator is
    intended primarily for stationary, reversible Markov chains.

    Parameters
    ----------
    x : array_like
        One-dimensional sample sequence.

    c : float, optional
        Minimum number of estimated autocorrelation times that should fit
        in the input series.  A warning is emitted when ``len(x)<c*tau``.

    Returns
    -------
    tau : float
        Estimated integrated autocorrelation time.
    """

    try:
        c = float(c)
    except (TypeError,ValueError):
        raise ValueError('c must be a positive finite number') from None
    if not np.isfinite(c) or c<=0.:
        raise ValueError('c must be a positive finite number')

    x = np.asarray(x)
    if np.iscomplexobj(x):
        raise ValueError('input must be real-valued')
    if x.ndim>1 and np.max(x.shape)==x.size:
        x = x.ravel()
    if x.ndim!=1:
        raise ValueError('input must be one-dimensional')
    if len(x)==0:
        raise ValueError('input must not be empty')

    x = np.asarray(x,dtype=float)
    if not np.all(np.isfinite(x)):
        raise ValueError('input must contain only finite values')
    if len(x)<2:
        return 1.0

    x = x-x.mean()
    variance = np.mean(x**2)
    if variance==0.:
        return 1.0

    # Use a common denominator at all lags.  Unlike unbiased lag-by-lag
    # normalization, this produces a positive-semidefinite autocovariance
    # sequence and is more stable in the noisy tail.
    n = len(x)
    nfft = 1 << (2*n-1).bit_length()
    transform = np.fft.rfft(x,nfft)
    acf = np.fft.irfft(transform*np.conjugate(transform),nfft)[:n]
    acf /= acf[0]

    # Geyer pairs are Gamma_k = rho_(2k) + rho_(2k+1), beginning with
    # rho_0 + rho_1.  Discard an unmatched final autocorrelation.
    npair = len(acf)//2
    gamma = acf[:2*npair:2]+acf[1:2*npair:2]

    # Initial positive sequence: later pairs are excluded as soon as the
    # first nonpositive pair is encountered.
    nonpositive = np.flatnonzero(gamma<=0.)
    if len(nonpositive)>0:
        gamma = gamma[:nonpositive[0]]
    if len(gamma)==0:
        return 0.0

    # Initial monotone sequence via a stack-based PAV implementation.
    # Each pooled block is represented by its mean and number of pairs.
    values = []
    weights = []
    for value in gamma:
        values.append(float(value))
        weights.append(1)
        while len(values)>1 and values[-2]<values[-1]:
            weight       = weights[-2]+weights[-1]
            value        = values[-2]*weights[-2]+values[-1]*weights[-1]/weight
            values[-2:]  = [value]
            weights[-2:] = [weight]

    gamma_mono = np.repeat(values,weights)
    tau = max(0.,float(-1.+2.*gamma_mono.sum()))

    if n<c*tau:
        import warnings
        warnings.warn(
            'the time series is shorter than c autocorrelation times; '
            'the estimate may be unreliable',
            RuntimeWarning,
            stacklevel=2)

    return tau
#end def integrated_autocorr_time



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
