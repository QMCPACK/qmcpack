
import numpy as np


def theil_sen(x,y):
    """Return the Theil--Sen slope and intercept for paired observations.

    Parameters
    ----------
    x : array_like
        Independent-variable observations.

    y : array_like
        Dependent-variable observations paired with ``x``.

    Returns
    -------
    slope : scalar
        Median of all pairwise slopes.

    intercept : scalar
        Median residual intercept for ``slope``.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    assert len(x)==x.size
    assert len(y)==y.size
    assert len(x)==len(y)
    x = x.ravel()
    y = y.ravel()
    n = len(x)
    npairs = n*(n-1)//2
    if npairs>0:
        # Fill the upper triangle a row at a time.  This keeps the quadratic
        # work in NumPy without allocating a pair of quadratic index arrays.
        first_slopes = (y[0]-y[1:])/(x[0]-x[1:])
        slopes = np.empty(npairs,dtype=first_slopes.dtype)
        slopes[:n-1] = first_slopes
        start = n-1
        for i in range(1,n-1):
            count = n-i-1
            slopes[start:start+count] = (y[i]-y[i+1:])/(x[i]-x[i+1:])
            start += count
    else:
        slopes = np.empty(0)
    m = np.median(slopes,overwrite_input=True)
    b = np.median(y-m*x)
    return m,b
#end def theil_sen


def theil_sen_stoch(x,y):
    """Estimate a Theil-Sen fit with stochastic pair sampling for large data.

    The number of sampled pairwise slopes is ``ceil(16000*sqrt(n))``.  This
    schedule was empirically calibrated on standardized linear-regression
    problems with Gaussian, heavy-tailed, heteroscedastic, and contaminated
    data; sampled slope angles stayed within one degree of their exact values.
    Pairwise slopes are sampled uniformly with replacement.  A fixed local
    random seed makes the estimate reproducible without altering NumPy's
    global random state.

    Parameters
    ----------
    x : array_like
        Independent-variable observations.

    y : array_like
        Dependent-variable observations paired with ``x``.

    Returns
    -------
    slope : scalar
        Median of the exact or sampled pairwise slopes.

    intercept : scalar
        Median residual intercept for ``slope``.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    assert len(x)==x.size
    assert len(y)==y.size
    assert len(x)==len(y)
    x = x.ravel()
    y = y.ravel()

    n = len(x)
    npairs = n*(n-1)//2
    sample_scale = 16000.
    nsampled = int(np.ceil(sample_scale*np.sqrt(n)))
    if npairs<=nsampled:
        return theil_sen(x,y)

    random_seed = 314159
    rng = np.random.Generator(np.random.PCG64(random_seed))
    i = rng.integers(0,n,size=nsampled)
    j = rng.integers(0,n-1,size=nsampled)
    j += j>=i
    slopes = (y[i]-y[j])/(x[i]-x[j])
    m = np.median(slopes,overwrite_input=True)
    b = np.median(y-m*x)
    return m,b
#end def theil_sen_stoch


def theil_sen_stoch_reblock(x,y):
    """Estimate a Theil-Sen fit using a reblocking-specific sample schedule.

    The number of sampled pairwise slopes is ``ceil(24*sqrt(n))``.  This
    schedule was empirically calibrated on reblocked IID, heavy-tailed,
    correlated, oscillatory, mixed-timescale, and nonstationary series; the
    sampled slope angle stayed within one degree of the exact slope angle.  A
    fixed local random seed makes the estimate reproducible without altering
    NumPy's global random state.

    Parameters
    ----------
    x : array_like
        Independent-variable observations.

    y : array_like
        Dependent-variable observations paired with ``x``.

    Returns
    -------
    slope : scalar
        Median of the exact or sampled pairwise slopes.

    intercept : scalar
        Median residual intercept for ``slope``.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    assert len(x)==x.size
    assert len(y)==y.size
    assert len(x)==len(y)
    x = x.ravel()
    y = y.ravel()

    n = len(x)
    npairs = n*(n-1)//2
    sample_scale = 24.
    nsampled = int(np.ceil(sample_scale*np.sqrt(n)))
    if npairs<=nsampled:
        return theil_sen(x,y)

    random_seed = 314159
    rng = np.random.Generator(np.random.PCG64(random_seed))
    i = rng.integers(0,n,size=nsampled)
    j = rng.integers(0,n-1,size=nsampled)
    j += j>=i
    slopes = (y[i]-y[j])/(x[i]-x[j])
    m = np.median(slopes,overwrite_input=True)
    b = np.median(y-m*x)
    return m,b
#end def theil_sen_stoch_reblock


def reblocked_autocorr_time(x,min_blocks=10,plot=False,show=False):
    """Estimate autocorrelation time from the growth of blocked errors.

    For every integer block length that leaves at least ``min_blocks`` blocks,
    contiguous block means and their standard error are computed.  Their
    ratio to the unblocked standard error is fitted as a function of block
    length with a Theil--Sen line (using reproducible stochastic pair sampling
    for large inputs).  The squared fitted ratio at the largest usable block
    length is used to obtain the auto-correlation time.

    Strengths: Independent block-based cross-check that can expose slow-mode
    uncertainty without requiring reversibility.

    Weaknesses: Its estimates have a broad 10--90% spread, and the calculated
    autocorrelation times can exhibit larger fluctutations toward overestimation.

    Parameters
    ----------
    x : numpy.ndarray
        Nonempty one-dimensional sample sequence.  Vector-shaped arrays are
        flattened.

    min_blocks : int, optional
        Minimum number of complete blocks retained at the largest block
        length.  Must be at least one.

    plot : bool, optional
        Plot normalized blocked errors, the robust fitted line, and the
        selected error estimate.

    show : bool, optional
        Display the plot immediately.  This has an effect only when ``plot``
        is true.

    Returns
    -------
    tau : float
        Estimated integrated autocorrelation time.  Constant and singleton
        sequences return one.
    """
    if x.ndim>1:
        assert np.max(x.shape)==x.size
        x = x.ravel()
    assert x.ndim==1,'data array must be 1-dimensional'
    assert len(x)>0,'data array must not be empty'
    assert min_blocks>=1,'minimum number of blocks must be at least one'
    t_auto = 1.
    if len(x)==1:
        return t_auto
    nblocks = len(x)
    nreblock_max = int(np.floor(nblocks/min_blocks))
    # length 1 "reblocking"
    data_errs1 = x.std()/np.sqrt(nblocks)
    if data_errs1==0:
        return t_auto

    block_lens = np.arange(1,max(1,nreblock_max)+1)
    data_errs = np.empty(len(block_lens),dtype=np.asarray(data_errs1).dtype)
    data_errs[0] = data_errs1

    if nreblock_max>=2:
        # A prefix sum gives every contiguous block sum with two indexed
        # reads.  Center first to limit cancellation in the differences.
        work_dtype = np.result_type(x.dtype,np.float64)
        centered = x.astype(work_dtype,copy=False)
        centered = centered-centered.mean()
        cumulative = np.empty(nblocks+1,dtype=work_dtype)
        cumulative[0] = 0.
        np.cumsum(centered,out=cumulative[1:])

        # Lay out the blocks for all reblocking lengths in one flat array so
        # their means and variances can be evaluated by grouped reductions.
        block_counts = nblocks//block_lens[1:]
        group_offsets = np.empty(len(block_counts),dtype=int)
        group_offsets[0] = 0
        np.cumsum(block_counts[:-1],out=group_offsets[1:])
        repeated_lens = np.repeat(block_lens[1:],block_counts)
        block_starts = (
            np.arange(block_counts.sum())
            -np.repeat(group_offsets,block_counts))
        block_starts *= repeated_lens
        block_means = (
            cumulative[block_starts+repeated_lens]-cumulative[block_starts]
            )/repeated_lens

        group_means = np.add.reduceat(block_means,group_offsets)/block_counts
        deviations = block_means-np.repeat(group_means,block_counts)
        group_variances = np.add.reduceat(
            np.abs(deviations)**2,group_offsets)/block_counts
        data_errs[1:] = np.sqrt(group_variances/block_counts)

    dem = data_errs/data_errs1
    des = np.zeros_like(dem)
    assert len(dem)==len(block_lens)
    if len(block_lens)>1:
        p = theil_sen_stoch_reblock(block_lens,dem)
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
        plt.title(f't_auto = {t_auto}')
        if show:
            plt.show()
    return t_auto
#end def reblocked_autocorr_time



def acf_autocorr_time(x,reliability=False):
    """Estimate autocorrelation time from a windowed sample ACF.

    The autocorrelation function is evaluated in ``O(N log N)`` time with an
    FFT and a common denominator at all lags.  A Bartlett noise estimate is
    used to locate the first sustained noise-dominated region.  Resolved lags
    are retained by a flat-top window, followed by a linear taper through the
    noisy boundary.  The returned value is the variance-inflation factor
    ``N*Var(mean)/Var(x)``.

    Strengths: Best long-chain accuracy, low variance, fast, and supports
    negative or oscillatory correlation.

    Weaknesses: Can truncate before detecting weak slow modes and generally
    underestimates the autocorrelation time for short time series.

    Parameters
    ----------
    x : array_like
        One-dimensional, finite, real-valued sample sequence.  Vector-shaped
        two-dimensional arrays are flattened.

    reliability : bool, optional
        If true, return ``(tau, not_reliable)`` instead of only ``tau``.

    Returns
    -------
    tau : float or (float, bool)
        Estimated integrated autocorrelation time.  A value of one denotes
        IID-like sampling; negative correlation can produce a value below
        one.  The optional Boolean is true when the ACF fails its reliability
        assessment.
    """
    x = np.asarray(x)
    if np.iscomplexobj(x):
        raise ValueError('data array must be real-valued')
    if x.ndim>1 and np.max(x.shape)==x.size:
        x = x.ravel()
    if x.ndim!=1:
        raise ValueError('data array must be 1-dimensional')
    if len(x)==0:
        raise ValueError('data array must not be empty')
    x = np.asarray(x,dtype=float)
    if not np.all(np.isfinite(x)):
        raise ValueError('data array must contain only finite values')
    not_reliable = False
    if len(x)==1:
        return (1.,not_reliable) if reliability else 1.

    x = x-x.mean()
    variance = np.mean(x**2)
    if variance==0.:
        return (1.,not_reliable) if reliability else 1.

    n = len(x)
    nfft = 1 << (2*n-1).bit_length()
    transform = np.fft.rfft(x,nfft)
    acf = np.fft.irfft(transform*np.conjugate(transform),nfft)[:n]
    acf /= acf[0]

    # Locate the noisy tail without treating a physical sign change as the
    # end of the correlation structure.  Bartlett's large-sample expression
    # supplies a lag-dependent noise scale; requiring several quiet lags in a
    # row avoids stopping at an isolated zero crossing.
    noise_scale = 1.5
    quiet_needed = 5
    max_lag = min(n-1,max(quiet_needed,n//2))
    rho2_sum = 0.
    quiet_count = 0
    quiet_start = None
    for lag in range(1,max_lag+1):
        noise = np.sqrt((1.+2.*rho2_sum)/n)
        if np.abs(acf[lag])<=noise_scale*noise:
            quiet_count += 1
        else:
            quiet_count = 0
        rho2_sum += acf[lag]**2
        if quiet_count>=quiet_needed:
            quiet_start = lag-quiet_needed+1
            break

    if quiet_start is None:
        # The ACF did not reach a noise-dominated region; the
        # autocorrelation-time estimate may be unreliable.
        not_reliable = True
        quiet_start = max_lag

    # If even the first nonzero lags are noise, the IID estimate is exact and
    # avoids adding pure-noise terms.  Otherwise, a flat-top lag window retains
    # the resolved ACF and smoothly damps the noisy region beyond it.
    if quiet_start==1:
        t_auto = 1.
    else:
        bandwidth = min(max_lag,max(1,2*quiet_start))
        lags = np.arange(1,bandwidth+1)
        fraction = lags/bandwidth
        window = np.where(fraction<=.5,1.,2.*(1.-fraction))
        t_auto = 1.+2.*np.sum(window*acf[1:bandwidth+1])
        t_auto = max(float(t_auto),np.finfo(float).eps)
    return (t_auto,not_reliable) if reliability else t_auto
#end def acf_autocorr_time



def geyer_ims_autocorr_time(x,c=5.0,reliability=False):
    """Estimate integrated autocorrelation time with Geyer's IMS method.

    Autocorrelations are computed with an FFT.  Geyer's initial positive
    sequence of adjacent autocorrelation pairs is then made non-increasing
    with a linear-time pool-adjacent-violators algorithm.  This estimator is
    intended primarily for stationary, reversible Markov chains.

    Strengths: Fast, stable noisy-tail treatment, and strong theoretical basis
    for reversible MCMC.

    Shortcomings: Relatively variable for negative correlation.

    Parameters
    ----------
    x : array_like
        One-dimensional sample sequence.

    c : float, optional
        Minimum number of estimated autocorrelation times that should fit
        in the input series.

    reliability : bool, optional
        If true, return ``(tau, not_reliable)`` instead of only ``tau``.

    Returns
    -------
    tau : float or (float, bool)
        Estimated integrated autocorrelation time.  The optional Boolean is
        true when fewer than ``c`` estimated autocorrelation times fit in the
        series.
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
    not_reliable = False
    if len(x)<2:
        return (1.,not_reliable) if reliability else 1.

    x = x-x.mean()
    variance = np.mean(x**2)
    if variance==0.:
        return (1.,not_reliable) if reliability else 1.

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
        return (0.,not_reliable) if reliability else 0.

    # Initial monotone sequence via a stack-based PAV implementation.
    # Each pooled block is represented by its mean and number of pairs.
    values = []
    weights = []
    for value in gamma:
        values.append(float(value))
        weights.append(1)
        while len(values)>1 and values[-2]<values[-1]:
            weight       = weights[-2]+weights[-1]
            value        = (
                values[-2]*weights[-2]+values[-1]*weights[-1]
                )/weight
            values[-2:]  = [value]
            weights[-2:] = [weight]

    gamma_mono = np.repeat(values,weights)
    tau = max(0.,float(-1.+2.*gamma_mono.sum()))

    if n<c*tau:
        # The time series is shorter than c autocorrelation times; the
        # estimate may be unreliable.
        not_reliable = True

    return (tau,not_reliable) if reliability else tau
#end def geyer_ims_autocorr_time



def autocorr_time(x,reliability=False):
    """Conservatively combine three autocorrelation-time estimates.

    The ACF and Geyer initial-monotone-sequence probe the correlation
    structure in different ways.  Since an underestimated autocorrelation
    time leads directly to an underestimated uncertainty, this function
    returns the largest of their estimates.

    Both ACF and Geyer show rapid convergence toward accurate estimates
    with time series lengths and low variability. Taking the max between
    them provides demonstrable stability against underestimation while
    retaining low bias for longer series.

    Parameters
    ----------
    x : array_like
        One-dimensional sample sequence.

    reliability : bool, optional
        If true, return ``(tau, not_reliable)`` instead of only ``tau``.  The
        flag combines the ACF and Geyer IMS reliability assessments.

    Returns
    -------
    tau : float or (float, bool)
        Maximum of ACF and Geyer autocorrelation times, optionally
        accompanied by the combined unreliability flag.
    """
    x = np.asarray(x)

    t_auto_acf,nr_acf     = acf_autocorr_time(x,reliability=True)
    t_auto_geyer,nr_geyer = geyer_ims_autocorr_time(x,reliability=True)

    t_auto = max(t_auto_acf,t_auto_geyer)
    not_reliable = nr_acf or nr_geyer

    return (t_auto,not_reliable) if reliability else t_auto
#end def autocorr_time



def series_stats(x,t_auto=None):
    if t_auto is None:
        t_auto = autocorr_time(x)
    N        = len(x)
    N_eff    = N/t_auto
    x_mean   = np.mean(x)
    x_stderr = np.std(x)/np.sqrt(N_eff)
    return x_mean,x_stderr,t_auto
#end def series_stats
