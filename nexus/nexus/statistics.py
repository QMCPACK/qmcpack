"""Robust statistics and autocorrelation analysis for simulation data."""

import numpy as np

############################################################################
#                                                                          #
#                Autocorrelation estimator stress testing                  #
#                ----------------------------------------                  #
#                                                                          #
# The following recommendations summarize empirical stress tests of the    #
# autocorrelation estimators for population-averaged QMC data.             #
#                                                                          #
# The estimators were judged with repeated equilibrium Markov-chain tests  #
# designed around population-averaged VMC and DMC data. Each test used the #
# mean of 32 walkers and an exactly known integrated autocorrelation time. #
# Equilibrium marginals included Gaussian, symmetric Laplace, centered     #
# exponential, Student-t(6), and Student-t(3) distributions; the last has  #
# the pessimistic 1/|x|**4 density tail considered in QMC.                 #
#                                                                          #
# Reversible retain-or-refresh chains covered IID data (tau=1), ordinary   #
# persistence (tau=9), and slow mixing (tau=39). Reversible sign-flip      #
# chains tested negative correlation (tau=1/9). Two-scale walker groups    #
# (tau=27) tested weak slow modes, while a slow population-wide common     #
# mode (tau=31.5) represented DMC-like branching or population-control     #
# correlations in data already averaged over walkers.                      #
#                                                                          #
# The main benchmark used 100 independent repetitions at lengths 64        #
# through 4096; a separate long-chain benchmark used 4096, 8192, and       #
# 16384. Comparisons included relative bias and RMSE, variability,         #
# empirical 10/50/90 percentiles relative to truth, failure rates,         #
# sustained convergence, weak-mode underestimation, and execution time.    #
# Every estimator received the same generated series within each trial.    #
#                                                                          #
# The prewhitened spectral method was removed because it provided no clear #
# long-chain accuracy advantage over ACF or Geyer, retained a broad upper  #
# tail, and developed a large runtime increase at the longest length. The  #
# Flyvbjerg--Petersen reblocking method was removed because it remained    #
# systematically low, particularly for weak slow modes, and converged in   #
# fewer cases than ACF or Geyer despite its relatively narrow spread.      #
#                                                                          #
# Deterministic test coverage is in tests/test_statistics.py.              #
#                                                                          #
############################################################################


def _paired_real_arrays(x,y):
    """Validate and flatten paired real-valued sample arrays."""
    x = np.asarray(x)
    y = np.asarray(y)
    if np.iscomplexobj(x) or np.iscomplexobj(y):
        msg = 'input arrays must be real-valued'
        raise ValueError(msg)
    if x.ndim>1 and np.max(x.shape)==x.size:
        x = x.ravel()
    if y.ndim>1 and np.max(y.shape)==y.size:
        y = y.ravel()
    if x.ndim!=1 or y.ndim!=1:
        msg = 'input arrays must be one-dimensional'
        raise ValueError(msg)
    if len(x)!=len(y):
        msg = 'input arrays must have equal lengths'
        raise ValueError(msg)
    try:
        x = np.asarray(x,dtype=float)
        y = np.asarray(y,dtype=float)
    except (TypeError,ValueError):
        msg = 'input arrays must be numeric'
        raise ValueError(msg) from None
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        msg = 'input arrays must contain only finite values'
        raise ValueError(msg)
    return x,y
#end def _paired_real_arrays


def theil_sen(x,y):
    """Return the Theil--Sen slope and intercept for paired observations.

    Parameters
    ----------
    x : array_like
        Finite, real-valued independent-variable observations containing at
        least two distinct values. Pairs with equal values are excluded from
        the slope sample.

    y : array_like
        Finite, real-valued dependent-variable observations paired with ``x``.

    Returns
    -------
    slope : scalar
        Median of all pairwise slopes.

    intercept : scalar
        Median residual intercept for ``slope``.
    """
    x,y = _paired_real_arrays(x,y)
    n = len(x)
    if n<2 or not np.any(x!=x[0]):
        msg = 'at least two distinct x values are required'
        raise ValueError(msg)

    npairs = n*(n-1)//2
    slope_dtype = np.result_type(x.dtype,y.dtype,float)
    slopes = np.empty(npairs,dtype=slope_dtype)
    start = 0
    # Fill the upper triangle a row at a time.  Pairs with equal independent
    # coordinates do not define a slope and are excluded.
    for i in range(n-1):
        dx = x[i]-x[i+1:]
        valid = dx!=0
        count = np.count_nonzero(valid)
        slopes[start:start+count] = (y[i]-y[i+1:])[valid]/dx[valid]
        start += count
    slopes = slopes[:start]
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
        Finite, real-valued independent-variable observations containing at
        least two distinct values. Pairs with equal values are excluded from
        the slope sample.

    y : array_like
        Finite, real-valued dependent-variable observations paired with ``x``.

    Returns
    -------
    slope : scalar
        Median of the exact or sampled pairwise slopes.

    intercept : scalar
        Median residual intercept for ``slope``.
    """
    x,y = _paired_real_arrays(x,y)

    n = len(x)
    npairs = n*(n-1)//2
    sample_scale = 16000.
    nsampled = int(np.ceil(sample_scale*np.sqrt(n)))
    if npairs<=nsampled:
        return theil_sen(x,y)

    random_seed = 314159
    rng = np.random.Generator(np.random.PCG64(random_seed))
    i   = rng.integers(0,n,size=nsampled)
    j   = rng.integers(0,n-1,size=nsampled)
    j += j>=i
    valid = x[i]!=x[j]
    if not np.any(valid):
        return theil_sen(x,y)
    i = i[valid]
    j = j[valid]
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
        Finite, real-valued independent-variable observations containing at
        least two distinct values. Pairs with equal values are excluded from
        the slope sample.

    y : array_like
        Finite, real-valued dependent-variable observations paired with ``x``.

    Returns
    -------
    slope : scalar
        Median of the exact or sampled pairwise slopes.

    intercept : scalar
        Median residual intercept for ``slope``.
    """
    x,y = _paired_real_arrays(x,y)

    n = len(x)
    npairs = n*(n-1)//2
    sample_scale = 24.
    nsampled = int(np.ceil(sample_scale*np.sqrt(n)))
    if npairs<=nsampled:
        return theil_sen(x,y)

    random_seed = 314159
    rng = np.random.Generator(np.random.PCG64(random_seed))
    i   = rng.integers(0,n,size=nsampled)
    j   = rng.integers(0,n-1,size=nsampled)
    j += j>=i
    valid = x[i]!=x[j]
    if not np.any(valid):
        return theil_sen(x,y)
    i = i[valid]
    j = j[valid]
    slopes = (y[i]-y[j])/(x[i]-x[j])
    m = np.median(slopes,overwrite_input=True)
    b = np.median(y-m*x)
    return m,b
#end def theil_sen_stoch_reblock


def reblocked_autocorr_time(x,min_blocks=10,plot=False,show=False):
    """Estimate autocorrelation time from the growth of blocked errors.

    This estimator currently overestimates the autocorrelation times in a 
    number of cases. Prefer the Geyer method.

    For MCMC data, just use the ``autocorr_time'' function.

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
    x : array_like
        Nonempty, finite, real-valued, one-dimensional sample sequence.
        Vector-shaped arrays are flattened.

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
    x = np.asarray(x)
    if np.iscomplexobj(x):
        msg = 'data array must be real-valued'
        raise ValueError(msg)
    if x.ndim>1 and np.max(x.shape)==x.size:
        x = x.ravel()
    if x.ndim!=1:
        msg = 'data array must be 1-dimensional'
        raise ValueError(msg)
    if len(x)==0:
        msg = 'data array must not be empty'
        raise ValueError(msg)
    try:
        x = np.asarray(x,dtype=float)
    except (TypeError,ValueError):
        msg = 'data array must be numeric'
        raise ValueError(msg) from None
    if not np.all(np.isfinite(x)):
        msg = 'data array must contain only finite values'
        raise ValueError(msg)
    if isinstance(min_blocks,(bool,np.bool_)) or not isinstance(
        min_blocks,(int,np.integer)
        ) or min_blocks<1:
        msg = 'minimum number of blocks must be a positive integer'
        raise ValueError(msg)
    t_auto = 1.
    if len(x)==1:
        return t_auto
    nblocks      = len(x)
    nreblock_max = int(np.floor(nblocks/min_blocks))
    # length 1 "reblocking"
    data_errs1 = x.std()/np.sqrt(nblocks)
    if data_errs1==0:
        return t_auto

    block_lens   = np.arange(1,max(1,nreblock_max)+1)
    data_errs    = np.empty(len(block_lens),dtype=np.asarray(data_errs1).dtype)
    data_errs[0] = data_errs1

    if nreblock_max>=2:
        # A prefix sum gives every contiguous block sum with two indexed
        # reads.  Center first to limit cancellation in the differences.
        work_dtype = np.result_type(x.dtype,np.float64)
        centered   = x.astype(work_dtype,copy=False)
        centered   = centered-centered.mean()
        cumulative = np.empty(nblocks+1,dtype=work_dtype)
        cumulative[0] = 0.
        np.cumsum(centered,out=cumulative[1:])

        # Lay out the blocks for all reblocking lengths in one flat array so
        # their means and variances can be evaluated by grouped reductions.
        block_counts  = nblocks//block_lens[1:]
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

        group_means     = np.add.reduceat(block_means,group_offsets)/block_counts
        deviations      = block_means-np.repeat(group_means,block_counts)
        group_variances = np.add.reduceat(
            np.abs(deviations)**2,group_offsets)/block_counts
        data_errs[1:] = np.sqrt(group_variances/block_counts)

    dem = data_errs/data_errs1
    des = np.zeros_like(dem)
    assert len(dem)==len(block_lens)
    if len(block_lens)>1:
        p = theil_sen_stoch_reblock(block_lens,dem)
        m,_ = p
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

    Best for long chains.  Generally prefer the Geyer method.

    For MCMC data, just use the ``autocorr_time'' function.

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
        msg = 'data array must be real-valued'
        raise ValueError(msg)
    if x.ndim>1 and np.max(x.shape)==x.size:
        x = x.ravel()
    if x.ndim!=1:
        msg = 'data array must be 1-dimensional'
        raise ValueError(msg)
    if len(x)==0:
        msg = 'data array must not be empty'
        raise ValueError(msg)
    x = np.asarray(x,dtype=float)
    if not np.all(np.isfinite(x)):
        msg = 'data array must contain only finite values'
        raise ValueError(msg)
    not_reliable = False
    if len(x)==1:
        return (1.,not_reliable) if reliability else 1.

    x = x-x.mean()
    variance = np.mean(x**2)
    if variance==0.:
        return (1.,not_reliable) if reliability else 1.

    n         = len(x)
    nfft      = 1 << (2*n-1).bit_length()
    transform = np.fft.rfft(x,nfft)
    acf       = np.fft.irfft(transform*np.conjugate(transform),nfft)[:n]
    acf /= acf[0]

    # Locate the noisy tail without treating a physical sign change as the
    # end of the correlation structure.  Bartlett's large-sample expression
    # supplies a lag-dependent noise scale; requiring several quiet lags in a
    # row avoids stopping at an isolated zero crossing.
    noise_scale  = 1.5
    quiet_needed = 5
    max_lag      = min(n-1,max(quiet_needed,n//2))
    rho2_sum     = 0.
    quiet_count  = 0
    quiet_start  = None
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
        quiet_start  = max_lag

    # If even the first nonzero lags are noise, the IID estimate is exact and
    # avoids adding pure-noise terms.  Otherwise, a flat-top lag window retains
    # the resolved ACF and smoothly damps the noisy region beyond it.
    if quiet_start==1:
        t_auto = 1.
    else:
        bandwidth = min(max_lag,max(1,2*quiet_start))
        lags      = np.arange(1,bandwidth+1)
        fraction  = lags/bandwidth
        window    = np.where(fraction<=.5,1.,2.*(1.-fraction))
        t_auto    = 1.+2.*np.sum(window*acf[1:bandwidth+1])
        t_auto    = max(float(t_auto),np.finfo(float).eps)
    return (t_auto,not_reliable) if reliability else t_auto
#end def acf_autocorr_time



def geyer_ims_autocorr_time(x,c=5.0,reliability=False,acf_fallback=True):
    """Estimate integrated autocorrelation time with Geyer's IMS method.

    This is the single best autocorrelation estimator.
    
    For MCMC data, just use the ``autocorr_time'' function.

    Autocorrelations are computed with an FFT.  Geyer's initial positive
    sequence of adjacent autocorrelation pairs is then made non-increasing
    with a linear-time pool-adjacent-violators algorithm.  This estimator is
    intended primarily for stationary, reversible Markov chains.

    By default, estimates below one are replaced by the result from
    :func:`acf_autocorr_time`.  Set ``acf_fallback=False`` to obtain the
    pure Geyer IMS estimate.

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

    acf_fallback : bool, optional
        If true, use :func:`acf_autocorr_time` when the Geyer IMS estimate
        would be less than one.

    Returns
    -------
    tau : float or (float, bool)
        Estimated integrated autocorrelation time.  The optional Boolean is
        the reliability assessment from the estimator that supplies the
        returned value.
    """

    msg = 'c must be a positive finite number'
    try:
        c = float(c)
    except (TypeError,ValueError):
        raise ValueError(msg) from None
    if not np.isfinite(c) or c<=0.:
        raise ValueError(msg)

    x = np.asarray(x)
    if np.iscomplexobj(x):
        msg = 'input must be real-valued'
        raise ValueError(msg)
    if x.ndim>1 and np.max(x.shape)==x.size:
        x = x.ravel()
    if x.ndim!=1:
        msg = 'input must be one-dimensional'
        raise ValueError(msg)
    if len(x)==0:
        msg = 'input must not be empty'
        raise ValueError(msg)

    x = np.asarray(x,dtype=float)
    if not np.all(np.isfinite(x)):
        msg = 'input must contain only finite values'
        raise ValueError(msg)
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
    n         = len(x)
    nfft      = 1 << (2*n-1).bit_length()
    transform = np.fft.rfft(x,nfft)
    acf       = np.fft.irfft(transform*np.conjugate(transform),nfft)[:n]
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
        tau = 0.
    else:
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
        tau        = max(0.,float(-1.+2.*gamma_mono.sum()))

    if acf_fallback and tau<1.:
        return acf_autocorr_time(x,reliability=reliability)

    if n<c*tau:
        # The time series is shorter than c autocorrelation times; the
        # estimate may be unreliable.
        not_reliable = True

    return (tau,not_reliable) if reliability else tau
#end def geyer_ims_autocorr_time



def autocorr_time(x,reliability=False):
    """Conservatively combine autocorrelation-time estimates.

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

    t_auto       =  max(t_auto_acf,t_auto_geyer)
    not_reliable = nr_acf or nr_geyer

    return (t_auto,not_reliable) if reliability else t_auto
#end def autocorr_time


def series_stats(x,t_auto=None):
    """Return the mean, autocorrelation-adjusted error, and correlation time.

    If ``t_auto`` is not supplied, it is estimated with
    :func:`autocorr_time`. A supplied value must be positive and finite.
    """
    x = np.asarray(x)
    if np.iscomplexobj(x):
        msg = 'data array must be real-valued'
        raise ValueError(msg)
    if x.ndim>1 and np.max(x.shape)==x.size:
        x = x.ravel()
    if x.ndim!=1:
        msg = 'data array must be 1-dimensional'
        raise ValueError(msg)
    if len(x)==0:
        msg = 'data array must not be empty'
        raise ValueError(msg)
    try:
        x = np.asarray(x,dtype=float)
    except (TypeError,ValueError):
        msg = 'data array must be numeric'
        raise ValueError(msg) from None
    if not np.all(np.isfinite(x)):
        msg = 'data array must contain only finite values'
        raise ValueError(msg)

    if t_auto is None:
        t_auto = autocorr_time(x)
    else:
        try:
            t_auto = float(t_auto)
        except (TypeError,ValueError):
            msg = 'autocorrelation time must be a positive finite number'
            raise ValueError(msg) from None
        if not np.isfinite(t_auto) or t_auto<=0.:
            msg = 'autocorrelation time must be a positive finite number'
            raise ValueError(msg)
    N        = len(x)
    N_eff    = N/t_auto
    x_mean   = np.mean(x)
    x_stderr = np.std(x)/np.sqrt(N_eff)
    return x_mean,x_stderr,t_auto
#end def series_stats
