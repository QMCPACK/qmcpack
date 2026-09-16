"""Robust statistics and autocorrelation analysis for simulation data."""

import numpy as np

from .developer_tools import DevBase, dotdict, obj

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


def _real_vector(x,name):
    """Return a real vector, flattening vector-shaped arrays."""
    x = np.asarray(x)
    if np.iscomplexobj(x):
        msg = f'{name} must be real-valued'
        raise ValueError(msg)
    if x.ndim>1 and np.max(x.shape)==x.size:
        x = x.ravel()
    if x.ndim!=1:
        msg = f'{name} must be one-dimensional'
        raise ValueError(msg)
    try:
        x = np.asarray(x,dtype=float)
    except (TypeError,ValueError):
        msg = f'{name} must be numeric'
        raise ValueError(msg) from None
    return x
#end def _real_vector


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


def reblocked_autocorr_time(
        x,
        min_blocks = 10,
        *,
        plot       = False,
        show       = False,
        ):
    """Estimate autocorrelation time from the growth of blocked errors.

    This estimator currently overestimates the autocorrelation times in a
    number of cases. Prefer the Geyer method.

    For MCMC data, just use the :func:`autocorr_time` function.

    For every integer block length that leaves at least ``min_blocks`` blocks,
    contiguous block means and their standard error are computed.  Their
    ratio to the unblocked standard error is fitted as a function of block
    length with a Theil--Sen line (using reproducible stochastic pair sampling
    for large inputs).  The squared fitted ratio at the largest usable block
    length is used to obtain the auto-correlation time.

    **Strengths**: Independent block-based cross-check that can expose slow-mode
    uncertainty without requiring reversibility.

    **Weaknesses**: Its estimates have a broad 10--90% spread, and the calculated
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



def acf_autocorr_time(x,*,reliability=False):
    """Estimate autocorrelation time from a windowed sample ACF.

    Best for long chains.  Generally prefer the Geyer method.

    For MCMC data, just use the :func:`autocorr_time` function.

    The autocorrelation function is evaluated in ``O(N log N)`` time with an
    FFT and a common denominator at all lags.  A Bartlett noise estimate is
    used to locate the first sustained noise-dominated region.  Resolved lags
    are retained by a flat-top window, followed by a linear taper through the
    noisy boundary.  The returned value is the variance-inflation factor
    ``N*Var(mean)/Var(x)``.

    **Strengths**: Best long-chain accuracy, low variance, fast, and supports
    negative or oscillatory correlation.

    **Weaknesses**: Can truncate before detecting weak slow modes and generally
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



def geyer_ims_autocorr_time(
        x,
        c            = 5.0,
        *,
        reliability  = False,
        acf_fallback = True,
        ):
    """Estimate integrated autocorrelation time with Geyer's IMS method.

    This is the single best autocorrelation estimator.

    For MCMC data, just use the :func:`autocorr_time` function.

    Autocorrelations are computed with an FFT.  Geyer's initial positive
    sequence of adjacent autocorrelation pairs is then made non-increasing
    with a linear-time pool-adjacent-violators algorithm.  This estimator is
    intended primarily for stationary, reversible Markov chains.

    By default, estimates below one are replaced by the result from
    :func:`acf_autocorr_time`.  Set ``acf_fallback=False`` to obtain the
    pure Geyer IMS estimate.

    **Strengths**: Fast, stable noisy-tail treatment, and strong theoretical basis
    for reversible MCMC.

    **Shortcomings**:: Relatively variable for negative correlation.

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



def autocorr_time(x,*,reliability=False):
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
    :func:`autocorr_time`.  The returned standard error is
    ``std(x) / sqrt(N / t_auto)``, where ``std`` uses NumPy's default ``ddof=0``
    and ``N / t_auto`` is the effective number of independent samples.  Thus,
    independently sampled data have ``t_auto`` near one, while positive
    serial correlation increases the reported uncertainty.

    Parameters
    ----------
    x : array_like
        Nonempty, finite, real-valued one-dimensional sample sequence.
        Vector-shaped arrays are flattened.

    t_auto : float, optional
        Positive, finite integrated autocorrelation time.  If omitted, it is
        estimated from ``x`` with :func:`autocorr_time`.

    Returns
    -------
    mean : float
        Arithmetic mean of the samples.

    error : float
        Autocorrelation-adjusted standard error of ``mean``.

    t_auto : float
        The supplied or estimated integrated autocorrelation time.
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


############################################################################
#                                                                          #
#              Line-crossing and interval-distribution analysis            #
#              ------------------------------------------------            #
#                                                                          #
# These functions represent a time series as intervals between neighboring #
# values and count their overlap along the value axis.  The resulting      #
# interval distribution is a line-crossing density: locations with many    #
# overlapping segments identify values persistently traversed by the       #
# series.                                                                  #
#                                                                          #
# The distribution and its peak provide robust center estimates that       #
# emphasize locally stable, equilibrium-like portions of a fluctuating     #
# series.  Rolling versions track this center over time, while related     #
# utilities support broader interval-distribution analysis.                #
############################################################################


def time_series_intervals(x,t=None):
    """Return ordered intervals between adjacent time-series values.

    Parameters
    ----------
    x : array_like
        Real one-dimensional series with at least two values. Vector-shaped
        arrays are flattened.

    t : array_like, optional
        Real times paired with ``x``. If supplied, must have the same length.

    Returns
    -------
    xi : ndarray
        ``(len(x)-1, 2)`` array of ordered adjacent-value intervals.

    ti : ndarray or None
        Adjacent-pair time midpoints, or ``None`` when ``t`` is omitted.
    """
    x = _real_vector(x,'data array')
    if len(x)<2:
        msg = 'data array must contain at least two values'
        raise ValueError(msg)
    if t is not None:
        t = _real_vector(t,'time array')
        if len(t)!=len(x):
            msg = 'time array must have the same length as data array'
            raise ValueError(msg)
    xi = np.empty((len(x)-1,2),dtype=x.dtype)
    for n in range(len(x)-1):
        xi[n,0] = x[n]
        xi[n,1] = x[n+1]
    xi = np.sort(xi,axis=1)
    if t is None:
        return xi,None
    else:
        ti = (t[:-1]+t[1:])/2
        return xi,ti
#end def time_series_intervals



def _int_dist_input(x1,x2=None):
    """Normalize one interval matrix or paired lower and upper endpoints.

    Returns ordered endpoint pairs and corresponding ``+1/-1`` edge signs.
    """
    if x2 is not None:
        x1 = _real_vector(x1,'lower endpoints')
        x2 = _real_vector(x2,'upper endpoints')
        if len(x1)!=len(x2):
            msg = 'interval endpoint arrays must have equal lengths'
            raise ValueError(msg)
        xi = np.vstack((x1,x2)).T
    else:
        xi = np.asarray(x1)
        if np.iscomplexobj(xi):
            msg = 'interval array must be real-valued'
            raise ValueError(msg)
    # xi is array of N intervals
    if xi.ndim!=2 or xi.shape[1]!=2:
        msg = 'interval array must have shape (n,2)'
        raise ValueError(msg)
    if len(xi)==0:
        msg = 'interval array must not be empty'
        raise ValueError(msg)
    try:
        xi = np.asarray(xi,dtype=float)
    except (TypeError,ValueError):
        msg = 'interval array must be numeric'
        raise ValueError(msg) from None
    # check endpoint ordering
    if np.any(xi[:,1]<xi[:,0]):
        msg = 'interval upper endpoints must not be less than lower endpoints'
        raise ValueError(msg)
    si = np.empty(xi.shape,dtype=int)
    si[:,0] =  1
    si[:,1] = -1
    return xi,si
#end def _int_dist_input



def _perturb_constant_intervals(xi,perturb_const):
    """Expand constant intervals by a fixed number of floating-point steps."""
    if isinstance(perturb_const,(bool,np.bool_)) or not isinstance(
        perturb_const,(int,np.integer)
        ) or perturb_const<0:
        msg = 'constant perturbation must be a nonnegative integer'
        raise ValueError(msg)
    if perturb_const==0:
        return xi
    constant = xi[:,0]==xi[:,1]
    if not constant.any():
        return xi
    xi = xi.copy()
    lower = xi[constant,0]
    upper = xi[constant,1]
    for _ in range(perturb_const):
        lower = np.nextafter(lower,-np.inf)
        upper = np.nextafter(upper,np.inf)
    xi[constant,0] = lower
    xi[constant,1] = upper
    return xi
#end def _perturb_constant_intervals



def interval_distribution(
        x1,
        x2             = None,
        *,
        perturb_const  = 1,
        ):
    """Return spans between interval edges and their overlap counts.

    Parameters
    ----------
    x1 : array_like
        ``(n,2)`` ordered interval array, or lower endpoints when ``x2`` is
        supplied.

    x2 : array_like, optional
        Upper endpoints paired with ``x1``.

    perturb_const : int, optional
        Number of floating-point steps used to expand each zero-width input
        interval by equal step counts toward negative and positive infinity.
        One gives a deterministic ULP-scale representation of constant
        intervals; zero leaves them unexpanded.

    Returns
    -------
    xi : ndarray
        Consecutive spans between sorted unique interval endpoints.

    ci : ndarray
        Number of input intervals overlapping each span in ``xi``.

    Notes
    -----
    Each row of ``xi`` denotes the open span between consecutive unique
    endpoints. Endpoint membership is not counted separately: touching
    intervals occupy adjacent spans.  By default, zero-width intervals are
    expanded by ``perturb_const`` representable floating-point values on each
    side before the distribution is constructed.  This preserves a narrow,
    deterministic LCD contribution for repeated adjacent time-series values.
    """
    xi,si = _int_dist_input(x1,x2)
    xi = _perturb_constant_intervals(xi,perturb_const)
    # organize by edge order
    edges = xi.ravel()
    signs = si.ravel()
    order = edges.argsort()
    edges = edges[order]
    signs = signs[order]

    # Combine all coincident edges before accumulating their net change.
    # This vectorized sweep avoids one Python dictionary entry and two list
    # appends per edge while preserving the span counts between unique edges.
    values,starts = np.unique(edges,return_index=True)
    counts = np.cumsum(np.add.reduceat(signs,starts))
    xi = np.empty((len(values)-1,2),dtype=values.dtype)
    xi[:,0] = values[:-1]
    xi[:,1] = values[1:]
    ci = counts[:-1]
    return xi,ci
#end def interval_distribution



def plot_interval_dist(
        xi,
        ci,
        style = 'b.-',
        ):
    """Plot an interval distribution as a piecewise-constant curve.

    Parameters
    ----------
    xi : array_like
        ``(n,2)`` interval-distribution spans.

    ci : array_like
        Counts paired with ``xi``.

    style : str, optional
        Matplotlib style specification for the distribution line.
    """
    import matplotlib.pyplot as plt
    xi,_ = _int_dist_input(xi)
    ci = _real_vector(ci,'interval counts')
    if len(ci)!=len(xi):
        msg = 'interval counts must have the same length as intervals'
        raise ValueError(msg)
    xif = xi.ravel()
    cif = np.zeros(xif.shape)
    cif[::2]  = ci
    cif[1::2] = ci
    plt.axhline(0,color='k')
    plt.plot(xif,cif,style)
#end def plot_interval_dist



def interval_dist_peak(
        xi,
        ci,
        method         = 'interval_mid',
        peak_frac      = 0.5,
        *,
        height         = False,
        quad_weighting = 'endpoint',
        perturb_const  = 1,
        ):
    """Return a representative location at the peak of an interval distribution.

    Parameters
    ----------
    xi : array_like
        ``(n,2)`` interval-distribution spans.

    ci : array_like
        Counts paired with ``xi``.

    method : {'interval_mid', 'interval_rand', 'quad_peak'}, optional
        Peak estimator. The first averages all maximum-count intervals, the
        second averages random interior samples of those intervals, and the
        third fits each separated high-count peak region quadratically.

    peak_frac : float, optional
        Fraction of the maximum count retained for each quadratic-fit region.
        Must lie in ``(0,1]``.

    height : bool, optional
        If true, return the peak location and its estimated height.

    quad_weighting : {'endpoint', 'width'}, optional
        Weighting used only by ``'quad_peak'``. ``'endpoint'`` gives every
        duplicated interval endpoint equal fit weight. ``'width'`` weights
        each endpoint by the square root of its interval width, making the
        least-squares objective proportional to interval width.

    perturb_const : int, optional
        Number of floating-point steps used to expand any zero-width spans
        before locating the peak.  One gives a deterministic ULP-scale
        representation; zero leaves the spans unexpanded.

    Returns
    -------
    peak : float or (float, float)
        Peak location, optionally followed by peak height. Separated
        equal-height modes are averaged.
    """
    xi,_ = _int_dist_input(xi)
    xi = _perturb_constant_intervals(xi,perturb_const)
    ci = _real_vector(ci,'interval counts')
    if len(ci)!=len(xi):
        msg = 'interval counts must have the same length as intervals'
        raise ValueError(msg)
    if not isinstance(method,str):
        msg = 'peak method must be a string'
        raise TypeError(msg)
    try:
        peak_frac = float(peak_frac)
    except (TypeError,ValueError):
        msg = 'peak fraction must be a finite number'
        raise ValueError(msg) from None
    if not np.isfinite(peak_frac) or not 0.<peak_frac<=1.:
        msg = 'peak fraction must be in the interval (0,1]'
        raise ValueError(msg)
    if not isinstance(quad_weighting,str) or quad_weighting not in (
        'endpoint','width'
        ):
        msg = 'quadratic weighting must be "endpoint" or "width"'
        raise ValueError(msg)
    if method=='interval_mid':
        cm = ci.max()
        xm = xi[ci==cm].mean()
    elif method=='interval_rand':
        cm    = ci.max()
        xi    = xi[ci==cm]
        u     = np.random.uniform(size=len(xi))
        x1,x2 = xi.T
        dx    = x2-x1
        xmid  = (x2+x1)/2
        x     = xmid + (u-0.5)*dx/2
        xm    = x.mean()
    elif method=='quad_peak':
        cm = ci.max()
        cf = peak_frac*cm
        high = ci>=cf
        edges = np.flatnonzero(np.diff(np.r_[False,high,False]))
        xpeaks = []
        cpeaks = []
        for i1,i2 in zip(edges[::2],edges[1::2]-1, strict=True):
            ci_region = ci[i1:i2+1]
            if ci_region.max()!=cm:
                continue
            xi_region = xi[i1:i2+1]
            peak_mean = xi_region[ci_region==cm].mean()
            xp = xi_region.ravel()
            cp = np.repeat(ci_region,2)
            weights = None
            if quad_weighting=='width':
                widths = xi_region[:,1]-xi_region[:,0]
                weights = np.repeat(np.sqrt(widths),2)
                nonzero = weights>0.
                xp = xp[nonzero]
                cp = cp[nonzero]
                weights = weights[nonzero]
            if len(np.unique(xp))<3:
                # A single usable span cannot determine a quadratic peak.
                xp = peak_mean
                cp = cm
            else:
                if weights is None:
                    p = np.polyfit(xp,cp,2)
                else:
                    p = np.polyfit(xp,cp,2,w=weights)
                if not np.isfinite(p[0]) or p[0]>=0.:
                    # A non-concave fit has no interior maximum.
                    xp = peak_mean
                    cp = cm
                else:
                    xp = -p[1]/(2*p[0])
                    xp = np.clip(xp,xi_region.min(),xi_region.max())
                    cp = np.polyval(p,xp)
            xpeaks.append(xp)
            cpeaks.append(cp)
        xm = np.mean(xpeaks)
        cm = np.mean(cpeaks)
    else:
        msg = f'unrecognized int. dist. max method: "{method}"'
        raise ValueError(msg)
    if not height:
        return xm
    else:
        return xm,cm
#end def interval_dist_peak



def rolling_interval_dist_peak(
        x1,
        x2             = None,
        window         = 10,
        step           = 5,
        method         = 'interval_mid',
        peak_frac      = 0.5,
        *,
        quad_weighting = 'endpoint',
        ret_height     = False,
        ret_windows    = False,
        ):
    """Return interval-distribution peaks for overlapping input windows.

    Parameters
    ----------
    x1 : array_like
        ``(n,2)`` interval array, or lower endpoints when ``x2`` is given.

    x2 : array_like, optional
        Upper endpoints paired with ``x1``.

    window : int, optional
        Number of input intervals in each rolling distribution.

    step : int, optional
        Number of intervals between successive window starts. It must not
        exceed ``window``.

    method : {'interval_mid', 'interval_rand', 'quad_peak'}, optional
        Peak estimator passed to :func:`interval_dist_peak`.

    peak_frac : float, optional
        Quadratic peak-region threshold passed to :func:`interval_dist_peak`.

    quad_weighting : {'endpoint', 'width'}, optional
        Quadratic-fit weighting passed to :func:`interval_dist_peak`.

    ret_height : bool, optional
        Include a peak-height array in the returned tuple.

    ret_windows : bool, optional
        Include ``(start, stop)`` bounds for each returned window.

    Returns
    -------
    result : tuple
        A tuple beginning with the peak-location array. When requested, it
        then contains the peak-height array and/or window-bound list, in that
        order. Windows advance by ``step``; a final window ending at the last
        input interval is appended when the regular sequence does not reach
        it exactly.
    """
    for value,name in ((window,'window'),(step,'step')):
        if isinstance(value,(bool,np.bool_)) or not isinstance(
            value,(int,np.integer)
            ) or value<1:
            msg = f'{name} must be a positive integer'
            raise ValueError(msg)
    if step>window:
        msg = 'step must not exceed window'
        raise ValueError(msg)
    interval_mid  = method=='interval_mid'
    interval_rand = method=='interval_rand'
    interval_quad = method=='quad_peak'
    if not interval_mid and not interval_rand and not interval_quad:
        msg = f'method "{method}" is unrecognized'
        raise ValueError(msg)
    # map inputs to intervals
    xia,_ = _int_dist_input(x1,x2)
    N = len(xia)
    if window > N:
        msg = 'window must not exceed the number of intervals'
        raise ValueError(msg)
    # find window segments
    starts = list(range(0,N-window+1,step))
    if starts[-1]!=N-window:
        starts.append(N-window)
    windows = [(i1,i1+window) for i1 in starts]
    # find interval dist peaks in each window
    xp = []
    cp = []
    for i1,i2 in windows:
        xi,ci = interval_distribution(xia[i1:i2],perturb_const=1)
        if len(ci)==0:
            msg = 'each rolling window must span a nonzero interval'
            raise ValueError(msg)
        xm,cm = interval_dist_peak(
            xi,
            ci,
            method         = method,
            peak_frac      = peak_frac,
            height         = True,
            quad_weighting = quad_weighting,
            perturb_const  = 1,
            )
        xp.append(xm)
        cp.append(cm)
    xp = np.array(xp)
    ret = [xp]
    if ret_height:
        cp = np.array(cp)
        ret.append(cp)
    if ret_windows:
        ret.append(windows)
    if len(ret)==0:
        return ret[0]
    else:
        return tuple(ret)
#end def rolling_interval_dist_peak



def _perturb_constant_series(x,perturb_const=1):
    """Alternately displace each repeated-value run by floating-point steps."""
    if isinstance(perturb_const,(bool,np.bool_)) or not isinstance(
        perturb_const,(int,np.integer)
        ) or perturb_const<0:
        msg = 'constant perturbation must be a nonnegative integer'
        raise ValueError(msg)
    if perturb_const==0:
        return x
    repeated = x[1:]==x[:-1]
    if not repeated.any():
        return x

    xp = x.copy()
    starts = np.r_[0,np.flatnonzero(x[1:]!=x[:-1])+1]
    stops = np.r_[starts[1:],len(x)]
    for i1,i2 in zip(starts,stops, strict=True):
        if i2-i1<2:
            continue
        lower = x[i1]
        upper = x[i1]
        for _ in range(perturb_const):
            lower = np.nextafter(lower,-np.inf)
            upper = np.nextafter(upper,np.inf)
        xp[i1:i2:2] = lower
        xp[i1+1:i2:2] = upper
    return xp
#end def _perturb_constant_series



def line_crossing_distribution(x,nperm=0, *, ret_x=False):
    """Return the line-crossing distribution of a series or its permutations.

    Parameters
    ----------
    x : array_like
        Real one-dimensional series with at least two values.

    nperm : int, optional
        Number of independently shuffled series to average. Zero evaluates
        the input series directly.

    ret_x : bool, optional
        Also return a copy of the input series in which each run of repeated
        values alternates between deterministic ULP-scale displacements about
        its original value.  This representation gives repeated values a
        nonzero side relative to an LCD peak while preserving their center.

    Returns
    -------
    xi : ndarray
        Line-crossing distribution spans.

    ci : ndarray
        Crossing counts, averaged over permutations when ``nperm`` is
        positive. Permutations are never connected to one another.

    xp : ndarray, optional
        Perturbed time series, returned only when ``ret_x`` is true.

    Notes
    -----
    Each adjacent pair defines an interval, and the distribution count at a
    value is the number of such intervals that span it.  For a continuous
    equilibrium series with independent samples ``X`` and ``Y`` drawn from
    CDF ``F``, the corresponding crossing probability is

    .. math::

       L(z) = P(\\min(X,Y) < z < \\max(X,Y)) = 2F(z)[1-F(z)].

    Thus, for a series with ``N`` samples, the expected count is
    :math:`(N - 1) L(z)`.  The distribution is maximized at a median of the
    sampled distribution, which motivates its use as a robust equilibrium
    location estimator.  It is a crossing-rate curve rather than a normalized
    probability density; when :math:`E[|X-Y|]` is finite, its normalized form is
    :math:`2 F(z) [1-F(z)] / E[|X-Y|]`.

    In the ideal continuous i.i.d. case, a probability-scale LCD can be
    inverted to obtain :math:`F(z) = (1 - \\sqrt{1 - 2 L(z)}) / 2` below a
    median and :math:`F(z) = (1 + \\sqrt{1 - 2 L(z)}) / 2` above one, followed
    by differentiation to obtain the density.  Empirical inversion is noisy,
    and the LCD does not uniquely determine distributions with atoms or gaps.
    """
    x = _real_vector(x,'data array')
    if len(x)<2:
        msg = 'data array must contain at least two values'
        raise ValueError(msg)
    if isinstance(nperm,(bool,np.bool_)) or not isinstance(
        nperm,(int,np.integer)
        ) or nperm<0:
        msg = 'number of permutations must be a nonnegative integer'
        raise ValueError(msg)
    if not isinstance(ret_x,(bool,np.bool_)):
        msg = 'ret_x must be a Boolean value'
        raise TypeError(msg)
    xp = _perturb_constant_series(x)

    # permutation-free (typical) case
    if nperm==0:
        xi,_ = time_series_intervals(x,t=None)
        xi,ci = interval_distribution(xi,perturb_const=1)
        if ret_x:
            return xi,ci,xp
        return xi,ci

    # use permutation shuffling
    permutation_intervals = []
    for _ in range(nperm):
        xp = x.copy()
        np.random.shuffle(xp)
        xi,_ = time_series_intervals(xp,t=None)
        permutation_intervals.append(xi)
    xi,ci = interval_distribution(
        np.vstack(permutation_intervals),perturb_const=1
        )
    ci = ci/nperm
    if ret_x:
        return xi,ci,xp
    return xi,ci
#end def line_crossing_distribution



def lcd_peak(
        x,
        method         = 'interval_mid',
        peak_frac      = 0.5,
        nperm          = 0,
        quad_weighting = 'endpoint',
        ):
    """Return a peak of a series line-crossing distribution.

    Parameters
    ----------
    x : array_like
        Time series supplied to :func:`line_crossing_distribution`.

    method, peak_frac, quad_weighting, nperm
        Options forwarded to :func:`interval_dist_peak` and
        :func:`line_crossing_distribution`.

    Returns
    -------
    float
        Estimated line-crossing-distribution peak.
    """
    xi,ci = line_crossing_distribution(x,nperm=nperm)
    xp = interval_dist_peak(
        xi,
        ci,
        method         = method,
        peak_frac      = peak_frac,
        quad_weighting = quad_weighting,
        perturb_const  = 1,
        )
    return xp
#end def lcd_peak



def lcd_smooth(
        x,
        t              = None,
        window         = 10,
        step           = 5,
        method         = 'interval_rand',
        peak_frac      = 0.5,
        quad_weighting = 'endpoint',
        ):
    """Return rolling line-crossing-distribution peaks for a time series.

    Parameters
    ----------
    x : array_like
        Time series to smooth with rolling line-crossing peaks.

    t : array_like, optional
        Times paired with ``x``.

    window, step, method, peak_frac, quad_weighting
        Options forwarded to :func:`rolling_interval_dist_peak`.

    Returns
    -------
    peaks : ndarray or (ndarray, ndarray)
        Rolling peak locations, optionally paired with their mean window
        times when ``t`` is supplied.
    """
    xi,ti = time_series_intervals(x,t)
    xp,windows = rolling_interval_dist_peak(
        xi,
        window         = window,
        step           = step,
        method         = method,
        peak_frac      = peak_frac,
        quad_weighting = quad_weighting,
        ret_windows    = True,
        )
    if t is None:
        return xp
    else:
        tp = np.array([ti[i1:i2].mean() for i1,i2 in windows])
        return xp,tp
#end def lcd_smooth


def pair_expand_ts_intervals(
        x,
        t      = None,
        expand = 10,
        ):
    """Return sorted pairs between samples in a bounded local neighborhood.

    Each sample is paired with up to ``expand/2`` earlier and later samples.
    Consequently, every unordered pair within that index separation occurs
    twice, once from each endpoint.  This supplies local multi-lag intervals
    for line-crossing analysis.

    Parameters
    ----------
    x : array_like
        Real vector-like time series.  Singleton dimensions are flattened.

    t : array_like, optional
        Real vector-like sample times paired with ``x``.

    expand : int, optional
        Positive even number of local neighbor positions.  It must be less
        than ``len(x)``.

    Returns
    -------
    xi : ndarray
        Sorted endpoint pairs for all local directed pairings.

    ti : ndarray or None
        Constructed sub-times paired with ``xi``, or ``None`` if ``t`` is not
        supplied.
    """
    x = _real_vector(x,'data array')
    if isinstance(expand,(bool,np.bool_)) or not isinstance(
        expand,(int,np.integer)
        ) or expand<1 or expand%2!=0:
        msg = 'expansion must be a positive even integer'
        raise ValueError(msg)
    if len(x)<=expand:
        msg = 'data array length must exceed expansion'
        raise ValueError(msg)
    if t is not None:
        t = _real_vector(t,'time array')
        if len(t)!=len(x):
            msg = 'time array must have the same length as data array'
            raise ValueError(msg)
    ne = expand//2
    N  = len(x)
    xi = []
    for i,x0 in enumerate(x):
        i1 = max(i-ne,0)
        i2 = min(i+ne,N-1)
        for j in range(i1,i2+1):
            if j==i:
                continue
            xi.append((x0,x[j]))
    xi = np.array(xi)
    xi = np.sort(xi,axis=1)
    ti = None
    if t is not None:
        ti = []
        for i,t0 in enumerate(t):
            if i==0:
                dt = t[i+1]-t0
                dtj = dt/2/(ne+1)
                ti.extend([t0+(j+1)*dtj for j in range(ne)])
            elif i==N-1:
                dt = t0-t[i-1]
                dtj = dt/2/(ne+1)
                ti.extend([t0-dt/2+(j+1)*dtj for j in range(ne)])
            else:
                n1 = min(ne,i)
                n2 = min(ne,N-1-i)
                dt = t0-t[i-1]
                dtj = dt/2/(n1+1)
                ti.extend([t0-dt/2+(j+1)*dtj for j in range(n1)])
                dt = t[i+1]-t0
                dtj = dt/2/(n2+1)
                ti.extend([t0+(j+1)*dtj for j in range(n2)])
        ti = np.array(ti)
        if len(ti)!=len(xi):
            msg = 'internal pair/time expansion length mismatch'
            raise RuntimeError(msg)
    return xi,ti
#end def pair_expand_ts_intervals



def _find_segments(
        x,
        mask,
        seg_min = 1,
        ):
    """Return contiguous true-mask index spans meeting a minimum length."""
    if isinstance(seg_min,(bool,np.bool_)) or not isinstance(
        seg_min,(int,np.integer)
        ) or seg_min<1:
        msg = 'minimum segment length must be a positive integer'
        raise ValueError(msg)
    if len(x)!=len(mask):
        msg = 'mask must have the same length as data array'
        raise ValueError(msg)
    seg = []
    n1 = None
    for n,in_seg in enumerate(mask):
        if n1 is None and in_seg:
            n1 = n
        if n1 is not None and not in_seg:
            if n-n1>=seg_min:
                seg.append((n1,n))
            n1 = None
    if n1 is not None and len(x)-n1>=seg_min:
        seg.append((n1,len(x)))
    return seg
#end def _find_segments



def _lcd_trim_input(x,niter):
    """Validate trimming inputs and obtain its LCD peak and perturbed series."""
    x = _real_vector(x,'data array')
    if len(x)<2:
        msg = 'data array must contain at least two values'
        raise ValueError(msg)
    if isinstance(niter,(bool,np.bool_)) or not isinstance(
        niter,(int,np.integer)
        ) or niter<1:
        msg = 'number of trim iterations must be a positive integer'
        raise ValueError(msg)
    xi,ci,x = line_crossing_distribution(x,ret_x=True)
    x_lcd = interval_dist_peak(xi,ci,perturb_const=1)
    return x,x_lcd,int(niter)
#end def _lcd_trim_input



def _lcd_trim_options(ret_seg,ret_mask):
    """Validate trim return selections."""
    for value,name in ((ret_seg,'ret_seg'),(ret_mask,'ret_mask')):
        if not isinstance(value,(bool,np.bool_)):
            msg = f'{name} must be a Boolean value'
            raise TypeError(msg)
#end def _lcd_trim_options



def _trim_run(
        x,
        x_lcd,
        start,
        stop,
        ):
    """Return the leading count before a crossing, retaining one endpoint."""
    if start>=stop:
        return 0
    val_sign = np.sign(x[start]-x_lcd)
    for n in range(start,stop):
        if np.sign(x[n]-x_lcd)*val_sign<0:
            return n-start
    return stop-start-1
#end def _trim_run



def lcd_trim_l(
        x,
        niter    = 3,
        *,
        ret_seg  = True,
        ret_mask = False,
        ):
    """Trim initial runs separated from the LCD peak by sign crossings.

    The LCD peak is calculated once from the full series.  Starting at the
    left endpoint, each iteration removes values up to the next crossing of
    that peak.  Repeated values receive deterministic ULP-scale perturbations
    for the crossing comparison, while returned masks still index the input
    series.

    Parameters
    ----------
    x : array_like
        Real vector-like time series with at least two values.

    niter : int, optional
        Positive number of consecutive left-side runs to remove.

    ret_seg : bool, optional
        Include clean and left-trim index spans in the result.

    ret_mask : bool, optional
        Include clean and left-trim Boolean masks in the result.

    Returns
    -------
    result : tuple
        Requested outputs in this order: ``clean_segments``,
        ``left_segment``, ``clean_mask``, and ``left_mask``.  Segment outputs
        are omitted when ``ret_seg`` is false and mask outputs are omitted
        when ``ret_mask`` is false.
    """
    x,x_lcd,niter = _lcd_trim_input(x,niter)
    _lcd_trim_options(ret_seg,ret_mask)
    mask_left  = np.zeros(len(x),dtype=bool)
    ntrim_l = 0
    # perform lcd trim
    for _ in range(niter):
        ntrim_l += _trim_run(x,x_lcd,ntrim_l,len(x))
        mask_left[:ntrim_l] = True
    mask_clean = ~mask_left
    seg_l = 0,ntrim_l
    seg_c = _find_segments(x,mask_clean,1)
    ret = []
    if ret_seg:
        ret.extend([seg_c,seg_l])
    if ret_mask:
        ret.extend([mask_clean,mask_left])
    return tuple(ret)
#end def lcd_trim_l


def lcd_trim_r(
        x,
        niter    = 3,
        *,
        ret_seg  = True,
        ret_mask = False,
        ):
    """Trim terminal runs separated from the LCD peak by sign crossings.

    This is the right-to-left counterpart of :func:`lcd_trim_l`: it starts at
    the final sample and removes up to ``niter`` consecutive runs before a
    crossing of the full-series LCD peak.  Constant-value runs are compared
    through the deterministic perturbed representation used by the LCD.

    Parameters
    ----------
    x : array_like
        Real vector-like time series with at least two values.

    niter : int, optional
        Positive number of consecutive right-side runs to remove.

    ret_seg : bool, optional
        Include clean and right-trim index spans in the result.

    ret_mask : bool, optional
        Include clean and right-trim Boolean masks in the result.

    Returns
    -------
    result : tuple
        Requested outputs in this order: ``clean_segments``,
        ``right_segment``, ``clean_mask``, and ``right_mask``.  Segment and
        mask groups are controlled by ``ret_seg`` and ``ret_mask``.
    """
    x,x_lcd,niter = _lcd_trim_input(x,niter)
    _lcd_trim_options(ret_seg,ret_mask)
    xr = np.flip(x)
    mask_right = np.zeros(len(x),dtype=bool)
    ntrim_r = 0
    # perform lcd trim
    for _ in range(niter):
        ntrim_r += _trim_run(xr,x_lcd,ntrim_r,len(x))
        mask_right[len(x)-ntrim_r:] = True
    mask_clean = ~mask_right
    seg_r = len(x)-ntrim_r,len(x)
    seg_c = _find_segments(x,mask_clean,1)
    ret = []
    if ret_seg:
        ret.extend([seg_c,seg_r])
    if ret_mask:
        ret.extend([mask_clean,mask_right])
    return tuple(ret)
#end def lcd_trim_r


def lcd_trim_lr(
        x,
        niter    = 3,
        *,
        ret_seg  = True,
        ret_mask = False,
        ):
    """Trim leading and trailing runs according to the full-series LCD peak.

    Left and right runs are removed independently on every iteration.  Each
    side stops at its next crossing of the LCD peak, and the opposite trim
    boundary limits the search so the masks remain disjoint and retain at
    least one clean sample.

    Parameters
    ----------
    x : array_like
        Real vector-like time series with at least two values.

    niter : int, optional
        Positive number of trim iterations on each endpoint.

    ret_seg : bool, optional
        Include clean, left-trim, and right-trim index spans.

    ret_mask : bool, optional
        Include corresponding Boolean masks.

    Returns
    -------
    result : tuple
        Requested outputs in this order: ``clean_segments``, ``left_segment``,
        ``right_segment``, ``clean_mask``, ``left_mask``, and ``right_mask``.
        Segment and mask groups are controlled by ``ret_seg`` and
        ``ret_mask``.
    """
    x,x_lcd,niter = _lcd_trim_input(x,niter)
    _lcd_trim_options(ret_seg,ret_mask)
    xr = np.flip(x)
    mask_left  = np.zeros(len(x),dtype=bool)
    mask_right = np.zeros(len(x),dtype=bool)
    ntrim_l = 0
    ntrim_r = 0
    # perform lcd trim
    for _ in range(niter):
        ntrim_l += _trim_run(x,x_lcd,ntrim_l,len(x)-ntrim_r)
        mask_left[:ntrim_l] = True
        ntrim_r += _trim_run(xr,x_lcd,ntrim_r,len(x)-ntrim_l)
        mask_right[len(x)-ntrim_r:] = True
    mask_clean = (~mask_left)&(~mask_right)
    seg_l = 0,ntrim_l
    seg_r = len(x)-ntrim_r,len(x)
    seg_c = _find_segments(x,mask_clean,1)
    ret = []
    if ret_seg:
        ret.extend([seg_c,seg_l,seg_r])
    if ret_mask:
        ret.extend([mask_clean,mask_left,mask_right])
    return tuple(ret)
#end def lcd_trim_lr


def lcd_trim_lrm(
        x,
        niter    = 3,
        low_scale = 2.,
        nseg_min = 4,
        *,
        ret_seg  = True,
        ret_mask = False,
        ):
    """Trim endpoint runs and sustained low-valued interior excursions.

    Endpoint removal follows :func:`lcd_trim_lr`.  Among the remaining
    samples, candidate middle regions below a cutoff derived from the LCD
    peak and retained maximum are retained only when they contain at least
    ``nseg_min`` samples.  Qualifying regions are extended until they reach
    values at or above the LCD peak.  Later iterations do not trim an endpoint
    through an identified middle region.

    Parameters
    ----------
    x : array_like
        Real vector-like time series with at least two values.

    niter : int, optional
        Positive number of endpoint and middle-trimming iterations.

    low_scale : float, optional
        Positive multiplier defining the low-value cutoff relative to the LCD
        peak and the maximum value retained after endpoint trimming.

    nseg_min : int, optional
        Minimum length of a candidate middle region before it is removed.

    ret_seg : bool, optional
        Include clean, left-trim, right-trim, and middle-trim index spans.

    ret_mask : bool, optional
        Include corresponding Boolean masks.

    Returns
    -------
    result : tuple
        Requested outputs in this order: ``clean_segments``, ``left_segment``,
        ``right_segment``, ``middle_segments``, ``clean_mask``, ``left_mask``,
        ``right_mask``, and ``middle_mask``.  Segment and mask groups are
        controlled by ``ret_seg`` and ``ret_mask``.
    """
    x,x_lcd,niter = _lcd_trim_input(x,niter)
    _lcd_trim_options(ret_seg,ret_mask)
    try:
        low_scale = float(low_scale)
    except (TypeError,ValueError):
        msg = 'low scale must be a positive finite number'
        raise ValueError(msg) from None
    if not np.isfinite(low_scale) or low_scale<=0.:
        msg = 'low scale must be a positive finite number'
        raise ValueError(msg)
    if isinstance(nseg_min,(bool,np.bool_)) or not isinstance(
        nseg_min,(int,np.integer)
        ) or nseg_min<1:
        msg = 'minimum segment length must be a positive integer'
        raise ValueError(msg)
    xr = np.flip(x)
    mask_left  = np.zeros(len(x),dtype=bool)
    mask_right = np.zeros(len(x),dtype=bool)
    ntrim_l = 0
    ntrim_r = 0
    seg_m = []
    # perform lcd trim
    for ni in range(niter):
        # left trim
        tleft = True
        if ni>0:
            for n1,n2 in seg_m:  # noqa: B007
                tleft &= ntrim_l < n1
        if tleft:
            ntrim_l += _trim_run(x,x_lcd,ntrim_l,len(x)-ntrim_r)
            mask_left[:ntrim_l] = True
        # right trim
        tright = True
        if ni>0:
            for n1,n2 in seg_m:  # noqa: B007
                tright &= len(x)-ntrim_r > n2
        if tright:
            ntrim_r += _trim_run(xr,x_lcd,ntrim_r,len(x)-ntrim_l)
            mask_right[len(x)-ntrim_r:] = True
        # mid trim
        #  initial low cut to identify mid trim segments
        mask_clean = (~mask_left)&(~mask_right)
        xc = x[mask_clean]
        xc_max = xc.max()
        xm_cut = x_lcd-low_scale*(xc_max-x_lcd)
        candidates = (x<xm_cut)&mask_clean
        seg_m = _find_segments(x,candidates,nseg_min)
        mask_mid = np.zeros(len(x),dtype=bool)
        #  extend the segments left and right
        for n1,n2 in seg_m:
            while n1>ntrim_l and x[n1]<x_lcd:
                n1 -= 1
            while n2<len(x)-ntrim_r and x[n2]<x_lcd:
                n2 += 1
            mask_mid[n1:n2] = True
        mask_mid = mask_mid&(~mask_left)&(~mask_right)
        #  find enlarged segments w/ >neg_min points
        seg_m = _find_segments(x,mask_mid,nseg_min)
        # find final clean points
        candidate_clean = (~mask_left)&(~mask_right)&(~mask_mid)
        if candidate_clean.any():
            mask_clean = candidate_clean
        else:
            mask_mid = np.zeros(len(x),dtype=bool)
            seg_m = []
            mask_clean = (~mask_left)&(~mask_right)
        seg_c = _find_segments(x,mask_clean,1)

    seg_l = 0,ntrim_l
    seg_r = len(x)-ntrim_r,len(x)
    ret = []
    if ret_seg:
        ret.extend([seg_c,seg_l,seg_r,seg_m])
    if ret_mask:
        ret.extend([mask_clean,mask_left,mask_right,mask_mid])
    return tuple(ret)
#end def lcd_trim_lrm


############################################################################
#                                                                          #
#                         Local series smoothers                           #
#                         ----------------------                           #
#                                                                          #
# The smoothers reduce short-scale variation while preserving the length   #
# and ordering of a series.  Their centered windows taper at endpoints,    #
# leaving the first and last values unchanged.                             #
#                                                                          #
# Mean smoothing provides simple local averaging.  Median smoothing is     #
# more resistant to isolated outliers and may be followed by mean          #
# smoothing.  Polynomial smoothing fits low-order local trends and may     #
# likewise receive a final mean pass.                                      #
#                                                                          #
# Local-median smoothing accepts one sample set per position.  It pools    #
# nearby sets while omitting the current one, takes a robust local median, #
# and then applies polynomial or mean smoothing to the resulting series.   #
#                                                                          #
############################################################################


def _smoothing_window_length(n,m,maximum=None):
    """Validate or select an odd smoothing-window length."""
    if m is None:
        if n==0:
            return None
        m = min(n//6,15)
        m = max(3,2*(m//2)+1)
        m = min(m,n if n%2 else n-1)
    elif not isinstance(m,(int,np.integer)) or isinstance(m,(bool,np.bool_)):
        msg = 'smoothing window length must be an integer'
        raise TypeError(msg)

    m = int(m)
    if m<1:
        msg = 'smoothing window length must be positive'
        raise ValueError(msg)
    if m%2==0:
        msg = 'smoothing window length must be odd'
        raise ValueError(msg)
    if m>n:
        msg = 'smoothing window length must not exceed the data length'
        raise ValueError(msg)
    if maximum is not None and m>maximum:
        msg = f'smoothing window length must not exceed {maximum}'
        raise ValueError(msg)
    return m
#end def _smoothing_window_length


def mean_smooth(x,m=None):
    """Smooth a sequence with tapered-endpoint moving averages.

    Each interior value is replaced by the mean in a centered, odd-length
    window of width ``m``.  Near either endpoint the window is shortened to
    remain symmetric about the current value, so the first and last samples
    are unchanged.  The result always has the same length as ``x``.

    Parameters
    ----------
    x : sequence
        Values to smooth.

    m : int, optional
        Positive odd window width no greater than ``len(x)``.  By default,
        use the largest applicable odd width no greater than
        ``min(len(x)/6, 15)``, nominally with a minimum of three.

    Returns
    -------
    ndarray
        Smoothed values.
    """
    N = len(x)
    m = _smoothing_window_length(N,m)
    if m is None:
        return np.array([])
    dm = m//2
    xs = []
    for n in range(N):
        if n<dm:
            n1 = 0
            n2 = 2*n+1
        elif N-1-n<dm:
            n1 = (N-1)-2*(N-1-n)
            n2 = N
        else:
            n1 = n-dm
            n2 = n+dm+1
        xsl = np.array(x[n1:n2])
        xsn = np.mean(xsl)
        xs.append(xsn)
    xs = np.array(xs)
    return xs
#end def mean_smooth


def median_smooth(x,m=None,*,post_mean=False):
    """Smooth a sequence with local medians, optionally followed by means.

    Local windows and endpoint treatment are the same as :func:`mean_smooth`,
    but each value is replaced by the window median.  This is less sensitive
    to isolated spikes.  When requested, a moving-mean pass is applied to the
    median-smoothed result.

    Parameters
    ----------
    x : sequence
        Values to smooth.

    m : int, optional
        Positive odd window width no greater than ``len(x)``.  The default is
        selected as in :func:`mean_smooth`.

    post_mean : bool, optional
        Apply :func:`mean_smooth` with the same width after the median pass.

    Returns
    -------
    ndarray
        Smoothed values.
    """
    N = len(x)
    m = _smoothing_window_length(N,m)
    if not isinstance(post_mean,(bool,np.bool_)):
        msg = 'post_mean must be a Boolean value'
        raise TypeError(msg)
    post_mean = bool(post_mean)
    if m is None:
        return np.array([])
    dm = m//2
    xs = []
    for n in range(N):
        if n<dm:
            n1 = 0
            n2 = 2*n+1
        elif N-1-n<dm:
            n1 = (N-1)-2*(N-1-n)
            n2 = N
        else:
            n1 = n-dm
            n2 = n+dm+1
        xsl = np.array(x[n1:n2])
        xsn = np.median(xsl)
        xs.append(xsn)
    xs = np.array(xs)
    if post_mean:
        xs = mean_smooth(xs,m=m)
    return xs
#end def median_smooth


def poly_smooth(x,m=None,*,post_mean=False):
    """Smooth a sequence by evaluating local polynomial fits.

    A polynomial is fitted in each centered window and evaluated at the
    current index.  Endpoint windows are shortened symmetrically as in
    :func:`mean_smooth`; a one-value window is returned unchanged.  Polynomial
    order increases gradually with window size, from constant for one sample
    to quartic for windows of 13--21 samples.  An optional moving-mean pass
    can further reduce residual variation.

    Parameters
    ----------
    x : sequence
        Values to smooth.

    m : int, optional
        Positive odd window width from 1 through 21 and no greater than
        ``len(x)``.  The default is selected as in :func:`mean_smooth` and is
        at most 15.

    post_mean : bool, optional
        Apply :func:`mean_smooth` with the same width after polynomial
        smoothing.

    Returns
    -------
    ndarray
        Smoothed values.
    """
    poly_order = {1:0,3:1,5:2,7:2,9:3,11:3,
                  13:4,15:4,17:4,19:4,21:4}
    N = len(x)
    m = _smoothing_window_length(N,m,maximum=21)
    if not isinstance(post_mean,(bool,np.bool_)):
        msg = 'post_mean must be a Boolean value'
        raise TypeError(msg)
    post_mean = bool(post_mean)
    if m is None:
        return np.array([])
    dm = m//2
    xs = []
    for n in range(N):
        if n<dm:
            n1 = 0
            n2 = 2*n+1
        elif N-1-n<dm:
            n1 = (N-1)-2*(N-1-n)
            n2 = N
        else:
            n1 = n-dm
            n2 = n+dm+1
        xsl = np.array(x[n1:n2])
        if len(xsl)>1:
            porder = poly_order[len(xsl)]
            p = np.polyfit(np.arange(n1,n2),xsl,porder)
            xsn = np.polyval(p,n)
        else:
            xsn = xsl[0]
        xs.append(xsn)
    xs = np.array(xs)
    if post_mean:
        xs = mean_smooth(xs,m=m)
    return xs
#end def poly_smooth
poly_smooth_ = poly_smooth


def local_median_smooth(x_list,m=None,*,poly_smooth=True,post_mean=False):
    """Smooth a sequence of sample sets through leave-one-out local medians.

    For each position, all neighboring sample sets in a centered window are
    pooled, excluding the sample set at that position, and their median is
    taken.  At an endpoint where the tapered window contains only that sample
    set, its own median is used.  The resulting median sequence is then
    polynomial-smoothed by default, or mean-smoothed when ``poly_smooth`` is
    false; either result can receive a final mean-smoothing pass.  This is
    not a batched version of :func:`median_smooth`.

    Parameters
    ----------
    x_list : sequence of array_like
        Per-position sample sets to pool locally.

    m : int, optional
        Positive odd window width no greater than ``len(x_list)``.  The
        default is selected as in :func:`mean_smooth` using ``len(x_list)``.

    poly_smooth : bool, optional
        Use :func:`poly_smooth` for the second pass.  If false, use
        :func:`mean_smooth` instead.

    post_mean : bool, optional
        Apply a final :func:`mean_smooth` pass with the same width.

    Returns
    -------
    ndarray
        One smoothed value for every input sample set.
    """
    # x_list: list of arrays containing trace/time-series data
    N = len(x_list)
    if not isinstance(poly_smooth,(bool,np.bool_)):
        msg = 'poly_smooth must be a Boolean value'
        raise TypeError(msg)
    poly_smooth = bool(poly_smooth)
    if not isinstance(post_mean,(bool,np.bool_)):
        msg = 'post_mean must be a Boolean value'
        raise TypeError(msg)
    post_mean = bool(post_mean)
    m = _smoothing_window_length(N,m,maximum=21 if poly_smooth else None)
    if m is None:
        return np.array([])
    dm = m//2
    # median smoother on data
    xs_list = [] # smoothed values
    for n in range(N):
        if n<dm:
            n1 = 0
            n2 = 2*n+1
        elif N-1-n<dm:
            n1 = (N-1)-2*(N-1-n)
            n2 = N
        else:
            n1 = n-dm
            n2 = n+dm+1
        nvals = list(range(n1,n2))
        if len(nvals)>1:
            nvals.remove(n)
        xl = [x_list[nv] for nv in nvals]
        xl = np.hstack(xl)
        xs = np.median(xl)
        xs_list.append(xs)
    xs_list = np.array(xs_list)
    # poly smoother on medians
    if poly_smooth:
        xs_list = poly_smooth_(xs_list,m=m)
    if post_mean:
        xs_list = mean_smooth(xs_list,m=m)
    return xs_list
#end def local_median_smooth




class TimeSeriesAnalyzer(DevBase):
    """Analyze a scalar time series sampled at uniform index intervals.

    LCD trimming can identify an initial, terminal, or interior region to
    remove before :func:`series_stats` estimates the retained mean and its
    autocorrelation-adjusted uncertainty.  The returned autocorrelation time
    is measured in sample-index units.  Because the clean region is selected
    from the observed data, the reported statistics describe that selected
    region and should be used as a diagnostic cleaning result rather than as
    a selection-free estimator.

    Parameters
    ----------
    arg0 : array_like or str, optional
        Real vector-like series, or path to a one-column text file containing
        one uniformly spaced sample per row.  If omitted, create an empty
        analyzer that can later receive data through :meth:`read`.

    clean_inp : {'simple', 'lcd_trim_l', 'lcd_trim_r', 'lcd_trim_lr', 'lcd_trim_lrm'} or mapping, optional
        Default analysis method.  A mapping must provide ``method`` and may
        provide supported options for that method, such as ``t_auto`` for
        ``'simple'`` or LCD-trim options such as ``niter``.

    label : str, optional
        Caller-supplied descriptive label retained without interpretation.

    analyze : bool, optional
        Analyze immediately when data are supplied.  If false, data are
        stored without calculating a clean partition or statistics.

    Attributes
    ----------
    x, ind : ndarray
        Original series and its uniform integer sample indices.

    xc, xl, xr, xm : ndarray or None
        Clean, left-trimmed, right-trimmed, and middle-trimmed values.

    x_mean, x_stderr, t_auto : float or None
        Mean, autocorrelation-adjusted standard error, and autocorrelation
        time for the current clean series.
    """
    def __init__(
            self,
            arg0      = None,
            clean_inp = 'lcd_trim_l',
            label     = '',
            *,
            analyze   = True,
            ):
        if not isinstance(analyze,(bool,np.bool_)):
            msg = 'analyze must be a Boolean value'
            raise TypeError(msg)
        self.filepath  = None
        self.clean_inp = clean_inp
        self.label     = label
        self.x         = None
        self.ind       = None
        self._reset()
        # process arg0
        if arg0 is None:
            self._check()
            return
        elif isinstance(arg0,str):
            self.read(arg0)
        else:
            self.x = _real_vector(arg0,'data array')
            self.ind = np.arange(len(self.x),dtype=int)
        # analyze time series
        if analyze:
            self.analyze()
        self._check()
    #end def __init_

    def _reset(self):
        """Clear all derived partition and statistical results."""
        #   results/outputs from analysis
        self.xc        = None # clean data
        self.indc      = None # indices of clean data
        self.xl        = None # left data removed
        self.indl      = None # indices of left data
        self.xr        = None # right data removed
        self.indr      = None # indices of right data
        self.xm        = None # middle data removed
        self.indm      = None # indices of middle data
        self.x_mean    = None # mean of clean data
        self.x_stderr  = None # errorbar of clean data
        self.t_auto    = None # autocorr time of clean data
    #end def _reset

    def _check(self):
        """Validate internal series, partition, and statistic consistency."""
        def check_x_ind(xk,indk):
            if self[indk] is None:
                if self[xk] is not None:
                    msg = f'{xk} requires matching {indk}'
                    raise RuntimeError(msg)
            elif self[xk] is None:
                msg = f'{indk} requires matching {xk}'
                raise RuntimeError(msg)
            else:
                if len(self[xk])==0 or len(self[indk])==0:
                    msg = f'{xk} and {indk} must not be empty'
                    raise RuntimeError(msg)
                if len(self[xk])!=len(self[indk]):
                    msg = f'{xk} and {indk} must have equal lengths'
                    raise RuntimeError(msg)
        check_x_ind('x','ind')
        check_x_ind('xc','indc')
        check_x_ind('xl','indl')
        check_x_ind('xr','indr')
        check_x_ind('xm','indm')
        if self.x_mean is not None and not np.isfinite(self.x_mean):
            msg = 'mean must be finite'
            raise RuntimeError(msg)
        if (
            self.x_stderr is not None
            and (not np.isfinite(self.x_stderr) or self.x_stderr<0.)
            ):
            msg = 'standard error must be finite and nonnegative'
            raise RuntimeError(msg)
        if (
            self.t_auto is not None
            and (not np.isfinite(self.t_auto) or self.t_auto<1.0-1e-12)
            ):
            msg = 'autocorrelation time must be finite and positive'
            raise RuntimeError(msg)
    #end def _check

    def read(self,filepath=None):
        """Load a one-dimensional uniformly sampled series from a text file.

        Parameters
        ----------
        filepath : str, optional
            File to load with :func:`numpy.loadtxt`.  By default, reload the
            previously recorded file path.

        Returns
        -------
        ndarray
            Loaded series.  Existing analysis results are cleared and integer
            sample indices are recreated.
        """
        if filepath is None:
            filepath = self.filepath
        if not isinstance(filepath,str):
            msg = 'filepath must be a string'
            raise TypeError(msg)
        x = _real_vector(np.loadtxt(filepath),'data array')
        self.filepath = filepath
        self.x   = x
        self.ind = np.arange(len(x),dtype=int)
        self._reset()
        self._check()
        return x
    #end def read

    def partition_from_timeseries(self,other):
        """Copy another analyzer's clean/removal index partition.

        Parameters
        ----------
        other : TimeSeriesAnalyzer
            Analyzer with a series of the same length.  Its clean and removed
            indices are applied to this analyzer's original values.
        """
        if not isinstance(other,TimeSeriesAnalyzer):
            msg = 'other must be a TimeSeriesAnalyzer'
            raise TypeError(msg)
        if len(other.x)!=len(self.x):
            msg = 'time series must have the same length'
            raise ValueError(msg)
        self._check()
        other._check()
        if other.indc is not None:
            self.xc   = self.x[other.indc]
            self.indc = self.ind[other.indc]
        if other.indl is not None:
            self.xl   = self.x[other.indl]
            self.indl = self.ind[other.indl]
        if other.indr is not None:
            self.xr   = self.x[other.indr]
            self.indr = self.ind[other.indr]
        if other.indm is not None:
            self.xm   = self.x[other.indm]
            self.indm = self.ind[other.indm]
        self._check()
        other._check()
    #end def partition_from_timeseries

    def clean_intersect(self,other):
        """Return values retained by both same-length analyzers.

        Parameters
        ----------
        other : TimeSeriesAnalyzer
            Analyzer defined on a series of the same length with a clean
            partition.

        Returns
        -------
        self_values, other_values, indices : ndarray
            Values from each analyzer and uniform indices at positions clean
            in both partitions.
        """
        if not isinstance(other,TimeSeriesAnalyzer):
            msg = 'other must be a TimeSeriesAnalyzer'
            raise TypeError(msg)
        if len(other.x)!=len(self.x):
            msg = 'time series must have the same length'
            raise ValueError(msg)
        self._check()
        other._check()
        count = np.zeros(len(self.x),dtype=int)
        count[self.indc]  += 1
        count[other.indc] += 1
        intersect = count==2
        ind = self.ind[intersect]
        xc1 = self.x[intersect]
        xc2 = other.x[intersect]
        self._check()
        other._check()
        return xc1,xc2,ind
    #end def clean_intersect

    def analyze(self,clean_inp=None):
        """Apply a cleaning method and calculate clean-series statistics.

        Parameters
        ----------
        clean_inp : {'simple', 'lcd_trim_l', 'lcd_trim_r', 'lcd_trim_lr', 'lcd_trim_lrm'} or mapping, optional
            Method and options for this analysis.  If omitted, use the stored
            default from construction or the preceding explicit call.  A
            mapping contains ``method`` plus options forwarded to that method.
            The analyzer manages trim return options internally.

        Returns
        -------
        None
            Results are stored in the clean/removal attributes and in
            ``x_mean``, ``x_stderr``, and ``t_auto``.
        """
        if self.x is None or self.ind is None:
            msg = 'a time series must be provided before analysis'
            raise ValueError(msg)
        self._check()
        if clean_inp is None:
            clean_inp = self.clean_inp
        else:
            self.clean_inp = clean_inp
        if clean_inp is None:
            method = 'simple'
            options = {}
        elif isinstance(clean_inp,str):
            method = clean_inp
            options = {}
        elif isinstance(clean_inp,(obj,dotdict,dict)):
            if 'method' not in clean_inp:
                msg = 'cleaning options must contain a method'
                raise ValueError(msg)
            method = clean_inp['method']
            options = dict(clean_inp.items())
            del options['method']
        else:
            msg = 'cleaning options must be a method string or mapping'
            raise TypeError(msg)
        x   = self.x
        ind = self.ind

        def calculate_stats():
            x_mean,x_stderr,t_auto = series_stats(self.xc)
            self.x_mean   = x_mean
            self.x_stderr = x_stderr
            self.t_auto   = t_auto
        #end def calculate_stats

        # no cleaning, straightforward data analysis
        if method=='simple':
            unknown = set(options)-{'t_auto'}
            if unknown:
                msg = f'unrecognized simple-analysis options: {sorted(unknown)}'
                raise ValueError(msg)
            t_auto = options.get('t_auto')
            xs   = x
            inds = ind
            x_mean,x_stderr,t_auto = series_stats(xs,t_auto=t_auto)
            self.x_mean   = x_mean
            self.x_stderr = x_stderr
            self.t_auto   = t_auto
            self.xc       = xs
            self.indc     = inds
            return
        self._reset()
        # clean the time series, then calculate stats
        #   (remove faulty data at beginning, middle, and/or end)
        if 'ret_seg' in options or 'ret_mask' in options:
            msg = 'trim return options are managed by TimeSeriesAnalyzer'
            raise ValueError(msg)
        options.update(ret_seg=False,ret_mask=True)
        # trim from left/right or both
        if method=='lcd_trim_lrm':
            mc,ml,mr,mm = lcd_trim_lrm(x,**options)
            self.xc   = self.x[mc]
            self.indc = self.ind[mc]
            if ml.sum()==0:
                self.xl   = None
                self.indl = None
            else:
                self.xl   = self.x[ml]
                self.indl = self.ind[ml]
            if mr.sum()==0:
                self.xr   = None
                self.indr = None
            else:
                self.xr   = self.x[mr]
                self.indr = self.ind[mr]
            if mm.sum()==0:
                self.xm   = None
                self.indm = None
            else:
                self.xm   = self.x[mm]
                self.indm = self.ind[mm]
            calculate_stats()
        elif method=='lcd_trim_lr':
            mc,ml,mr = lcd_trim_lr(x,**options)
            self.xc   = self.x[mc]
            self.indc = self.ind[mc]
            if ml.sum()==0:
                self.xl   = None
                self.indl = None
            else:
                self.xl   = self.x[ml]
                self.indl = self.ind[ml]
            if mr.sum()==0:
                self.xr   = None
                self.indr = None
            else:
                self.xr   = self.x[mr]
                self.indr = self.ind[mr]
            self.xm   = None
            self.indm = None
            calculate_stats()
        elif method=='lcd_trim_l':
            mc,ml = lcd_trim_l(x,**options)
            self.xc   = self.x[mc]
            self.indc = self.ind[mc]
            if ml.sum()==0:
                self.xl   = None
                self.indl = None
            else:
                self.xl   = self.x[ml]
                self.indl = self.ind[ml]
            self.xr   = None
            self.indr = None
            self.xm   = None
            self.indm = None
            calculate_stats()
        elif method=='lcd_trim_r':
            mc,mr = lcd_trim_r(x,**options)
            self.xc   = self.x[mc]
            self.indc = self.ind[mc]
            if mr.sum()==0:
                self.xr   = None
                self.indr = None
            else:
                self.xr   = self.x[mr]
                self.indr = self.ind[mr]
            self.xl   = None
            self.indl = None
            self.xm   = None
            self.indm = None
            calculate_stats()
        else:
            msg = f'unrecognized data cleaning method "{method}"'
            raise ValueError(msg)
        self._check()
    #end def analyze

    def plot(
            self,
            *,
            fig    = False,
            show   = False,
            ishift = 0,
            legend = True,
            ):
        """Plot the trace, cleaning partition, and clean-data reference lines.

        The full series is gray, clean segments are black, left and right
        trims are red and blue, and middle trims are magenta.  Green lines
        show the clean mean and one clean-series standard deviation on either
        side.  The existing pyplot axes are used unless ``fig`` is true.

        Parameters
        ----------
        fig : bool, optional
            Create a new tight-layout matplotlib figure before plotting.

        show : bool, optional
            Display the current figure after plotting.

        ishift : float, optional
            Add this offset to each integer sample index on the horizontal
            axis.

        legend : bool, optional
            Add a legend describing the plotted series and reference lines.

        Returns
        -------
        None
            Lines are added to the current matplotlib axes.
        """
        import matplotlib.pyplot as plt
        self._check()
        if self.xc is None or self.x_mean is None or self.x_stderr is None:
            msg = 'analysis must be completed before plotting'
            raise ValueError(msg)
        if fig:
            plt.figure(tight_layout=True)

        plt.plot(self.ind+ishift,self.x,color='Grey',label='full series')

        imin = self.ind[0]+ishift
        imax = self.ind[-1]+ishift
        plt.plot([imin,imax],2*[self.x_mean],'g-',label='clean mean')
        plt.plot(
            [imin,imax],
            2*[self.x_mean+np.std(self.xc)],
            'g-.',
            label='clean mean ± std. dev.',
            )
        plt.plot([imin,imax],2*[self.x_mean-np.std(self.xc)],'g-.')

        mask = np.zeros(len(self.x),dtype=bool)
        mask[self.indc]=True
        segs = _find_segments(self.x,mask,1)
        for iseg,(n1,n2) in enumerate(segs):
            label = 'clean series' if iseg==0 else None
            plt.plot(self.ind[n1:n2]+ishift,self.x[n1:n2],'k-',label=label)

        if self.xl is not None:
            plt.plot(self.indl+ishift,self.xl,'r-',label='left trim')

        if self.xr is not None:
            plt.plot(self.indr+ishift,self.xr,'b-',label='right trim')

        if self.xm is not None:
            mask = np.zeros(len(self.x),dtype=bool)
            mask[self.indm]=True
            segs = _find_segments(self.x,mask,1)
            for iseg,(n1,n2) in enumerate(segs):
                label = 'middle trim' if iseg==0 else None
                plt.plot(self.ind[n1:n2]+ishift,self.x[n1:n2],'m-',label=label)

        if legend:
            plt.legend()
        plt.xlabel('time index')
        ylabel = self.label if self.label is not None else 'time series'
        plt.ylabel(ylabel)
        if show:
            plt.show()
        self._check()
    #end def plot
#end class TimeSeriesAnalyzer
