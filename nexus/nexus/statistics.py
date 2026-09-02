
import numpy as np


def theil_sen(x,y):
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

    The exact estimator is used through approximately 700,000 pairwise
    slopes.  Above that crossover, the number of sampled slopes grows
    linearly with the number of input points, with a rate calibrated to
    100,000 samples at 448 points.  Pairwise slopes are sampled uniformly with
    replacement.  A fixed local random seed makes the stochastic estimate
    reproducible without altering NumPy's global random state.
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
    exact_pair_limit = 700_000
    if npairs<=exact_pair_limit:
        return theil_sen(x,y)

    sample_pair_reference = 100_000
    sample_n_reference = 448
    nsampled = (
        sample_pair_reference*n+sample_n_reference//2)//sample_n_reference
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


def acf_autocorr_time(x):
    """Estimate IAT from a noise-truncated, flat-top-windowed ACF."""
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
    if len(x)==1:
        return 1.

    x = x-x.mean()
    variance = np.mean(x**2)
    if variance==0.:
        return 1.

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
        import warnings
        quiet_start = max_lag
        warnings.warn(
            'the ACF did not reach a noise-dominated region; '
            'the autocorrelation-time estimate may be unreliable',
            RuntimeWarning,
            stacklevel=2)

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
    # length 1 "reblocking"
    data_errs1 = x.std()/np.sqrt(nblocks)
    if data_errs1==0:
        return 1.

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
        plt.title('t_auto = {}'.format(t_auto))
        if show:
            plt.show()
    return t_auto
#end def reblocked_autocorr_time


def reblocked_autocorr_time2(
    x,confidence=0.99,plot=False,show=False,min_blocks=16,block_scale=5.):
    """Estimate the integrated autocorrelation time by data blocking.

    This implementation follows the blocking transformation and automatic
    stopping test of Flyvbjerg and Petersen.  To avoid selecting a level before
    the blocking curve has equilibrated, the selected block length must also
    be at least ``block_scale`` times its IAT estimate while retaining at least
    ``min_blocks`` blocks.  The returned value is the variance-inflation factor
    ``N*Var(mean)/Var(data)``, consistent with an effective sample size of
    ``N/t_auto``.
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
    try:
        confidence = float(confidence)
        block_scale = float(block_scale)
    except (TypeError,ValueError):
        raise ValueError(
            'confidence and block scale must be finite numbers') from None
    if not np.isfinite(confidence) or not 0.<confidence<1.:
        raise ValueError('confidence must be between zero and one')
    if not isinstance(min_blocks,(int,np.integer)) or min_blocks<2:
        raise ValueError('minimum number of blocks must be an integer >= 2')
    if not np.isfinite(block_scale) or block_scale<=0.:
        raise ValueError('block scale must be a positive finite number')
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
    warning_reasons = []
    if len(passing_levels)>0:
        test_level = int(passing_levels[0])
    else:
        test_level = nlevel-1
        warning_reasons.append('the blocking test did not reach an uncorrelated level')

    tau_levels = block_sizes*variances/variances[0]
    reliable_levels = np.flatnonzero(block_counts>=min_blocks)
    if len(reliable_levels)>0:
        max_reliable_level = int(reliable_levels[-1])
    else:
        max_reliable_level = 0
        warning_reasons.append(
            f'the series contains fewer than {min_blocks} samples')
    if test_level>max_reliable_level:
        warning_reasons.append(
            f'the blocking test requires fewer than {min_blocks} blocks')
    test_level = min(test_level,max_reliable_level)

    if test_level==0:
        selected_level = 0
    else:
        levels = np.arange(nlevel)
        equilibrated = np.flatnonzero(
            (levels>=test_level)
            &(levels<=max_reliable_level)
            &(block_sizes>=block_scale*tau_levels))
        if len(equilibrated)>0:
            selected_level = int(equilibrated[0])
        else:
            selected_level = max_reliable_level
            warning_reasons.append(
                'no reliable block length is sufficiently long relative to the IAT')

    t_auto = max(float(tau_levels[selected_level]),np.finfo(float).eps)

    if warning_reasons:
        import warnings
        warnings.warn(
            '; '.join(warning_reasons)+'; the estimate may be unreliable',
            RuntimeWarning,
            stacklevel=2)

    if plot:
        import matplotlib.pyplot as plt
        tau_errors = tau_levels*np.sqrt(2./(block_counts-1))
        plt.figure(tight_layout=True)
        plt.errorbar(block_sizes,tau_levels,yerr=tau_errors,fmt='b.-')
        plt.axvline(block_sizes[test_level],color='grey',linestyle=':')
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
            value        = (
                values[-2]*weights[-2]+values[-1]*weights[-1]
                )/weight
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


def autocorr_time(x, disagreement_threshold=1.5):
    """Conservatively combine three autocorrelation-time estimates.

    The ACF, Geyer initial-monotone-sequence, and legacy reblocking
    estimators probe the correlation structure in different ways.  Since an
    underestimated autocorrelation time leads directly to an underestimated
    uncertainty, this function returns the largest of their estimates.  A
    warning is emitted when the largest and smallest estimates differ by more
    than ``disagreement_threshold``; such disagreement can indicate that the
    series is too short to resolve a slow correlation mode.

    Parameters
    ----------
    x : array_like
        One-dimensional sample sequence.

    disagreement_threshold : float, optional
        Warn when ``max(estimates)/min(estimates)`` exceeds this value.  Must
        be a finite number greater than one.

    Returns
    -------
    tau : float
        Maximum of ``acf_autocorr_time(x)``,
        ``integrated_autocorr_time(x)``, and
        ``reblocked_autocorr_time(x)``.
    """

    try:
        disagreement_threshold = float(disagreement_threshold)
    except (TypeError,ValueError):
        raise ValueError(
            'disagreement threshold must be a finite number greater than one'
            ) from None
    if (
        not np.isfinite(disagreement_threshold)
        or disagreement_threshold<=1.
    ):
        raise ValueError(
            'disagreement threshold must be a finite number greater than one')

    x = np.asarray(x)
    names = ('ACF','integrated','reblocking')
    estimates = np.array([
        acf_autocorr_time(x),
        integrated_autocorr_time(x),
        reblocked_autocorr_time(x),
        ],dtype=float)
    if not np.all(np.isfinite(estimates)):
        raise ValueError('component autocorrelation-time estimate is nonfinite')

    estimate_min = float(estimates.min())
    estimate_max = float(estimates.max())
    if estimate_min>0.:
        disagreement = estimate_max/estimate_min
    elif estimate_max>0.:
        disagreement = np.inf
    else:
        disagreement = 1.

    if disagreement>disagreement_threshold:
        import warnings
        detail = ', '.join(
            f'{name}={estimate:.6g}'
            for name,estimate in zip(names,estimates))
        warnings.warn(
            'autocorrelation-time estimators disagree by a factor of '
            f'{disagreement:.3g} ({detail}); the estimate may be unreliable',
            RuntimeWarning,
            stacklevel=2)

    return estimate_max
#end def autocorr_time



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
