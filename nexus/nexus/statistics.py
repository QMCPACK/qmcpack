
import numpy as np


def theil_sen(x,y):
    """Return the Theil--Sen slope and intercept for paired observations."""
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
    """Estimate autocorrelation time from a windowed sample ACF.

    The autocorrelation function is evaluated in ``O(N log N)`` time with an
    FFT and a common denominator at all lags.  A Bartlett noise estimate is
    used to locate the first sustained noise-dominated region.  Resolved lags
    are retained by a flat-top window, followed by a linear taper through the
    noisy boundary.  The returned value is the variance-inflation factor
    ``N*Var(mean)/Var(x)``.

    Parameters
    ----------
    x : array_like
        One-dimensional, finite, real-valued sample sequence.  Vector-shaped
        two-dimensional arrays are flattened.

    Returns
    -------
    tau : float
        Estimated integrated autocorrelation time.  A value of one denotes
        IID-like sampling; negative correlation can produce a value below
        one.

    Warns
    -----
    RuntimeWarning
        If the ACF does not reach a noise-dominated region within the inspected
        lag range.
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
    """Estimate autocorrelation time from the growth of blocked errors.

    For every integer block length that leaves at least ``min_blocks`` blocks,
    contiguous block means and their standard error are computed.  Their
    ratio to the unblocked standard error is fitted as a function of block
    length with a Theil--Sen line (using reproducible stochastic pair sampling
    for large inputs).  The squared fitted ratio at the largest usable block
    length is used to obtain the auto-correlation time.

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


def reblocked_autocorr_time_fp(
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
#end def reblocked_autocorr_time_fp


def geyer_ims_autocorr_time(x, c=5.0):
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
#end def geyer_ims_autocorr_time


def prewhitened_spectral_autocorr_time(
    x,max_order=10,stability_threshold=2.):
    """Estimate autocorrelation time with a prewhitened QS spectrum.

    A low-order autoregressive model is used only as a prewhitening filter;
    its order is selected by an AICc comparison of winsorized Yule-Walker
    pilots.  The selected order is refitted by least squares to preserve the
    second-moment target.  The zero-frequency spectrum of the filtered series
    is estimated with a quadratic-spectral (QS) lag window and Andrews'
    automatic plug-in bandwidth, then recolored through the fitted filter.

    Prewhitening is rejected when the zero-frequency filter gain is not
    reliably separated from zero.  A geometric bandwidth diagnostic also
    expands beyond the plug-in value when the residual spectrum contains a
    resolved slow mode.  Dividing the resulting long-run variance by the
    marginal variance gives the integrated autocorrelation time.

    Unlike estimators designed specifically to guard against underestimated
    error bars, this estimator does not impose ``tau >= 1`` or select the
    largest of several estimates.  It can therefore represent variance
    reduction due to negative or oscillatory correlation.

    Parameters
    ----------
    x : array_like
        One-dimensional, real-valued sample sequence.

    max_order : int, optional
        Largest autoregressive pilot order to consider.  The effective upper
        bound is also limited by the series length.

    stability_threshold : float, optional
        Warn when positive long-run-variance estimates obtained at half and
        twice the selected bandwidth span a factor larger than this value.
        Must be a finite number greater than one.

    Returns
    -------
    tau : float
        Estimated long-run variance divided by the marginal variance.
    """

    import warnings

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
    if (
        isinstance(max_order,(bool,np.bool_))
        or not isinstance(max_order,(int,np.integer))
        or max_order<0
    ):
        raise ValueError('maximum order must be a nonnegative integer')
    try:
        stability_threshold = float(stability_threshold)
    except (TypeError,ValueError):
        raise ValueError(
            'stability threshold must be a finite number greater than one'
            ) from None
    if not np.isfinite(stability_threshold) or stability_threshold<=1.:
        raise ValueError(
            'stability threshold must be a finite number greater than one')
    if len(x)<2:
        return 1.

    centered = x-x.mean()
    marginal_variance = np.mean(centered**2)
    if marginal_variance==0.:
        return 1.

    tiny = np.finfo(float).eps
    normal_mad = 0.6744897501960817

    def robust_scale(values):
        """MAD scale with an RMS fallback for discrete/small samples."""
        residual = values-np.median(values)
        scale = np.median(np.abs(residual))/normal_mad
        if not np.isfinite(scale) or scale<=tiny:
            scale = np.sqrt(np.mean(residual**2))
        if not np.isfinite(scale) or scale<=tiny:
            scale = np.sqrt(np.mean(values**2))
        return max(float(scale),tiny)

    def stable_coefficients(coefficients):
        """Move fitted AR roots just inside the stationary region."""
        if len(coefficients)==0:
            return coefficients,False
        roots = np.roots(np.r_[1.,-coefficients])
        root_limit = .995
        root_sizes = np.abs(roots)
        clipped = root_sizes>=root_limit
        if np.any(clipped):
            roots[clipped] *= root_limit/root_sizes[clipped]
        polynomial = np.real_if_close(np.poly(roots),tol=1000).real
        return -polynomial[1:],bool(np.any(clipped))

    # Estimate all candidate orders from one winsorized autocovariance
    # sequence.  Levinson recursion makes this O(n*max_order+max_order**2),
    # instead of performing iterative robust regressions at every order.
    pilot_location = np.median(x)
    pilot_scale = robust_scale(x)
    pilot = np.clip((x-pilot_location)/pilot_scale,-6.,6.)
    pilot -= pilot.mean()
    n = len(pilot)
    order_limit = min(
        int(max_order),int(np.floor(n**(1./3.))),max(0,(n-5)//3))
    pilot_autocovariance = np.array([
        np.dot(pilot[:n-lag],pilot[lag:])/n
        for lag in range(order_limit+1)])
    variance_floor = tiny*max(float(pilot_autocovariance[0]),1.)

    def aicc_score(order,prediction_variance):
        nscore = n-order
        nparameter = order+1
        score = nscore*np.log(max(prediction_variance,variance_floor))
        score += 2.*nparameter
        if nscore>nparameter+1:
            score += (
                2.*nparameter*(nparameter+1)
                /(nscore-nparameter-1))
        else:
            score = np.inf
        return score

    coefficients = np.empty(0,dtype=float)
    prediction_variance = max(
        float(pilot_autocovariance[0]),variance_floor)
    candidate_coefficients = [coefficients]
    candidate_clipped = [False]
    candidate_scores = [aicc_score(0,prediction_variance)]
    recursion_clipped = False
    for order in range(1,order_limit+1):
        numerator = pilot_autocovariance[order]
        if order>1:
            numerator -= np.dot(
                coefficients,pilot_autocovariance[1:order][::-1])
        reflection = numerator/prediction_variance
        if abs(reflection)>=.995:
            reflection = np.copysign(.995,reflection)
            recursion_clipped = True
        updated = np.empty(order,dtype=float)
        if order>1:
            updated[:-1] = coefficients-reflection*coefficients[::-1]
        updated[-1] = reflection
        coefficients = updated
        prediction_variance = max(
            prediction_variance*(1.-reflection**2),variance_floor)
        candidate_coefficients.append(coefficients.copy())
        candidate_clipped.append(recursion_clipped)
        candidate_scores.append(aicc_score(order,prediction_variance))

    candidate_scores = np.asarray(candidate_scores)
    best_score = np.min(candidate_scores)
    # Select the AICc minimum; ties are resolved in favor of the lower order.
    selected_order = int(np.flatnonzero(candidate_scores==best_score)[0])
    coefficients = candidate_coefficients[selected_order]
    warning_reasons = []
    coefficients_clipped = candidate_clipped[selected_order]
    selected_condition = 1.
    filter_gain_se = 0.
    if selected_order>0:
        # Refit only the selected order on the original observations.  This
        # preserves the second-moment target for nonlinear Markov transitions.
        response = centered[selected_order:]
        predictors = np.column_stack([
            centered[selected_order-lag:n-lag]
            for lag in range(1,selected_order+1)])
        coefficients = np.linalg.lstsq(
            predictors,response,rcond=None)[0]
        selected_condition = np.linalg.cond(predictors)
        fit_residual = response-predictors@coefficients
        residual_degrees_of_freedom = max(
            len(response)-selected_order,1)
        innovation_variance = (
            np.dot(fit_residual,fit_residual)
            /residual_degrees_of_freedom)
        coefficient_covariance = innovation_variance*np.linalg.pinv(
            predictors.T@predictors)
        sum_direction = np.ones(selected_order,dtype=float)
        filter_gain_variance = sum_direction@coefficient_covariance@sum_direction
        filter_gain_se = np.sqrt(max(float(filter_gain_variance),0.))
        coefficients,refit_clipped = stable_coefficients(coefficients)
        coefficients_clipped = coefficients_clipped or refit_clipped
    if coefficients_clipped:
        warning_reasons.append(
            'the autoregressive pilot required a stationarity correction')
    if selected_condition>1.e10:
        warning_reasons.append(
            'the autoregressive pilot is ill-conditioned')

    # Recoloring divides by A(1)**2, where A(1)=1-sum(phi).  Require that this
    # gain be separated from zero by three estimated standard errors.  When it
    # is not, use an unprewhitened spectral estimate rather than amplify an
    # uncertain fitted quantity.
    filter_gain = 1.-np.sum(coefficients)
    gain_unresolved = (
        selected_order>0
        and (
            not np.isfinite(filter_gain)
            or abs(filter_gain)<=3.*filter_gain_se
            or abs(filter_gain)<1.e-6
            or coefficients_clipped
            or selected_condition>1.e10
        )
    )
    if gain_unresolved:
        warning_reasons.append(
            'the autoregressive recoloring gain is not reliably resolved; '
            'using an unprewhitened spectral estimate')
        coefficients = np.empty(0,dtype=float)
        selected_order = 0
        filter_gain = 1.
    elif abs(filter_gain)<.05:
        warning_reasons.append(
            'the recoloring factor is sensitive to the autoregressive pilot')

    # Apply the accepted filter to the original, mean-centered observations.
    residual = centered[selected_order:].copy()
    for lag,coefficient in enumerate(coefficients,start=1):
        residual -= coefficient*centered[
            selected_order-lag:n-lag]
    residual -= residual.mean()
    residual_variance = np.mean(residual**2)
    if residual_variance<=tiny*marginal_variance:
        warning_reasons.append(
            'the prewhitened residual variance is numerically zero')
        residual_variance = max(residual_variance,tiny*marginal_variance)

    # Biased (common-denominator) autocovariances yield a stable finite-sample
    # sequence.  The QS window is evaluated over all available residual lags.
    nresidual = len(residual)
    nfft = 1 << (2*nresidual-1).bit_length()
    transform = np.fft.rfft(residual,nfft)
    autocovariance = np.fft.irfft(
        transform*np.conjugate(transform),nfft)[:nresidual]/nresidual

    # Andrews' plug-in rule for the QS bandwidth.  For an AR(p) pilot, the
    # required normalized spectral curvature follows analytically from its
    # zero-frequency filter gain.  The p=0 fallback uses an AR(1) curvature
    # estimate directly from the unfiltered observations.
    if len(coefficients)>0:
        lag_numbers = np.arange(1,len(coefficients)+1,dtype=float)
        first_moment = np.dot(lag_numbers,coefficients)
        second_moment = np.dot(lag_numbers**2,coefficients)
        spectral_curvature = -2.*(
            filter_gain*second_moment+first_moment**2
            )/filter_gain**2
        alpha = spectral_curvature**2
    else:
        if nresidual>1:
            denominator = np.dot(residual[:-1],residual[:-1])
            if denominator>0.:
                rho = np.dot(residual[:-1],residual[1:])/denominator
            else:
                rho = 0.
        else:
            rho = 0.
        rho = float(np.clip(rho,-.99,.99))
        alpha = 4.*rho**2/(1.-rho)**4
    bandwidth = 1.3221*(alpha*nresidual)**.2 if alpha>0. else 1.
    maximum_bandwidth = max(1.,nresidual/4.)
    if bandwidth>maximum_bandwidth:
        bandwidth = maximum_bandwidth
        warning_reasons.append(
            'the automatic spectral bandwidth reached its sample-size limit')
    bandwidth = max(1.,float(bandwidth))

    def quadratic_spectral_lrv(selected_bandwidth):
        if nresidual==1:
            return float(autocovariance[0])
        lags = np.arange(1,nresidual,dtype=float)
        ratio = lags/selected_bandwidth
        argument = 6.*np.pi*ratio/5.
        weights = 25./(12.*np.pi**2*ratio**2)*(
            np.sin(argument)/argument-np.cos(argument))
        return float(
            autocovariance[0]
            +2.*np.dot(weights,autocovariance[1:]))

    # Search geometrically beyond the plug-in bandwidth.  A low-amplitude
    # slow mode can be nearly irrelevant to one-step AR prediction while still
    # dominating the spectrum at zero.  Expand only when the estimated
    # spectrum grows by more than its bandwidth-dependent noise allowance,
    # and stop at the first plateau or downturn.  Retaining at least 16
    # bandwidth-widths in the record limits the variance of this diagnostic.
    diagnostic_maximum = max(1.,nresidual/16.)
    bandwidth_grid = [bandwidth]
    while bandwidth_grid[-1]<diagnostic_maximum:
        next_bandwidth = min(2.*bandwidth_grid[-1],diagnostic_maximum)
        if next_bandwidth<=bandwidth_grid[-1]*(1.+tiny):
            break
        bandwidth_grid.append(next_bandwidth)
    grid_lrvs = np.array([
        quadratic_spectral_lrv(value)
        for value in bandwidth_grid])

    selected_index = 0
    peak_index = 0
    low_frequency_growth = False
    for index in range(1,len(bandwidth_grid)):
        previous = grid_lrvs[index-1]
        current = grid_lrvs[index]
        if not np.isfinite(current) or current<=0. or previous<=0.:
            break
        noise_allowance = min(
            .35,max(.12,1.5*np.sqrt(bandwidth_grid[index]/nresidual)))
        if not low_frequency_growth:
            if current>grid_lrvs[0]*(1.+noise_allowance):
                low_frequency_growth = True
                selected_index = index
                peak_index = index
        else:
            if current>grid_lrvs[peak_index]:
                peak_index = index
            if abs(current/previous-1.)<=noise_allowance:
                selected_index = index
                break
            if current<previous*(1.-noise_allowance):
                selected_index = peak_index
                break
            selected_index = index

    if low_frequency_growth:
        bandwidth = bandwidth_grid[selected_index]
        residual_lrv = grid_lrvs[selected_index]
        warning_reasons.append(
            'residual low-frequency growth required a larger spectral bandwidth')
    else:
        residual_lrv = grid_lrvs[0]

    comparison_bandwidths = np.unique(np.clip(
        [bandwidth/2.,bandwidth,2.*bandwidth],1.,maximum_bandwidth))
    comparison_lrvs = np.array([
        quadratic_spectral_lrv(value)
        for value in comparison_bandwidths])
    positive_lrvs = comparison_lrvs[comparison_lrvs>0.]
    if len(positive_lrvs)!=len(comparison_lrvs):
        warning_reasons.append(
            'the spectral estimate is not positive at nearby bandwidths')
    elif (
        len(positive_lrvs)>1
        and positive_lrvs.max()/positive_lrvs.min()>stability_threshold
    ):
        warning_reasons.append(
            'the spectral estimate is sensitive to the bandwidth')

    if not np.isfinite(residual_lrv) or residual_lrv<=0.:
        warning_reasons.append(
            'the estimated residual spectrum is not numerically positive')
        residual_lrv = tiny*residual_variance
    long_run_variance = residual_lrv/filter_gain**2
    tau = max(float(long_run_variance/marginal_variance),tiny)
    if n<20.*tau:
        warning_reasons.append(
            'the series contains fewer than 20 estimated correlation times')

    if warning_reasons:
        warnings.warn(
            '; '.join(dict.fromkeys(warning_reasons))
            +'; the estimate may be unreliable',
            RuntimeWarning,
            stacklevel=2)

    return tau
#end def prewhitened_spectral_autocorr_time


def autocorr_time(x):
    """Conservatively combine three autocorrelation-time estimates.

    The ACF, Geyer initial-monotone-sequence, and reblocking estimators 
    probe the correlation structure in different ways.  Since an
    underestimated autocorrelation time leads directly to an underestimated
    uncertainty, this function returns the largest of their estimates. 

    Parameters
    ----------
    x : array_like
        One-dimensional sample sequence.

    Returns
    -------
    tau : float
        Maximum of ACF, Geyer, and reblocked auto-correlation times.
    """

    x = np.asarray(x)
    t_auto_acf     = acf_autocorr_time(x)
    t_auto_reblock = reblocked_autocorr_time(x)
    t_auto_geyer   = geyer_ims_autocorr_time(x)

    t_auto = max(t_auto_acf,t_auto_reblock,t_auto_geyer)

    return t_auto
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
