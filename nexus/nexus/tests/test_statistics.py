import numpy as np
import pytest

from .. import statistics
from . import NexusTestOrder

pytestmark = pytest.mark.order(NexusTestOrder.STATISTICS)



def test_theil_sen():
    """Verify exact robust fitting with an extreme outlier.

    Also check vector flattening and mismatched-length rejection.
    """
    x = np.arange(5,dtype=float)
    y = np.array([1.,3.,5.,7.,101.])

    slope,intercept = statistics.theil_sen(x,y)

    assert(slope==pytest.approx(2.))
    assert(intercept==pytest.approx(1.))

    slope_column,intercept_column = statistics.theil_sen(
        x.reshape(-1,1),
        y.reshape(-1,1),
        )
    assert(slope_column==pytest.approx(slope))
    assert(intercept_column==pytest.approx(intercept))

    x_duplicate = np.array([0.,0.,1.,2.])
    y_duplicate = 2.*x_duplicate+1.
    slope_duplicate,intercept_duplicate = statistics.theil_sen(
        x_duplicate,
        y_duplicate,
        )
    assert(slope_duplicate==pytest.approx(2.))
    assert(intercept_duplicate==pytest.approx(1.))

    for x_degenerate in (np.array([]),np.array([1.]),np.ones(3)):
        with pytest.raises(
            ValueError,
            match=r'at least two distinct x values are required',
            ):
            statistics.theil_sen(x_degenerate,x_degenerate)

    with pytest.raises(ValueError,match=r'input arrays must have equal lengths'):
        statistics.theil_sen(x,y[:-1])

    invalid_cases = [
        (1.,1.,r'input arrays must be one-dimensional'),
        (np.ones((2,2)),np.ones((2,2)),r'input arrays must be one-dimensional'),
        (x,np.full(len(x),np.nan),r'input arrays must contain only finite values'),
        (x.astype(complex),y,r'input arrays must be real-valued'),
        ]
    for x_invalid,y_invalid,message in invalid_cases:
        with pytest.raises(ValueError,match=message):
            statistics.theil_sen(x_invalid,y_invalid)
#end def test_theil_sen



def test_theil_sen_stochastic_exact_path(monkeypatch):
    """Verify both stochastic estimators delegate below crossover.

    The delegated inputs and returned values should remain unchanged.
    """
    x        = np.arange(10,dtype=float)
    y        = 1.5*x-4.
    expected = (7.,-3.)
    calls    = []

    def exact_theil_sen(x_arg,y_arg):
        calls.append((x_arg.copy(),y_arg.copy()))
        return expected
    #end def exact_theil_sen

    monkeypatch.setattr(statistics,'theil_sen',exact_theil_sen)

    assert(statistics.theil_sen_stoch(x,y)==expected)
    assert(statistics.theil_sen_stoch_reblock(x,y)==expected)
    assert(len(calls)==2)
    for x_arg,y_arg in calls:
        np.testing.assert_array_equal(x_arg,x)
        np.testing.assert_array_equal(y_arg,y)
#end def test_theil_sen_stochastic_exact_path



def test_theil_sen_stochastic_sampled_path(monkeypatch):
    """Exercise the sampled paths without calling the exact estimator.

    Check linear fits and reproducibility from the fixed random seed.
    """
    def reject_exact_path(x,y):
        pytest.fail('sampled path unexpectedly called the exact estimator')
    #end def reject_exact_path

    monkeypatch.setattr(statistics,'theil_sen',reject_exact_path)

    x = np.arange(1024,dtype=float)
    x[1] = x[0]
    y = .75*x-2.
    slope,intercept = statistics.theil_sen_stoch(x,y)
    assert(slope==pytest.approx(.75))
    assert(intercept==pytest.approx(-2.))

    x_reblock = np.arange(16,dtype=float)
    y_reblock = -.5*x_reblock+3.
    slope,intercept = statistics.theil_sen_stoch_reblock(
        x_reblock,
        y_reblock,
        )
    assert(slope==pytest.approx(-.5))
    assert(intercept==pytest.approx(3.))

    y_noisy = .2*x+np.sin(x/7.)
    first   = statistics.theil_sen_stoch(x,y_noisy)
    second  = statistics.theil_sen_stoch(x,y_noisy)
    assert(first==second)
#end def test_theil_sen_stochastic_sampled_path



def test_reblocked_autocorr_time():
    """Check reblocking for trivial, IID, and correlated series.

    Correlation should increase the estimate for flat and column input.
    """
    singleton = np.array([3.])
    constant  = np.ones(32)
    assert(statistics.reblocked_autocorr_time(singleton)==1.)
    assert(statistics.reblocked_autocorr_time(constant)==1.)

    rng = np.random.default_rng(90210)
    iid = rng.normal(size=256)
    correlated = iid.copy()
    for index in range(1,len(correlated)):
        correlated[index] += .8*correlated[index-1]

    tau_iid        = statistics.reblocked_autocorr_time(iid)
    tau_correlated = statistics.reblocked_autocorr_time(correlated)
    tau_column = statistics.reblocked_autocorr_time(correlated.reshape(-1,1))
    tau_list = statistics.reblocked_autocorr_time(correlated.tolist())

    assert(np.isfinite(tau_iid))
    assert(tau_iid>0.)
    assert(tau_correlated>3.*tau_iid)
    assert(tau_column==pytest.approx(tau_correlated))
    assert(tau_list==pytest.approx(tau_correlated))
#end def test_reblocked_autocorr_time



def test_reblocked_autocorr_time_invalid_input():
    """Check reblocking validation for invalid shapes and limits.

    Each diagnostic-bearing assertion must report its expected message.
    """
    test_cases = [
        (np.array([]),10,r'data array must not be empty'),
        (np.arange(8),0,r'minimum number of blocks must be a positive integer'),
        (np.arange(8),1.5,r'minimum number of blocks must be a positive integer'),
        (np.arange(8),np.inf,r'minimum number of blocks must be a positive integer'),
        (np.array([1.,np.nan]),10,r'data array must contain only finite values'),
        (np.array([1.+1.j]),10,r'data array must be real-valued'),
        ]
    for x,min_blocks,message in test_cases:
        with pytest.raises(ValueError,match=message):
            statistics.reblocked_autocorr_time(x,min_blocks=min_blocks)

    with pytest.raises(ValueError,match=r'data array must be 1-dimensional'):
        statistics.reblocked_autocorr_time(np.ones((2,2)))
#end def test_reblocked_autocorr_time_invalid_input



def test_acf_autocorr_time():
    """Check ACF estimates for trivial, IID, and correlated series.

    Also verify reliability results and vector-shaped input handling.
    """
    assert(statistics.acf_autocorr_time(np.array([2.]))==1.)
    assert(statistics.acf_autocorr_time(np.ones(32))==1.)
    assert(
        statistics.acf_autocorr_time(np.ones(32),reliability=True)
        ==(1.,False)
        )

    rng        = np.random.default_rng(90210)
    iid        = rng.normal(size=256)
    correlated = iid.copy()
    for index in range(1,len(correlated)):
        correlated[index] += .8*correlated[index-1]

    tau_iid,unreliable_iid = statistics.acf_autocorr_time(
        iid,
        reliability=True,
        )
    tau_correlated,unreliable_correlated = statistics.acf_autocorr_time(
        correlated,
        reliability=True,
        )
    tau_column = statistics.acf_autocorr_time(correlated.reshape(-1,1))

    assert(.5<tau_iid<2.)
    assert(tau_correlated>3.*tau_iid)
    assert(tau_column==pytest.approx(tau_correlated))
    assert(not unreliable_iid)
    assert(not unreliable_correlated)
#end def test_acf_autocorr_time



def test_acf_autocorr_time_invalid_input():
    """Verify ACF validation rejects unsupported input arrays.

    Error messages should identify the violated input requirement.
    """
    test_cases = [
        (np.array([1.+1.j]),r'data array must be real-valued'),
        (np.ones((2,2)),r'data array must be 1-dimensional'),
        (np.array([]),r'data array must not be empty'),
        (np.array([1.,np.nan]),r'data array must contain only finite values'),
        ]
    for x,message in test_cases:
        with pytest.raises(ValueError,match=message):
            statistics.acf_autocorr_time(x)
#end def test_acf_autocorr_time_invalid_input



def test_geyer_ims_autocorr_time():
    """Check Geyer IMS edge cases, fallback, and reliability.

    Alternating data exercises the pure estimate and ACF fallback.
    """
    assert(statistics.geyer_ims_autocorr_time(np.array([2.]))==1.)
    assert(statistics.geyer_ims_autocorr_time(np.ones(32))==1.)

    alternating = np.tile([-1.,1.],64)
    tau_acf,reliability_acf = statistics.acf_autocorr_time(
        alternating,
        reliability=True,
        )
    tau_fallback,reliability_fallback = statistics.geyer_ims_autocorr_time(
        alternating,
        reliability=True,
        )
    tau_pure = statistics.geyer_ims_autocorr_time(
        alternating,
        acf_fallback=False,
        )

    assert(tau_fallback==pytest.approx(tau_acf))
    assert(reliability_fallback==reliability_acf)
    assert(tau_pure==0.)

    trend = np.arange(64,dtype=float)
    tau_trend,unreliable_trend = statistics.geyer_ims_autocorr_time(
        trend,
        reliability=True,
        )
    tau_column = statistics.geyer_ims_autocorr_time(trend.reshape(-1,1))
    assert(tau_trend>1.)
    assert(tau_column==pytest.approx(tau_trend))
    assert(unreliable_trend)
#end def test_geyer_ims_autocorr_time



def test_geyer_ims_autocorr_time_invalid_input():
    """Verify Geyer IMS validation for options and input arrays.

    Each invalid value should produce its documented diagnostic.
    """
    invalid_c_values = [0.,-1.,np.inf,'invalid']
    for c in invalid_c_values:
        with pytest.raises(
            ValueError,
            match=r'c must be a positive finite number',
            ):
            statistics.geyer_ims_autocorr_time(np.arange(8),c=c)

    test_cases = [
        (np.array([1.+1.j]),r'input must be real-valued'),
        (np.ones((2,2)),r'input must be one-dimensional'),
        (np.array([]),r'input must not be empty'),
        (np.array([1.,np.nan]),r'input must contain only finite values'),
        ]
    for x,message in test_cases:
        with pytest.raises(ValueError,match=message):
            statistics.geyer_ims_autocorr_time(x)
#end def test_geyer_ims_autocorr_time_invalid_input



def test_autocorr_time(monkeypatch):
    """Verify combined estimation selects the conservative maximum.

    Reliability requests and component flags must also be propagated.
    """
    calls = []

    def fake_acf(x,reliability=False):
        calls.append(('acf',x.copy(),reliability))
        return 2.,False
    #end def fake_acf

    def fake_geyer(x,reliability=False):
        calls.append(('geyer',x.copy(),reliability))
        return 3.,True
    #end def fake_geyer

    monkeypatch.setattr(statistics,'acf_autocorr_time',fake_acf)
    monkeypatch.setattr(statistics,'geyer_ims_autocorr_time',fake_geyer)

    x = [1.,2.,3.]
    assert(statistics.autocorr_time(x)==3.)
    assert(statistics.autocorr_time(x,reliability=True)==(3.,True))
    assert(len(calls)==4)
    for name,x_arg,reliability in calls:
        assert(name in {'acf','geyer'})
        np.testing.assert_array_equal(x_arg,x)
        assert(reliability)
#end def test_autocorr_time



def test_series_stats(monkeypatch):
    """Check series statistics with supplied and estimated timing.

    The standard error must use the autocorrelation-adjusted sample size.
    """
    x = np.array([1.,2.,3.,4.])
    mean,error,tau = statistics.series_stats(x,t_auto=4.)
    assert(mean==pytest.approx(np.mean(x)))
    assert(error==pytest.approx(np.std(x)))
    assert(tau==4.)

    def fake_autocorr_time(x_arg):
        np.testing.assert_array_equal(x_arg,x)
        return 2.
    #end def fake_autocorr_time

    monkeypatch.setattr(statistics,'autocorr_time',fake_autocorr_time)
    mean,error,tau = statistics.series_stats(x)
    expected_error = np.std(x)/np.sqrt(len(x)/tau)

    assert(mean==pytest.approx(np.mean(x)))
    assert(error==pytest.approx(expected_error))
    assert(tau==2.)

    invalid_x = [
        np.array([]),
        np.array([1.,np.nan]),
        np.array([1.+1.j]),
        np.ones((2,2)),
        ]
    for x_invalid in invalid_x:
        with pytest.raises(ValueError):
            statistics.series_stats(x_invalid,t_auto=1.)

    for t_auto_invalid in (0.,-1.,np.nan,np.inf,'invalid'):
        with pytest.raises(
            ValueError,
            match=r'autocorrelation time must be a positive finite number',
            ):
            statistics.series_stats(x,t_auto=t_auto_invalid)
#end def test_series_stats



def test_time_series_intervals():
    """Check adjacent-value intervals and their associated midpoint times."""
    x = np.array([3.,1.,2.])
    t = np.array([0.,2.,5.])

    intervals,times = statistics.time_series_intervals(x,t)
    np.testing.assert_array_equal(intervals,[[1.,3.],[1.,2.]])
    np.testing.assert_array_equal(times,[1.,3.5])

    intervals,no_times = statistics.time_series_intervals(x)
    np.testing.assert_array_equal(intervals,[[1.,3.],[1.,2.]])
    assert(no_times is None)

    column_intervals,column_times = statistics.time_series_intervals(
        x.reshape(-1,1),
        t.reshape(-1,1),
        )
    np.testing.assert_array_equal(column_intervals,intervals)
    np.testing.assert_array_equal(column_times,times)
#end def test_time_series_intervals



def test_interval_and_lcd_input_validation():
    """Check diagnostics for malformed interval and LCD inputs."""
    with pytest.raises(ValueError,match=r'data array must contain at least two values'):
        statistics.time_series_intervals([1.])
    with pytest.raises(ValueError,match=r'time array must have the same length'):
        statistics.time_series_intervals([1.,2.],[0.])
    with pytest.raises(ValueError,match=r'time array must be real-valued'):
        statistics.time_series_intervals([1.,2.],[0.+1.j,1.+1.j])

    invalid_intervals = [
        (np.array([1.,2.]),None,r'interval array must have shape'),
        (np.empty((0,2)),None,r'interval array must not be empty'),
        (np.array([[2.,1.]]),None,r'upper endpoints must not be less'),
        (np.array([1.,2.]),np.array([3.]),r'endpoint arrays must have equal lengths'),
        ]
    for lower,upper,message in invalid_intervals:
        with pytest.raises(ValueError,match=message):
            statistics.interval_distribution(lower,upper)

    intervals = np.array([[0.,1.],[1.,2.]])
    for counts,message in [
        ([1.],r'counts must have the same length'),
        ]:
        with pytest.raises(ValueError,match=message):
            statistics.interval_dist_peak(intervals,counts)
    with pytest.raises(ValueError,match=r'peak method must be a string'):
        statistics.interval_dist_peak(intervals,[1.,2.],method=1)
    for peak_frac in (np.nan,0.,-1.,1.1):
        with pytest.raises(ValueError,match=r'peak fraction must be in the interval'):
            statistics.interval_dist_peak(intervals,[1.,2.],peak_frac=peak_frac)

    for window,step,message in [
        (0,1,r'window must be a positive integer'),
        (1,0,r'step must be a positive integer'),
        (1,2,r'step must not exceed window'),
        (3,1,r'window must not exceed the number of intervals'),
        ]:
        with pytest.raises(ValueError,match=message):
            statistics.rolling_interval_dist_peak(intervals,window=window,step=step)
    with pytest.raises(ValueError,match=r'method "invalid" is unrecognized'):
        statistics.rolling_interval_dist_peak(
            intervals,window=1,step=1,method='invalid'
            )

    for x,nperm,message in [
        ([1.],0,r'data array must contain at least two values'),
        ([1.,2.],-1,r'number of permutations must be a nonnegative integer'),
        ([1.,2.],True,r'number of permutations must be a nonnegative integer'),
        ]:
        with pytest.raises(ValueError,match=message):
            statistics.line_crossing_distribution(x,nperm=nperm)

    with pytest.raises(ValueError,match=r'window must not exceed the number of intervals'):
        statistics.lcd_smooth([1.,2.],window=2,step=1)
#end def test_interval_and_lcd_input_validation



def test_interval_distribution_and_peak(monkeypatch):
    """Check interval-overlap counts and the supported peak selections."""
    endpoints = np.array([1.,2.,3.])
    upper     = np.array([4.,5.,6.])
    intervals,counts = statistics.interval_distribution(endpoints,upper)
    expected_intervals = np.array(
        [[1.,2.],[2.,3.],[3.,4.],[4.,5.],[5.,6.]]
        )
    np.testing.assert_array_equal(intervals,expected_intervals)
    np.testing.assert_array_equal(counts,[1,2,3,2,1])

    matrix_intervals,matrix_counts = statistics.interval_distribution(
        np.column_stack((endpoints,upper))
        )
    np.testing.assert_array_equal(matrix_intervals,expected_intervals)
    np.testing.assert_array_equal(matrix_counts,counts)

    column_intervals,column_counts = statistics.interval_distribution(
        endpoints.reshape(-1,1),
        upper.reshape(1,-1),
        )
    np.testing.assert_array_equal(column_intervals,expected_intervals)
    np.testing.assert_array_equal(column_counts,counts)

    touching_intervals,touching_counts = statistics.interval_distribution(
        np.array([[1.,2.],[2.,3.],[3.,4.]])
        )
    np.testing.assert_array_equal(
        touching_intervals,
        [[1.,2.],[2.,3.],[3.,4.]],
        )
    np.testing.assert_array_equal(touching_counts,[1,1,1])

    repeated_intervals,repeated_counts = statistics.interval_distribution(
        np.array([[1.,1.],[1.,2.],[2.,2.]]),perturb_const=0
        )
    np.testing.assert_array_equal(repeated_intervals,[[1.,2.]])
    np.testing.assert_array_equal(repeated_counts,[1])

    constant_intervals = np.array([[1.,1.]])
    perturbed_intervals,perturbed_counts = statistics.interval_distribution(
        constant_intervals
        )
    np.testing.assert_array_equal(
        perturbed_intervals,
        [[np.nextafter(1.,-np.inf),np.nextafter(1.,np.inf)]],
        )
    np.testing.assert_array_equal(perturbed_counts,[1])
    constant_peak,constant_height = statistics.interval_dist_peak(
        constant_intervals,
        [3.],
        height=True,
        )
    assert(constant_peak==pytest.approx(1.))
    assert(constant_height==3.)

    irregular_intervals,irregular_counts = statistics.interval_distribution(
        np.array([[0.,10.],[.5,.75]])
        )
    np.testing.assert_array_equal(
        irregular_intervals,
        [[0.,.5],[.5,.75],[.75,10.]],
        )
    np.testing.assert_array_equal(irregular_counts,[1,2,1])

    peak_intervals = np.array([[0.,1.],[1.,2.],[2.,3.]])
    peak_counts    = np.array([1,3,3])
    peak,height = statistics.interval_dist_peak(
        peak_intervals,
        peak_counts,
        height=True,
        )
    assert(peak==pytest.approx(2.))
    assert(height==3)

    column_peak = statistics.interval_dist_peak(
        peak_intervals,
        peak_counts.reshape(-1,1),
        )
    assert(column_peak==pytest.approx(peak))

    quadratic_intervals = np.array(
        [[0.,1.],[1.,2.],[2.,3.],[3.,4.],[4.,5.]]
        )
    quadratic_counts = np.array([1.,2.,3.,2.,1.])
    quadratic_peak,quadratic_height = statistics.interval_dist_peak(
        quadratic_intervals,
        quadratic_counts,
        method='quad_peak',
        height=True,
    )
    assert(quadratic_peak==pytest.approx(2.5))
    assert(np.isfinite(quadratic_height))
    assert(quadratic_height>0.)

    fallback_peak,fallback_height = statistics.interval_dist_peak(
        np.array([[0.,2.]]),
        np.array([4.]),
        method='quad_peak',
        height=True,
        )
    assert(fallback_peak==pytest.approx(1.))
    assert(fallback_height==4.)

    multimodal_intervals = np.array([[0.,1.],[1.,2.],[2.,3.],[3.,4.]])
    multimodal_counts    = np.array([4.,1.,1.,4.])
    assert(
        statistics.interval_dist_peak(
            multimodal_intervals,multimodal_counts,method='quad_peak'
            )
        ==pytest.approx(2.)
        )

    monkeypatch.setattr(
        statistics.np,
        'polyfit',
        lambda x,y,degree: np.array([1.,0.,0.]),
        )
    assert(
        statistics.interval_dist_peak(
            quadratic_intervals,quadratic_counts,method='quad_peak'
            )
        ==pytest.approx(2.5)
        )

    monkeypatch.setattr(
        statistics.np.random,
        'uniform',
        lambda size: np.full(size,.5),
        )
    assert(
        statistics.interval_dist_peak(
            peak_intervals,
            peak_counts,
            method='interval_rand',
            )
        ==pytest.approx(2.)
        )

    with pytest.raises(ValueError,match=r'unrecognized int. dist. max method'):
        statistics.interval_dist_peak(peak_intervals,peak_counts,method='invalid')
    with pytest.raises(ValueError,match=r'quadratic weighting'):
        statistics.interval_dist_peak(
            peak_intervals,peak_counts,quad_weighting='invalid'
            )
    for perturb_const in (-1,1.5,True):
        with pytest.raises(
            ValueError,
            match=r'constant perturbation must be a nonnegative integer',
            ):
            statistics.interval_distribution(
                constant_intervals,perturb_const=perturb_const
                )
        with pytest.raises(
            ValueError,
            match=r'constant perturbation must be a nonnegative integer',
            ):
            statistics.interval_dist_peak(
                constant_intervals,[1.],perturb_const=perturb_const
                )
#end def test_interval_distribution_and_peak



def test_quad_peak_width_weighting(monkeypatch):
    """Check width-based quadratic weights and convenience-API forwarding."""
    intervals = np.array(
        [[0.,1.],[1.,3.],[3.,6.],[6.,10.],[10.,15.]]
        )
    counts = np.array([1.,2.,3.,2.,1.])
    polyfit = statistics.np.polyfit
    weights = []

    def capture_polyfit(x,y,degree,**kwargs):
        weights.append(kwargs.get('w'))
        return polyfit(x,y,degree,**kwargs)
    #end def capture_polyfit

    monkeypatch.setattr(statistics.np,'polyfit',capture_polyfit)
    peak = statistics.interval_dist_peak(
        intervals,counts,method='quad_peak',quad_weighting='width'
        )
    assert(np.isfinite(peak))
    np.testing.assert_allclose(weights[0],np.sqrt([2.,2.,3.,3.,4.,4.]))

    rolling_peak, = statistics.rolling_interval_dist_peak(
        intervals,
        window=5,
        step=1,
        method='quad_peak',
        quad_weighting='width',
        )
    assert(np.isfinite(rolling_peak[0]))
#end def test_quad_peak_width_weighting



def test_rolling_interval_dist_peak_and_lcd_smooth():
    """Check rolling peak locations, heights, window bounds, and times."""
    intervals = np.array([[0.,2.],[1.,3.],[2.,4.],[3.,5.]])
    peaks,heights,windows = statistics.rolling_interval_dist_peak(
        intervals,
        window=2,
        step=2,
        ret_height=True,
        ret_windows=True,
        )
    np.testing.assert_allclose(peaks,[1.5,3.5])
    np.testing.assert_array_equal(heights,[2,2])
    assert(windows==[(0,2),(2,4)])

    quadratic_peaks,quadratic_heights = statistics.rolling_interval_dist_peak(
        np.array([[1.,4.],[2.,5.],[3.,6.]]),
        window=3,
        step=1,
        method='quad_peak',
        ret_height=True,
        )
    np.testing.assert_allclose(quadratic_peaks,[3.5])
    assert(quadratic_heights[0]>0.)

    x = np.array([0.,2.,1.,3.])
    t = np.array([0.,1.,3.,6.])
    smooth,times = statistics.lcd_smooth(x,t,window=2,step=1,method='interval_mid')
    np.testing.assert_allclose(smooth,[1.5,1.5])
    np.testing.assert_allclose(times,[1.25,3.25])
    np.testing.assert_allclose(
        statistics.lcd_smooth(x,window=2,step=1,method='interval_mid'),
        smooth,
        )
#end def test_rolling_interval_dist_peak_and_lcd_smooth



def test_line_crossing_distribution_and_lcd_peak(monkeypatch):
    """Check LCD counts, peaks, and independently accumulated permutations."""
    x = np.array([0.,2.,1.])

    def fail_shuffle(values):
        pytest.fail('the nperm=0 path must not shuffle data')
    #end def fail_shuffle

    monkeypatch.setattr(statistics.np.random,'shuffle',fail_shuffle)
    intervals,counts = statistics.line_crossing_distribution(x)
    np.testing.assert_array_equal(intervals,[[0.,1.],[1.,2.]])
    np.testing.assert_array_equal(counts,[1,2])
    assert(statistics.lcd_peak(x)==pytest.approx(1.5))

    constant = np.full(4,5.)
    constant_intervals,constant_counts = statistics.line_crossing_distribution(
        constant
        )
    np.testing.assert_array_equal(
        constant_intervals,
        [[np.nextafter(5.,-np.inf),np.nextafter(5.,np.inf)]],
        )
    np.testing.assert_array_equal(constant_counts,[3])
    assert(statistics.lcd_peak(constant)==pytest.approx(5.))

    _,_,perturbed_constant = statistics.line_crossing_distribution(
        constant,
        ret_x=True,
        )
    np.testing.assert_array_equal(
        perturbed_constant,
        [
            np.nextafter(5.,-np.inf),np.nextafter(5.,np.inf),
            np.nextafter(5.,-np.inf),np.nextafter(5.,np.inf),
            ],
        )

    def reverse(values):
        values[:] = values[::-1]
    #end def reverse

    monkeypatch.setattr(statistics.np.random,'shuffle',reverse)
    intervals,counts = statistics.line_crossing_distribution(x,nperm=2)
    np.testing.assert_array_equal(intervals,[[0.,1.],[1.,2.]])
    np.testing.assert_array_equal(counts,[1.,2.])
    assert(statistics.lcd_peak(x,nperm=2)==pytest.approx(1.5))

    permutations = [
        np.array([0.,1.,3.,6.]),
        np.array([0.,3.,1.,6.]),
        ]
    def set_permutation(values):
        values[:] = permutations.pop(0)
    #end def set_permutation

    monkeypatch.setattr(statistics.np.random,'shuffle',set_permutation)
    intervals,counts = statistics.line_crossing_distribution(
        np.array([0.,1.,3.,6.]),nperm=2
        )
    np.testing.assert_array_equal(intervals,[[0.,1.],[1.,3.],[3.,6.]])
    np.testing.assert_array_equal(counts,[1.,2.,1.])
#end def test_line_crossing_distribution_and_lcd_peak



def test_lcd_trim_find_segments():
    """Identify contiguous masked regions, including a final region."""
    x = np.arange(9.)
    mask = np.array([False,True,True,False,True,False,True,True,True])

    assert(statistics.find_segments(x,mask)==[(1,3),(4,5),(6,9)])
    assert(statistics.find_segments(x,mask,seg_min=2)==[(1,3),(6,9)])
    assert(statistics.find_segments(x,[False]*8+[True],seg_min=2)==[])

    with pytest.raises(ValueError,match=r'mask must have the same length'):
        statistics.find_segments(x,mask[:-1])
    with pytest.raises(ValueError,match=r'minimum segment length'):
        statistics.find_segments(x,mask,seg_min=0)
#end def test_lcd_trim_find_segments



@pytest.mark.parametrize(
    'trim_function,nmasks',
    [
        (statistics.lcd_trim_l,2),
        (statistics.lcd_trim_r,2),
        (statistics.lcd_trim_lr,3),
        (statistics.lcd_trim_lrm,4),
        ],
    )
def test_lcd_trim_masks_and_segments(trim_function,nmasks):
    """Return disjoint masks and clean segments consistent with those masks."""
    x = np.sin(np.linspace(0.,8.*np.pi,128))

    result = trim_function(x,niter=1,ret_seg=True,ret_mask=True)
    segments = result[:nmasks]
    masks = result[nmasks:]

    assert(len(masks)==nmasks)
    assert(all(mask.dtype==bool and mask.shape==x.shape for mask in masks))
    np.testing.assert_array_equal(
        sum(mask.astype(int) for mask in masks),
        np.ones(len(x),dtype=int),
        )
    assert(masks[0].any())
    assert(segments[0]==statistics.find_segments(x,masks[0]))
#end def test_lcd_trim_masks_and_segments



@pytest.mark.parametrize(
    'trim_function,nmasks',
    [
        (statistics.lcd_trim_l,2),
        (statistics.lcd_trim_r,2),
        (statistics.lcd_trim_lr,3),
        (statistics.lcd_trim_lrm,4),
        ],
    )
def test_lcd_trim_handles_constant_and_vector_shaped_data(
        trim_function,
        nmasks,
        ):
    """Preserve a clean point and a mask partition for constant input."""
    x = np.full(4,5.)
    masks = trim_function(x,niter=3,ret_seg=False,ret_mask=True)
    shaped_masks = trim_function(
        x.reshape(1,len(x),1),
        niter=3,
        ret_seg=False,
        ret_mask=True,
        )

    assert(len(masks)==nmasks)
    np.testing.assert_array_equal(
        sum(mask.astype(int) for mask in masks),
        np.ones(len(x),dtype=int),
        )
    assert(masks[0].any())
    for mask,shaped_mask in zip(masks,shaped_masks):
        np.testing.assert_array_equal(shaped_mask,mask)

    with pytest.raises(ValueError,match=r'data array must be one-dimensional'):
        trim_function(np.ones((2,2)))
    for invalid_niter in (0,1.5,True):
        with pytest.raises(ValueError,match=r'number of trim iterations'):
            trim_function(x,niter=invalid_niter)
#end def test_lcd_trim_handles_constant_and_vector_shaped_data



def test_lcd_trim_lrm_respects_minimum_middle_segment_length():
    """Do not remove an isolated low point below a multi-point threshold."""
    x = np.array([0.,1.,0.,1.,0.,-100.,0.,1.,0.,1.,0.])

    clean,left,right,middle = statistics.lcd_trim_lrm(
        x,
        niter=1,
        nseg_min=2,
        ret_seg=False,
        ret_mask=True,
        )

    assert(middle.sum()==0)
    np.testing.assert_array_equal(
        clean|left|right,
        np.ones(len(x),dtype=bool),
        )
    for low_scale in (0.,-1.,np.nan):
        with pytest.raises(ValueError,match=r'low scale must be a positive finite'):
            statistics.lcd_trim_lrm(x,low_scale=low_scale)
    with pytest.raises(ValueError,match=r'minimum segment length'):
        statistics.lcd_trim_lrm(x,nseg_min=0)
#end def test_lcd_trim_lrm_respects_minimum_middle_segment_length



def test_pair_expand_ts_intervals():
    """Expand each sample into duplicated, bounded-lag sorted intervals."""
    x = np.arange(5.)
    t = 2.*x

    intervals,times = statistics.pair_expand_ts_intervals(x,t,expand=2)
    expected_intervals = np.array(
        [
            [0.,1.],[0.,1.],[1.,2.],[1.,2.],
            [2.,3.],[2.,3.],[3.,4.],[3.,4.],
            ],
        )
    np.testing.assert_array_equal(intervals,expected_intervals)
    np.testing.assert_allclose(times,[.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5])

    intervals_without_times,times_without_times = (
        statistics.pair_expand_ts_intervals(x.reshape(1,len(x),1),expand=2)
        )
    np.testing.assert_array_equal(intervals_without_times,expected_intervals)
    assert(times_without_times is None)

    invalid_inputs = [
        (np.arange(2.),None,2,r'data array length must exceed expansion'),
        (x,None,0,r'expansion must be a positive even integer'),
        (x,None,3,r'expansion must be a positive even integer'),
        (x,None,True,r'expansion must be a positive even integer'),
        (np.ones((2,2)),None,2,r'data array must be one-dimensional'),
        (x,[0.,1.],2,r'time array must have the same length'),
        ]
    for x_invalid,t_invalid,expand,message in invalid_inputs:
        with pytest.raises(ValueError,match=message):
            statistics.pair_expand_ts_intervals(x_invalid,t_invalid,expand=expand)
#end def test_pair_expand_ts_intervals



@pytest.mark.parametrize(
    'method',
    ['simple','lcd_trim_l','lcd_trim_r','lcd_trim_lr','lcd_trim_lrm'],
    )
def test_time_series_analyzer_modes(method):
    """Analyze a representative uniformly sampled trace in each mode."""
    x = np.sin(np.linspace(0.,8.*np.pi,128))
    analyzer = statistics.TimeSeriesAnalyzer(x,clean_inp=method)

    assert(analyzer.clean_inp==method)
    assert(len(analyzer.xc)==len(analyzer.indc))
    assert(len(analyzer.xc)>=2)
    assert(np.isfinite(analyzer.x_mean))
    assert(np.isfinite(analyzer.x_stderr))
    assert(np.isfinite(analyzer.t_auto))
#end def test_time_series_analyzer_modes



def test_time_series_analyzer_read_partition_and_intersect(tmp_path):
    """Initialize from a file and transfer a clean partition by index."""
    x = np.sin(np.linspace(0.,8.*np.pi,128))
    filepath = tmp_path/'series.dat'
    np.savetxt(filepath,x)

    source = statistics.TimeSeriesAnalyzer(str(filepath),clean_inp='lcd_trim_l')
    target = statistics.TimeSeriesAnalyzer(analyze=False)
    target.read(str(filepath))
    target.partition_from_timeseries(source)
    xc_source,xc_target,indices = source.clean_intersect(target)

    np.testing.assert_array_equal(source.x,x)
    np.testing.assert_array_equal(target.xc,source.xc)
    np.testing.assert_array_equal(target.indc,source.indc)
    np.testing.assert_array_equal(indices,source.indc)
    np.testing.assert_array_equal(xc_source,source.xc)
    np.testing.assert_array_equal(xc_target,target.xc)
#end def test_time_series_analyzer_read_partition_and_intersect



def test_time_series_analyzer_forwards_trim_options(monkeypatch):
    """Pass LCD trim settings without mutating the caller configuration."""
    captured = {}

    def trim(x,**kwargs):
        captured.update(kwargs)
        return np.ones(len(x),dtype=bool),np.zeros(len(x),dtype=bool)
    #end def trim

    monkeypatch.setattr(statistics,'lcd_trim_l',trim)
    clean_inp = statistics.obj(method='lcd_trim_l',niter=1)
    statistics.TimeSeriesAnalyzer(np.arange(8.),clean_inp=clean_inp)

    assert(captured=={'niter':1,'ret_seg':False,'ret_mask':True})
    assert(clean_inp.method=='lcd_trim_l')

    single = statistics.TimeSeriesAnalyzer([1.],clean_inp='simple')
    assert(single.x_mean==1.)
    assert(single.x_stderr==0.)
    assert(single.t_auto==1.)
    with pytest.raises(ValueError,match=r'unrecognized data cleaning method'):
        statistics.TimeSeriesAnalyzer(np.arange(8.),clean_inp='invalid')
#end def test_time_series_analyzer_forwards_trim_options



def test_time_series_analyzer_reanalysis_is_configuration_independent():
    """Repeated and explicit analyses depend only on the requested method."""
    x = np.sin(np.linspace(0.,8.*np.pi,128))
    analyzer = statistics.TimeSeriesAnalyzer(x,clean_inp='lcd_trim_l')
    first = (analyzer.xc.copy(),analyzer.indc.copy(),analyzer.x_mean,
             analyzer.x_stderr,analyzer.t_auto)

    analyzer.analyze()
    second = (analyzer.xc,analyzer.indc,analyzer.x_mean,
              analyzer.x_stderr,analyzer.t_auto)
    for value,expected in zip(second,first):
        if isinstance(value,np.ndarray):
            np.testing.assert_array_equal(value,expected)
        else:
            assert(value==expected)

    analyzer.analyze(clean_inp='simple')
    simple = statistics.TimeSeriesAnalyzer(x,clean_inp='simple')
    np.testing.assert_array_equal(analyzer.xc,simple.xc)
    np.testing.assert_array_equal(analyzer.indc,simple.indc)
    assert(analyzer.x_mean==simple.x_mean)
    assert(analyzer.x_stderr==simple.x_stderr)
    assert(analyzer.t_auto==simple.t_auto)
#end def test_time_series_analyzer_reanalysis_is_configuration_independent



def test_time_series_analyzer_plot():
    """Plot labeled clean and full traces only after analysis has run."""
    matplotlib = pytest.importorskip('matplotlib')
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    x = np.sin(np.linspace(0.,4.*np.pi,64))
    uninitialized = statistics.TimeSeriesAnalyzer(x,analyze=False)
    with pytest.raises(ValueError,match=r'analysis must be completed'):
        uninitialized.plot()

    analyzer = statistics.TimeSeriesAnalyzer(x,clean_inp='lcd_trim_l')
    figure,axis = plt.subplots()
    plt.sca(axis)
    analyzer.plot(legend=True)
    labels = {text.get_text() for text in axis.get_legend().get_texts()}
    assert({'full series','clean mean','clean mean ± std. dev.','clean series'}<=labels)
    plt.close(figure)
#end def test_time_series_analyzer_plot



def test_plot_interval_dist():
    """Check that interval-distribution plotting adds the expected lines."""
    matplotlib = pytest.importorskip('matplotlib')
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    figure,axis = plt.subplots()
    plt.sca(axis)
    statistics.plot_interval_dist(
        np.array([[0.,1.],[1.,2.]]),
        np.array([1,2]),
        )
    assert(len(axis.lines)==2)
    plt.close(figure)
#end def test_plot_interval_dist
def test_mean_smooth():
    """Check moving means, tapered endpoints, and automatic width."""
    x = np.array([0.,10.,0.,10.,0.])
    expected = np.array([0.,10./3.,20./3.,10./3.,0.])

    np.testing.assert_allclose(statistics.mean_smooth(x,m=3),expected)
    np.testing.assert_allclose(statistics.mean_smooth(x),expected)
#end def test_mean_smooth



def test_median_smooth():
    """Check median smoothing suppresses spikes and supports a mean pass."""
    x = np.array([0.,10.,0.,10.,0.])
    median = np.array([0.,0.,10.,0.,0.])
    post_mean = np.array([0.,10./3.,10./3.,10./3.,0.])

    np.testing.assert_allclose(statistics.median_smooth(x,m=3),median)
    np.testing.assert_allclose(statistics.median_smooth(x),median)
    np.testing.assert_allclose(
        statistics.median_smooth(x,m=3,post_mean=True),
        post_mean,
        )
#end def test_median_smooth



def test_poly_smooth(capsys):
    """Check local linear fits and the optional mean post-processing pass."""
    x = np.array([0.,10.,0.,10.,0.])
    polynomial = np.array([0.,10./3.,20./3.,10./3.,0.])
    post_mean = np.array([0.,10./3.,40./9.,10./3.,0.])

    np.testing.assert_allclose(statistics.poly_smooth(x,m=3),polynomial)
    np.testing.assert_allclose(
        statistics.poly_smooth(x,m=3,post_mean=True),
        post_mean,
        )
    np.testing.assert_allclose(statistics.poly_smooth(np.arange(24.)),np.arange(24.))
    assert(capsys.readouterr().out=='')
#end def test_poly_smooth



def test_local_median_smooth():
    """Check leave-one-out pooling and each selectable second pass."""
    x_list = [np.array([value]) for value in (0.,10.,0.,10.,0.)]
    median = np.array([0.,0.,10.,0.,0.])
    mean = np.array([0.,10./3.,10./3.,10./3.,0.])

    np.testing.assert_allclose(
        statistics.local_median_smooth(x_list,m=3,poly_smooth=False),
        median,
        )
    np.testing.assert_allclose(
        statistics.local_median_smooth(x_list,m=3),
        mean,
        )
    np.testing.assert_allclose(statistics.local_median_smooth(x_list),mean)
    np.testing.assert_allclose(
        statistics.local_median_smooth(
            x_list,
            m=3,
            poly_smooth=False,
            post_mean=True,
            ),
        mean,
        )
#end def test_local_median_smooth



@pytest.mark.parametrize(
    'smoother,kwargs',
    [
        (statistics.mean_smooth,{}),
        (statistics.median_smooth,{}),
        (statistics.poly_smooth,{}),
        (statistics.local_median_smooth,{'poly_smooth':False}),
        ],
    )
def test_smoothers_validate_window_length(smoother,kwargs):
    """Require an in-range, positive odd integer smoothing window."""
    x = np.arange(5.)
    if smoother is statistics.local_median_smooth:
        x = [np.array([value]) for value in x]

    result = smoother(x,m=np.int64(3),**kwargs)
    assert(len(result)==len(x))

    invalid_windows = [
        (3.,TypeError,r'smoothing window length must be an integer'),
        (True,TypeError,r'smoothing window length must be an integer'),
        (0,ValueError,r'smoothing window length must be positive'),
        (-1,ValueError,r'smoothing window length must be positive'),
        (2,ValueError,r'smoothing window length must be odd'),
        (7,ValueError,r'smoothing window length must not exceed the data length'),
        ]
    for m,error_type,message in invalid_windows:
        with pytest.raises(error_type,match=message):
            smoother(x,m=m,**kwargs)
#end def test_smoothers_validate_window_length



def test_poly_smooth_validates_maximum_window_length():
    """Reject polynomial widths beyond the supported order table."""
    with pytest.raises(
        ValueError,
        match=r'smoothing window length must not exceed 21',
        ):
        statistics.poly_smooth(np.arange(23.),m=23)
#end def test_poly_smooth_validates_maximum_window_length



def test_smoothers_validate_boolean_options():
    """Require explicit Boolean values for smoothing options."""
    x = np.arange(5.)
    x_list = [np.array([value]) for value in x]

    for smoother,kwargs,name in (
        (statistics.median_smooth,{'post_mean':1},'post_mean'),
        (statistics.poly_smooth,{'post_mean':'yes'},'post_mean'),
        (statistics.local_median_smooth,{'poly_smooth':0},'poly_smooth'),
        (statistics.local_median_smooth,{'post_mean':None},'post_mean'),
        ):
        data = x_list if smoother is statistics.local_median_smooth else x
        with pytest.raises(TypeError,match=rf'{name} must be a Boolean value'):
            smoother(data,m=3,**kwargs)

    assert(len(statistics.median_smooth(x,m=3,post_mean=np.bool_(True)))==len(x))
    assert(len(statistics.local_median_smooth(
        x_list,
        m=3,
        poly_smooth=np.bool_(False),
        post_mean=np.bool_(True),
        ))==len(x_list))
#end def test_smoothers_validate_boolean_options
