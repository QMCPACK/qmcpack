import numpy as np
import pytest

from . import NexusTestOrder
from .. import statistics

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

    with pytest.raises(AssertionError):
        statistics.theil_sen(x,y[:-1])
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

    assert(np.isfinite(tau_iid))
    assert(tau_iid>0.)
    assert(tau_correlated>3.*tau_iid)
    assert(tau_column==pytest.approx(tau_correlated))
#end def test_reblocked_autocorr_time



def test_reblocked_autocorr_time_invalid_input():
    """Check reblocking validation for invalid shapes and limits.

    Each diagnostic-bearing assertion must report its expected message.
    """
    test_cases = [
        (np.array([]),10,r'data array must not be empty'),
        (np.arange(8),0,r'minimum number of blocks must be at least one'),
        ]
    for x,min_blocks,message in test_cases:
        with pytest.raises(AssertionError,match=message):
            statistics.reblocked_autocorr_time(x,min_blocks=min_blocks)

    with pytest.raises(AssertionError):
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
    assert(tau_trend>1.)
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
#end def test_series_stats
