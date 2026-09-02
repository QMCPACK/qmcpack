


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

