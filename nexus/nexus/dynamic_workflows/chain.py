##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  chain.py                                                          #
#    Chain mode: sequential parameter chain (initial / propose / stop).#
#                                                                    #
#  Content summary:                                                  #
#    SuccessiveChange(Target)                                        #
#      Successive |Δ| of one or more products within atol/rtol.      #
#    ChainDecision(DynamicDecision)                                  #
#      Outcome of one chain observation.                             #
#    DynamicChain(DynamicMode)                                       #
#      Sequential chain.  stop is Target.reached; observe returns    #
#      ChainDecision.  drive polls until completed / max_runs /      #
#      failed.  Subclasses implement initial / propose.              #
#====================================================================#


from ..developer import error
from .base import DynamicDecision, DynamicMode, Target, _as_products, _iter_work


def _as_params(params):
    if not isinstance(params, dict) or len(params) == 0:
        error(
            f'chain parameters must be a non-empty dict, received {params}',
            header='DynamicChain',
            )
    #end if
    return dict(params)
#end def _as_params


class ChainDecision(DynamicDecision):
    """Outcome of ``DynamicChain.observe`` or ``drive``.

    status
      ``continue``   propose a new parameter point
      ``completed``  target reached
      ``max_runs``   history length reached max_runs
      ``failed``     no next work (propose returned None), or a sim failed
    """

    def __init__(self, status, params, next_params=None, products=None):
        DynamicDecision.__init__(self, status, products=products)
        self.params      = None if params is None else dict(params)
        self.next_params = None if next_params is None else dict(next_params)
    #end def __init__
#end class ChainDecision


class SuccessiveChange(Target):
    """Target: successive values of one or more products have settled.

    A pair ``(old, new)`` of a product is within tolerance if::

        abs(new - old) <= atol + rtol * abs(old)

    (same combination as ``numpy.isclose``).  Every listed key must pass
    for the last ``consecutive`` successive pairs.

    Parameters
    ----------
    keys : str or sequence of str
        Product names, e.g. ``'energy'`` or ``('energy', 'variance')``.
    atol : float or dict, optional
        Absolute tolerance.  A scalar applies to every key; a dict maps
        key → atol.  Default is ``0.0`` when ``rtol`` is set.
    rtol : float or dict, optional
        Relative tolerance against the previous value.  Default is ``0.0``
        when ``atol`` is set.
    consecutive : int, optional
        Number of successive pairs that must all pass.  Default is ``1``.

    Examples
    --------
    abs_energy = SuccessiveChange('energy', atol=1e-4)
    rel_energy = SuccessiveChange('energy', rtol=1e-3)
    mixed_target = SuccessiveChange(
        ('energy', 'variance'),
        atol={'energy': 1e-4, 'variance': 0.01},
        )
    """

    def __init__(self, keys='energy', atol=None, rtol=None, consecutive=1):
        if isinstance(keys, str):
            keys = (keys,)
        else:
            keys = tuple(keys)
        #end if
        if len(keys) == 0:
            self.error('SuccessiveChange requires at least one product key')
        #end if
        if atol is None and rtol is None:
            self.error('SuccessiveChange requires atol or rtol')
        #end if
        consecutive = int(consecutive)
        if consecutive < 1:
            self.error(f'consecutive must be >= 1, received {consecutive}')
        #end if
        self.keys        = keys
        self.atol        = _tol_map(atol, keys, 0.0)
        self.rtol        = _tol_map(rtol, keys, 0.0)
        self.consecutive = consecutive
    #end def __init__

    def _value(self, record, key):
        products = record.get('products') or {}
        if key not in products:
            self.error(f'SuccessiveChange target requires products[{key!r}]')
        #end if
        value = float(products[key])
        return value
    #end def _value

    def _within(self, key, old, new):
        cond_atol = abs(new - old) <= self.atol[key]
        cond_rtol = abs(new - old) <= self.rtol[key] * abs(old)
        return cond_atol and cond_rtol
    #end def _within

    def reached(self, history):
        needed = self.consecutive + 1
        if len(history) < needed:
            return False
        #end if
        for offset in range(self.consecutive):
            newer = history[-(offset + 1)]
            older = history[-(offset + 2)]
            for key in self.keys:
                if not self._within(key, self._value(older, key), self._value(newer, key)):
                    return False
                #end if
            #end for
        #end for
        return True
    #end def reached
#end class SuccessiveChange


available_targets = {
    'successive_change': SuccessiveChange,
    }


def _tol_map(spec, keys, default):
    if spec is None:
        return {key: float(default) for key in keys}
    #end if
    if isinstance(spec, dict):
        missing = [key for key in keys if key not in spec]
        if missing and default is None:
            error(
                f'missing tolerances for product keys: {missing}',
                header='SuccessiveChange',
                )
        #end if
        return {key: float(spec[key]) if key in spec else float(default) for key in keys}
    #end if
    val = float(spec)
    return {key: val for key in keys}
#end def _tol_map


class DynamicChain(DynamicMode):
    """Sequential parameter chain (chain mode).

    ``stop`` is ``self.target.reached(history)``.  ``observe`` returns
    ``ChainDecision``.  ``drive`` polls until the chain is not continue.
    """

    def __init__(self, target=None, max_runs=10):
        DynamicMode.__init__(self, max_runs=max_runs)
        if target is not None and not isinstance(target, Target):
            tname = type(target).__name__
            self.error(
                f'target must be a Target instance with reached(history), received {tname}'
                )
        #end if
        self.target = target
    #end def __init__

    def drive(self, sim_generator, wm, products, poll=1):
        """Poll until the chain completes, hits max_runs, or a sim fails."""
        sims = []
        for params in _iter_work(self.initial()):
            sims.append([params, sim_generator(params)])
        #end for
        if len(sims) == 0:
            self.error('initial() produced no work to launch')
        #end if
        while True:
            for sim_info in list(sims):
                params, sim = sim_info
                if sim.succ:
                    decision = self.observe(params, _as_products(products, sim))
                    sims.remove(sim_info)
                    if decision.status != 'continue':
                        return decision
                    #end if
                    nxt = list(_iter_work(decision.next_params))
                    if len(nxt) == 0:
                        self.error('continue decision has no next_params')
                    #end if
                    for p in nxt:
                        sims.append([p, sim_generator(p)])
                    #end for
                elif sim.fail:
                    return ChainDecision(
                        'failed',
                        params = params,
                        )
                #end if
            #end for
            if len(sims) == 0:
                self.error('no remaining jobs')
            #end if
            wm.poll(poll)
        #end while
    #end def drive

    def stop(self, history):
        if self.target is None:
            return DynamicMode.stop(self, history)
        #end if
        return self.target.reached(history)
    #end def stop

    def observe(self, params, products):
        params   = _as_params(params)
        products = dict(products)
        self.history.append(dict(params=params, products=products))
        if self.stop(self.history):
            return ChainDecision(
                'completed',
                params   = params,
                products = products,
                )
        #end if
        nxt = None
        if len(self.history) < self.max_runs:
            nxt = self.propose(self.history)
        #end if
        if nxt is None:
            if len(self.history) >= self.max_runs:
                status = 'max_runs'
            else:
                status = 'failed'
            #end if
            return ChainDecision(
                status,
                params   = params,
                products = products,
                )
        #end if
        return ChainDecision(
            'continue',
            params      = params,
            next_params = _as_params(nxt),
            products    = products,
            )
    #end def observe
#end class DynamicChain
