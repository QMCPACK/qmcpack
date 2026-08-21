##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  walk.py                                                           #
#    Walk mode: sequential parameter walk (initial / propose / stop).#
#                                                                    #
#  Content summary:                                                  #
#    SuccessiveChange(Target)                                        #
#      Successive |Δ| of one or more products within atol/rtol.      #
#    WalkDecision(DynamicDecision)                                   #
#      Outcome of one walk observation.                              #
#    DynamicWalk(DynamicMode)                                        #
#      Sequential walk.  stop is Target.reached; observe returns     #
#      WalkDecision.  Subclasses implement initial / propose.        #
#====================================================================#


from ..developer import error
from .base import DynamicDecision, DynamicMode, Target


def _as_params(params):
    if not isinstance(params, dict) or len(params) == 0:
        error(
            f'walk parameters must be a non-empty dict, received {params}',
            header='DynamicWalk',
            )
    #end if
    return dict(params)
#end def _as_params


class WalkDecision(DynamicDecision):
    """Outcome of ``DynamicWalk.observe``.

    status
      ``continue``   propose a new parameter point
      ``completed``  target reached
      ``max_runs``   history length reached max_runs, or propose() returned None
    """

    def __init__(self, status, params, next_params=None, products=None):
        DynamicDecision.__init__(self, status, products=products)
        self.params      = None if params is None else dict(params)
        self.next_params = None if next_params is None else dict(next_params)
    #end def __init__
#end class WalkDecision


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
            self.error('SuccessiveChange requires atol and/or rtol')
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
        return abs(new - old) <= self.atol[key] + self.rtol[key] * abs(old)
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


class DynamicWalk(DynamicMode):
    """Sequential parameter walk (walk mode).

    ``stop`` is ``self.target.reached(history)``.  ``observe`` returns
    ``WalkDecision``.
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
            return WalkDecision(
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
            return WalkDecision(
                'max_runs',
                params   = params,
                products = products,
                )
        #end if
        return WalkDecision(
            'continue',
            params      = params,
            next_params = _as_params(nxt),
            products    = products,
            )
    #end def observe
#end class DynamicWalk
