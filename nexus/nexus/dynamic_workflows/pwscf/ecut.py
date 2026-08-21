##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  ecut.py                                                           #
#    Walk-mode PWscf planewave-cutoff (example 02 repeat).           #
#====================================================================#


from ..walk import DynamicWalk, SuccessiveChange
from ...developer import error


def next_ecutwfc(ecut, factor=1.5, round_ry=10):
    """Increase cutoff as in the diamond energy-convergence example.

    ``int(((factor * ecut) // round_ry) * round_ry)``, then bump by
    ``round_ry`` if the rounded value did not increase.
    """
    ecut     = float(ecut)
    factor   = float(factor)
    round_ry = int(round_ry)
    if ecut <= 0.0:
        error(
            f'ecutwfc must be positive, received {ecut}',
            header='Ecut',
            )
    #end if
    if factor <= 1.0:
        error(
            f'ecut factor must be > 1, received {factor}',
            header='Ecut',
            )
    #end if
    if round_ry < 1:
        error(
            f'round_ry must be >= 1, received {round_ry}',
            header='Ecut',
            )
    #end if
    nxt = int(((factor * ecut) // round_ry) * round_ry)
    if nxt <= int(ecut):
        nxt = int(ecut) + round_ry
    #end if
    return nxt
#end def next_ecutwfc


class Ecut(DynamicWalk):
    """Walk ``ecutwfc`` until a successive-change target is reached.

    Parameter dict: ``{'ecutwfc': int}``.  Default stepping matches
    ``examples/dynamic_workflows/02_energy_convergence`` (×1.5, rounded
    down to 10 Ry).  The poll loop and ``generate_pwscf`` stay in the
    user script.

    ``target`` may be a target object (``reached(history)``), or the
    name ``'consecutive'`` (default), which builds
    ``SuccessiveChange('energy', atol=1e-4, ...)``.  ``tol=`` is accepted
    as an alias for ``atol=``.
    """

    def __init__(
        self,
        start=50,
        end=None,
        factor=1.5,
        round_ry=10,
        max_runs=10,
        target='consecutive',
        **target_kw
        ):
        start = int(start)
        if start < 1:
            self.error(f'start ecutwfc must be >= 1, received {start}')
        #end if
        if end is not None:
            end = int(end)
            if end < start:
                self.error(f'end must be >= start, received {end} < {start}')
            #end if
        #end if
        if target is None or target == 'consecutive':
            if 'tol' in target_kw and 'atol' not in target_kw:
                target_kw['atol'] = target_kw.pop('tol')
            #end if
            if 'atol' not in target_kw and 'rtol' not in target_kw:
                target_kw['atol'] = 1e-4
            #end if
            target = SuccessiveChange('energy', **target_kw)
        elif not hasattr(target, 'reached'):
            self.error(f'invalid target: {target}')
        #end if
        DynamicWalk.__init__(self, target=target, max_runs=max_runs)
        self.start    = start
        self.end      = end
        self.factor   = float(factor)
        self.round_ry = int(round_ry)
    #end def __init__

    def initial(self):
        return dict(ecutwfc=self.start)
    #end def initial

    def propose(self, history):
        ecut = int(history[-1]['params']['ecutwfc'])
        nxt  = next_ecutwfc(ecut, factor=self.factor, round_ry=self.round_ry)
        if self.end is not None and nxt > self.end:
            return None
        #end if
        return dict(ecutwfc=nxt)
    #end def propose
#end class Ecut
