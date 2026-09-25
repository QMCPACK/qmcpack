##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  ecut.py                                                           #
#    Chain-mode PWscf planewave-cutoff (example 02 repeat).          #
#====================================================================#


from ..chain import DynamicChain, SuccessiveChange
from ...developer import error


def next_ecutwfc(ecut, factor=1.5, round_ry=10):
    """Increase cutoff for convergence example. 
    Uses ecut_nxt = factor * ecut, but rounded using round_ry.
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


class Ecut(DynamicChain):
    """1D Ecut parameter scan until target is reached. 
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
        DynamicChain.__init__(self, target=target, max_runs=max_runs)
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
