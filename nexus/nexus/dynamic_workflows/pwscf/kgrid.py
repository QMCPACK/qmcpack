##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  kgrid.py                                                          #
#    Walk-mode PWscf k-point mesh (example 02 policy).               #
#                                                                    #
#    increment  add to the integer Monkhorst-Pack grid each step.    #
#    spacing    initial mesh from Structure.kgrid_from_kspacing.     #
#====================================================================#


from ..walk import DynamicWalk, SuccessiveChange
from ...developer import error


def _as_kgrid(kgrid, name='kgrid'):
    if isinstance(kgrid, (int, float)):
        n = int(kgrid)
        kgrid = (n, n, n)
    else:
        kgrid = tuple(int(n) for n in kgrid)
    #end if
    if len(kgrid) != 3 or min(kgrid) < 1:
        error(
            f'{name} must be a positive int or 3-tuple, received {kgrid}',
            header='Kgrid',
            )
    #end if
    return kgrid
#end def _as_kgrid


def _as_increment(increment):
    if isinstance(increment, (int, float)):
        d = int(increment)
        increment = (d, d, d)
    else:
        increment = tuple(int(n) for n in increment)
    #end if
    if len(increment) != 3 or min(increment) < 0 or max(increment) < 1:
        error(
            f'increment must be a positive int or 3-tuple with at least one > 0, received {increment}',
            header='Kgrid',
            )
    #end if
    return increment
#end def _as_increment


def _structure(system):
    if system is None:
        error(
            'structure is required when spacing is set',
            header='Kgrid',
            )
    #end if
    struct = system.structure if hasattr(system, 'structure') else system
    if not hasattr(struct, 'kgrid_from_kspacing'):
        tname = type(system).__name__
        error(
            f'structure must provide kgrid_from_kspacing, received {tname}',
            header='Kgrid',
            )
    #end if
    return struct
#end def _structure


def kgrid_from_kspacing(structure, kspacing):
    """Integer mesh from ``Structure.kgrid_from_kspacing``."""
    kspacing = float(kspacing)
    if kspacing <= 0.0:
        error(
            f'kspacing must be positive, received {kspacing}',
            header='Kgrid',
            )
    #end if
    return tuple(int(n) for n in _structure(structure).kgrid_from_kspacing(kspacing))
#end def kgrid_from_kspacing


def next_kgrid(kgrid, increment=1):
    """Increase a Monkhorst-Pack grid by ``increment`` (example 02: +1).

    ``kgrid`` and ``increment`` may be a positive int (applied on all
    three axes) or a 3-tuple.
    """
    kgrid     = _as_kgrid(kgrid)
    increment = _as_increment(increment)
    return tuple(n + d for n, d in zip(kgrid, increment))
#end def next_kgrid


class Kgrid(DynamicWalk):
    """Walk ``kgrid`` until a successive-change target is reached.

    Parameter dict: ``{'kgrid': (nx, ny, nz)}``.  Default stepping
    matches ``examples/dynamic_workflows/02_energy_convergence``
    (start ``(1,1,1)``, add 1 to each axis).  The poll loop and
    ``generate_pwscf`` stay in the user script.

    ``spacing``
      If set, the first mesh is ``structure.kgrid_from_kspacing(spacing)``.
      Later points still use ``increment``.  ``structure`` may be a
      ``Structure`` or a ``PhysicalSystem``.
    ``increment``
      Added to the integer grid each step (int or 3-tuple).
    ``end``
      Optional cap (int or 3-tuple).  ``propose`` returns ``None`` if
      the next mesh would exceed it on any axis.
    """

    def __init__(
        self,
        start=(1, 1, 1),
        end=None,
        increment=1,
        spacing=None,
        structure=None,
        max_runs=10,
        target='successive_change',
        **target_kw
        ):
        increment = _as_increment(increment)
        if spacing is not None:
            start = kgrid_from_kspacing(structure, spacing)
        else:
            start = _as_kgrid(start, 'start')
        #end if
        if end is not None:
            end = _as_kgrid(end, 'end')
            if any(s > e for s, e in zip(start, end)):
                self.error(f'end must be >= start, received {end} < {start}')
            #end if
        #end if
        if target is None or target == 'successive_change':
            if 'tol' in target_kw and 'atol' not in target_kw:
                target_kw['atol'] = target_kw.pop('tol')
            #end if
            if 'atol' not in target_kw and 'rtol' not in target_kw:
                target_kw['atol'] = 1e-3
            #end if
            target = SuccessiveChange('energy', **target_kw)
        elif not hasattr(target, 'reached'):
            self.error(f'invalid target: {target}')
        #end if
        DynamicWalk.__init__(self, target=target, max_runs=max_runs)
        self.start     = start
        self.end       = end
        self.increment = increment
        self.spacing   = None if spacing is None else float(spacing)
    #end def __init__

    def initial(self):
        return dict(kgrid=self.start)
    #end def initial

    def propose(self, history):
        kgrid = _as_kgrid(history[-1]['params']['kgrid'])
        nxt   = next_kgrid(kgrid, increment=self.increment)
        if self.end is not None and any(n > e for n, e in zip(nxt, self.end)):
            return None
        #end if
        return dict(kgrid=nxt)
    #end def propose
#end class Kgrid
