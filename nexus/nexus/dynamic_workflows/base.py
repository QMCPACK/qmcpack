##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  base.py                                                           #
#    Shared types for dynamic-workflow modes (walk, later spawn).    #
#                                                                    #
#  Content summary:                                                  #
#    Target                                                          #
#      Sequential stop rule: reached(history).  Used by walk.        #
#      Spawn pick-rules should not inherit Target.                   #
#    DynamicDecision                                                 #
#      status / products / completed / max_runs.                     #
#      WalkDecision and spawn decisions subclass this.               #
#    DynamicMode                                                     #
#      Shared protocol: initial / propose / stop / observe / drive.  #
#      history, max_runs.                                            #
#====================================================================#


from ..developer import DevBase


def _iter_work(params):
    """Yield parameter dicts.  Walk: one dict.  Spawn: a sequence of dicts."""
    if params is None:
        return
    #end if
    if isinstance(params, dict):
        yield params
        return
    #end if
    for item in params:
        yield item
    #end for
#end def _iter_work


def _as_products(spec, sim):
    if callable(spec):
        return dict(spec(sim))
    #end if
    if isinstance(spec, str):
        spec = (spec,)
    #end if
    products = {}
    for key in spec:
        products[key] = getattr(sim.products, key)
    #end for
    return products
#end def _as_products


class Target(DevBase):
    """Sequential stop rule.  Subclasses implement ``reached(history)``.

    ``history`` is a list of ``{params, products}`` records in time order.
    This protocol is sequential; spawn pick-rules do not inherit Target.
    """

    def reached(self, history):
        self.error('Target.reached() must be implemented in a subclass')
    #end def reached
#end class Target


class DynamicDecision(DevBase):
    """Outcome of one dynamic-mode observation.

    ``status`` and ``products`` are shared with future spawn decisions.

    status
      ``continue``   keep polling this step
      ``completed``  finished successfully (walk: target reached)
      ``max_runs``   stopped without completing (cap, or propose() is None)
    """

    def __init__(self, status, products=None):
        self.status   = status
        self.products = None if products is None else dict(products)
    #end def __init__

    @property
    def completed(self):
        return self.status == 'completed'
    #end def completed

    @property
    def max_runs(self):
        return self.status == 'max_runs'
    #end def max_runs
#end class DynamicDecision


class DynamicMode(DevBase):
    """Shared protocol for dynamic workflows: initial/propose/drive/observe/stop.

    initial: First work to launch.
    propose(history): Next work(s) to launch, or None if there is no further work.
    drive(make, wm, products): Poll until observe is not continue.
    observe(params, products): Record a finished job and return a ``DynamicDecision``.
    stop(history): True when the mode has finished successfully.  Default is False.
    """

    def __init__(self, max_runs=10):
        max_runs = int(max_runs)
        if max_runs < 1:
            self.error(f'max_runs must be >= 1, received {max_runs}')
        #end if
        self.max_runs = max_runs
        self.history  = []
    #end def __init__

    def initial(self):
        self.error('initial() must be implemented in a subclass')
    #end def initial

    def propose(self, history):
        self.error('propose() must be implemented in a subclass')
    #end def propose

    def drive(self, make, wm, products, poll=1):
        """Return a DynamicDecision after polling until the mode is finished."""
        jobs = []
        for params in _iter_work(self.initial()):
            jobs.append([params, make(params)])
        #end for
        if len(jobs) == 0:
            self.error('initial() produced no work to launch')
        #end if
        while True:
            for job in list(jobs):
                params, sim = job
                if sim.succ:
                    decision = self.observe(params, _as_products(products, sim))
                    jobs.remove(job)
                    if decision.status != 'continue':
                        return decision
                    #end if
                    nxt = list(_iter_work(decision.next_params))
                    if len(nxt) == 0:
                        self.error('continue decision has no next_params')
                    #end if
                    for p in nxt:
                        jobs.append([p, make(p)])
                    #end for
                elif sim.fail:
                    self.error('simulation failed')
                #end if
            #end for
            if len(jobs) == 0:
                self.error('no remaining jobs')
            #end if
            wm.poll(poll)
        #end while
    #end def drive

    def observe(self, params, products):
        self.error('observe() must be implemented in a subclass')
    #end def observe

    def stop(self, history):
        return False
    #end def stop
#end class DynamicMode
