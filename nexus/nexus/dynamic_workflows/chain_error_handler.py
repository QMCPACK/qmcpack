##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  chain_error_handler.py                                            #
#    Chain recovery for a failed simulation.                         #
#                                                                    #
#  Content summary:                                                  #
#    ChainErrorHandler(DynamicChain)                                 #
#      Parse a failed job, propose one input patch per attempt.      #
#      First matching handler wins.                                  #
#      Subclasses implement parse / apply_patch / params_from_input. #
#      Spawn can later apply every matching handler in parallel.     #
#====================================================================#


import os

from .base import DynamicMode
from .chain import DynamicChain


BOOKKEEPING = frozenset(('attempt',))
CHUNK_BYTES = 64 * 1024


def as_simulation(sim):
    inner = getattr(sim, 'sim', None)
    if inner is not None and hasattr(sim, 'dpid'):
        return inner
    #end if
    return sim
#end def as_simulation


def scan_file(path, markers, extra_re=None):
    """Return (matched marker strings, extra_re match or None).  Chunked."""
    found = []
    extra = None
    if not path or not os.path.isfile(path):
        return found, extra
    #end if
    remaining = []
    overlap = 80
    for marker in markers:
        raw = marker.encode('utf-8') if isinstance(marker, str) else marker
        remaining.append((marker, raw))
        overlap = max(overlap, len(raw))
    #end for
    try:
        with open(path, 'rb') as handle:
            prev = b''
            while True:
                chunk = handle.read(CHUNK_BYTES)
                if not chunk:
                    break
                #end if
                window = prev + chunk
                if extra_re is not None and extra is None:
                    match = extra_re.search(window)
                    if match:
                        extra = match
                    #end if
                #end if
                still = []
                for marker, raw in remaining:
                    if raw in window:
                        if isinstance(marker, str):
                            found.append(marker)
                        else:
                            found.append(marker.decode('utf-8', 'replace'))
                        #end if
                    else:
                        still.append((marker, raw))
                    #end if
                #end for
                remaining = still
                if not remaining and (extra_re is None or extra is not None):
                    break
                #end if
                prev = window[-overlap:]
            #end while
    except (OSError, MemoryError):
        return found, extra
    #end try
    return found, extra
#end def scan_file


def drop_bookkeeping(params):
    out = {}
    for key, value in params.items():
        if key in BOOKKEEPING or str(key).startswith('_'):
            continue
        #end if
        out[key] = value
    #end for
    return out
#end def drop_bookkeeping


class ChainErrorHandler(DynamicChain):
    """Chain recovery: parse a failure and propose one input patch.
    """

    unrecoverable = frozenset()

    def __init__(self, start=None, max_runs=3, handlers=None):
        DynamicChain.__init__(self, max_runs=max_runs)
        self.start    = dict(start or {})
        self.handlers = tuple(self.default_handlers() if handlers is None else handlers)
    #end def __init__

    def default_handlers(self):
        return ()
    #end def default_handlers

    def parse(self, sim):
        self.error('parse() must be implemented in a subclass')
    #end def parse

    def apply_patch(self, inp, patch):
        self.error('apply_patch() must be implemented in a subclass')
    #end def apply_patch

    def params_from_input(self, inp):
        return {}
    #end def params_from_input

    def recover_failed(self, sim):
        """If a patch exists, archive the attempt, apply it, and resubmit.

        Returns True when the simulation should be retried.
        """
        sim.log('recovering failed run'+sim.idstr(), n=3)
        products = self.parse(sim)
        sim.log('parsed products'+str(products), n=3)
        params = self.params_from_input(sim.input)
        params['attempt'] = len(self.history)
        decision = self.observe(params, products)
        sim.log('observed decision'+str(decision), n=3)
        if decision.status != 'continue':
            return False
        #end if
        sim.save_attempt()
        sim.log('saved attempt', n=3)
        self.apply_patch(sim.input, decision.next_params)
        sim.log('applied patch', n=3)
        sim.input.write(os.path.join(sim.locdir, sim.infile))
        sim.log('wrote input files'+sim.idstr(), n=3)
        sim.reset_indicators()
        sim.log('reset indicators', n=3)
        return True
    #end def recover_failed

    def reset(self):
        DynamicMode.reset(self)
    #end def reset

    def initial(self):
        if not self.start:
            self.error(f'{type(self).__name__}.initial needs start params')
        #end if
        params = dict(self.start)
        params.setdefault('attempt', 0)
        return params
    #end def initial

    def propose(self, history):
        last     = history[-1]
        params   = dict(last['params'])
        products = last.get('products') or {}
        errors   = list(products.get('errors') or [])
        if any(tag in self.unrecoverable for tag in errors):
            return None
        #end if
        patch = self._patch(params, products)
        if patch is None:
            return None
        #end if
        nxt = dict(params)
        nxt.update(patch)
        nxt['attempt'] = int(params.get('attempt', 0)) + 1
        return nxt
    #end def propose

    def _patch(self, params, products):
        """Chain: first matching handler."""
        for handler in self.handlers:
            patch = handler(params, products)
            if patch is not None:
                return patch
            #end if
        #end for
        return None
    #end def _patch
#end class ChainErrorHandler
