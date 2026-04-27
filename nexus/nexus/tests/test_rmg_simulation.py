try:
    import pytest
    from . import NexusTestOrder
    pytestmark = pytest.mark.order(NexusTestOrder.RMG_SIMULATION)
    from ..generic import generic_settings
    generic_settings.raise_error = True
except ImportError:
    pass

from ..testing import clear_all_sims



def test_minimal_init():
    from ..machines import job
    from ..rmg import Rmg,generate_rmg

    sim = generate_rmg(
        job    = job(machine='ws1',cores=1),
        )

    assert(isinstance(sim,Rmg))

    clear_all_sims()
#end def test_minimal_init

