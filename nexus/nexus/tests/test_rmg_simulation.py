import pytest

from ..machines import job
from ..rmg import Rmg, generate_rmg
from ..testing import clear_all_sims
from . import NexusTestOrder

pytestmark = pytest.mark.order(NexusTestOrder.RMG_SIMULATION)

def test_minimal_init():
    sim = generate_rmg(
        job    = job(machine='ws1',cores=1),
        )

    assert(isinstance(sim,Rmg))

    clear_all_sims()
#end def test_minimal_init

