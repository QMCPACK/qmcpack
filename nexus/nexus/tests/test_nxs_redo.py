import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.NXS_REDO)

from . import TEST_DIR
from ..testing import execute


def test_redo(tmp_path):

    exe = TEST_DIR.parent / "bin/nxs-redo"

    command = f'{exe} {tmp_path}'


    # empty directory
    assert(list(tmp_path.iterdir())==[])

    with pytest.raises(
        AssertionError,
        match="no simulation directories found",
        ):
        out,err,rc = execute(command)


    # directory w/ files, but not nexus simulation directory
    (tmp_path / "qmc.in.xml").touch()

    with pytest.raises(
        AssertionError,
        match="no simulation directories found",
        ):
        out,err,rc = execute(command)

    assert(set(tmp_path.iterdir())==set([tmp_path / 'qmc.in.xml']))


    # nexus simulation directory
    (tmp_path / "sim_qmc").mkdir()

    assert(set(tmp_path.iterdir())==set([tmp_path / 'qmc.in.xml', tmp_path / 'sim_qmc']))

    out,err,rc = execute(command)

    assert(set(tmp_path.iterdir())==set([tmp_path / 'attempt1']))


    # nexus simulation directory w/ prior attempt
    (tmp_path / "qmc.in.xml").touch()
    (tmp_path / "sim_qmc").mkdir()

    assert(set(tmp_path.iterdir())==set([tmp_path / 'attempt1',tmp_path / 'qmc.in.xml',tmp_path / 'sim_qmc']))

    out,err,rc = execute(command)

    assert(set(tmp_path.iterdir())==set([tmp_path / 'attempt1',tmp_path / 'attempt2']))
#end def test_redo
