import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.NXS_REDO)

from ..generic import generic_settings
generic_settings.raise_error = True

from . import TEST_DIR
from ..testing import execute


def test_redo(tmp_path):

    exe = TEST_DIR.parent / "bin/nxs-redo"

    command = '{} {}'.format(exe,tmp_path)


    # empty directory
    assert(list(tmp_path.iterdir())==[])

    out,err,rc = execute(command)

    assert('no simulation directories found' in out)


    # directory w/ files, but not nexus simulation directory
    (tmp_path / "qmc.in.xml").touch()

    out,err,rc = execute(command)

    assert('no simulation directories found' in out)

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
