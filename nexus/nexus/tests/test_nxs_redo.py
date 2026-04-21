try:
    import pytest
    from . import NexusTestOrder
    pytestmark = pytest.mark.order(NexusTestOrder.NXS_REDO)
    from ..generic import generic_settings
    generic_settings.raise_error = True
except ImportError:
    pass

from .. import testing
from ..testing import execute


def test_redo(tmp_path):

    tmp_dir = tmp_path / "test_nxs_redo_output"
    tmp_dir.mkdir(exist_ok=True)

    exe = testing.executable_path('nxs-redo')

    command = '{} {}'.format(exe,tmp_dir)


    # empty directory
    assert(list(tmp_dir.iterdir())==[])

    out,err,rc = execute(command)

    assert('no simulation directories found' in out)


    # directory w/ files, but not nexus simulation directory
    (tmp_dir / "qmc.in.xml").touch()

    out,err,rc = execute(command)

    assert('no simulation directories found' in out)

    assert(set(tmp_dir.iterdir())==set([tmp_dir / 'qmc.in.xml']))


    # nexus simulation directory
    (tmp_dir / "sim_qmc").mkdir()

    assert(set(tmp_dir.iterdir())==set([tmp_dir / 'qmc.in.xml', tmp_dir / 'sim_qmc']))

    out,err,rc = execute(command)

    assert(set(tmp_dir.iterdir())==set([tmp_dir / 'attempt1']))


    # nexus simulation directory w/ prior attempt
    (tmp_dir / "qmc.in.xml").touch()
    (tmp_dir / "sim_qmc").mkdir()

    assert(set(tmp_dir.iterdir())==set([tmp_dir / 'attempt1',tmp_dir / 'qmc.in.xml',tmp_dir / 'sim_qmc']))

    out,err,rc = execute(command)

    assert(set(tmp_dir.iterdir())==set([tmp_dir / 'attempt1',tmp_dir / 'attempt2']))
#end def test_redo
