import numpy as np
import pytest

from ..testing import execute
from . import TEST_DIR, NexusTestOrder

_ = pytest.importorskip("h5py")

pytestmark = pytest.mark.order(NexusTestOrder.ESHDF)

TEST_FILES = {
    "small_archive.h5":      TEST_DIR / "reference/eshdf/small_eshdf.h5",
    "eshdf_write_nk_ref.h5": TEST_DIR / "reference/eshdf/eshdf_write_nk_ref.h5",
    }

ESHDF_EXECUTABLE = TEST_DIR.parent / "bin/eshdf"
E_FERMI = 19.1114

for file in TEST_FILES.values():
    assert(file.exists()), f"Test file {file} does not exist!"

def test_kinetic():
    _ = pytest.importorskip("h5py")

    command = f"{ESHDF_EXECUTABLE} kinetic {TEST_FILES['small_archive.h5']} --Ef={E_FERMI}"
    out, _err, rc = execute(command)

    # Assert that return code is 0
    assert(rc==0)

    ref_output = """\
Number of spins              : 1
Number kpoints               : 4
Number of electrons per spin : 16
Summed orbital norm per spin : 16.
Total kinetic energy         : 13.315652678318026 Ha
Kinetic energy per spin      : 13.31565268 Ha
"""

    assert(out.strip() == ref_output.strip())
#end def test_kinetic


def test_kinetic_orb():
    command = f"{ESHDF_EXECUTABLE} kinetic {TEST_FILES['small_archive.h5']} --Ef={E_FERMI} --orb"
    out, _err, rc = execute(command)

    # Assert that return code is 0
    assert(rc==0)

    ref_output = """\
Number of spins              : 1
Number kpoints               : 4
Number of electrons per spin : 16
Summed orbital norm per spin : 16.
Total kinetic energy         : 13.315652678318026 Ha
Kinetic energy per spin      : 13.31565268 Ha

Per orbital kinetic energies
  Spin up energies
    index kpoint_index  KS eig (eV)  kinetic (Ha)
      0        0         -8.582273     0.079274
      1        2          1.940267     0.507421
      2        1          1.940267     0.507421
      3        3          1.940267     0.507421
      4        3          1.940267     0.507421
      5        2          1.940267     0.507421
      6        1          1.940267     0.507421
      7        1         11.328221     0.980532
      8        3         11.328221     0.980532
      9        2         11.328221     0.980532
     10        3         11.328221     0.980532
     11        1         11.328221     0.980532
     12        2         11.328221     0.980532
     13        0         19.111379     1.436222
     14        0         19.111379     1.436222
     15        0         19.111379     1.436222
"""

    assert(out.strip() == ref_output.strip())
#end def test_kinetic_orb


def test_write_nk(tmp_path):
    import h5py

    outfile = tmp_path / "eshdf_write_nk.h5"
    command = f"{ESHDF_EXECUTABLE} write_nk {TEST_FILES['small_archive.h5']} --Ef={E_FERMI} --outfile={outfile}"
    print(command)
    out, _, rc = execute(command)

    # Assert that return code is 0
    assert(rc == 0)
    assert(f"Writing n(k) to HDF5 {outfile}" in out)
    assert(outfile.exists())

    ref = h5py.File(TEST_FILES["eshdf_write_nk_ref.h5"], mode="r")
    calc = h5py.File(outfile, mode="r")

    assert(calc.keys() == {"data"})

    ref_data  = np.asarray(ref.get("data"), dtype=float)
    calc_data = np.asarray(calc.get("data"), dtype=float)

    np.testing.assert_allclose(ref_data, calc_data)
