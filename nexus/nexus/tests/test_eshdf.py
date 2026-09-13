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
E_FERMI = 13.9291

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
Total kinetic energy         : 17.62702973479747 Ha
Kinetic energy per spin      : 17.62702973 Ha
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
Total kinetic energy         : 17.62702973479747 Ha
Kinetic energy per spin      : 17.62702973 Ha

Per orbital kinetic energies
  Spin up energies
    index kpoint_index  KS eig (eV)  kinetic (Ha)
      0        0         -7.937200     0.088413
      1        1          0.831948     0.654364
      2        3          0.831948     0.654364
      3        2          0.831948     0.654364
      4        1          0.831948     0.654364
      5        2          0.831948     0.654364
      6        3          0.831948     0.654364
      7        1          7.392005     1.333446
      8        2          7.392005     1.333446
      9        3          7.392005     1.333446
     10        3          7.392006     1.333446
     11        2          7.392006     1.333446
     12        1          7.392006     1.333446
     13        0         13.929089     1.870585
     14        0         13.929089     1.870585
     15        0         13.929089     1.870585
"""

    assert(out.strip() == ref_output.strip())
#end def test_kinetic_orb


def test_write_nk(tmp_path):
    import h5py
    _ = pytest.importorskip("tables")

    outfile = tmp_path / "eshdf_write_nk.h5"
    command = f"{ESHDF_EXECUTABLE} write_nk {TEST_FILES['small_archive.h5']} --Ef={E_FERMI} --outfile={outfile}"
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
