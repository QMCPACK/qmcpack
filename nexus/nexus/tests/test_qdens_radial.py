import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.QDENS_RADIAL)

from ..generic import generic_settings
generic_settings.raise_error = True

from . import TEST_DIR
from ..testing import execute,text_eq,check_value_eq


def test_radial_density():
    _ = pytest.importorskip("spglib")
    import os

    exe = TEST_DIR.parent / "bin/qdens-radial"

    vmc_path = TEST_DIR / "test_qdens_radial_files/diamond_twist/vmc"

    vmc_files_bef = (
        vmc_path / "vmc.avg.s000.SpinDensity_u+d+err.xsf",
        vmc_path / "vmc.avg.s000.SpinDensity_u+d-err.xsf",
        vmc_path / "vmc.avg.s000.SpinDensity_u+d.xsf",
        vmc_path / "vmc.g000.s000.scalar.dat",
        vmc_path / "vmc.g001.s000.scalar.dat",
        vmc_path / "vmc.g002.s000.scalar.dat",
        vmc_path / "vmc.g003.s000.scalar.dat",
        vmc_path / "vmc.g000.twistnum_0.in.g000.qmc",
        vmc_path / "vmc.g001.twistnum_1.in.g001.qmc",
        vmc_path / "vmc.g002.twistnum_2.in.g002.qmc",
        vmc_path / "vmc.g003.twistnum_3.in.g003.qmc",
        vmc_path / "vmc.g000.twistnum_0.in.xml",
        vmc_path / "vmc.g001.twistnum_1.in.xml",
        vmc_path / "vmc.g002.twistnum_2.in.xml",
        vmc_path / "vmc.g003.twistnum_3.in.xml",
        vmc_path / "vmc.in",
        vmc_path / "vmc.out",
        vmc_path / "vmc.info.xml",
    )

    assert(check_value_eq(set(vmc_path.iterdir()),set(vmc_files_bef)))

    dmc_path = TEST_DIR / "test_qdens_radial_files/diamond_twist/dmc"

    dmc_files_bef = (
        dmc_path / "dmc.avg.s001.SpinDensity_u+d+err.xsf",
        dmc_path / "dmc.avg.s001.SpinDensity_u+d-err.xsf",
        dmc_path / "dmc.avg.s001.SpinDensity_u+d.xsf",
        dmc_path / "dmc.g000.s000.scalar.dat",
        dmc_path / "dmc.g000.s001.scalar.dat",
        dmc_path / "dmc.g001.s000.scalar.dat",
        dmc_path / "dmc.g001.s001.scalar.dat",
        dmc_path / "dmc.g002.s000.scalar.dat",
        dmc_path / "dmc.g002.s001.scalar.dat",
        dmc_path / "dmc.g003.s000.scalar.dat",
        dmc_path / "dmc.g003.s001.scalar.dat",
        dmc_path / "dmc.g000.twistnum_0.in.g000.qmc",
        dmc_path / "dmc.g001.twistnum_1.in.g001.qmc",
        dmc_path / "dmc.g002.twistnum_2.in.g002.qmc",
        dmc_path / "dmc.g003.twistnum_3.in.g003.qmc",
        dmc_path / "dmc.g000.twistnum_0.in.xml",
        dmc_path / "dmc.g001.twistnum_1.in.xml",
        dmc_path / "dmc.g002.twistnum_2.in.xml",
        dmc_path / "dmc.g003.twistnum_3.in.xml",
        dmc_path / "dmc.in",
        dmc_path / "dmc.out",
    )

    assert(check_value_eq(set(dmc_path.iterdir()),set(dmc_files_bef)))

        # VMC non-cumulative
    command = f'{exe} -s C -r 0.6 {vmc_path}/vmc.avg.s000.SpinDensity_u+d.xsf'
    out,err,rc = execute(command)

    # Assert that return code is 0
    assert(rc==0)

    # Assert that output is consistent with reference
    out_ref = '''
Norm:   tot             = 8.00000003

Non-Cumulative Value of C Species at Cutoff 0.6 is: 4.77326491
'''
    assert(text_eq(out,out_ref))

        # VMC cumulative
    command = f'{exe} -s C -r 0.6 -c {vmc_path}/vmc.avg.s000.SpinDensity_u+d.xsf'
    out,err,rc = execute(command)

    # Assert that return code is 0
    assert(rc==0)

    # Assert that output is consistent with reference
    out_ref = '''
Norm:   tot             = 8.00000003

Cumulative Value of C Species at Cutoff 0.6 is: 1.06842486
'''

        # DMC extrapolated non-cumulative
    command = f'{exe} -s C -r 0.6 --vmc={vmc_path}/vmc.avg.s000.SpinDensity_u+d.xsf {dmc_path}/dmc.avg.s001.SpinDensity_u+d.xsf'
    out,err,rc = execute(command)

    # Assert that return code is 0
    assert(rc==0)

    # Assert that output is consistent with reference
    out_ref = '''
Extrapolating from VMC and DMC densities...

Norm:   tot             = 7.999999969999998

Non-Cumulative Value of C Species at Cutoff 0.6 is: 4.91009396
'''
    assert(text_eq(out,out_ref))

        # DMC extrapolated cumulative
    command = f'{exe} -s C -r 0.6 -c --vmc={vmc_path}/vmc.avg.s000.SpinDensity_u+d.xsf {dmc_path}/dmc.avg.s001.SpinDensity_u+d.xsf'
    out,err,rc = execute(command)

    # Assert that return code is 0
    assert(rc==0)

    # Assert that output is consistent with reference
    out_ref = '''
Extrapolating from VMC and DMC densities...

Norm:   tot             = 7.999999969999998

Cumulative Value of C Species at Cutoff 0.6 is: 1.10782675
'''
    assert(text_eq(out,out_ref))


        # DMC extrapolated cumulative with error bar
    command = f'{exe} -s C -r 0.6 -c -n 3 --seed=0 --vmc={vmc_path}/vmc.avg.s000.SpinDensity_u+d.xsf --vmcerr={vmc_path}/vmc.avg.s000.SpinDensity_u+d+err.xsf --dmcerr={dmc_path}/dmc.avg.s001.SpinDensity_u+d+err.xsf {dmc_path}/dmc.avg.s001.SpinDensity_u+d.xsf'
    out,err,rc = execute(command)

    # Assert that return code is 0
    assert(rc==0)

    # Assert that output is consistent with reference
    out_ref = '''
Extrapolating from VMC and DMC densities...
Resampling to obtain error bar (NOTE: This can be slow)...
Will compute 3 samples...
sample: 0
sample: 1
sample: 2

Norm:   tot             = 7.999999969999998

Cumulative Value of C Species at Cutoff 0.6 is: 1.10782675 +/- 0.00160665
'''
    assert(text_eq(out,out_ref,atol=1e-8))
#end def test_radial_density
