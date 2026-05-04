import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.QDENS)

from ..generic import generic_settings
generic_settings.raise_error = True


from . import TEST_DIR
import shutil
from .. import versions
from ..testing import execute, text_eq


if versions.h5py_available:
    def test_density(tmp_path):

        exe = TEST_DIR.parent / "bin/qdens"

        qa_files_path = TEST_DIR / "test_qmcpack_analyzer_files/diamond_gamma"
        shutil.copytree(qa_files_path, tmp_path, dirs_exist_ok=True)

        dmc_path = tmp_path / 'dmc'

        files_bef = (
            dmc_path / "dmc.s000.stat.h5",
            dmc_path / "dmc.s000.scalar.dat",
            dmc_path / "dmc.s001.stat.h5",
            dmc_path / "dmc.s001.scalar.dat",
            dmc_path / "dmc.s002.stat.h5",
            dmc_path / "dmc.s002.scalar.dat",
            dmc_path / "dmc.s003.stat.h5",
            dmc_path / "dmc.s003.scalar.dat",
            dmc_path / "dmc.in.xml",
            dmc_path / "dmc.out",
            dmc_path / "dmc.err",
        )

        assert(set(dmc_path.iterdir()) == set(files_bef))

        command = f'{exe} -v -e 4 -f xsf -i {dmc_path}/dmc.in.xml {dmc_path}/*stat.h5'

        out,err,rc = execute(command)

        files_aft = (
            tmp_path / "dmc/dmc.err",
            tmp_path / "dmc/dmc.in.xml",
            tmp_path / "dmc/dmc.out",
            tmp_path / "dmc/dmc.s000.scalar.dat",
            tmp_path / "dmc/dmc.s000.SpinDensity_d+err.xsf",
            tmp_path / "dmc/dmc.s000.SpinDensity_d-err.xsf",
            tmp_path / "dmc/dmc.s000.SpinDensity_d.xsf",
            tmp_path / "dmc/dmc.s000.SpinDensity_u+d+err.xsf",
            tmp_path / "dmc/dmc.s000.SpinDensity_u+d-err.xsf",
            tmp_path / "dmc/dmc.s000.SpinDensity_u-d+err.xsf",
            tmp_path / "dmc/dmc.s000.SpinDensity_u-d-err.xsf",
            tmp_path / "dmc/dmc.s000.SpinDensity_u+d.xsf",
            tmp_path / "dmc/dmc.s000.SpinDensity_u-d.xsf",
            tmp_path / "dmc/dmc.s000.SpinDensity_u+err.xsf",
            tmp_path / "dmc/dmc.s000.SpinDensity_u-err.xsf",
            tmp_path / "dmc/dmc.s000.SpinDensity_u.xsf",
            tmp_path / "dmc/dmc.s000.stat.h5",
            tmp_path / "dmc/dmc.s001.scalar.dat",
            tmp_path / "dmc/dmc.s001.SpinDensity_d+err.xsf",
            tmp_path / "dmc/dmc.s001.SpinDensity_d-err.xsf",
            tmp_path / "dmc/dmc.s001.SpinDensity_d.xsf",
            tmp_path / "dmc/dmc.s001.SpinDensity_u+d+err.xsf",
            tmp_path / "dmc/dmc.s001.SpinDensity_u+d-err.xsf",
            tmp_path / "dmc/dmc.s001.SpinDensity_u-d+err.xsf",
            tmp_path / "dmc/dmc.s001.SpinDensity_u-d-err.xsf",
            tmp_path / "dmc/dmc.s001.SpinDensity_u+d.xsf",
            tmp_path / "dmc/dmc.s001.SpinDensity_u-d.xsf",
            tmp_path / "dmc/dmc.s001.SpinDensity_u+err.xsf",
            tmp_path / "dmc/dmc.s001.SpinDensity_u-err.xsf",
            tmp_path / "dmc/dmc.s001.SpinDensity_u.xsf",
            tmp_path / "dmc/dmc.s001.stat.h5",
            tmp_path / "dmc/dmc.s002.scalar.dat",
            tmp_path / "dmc/dmc.s002.SpinDensity_d+err.xsf",
            tmp_path / "dmc/dmc.s002.SpinDensity_d-err.xsf",
            tmp_path / "dmc/dmc.s002.SpinDensity_d.xsf",
            tmp_path / "dmc/dmc.s002.SpinDensity_u+d+err.xsf",
            tmp_path / "dmc/dmc.s002.SpinDensity_u+d-err.xsf",
            tmp_path / "dmc/dmc.s002.SpinDensity_u-d+err.xsf",
            tmp_path / "dmc/dmc.s002.SpinDensity_u-d-err.xsf",
            tmp_path / "dmc/dmc.s002.SpinDensity_u+d.xsf",
            tmp_path / "dmc/dmc.s002.SpinDensity_u-d.xsf",
            tmp_path / "dmc/dmc.s002.SpinDensity_u+err.xsf",
            tmp_path / "dmc/dmc.s002.SpinDensity_u-err.xsf",
            tmp_path / "dmc/dmc.s002.SpinDensity_u.xsf",
            tmp_path / "dmc/dmc.s002.stat.h5",
            tmp_path / "dmc/dmc.s003.scalar.dat",
            tmp_path / "dmc/dmc.s003.SpinDensity_d+err.xsf",
            tmp_path / "dmc/dmc.s003.SpinDensity_d-err.xsf",
            tmp_path / "dmc/dmc.s003.SpinDensity_d.xsf",
            tmp_path / "dmc/dmc.s003.SpinDensity_u+d+err.xsf",
            tmp_path / "dmc/dmc.s003.SpinDensity_u+d-err.xsf",
            tmp_path / "dmc/dmc.s003.SpinDensity_u-d+err.xsf",
            tmp_path / "dmc/dmc.s003.SpinDensity_u-d-err.xsf",
            tmp_path / "dmc/dmc.s003.SpinDensity_u+d.xsf",
            tmp_path / "dmc/dmc.s003.SpinDensity_u-d.xsf",
            tmp_path / "dmc/dmc.s003.SpinDensity_u+err.xsf",
            tmp_path / "dmc/dmc.s003.SpinDensity_u-err.xsf",
            tmp_path / "dmc/dmc.s003.SpinDensity_u.xsf",
            tmp_path / "dmc/dmc.s003.stat.h5",
        )

        assert(set(dmc_path.iterdir()) == set(files_aft))

        tot_file = dmc_path / 'dmc.s003.SpinDensity_u+d.xsf'
        pol_file = dmc_path / 'dmc.s003.SpinDensity_u-d.xsf'

        tot = tot_file.read_text()
        pol = pol_file.read_text()

        tot_ref = '''
             CRYSTAL
             PRIMVEC 
                 1.78500000   1.78500000   0.00000000
                -0.00000000   1.78500000   1.78500000
                 1.78500000  -0.00000000   1.78500000
             PRIMCOORD 
               2 1
                 6   0.00000000   0.00000000   0.00000000
                 6   0.89250000   0.89250000   0.89250000
             BEGIN_BLOCK_DATAGRID_3D
               density
               BEGIN_DATAGRID_3D_density
                 4 4 4
                 0.59500000   0.59500000   0.59500000
                 1.78500000   1.78500000   0.00000000
                -0.00000000   1.78500000   1.78500000
                 1.78500000  -0.00000000   1.78500000
                   0.73126076   0.62407496   0.51676366   0.73126076
                   0.62575089   0.19225114   0.18686389   0.62575089
                   0.51847569   0.18457799   0.42203355   0.51847569
                   0.73126076   0.62407496   0.51676366   0.73126076
                   0.62659840   0.19325900   0.18422995   0.62659840
                   0.19219866   0.04873728   0.13184395   0.19219866
                   0.18474638   0.13013188   0.10227670   0.18474638
                   0.62659840   0.19325900   0.18422995   0.62659840
                   0.51793019   0.18615766   0.41806405   0.51793019
                   0.18425005   0.13092538   0.10088238   0.18425005
                   0.41967003   0.10133434   0.14471118   0.41967003
                   0.51793019   0.18615766   0.41806405   0.51793019
                   0.73126076   0.62407496   0.51676366   0.73126076
                   0.62575089   0.19225114   0.18686389   0.62575089
                   0.51847569   0.18457799   0.42203355   0.51847569
                   0.73126076   0.62407496   0.51676366   0.73126076
               END_DATAGRID_3D_density
             END_BLOCK_DATAGRID_3D
            '''

        pol_ref = '''
             CRYSTAL
             PRIMVEC 
                 1.78500000   1.78500000   0.00000000
                -0.00000000   1.78500000   1.78500000
                 1.78500000  -0.00000000   1.78500000
             PRIMCOORD 
               2 1
                 6   0.00000000   0.00000000   0.00000000
                 6   0.89250000   0.89250000   0.89250000
             BEGIN_BLOCK_DATAGRID_3D
               density
               BEGIN_DATAGRID_3D_density
                 4 4 4
                 0.59500000   0.59500000   0.59500000
                 1.78500000   1.78500000   0.00000000
                -0.00000000   1.78500000   1.78500000
                 1.78500000  -0.00000000   1.78500000
                   0.00106753   0.00015792  -0.00122859   0.00106753
                  -0.00003402   0.00018762  -0.00051347  -0.00003402
                   0.00154254   0.00067654   0.00073434   0.00154254
                   0.00106753   0.00015792  -0.00122859   0.00106753
                   0.00263956   0.00079744  -0.00118289   0.00263956
                  -0.00039348  -0.00026396  -0.00069392  -0.00039348
                   0.00087929   0.00000719   0.00113934   0.00087929
                   0.00263956   0.00079744  -0.00118289   0.00263956
                  -0.00013655  -0.00041508  -0.00235212  -0.00013655
                   0.00003805  -0.00025962  -0.00133495   0.00003805
                   0.00040692   0.00051699  -0.00198263   0.00040692
                  -0.00013655  -0.00041508  -0.00235212  -0.00013655
                   0.00106753   0.00015792  -0.00122859   0.00106753
                  -0.00003402   0.00018762  -0.00051347  -0.00003402
                   0.00154254   0.00067654   0.00073434   0.00154254
                   0.00106753   0.00015792  -0.00122859   0.00106753
               END_DATAGRID_3D_density
             END_BLOCK_DATAGRID_3D
            '''

        assert(text_eq(tot,tot_ref,atol=1e-7))
        assert(text_eq(pol,pol_ref,atol=1e-7))
    #end def test_density
#end if
