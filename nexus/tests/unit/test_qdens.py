
import versions
import testing
from testing import execute,text_eq,check_value_eq


if versions.h5py_available:
    def test_density():
        import os

        tpath = testing.setup_unit_test_output_directory('qdens','test_density')

        exe = testing.executable_path('qdens')

        qa_files_path = testing.unit_test_file_path('qmcpack_analyzer','diamond_gamma/dmc')
        command = 'rsync -a {} {}'.format(qa_files_path,tpath)
        out,err,rc = execute(command)
        assert(rc==0)
        dmc_path = os.path.join(tpath,'dmc')
        dmc_infile = os.path.join(dmc_path,'dmc.in.xml')
        assert(os.path.exists(dmc_infile))

        files_bef = '''
    dmc.err     dmc.s000.scalar.dat  dmc.s001.stat.h5     dmc.s003.scalar.dat
    dmc.in.xml  dmc.s000.stat.h5     dmc.s002.scalar.dat  dmc.s003.stat.h5
    dmc.out     dmc.s001.scalar.dat  dmc.s002.stat.h5
        '''.split()

        assert(check_value_eq(set(os.listdir(dmc_path)),set(files_bef)))

        command = '{0} -v -e 4 -f xsf -i {1}/dmc.in.xml {1}/*stat.h5'.format(exe,dmc_path)

        out,err,rc = execute(command)

        files_aft = '''
            dmc.err                           dmc.s001.stat.h5
            dmc.in.xml                        dmc.s002.scalar.dat
            dmc.out                           dmc.s002.SpinDensity_d-err.xsf
            dmc.s000.scalar.dat               dmc.s002.SpinDensity_d+err.xsf
            dmc.s000.SpinDensity_d-err.xsf    dmc.s002.SpinDensity_d.xsf
            dmc.s000.SpinDensity_d+err.xsf    dmc.s002.SpinDensity_u-d-err.xsf
            dmc.s000.SpinDensity_d.xsf        dmc.s002.SpinDensity_u-d+err.xsf
            dmc.s000.SpinDensity_u-d-err.xsf  dmc.s002.SpinDensity_u+d-err.xsf
            dmc.s000.SpinDensity_u-d+err.xsf  dmc.s002.SpinDensity_u+d+err.xsf
            dmc.s000.SpinDensity_u+d-err.xsf  dmc.s002.SpinDensity_u-d.xsf
            dmc.s000.SpinDensity_u+d+err.xsf  dmc.s002.SpinDensity_u+d.xsf
            dmc.s000.SpinDensity_u-d.xsf      dmc.s002.SpinDensity_u-err.xsf
            dmc.s000.SpinDensity_u+d.xsf      dmc.s002.SpinDensity_u+err.xsf
            dmc.s000.SpinDensity_u-err.xsf    dmc.s002.SpinDensity_u.xsf
            dmc.s000.SpinDensity_u+err.xsf    dmc.s002.stat.h5
            dmc.s000.SpinDensity_u.xsf        dmc.s003.scalar.dat
            dmc.s000.stat.h5                  dmc.s003.SpinDensity_d-err.xsf
            dmc.s001.scalar.dat               dmc.s003.SpinDensity_d+err.xsf
            dmc.s001.SpinDensity_d-err.xsf    dmc.s003.SpinDensity_d.xsf
            dmc.s001.SpinDensity_d+err.xsf    dmc.s003.SpinDensity_u-d-err.xsf
            dmc.s001.SpinDensity_d.xsf        dmc.s003.SpinDensity_u-d+err.xsf
            dmc.s001.SpinDensity_u-d-err.xsf  dmc.s003.SpinDensity_u+d-err.xsf
            dmc.s001.SpinDensity_u-d+err.xsf  dmc.s003.SpinDensity_u+d+err.xsf
            dmc.s001.SpinDensity_u+d-err.xsf  dmc.s003.SpinDensity_u-d.xsf
            dmc.s001.SpinDensity_u+d+err.xsf  dmc.s003.SpinDensity_u+d.xsf
            dmc.s001.SpinDensity_u-d.xsf      dmc.s003.SpinDensity_u-err.xsf
            dmc.s001.SpinDensity_u+d.xsf      dmc.s003.SpinDensity_u+err.xsf
            dmc.s001.SpinDensity_u-err.xsf    dmc.s003.SpinDensity_u.xsf
            dmc.s001.SpinDensity_u+err.xsf    dmc.s003.stat.h5
            dmc.s001.SpinDensity_u.xsf
            '''.split()

        assert(check_value_eq(set(os.listdir(dmc_path)),set(files_aft),verbose=True))

        tot_file = os.path.join(dmc_path,'dmc.s003.SpinDensity_u+d.xsf')
        pol_file = os.path.join(dmc_path,'dmc.s003.SpinDensity_u-d.xsf')

        tot = open(tot_file,'r').read()
        pol = open(pol_file,'r').read()

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
