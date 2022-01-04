
import versions
import testing
from testing import execute,text_eq,check_value_eq


if versions.spglib_available:
    def test_radial_density():
        import os

        tpath = testing.setup_unit_test_output_directory('qdens_radial','test_radial_density')

        exe = testing.executable_path('qdens-radial')

        qr_vmc_files_path = testing.unit_test_file_path('qdens_radial','diamond_twist/vmc')
        command = 'rsync -a {} {}'.format(qr_vmc_files_path,tpath)
        out,err,rc = execute(command)
        assert(rc==0)
        vmc_path = os.path.join(tpath,'vmc')
        vmc_infile = os.path.join(vmc_path,'vmc.g000.twistnum_0.in.xml')
        assert(os.path.exists(vmc_infile))
        vmc_infile = os.path.join(vmc_path,'vmc.g001.twistnum_1.in.xml')
        assert(os.path.exists(vmc_infile))
        vmc_infile = os.path.join(vmc_path,'vmc.g002.twistnum_2.in.xml')
        assert(os.path.exists(vmc_infile))
        vmc_infile = os.path.join(vmc_path,'vmc.g003.twistnum_3.in.xml')
        assert(os.path.exists(vmc_infile))

        files_bef = '''
        vmc.avg.s000.SpinDensity_u+d+err.xsf  vmc.avg.s000.SpinDensity_u+d-err.xsf
        vmc.avg.s000.SpinDensity_u+d.xsf      vmc.g000.twistnum_0.in.xml
        vmc.g000.s000.scalar.dat              vmc.g001.s000.scalar.dat
        vmc.g000.twistnum_0.in.g000.qmc       vmc.g001.twistnum_1.in.g001.qmc
        vmc.g001.twistnum_1.in.xml            vmc.g002.twistnum_2.in.xml
        vmc.g002.s000.scalar.dat              vmc.g003.s000.scalar.dat
        vmc.g002.twistnum_2.in.g002.qmc       vmc.g003.twistnum_3.in.g003.qmc
        vmc.g003.twistnum_3.in.xml            vmc.out
        vmc.in                                vmc.info.xml
        '''.split()

        qr_dmc_files_path = testing.unit_test_file_path('qdens_radial','diamond_twist/dmc')
        command = 'rsync -a {} {}'.format(qr_dmc_files_path,tpath)
        out,err,rc = execute(command)
        assert(rc==0)
        dmc_path = os.path.join(tpath,'dmc')
        dmc_infile = os.path.join(dmc_path,'dmc.g000.twistnum_0.in.xml')
        assert(os.path.exists(dmc_infile))
        dmc_infile = os.path.join(dmc_path,'dmc.g001.twistnum_1.in.xml')
        assert(os.path.exists(dmc_infile))
        dmc_infile = os.path.join(dmc_path,'dmc.g002.twistnum_2.in.xml')
        assert(os.path.exists(dmc_infile))
        dmc_infile = os.path.join(dmc_path,'dmc.g003.twistnum_3.in.xml')
        assert(os.path.exists(dmc_infile))

        files_bef = '''
        dmc.avg.s001.SpinDensity_u+d+err.xsf  dmc.avg.s001.SpinDensity_u+d-err.xsf
        dmc.avg.s001.SpinDensity_u+d.xsf      dmc.g001.s001.scalar.dat                               
        dmc.g000.s000.scalar.dat              dmc.g000.s001.scalar.dat
        dmc.g000.twistnum_0.in.g000.qmc       dmc.g001.twistnum_1.in.g001.qmc
        dmc.g000.twistnum_0.in.xml            dmc.g001.twistnum_1.in.xml     
        dmc.g001.s000.scalar.dat              dmc.g002.s000.scalar.dat       
        dmc.g002.s001.scalar.dat              dmc.g003.s001.scalar.dat       
        dmc.g002.twistnum_2.in.g002.qmc       dmc.g003.twistnum_3.in.g003.qmc
        dmc.g002.twistnum_2.in.xml            dmc.g003.twistnum_3.in.xml
        dmc.g003.s000.scalar.dat              dmc.in
        dmc.out
        '''.split()

        assert(check_value_eq(set(os.listdir(dmc_path)),set(files_bef)))

        # VMC non-cumulative
        command = '{0} -s C -r 0.6 {1}/vmc.avg.s000.SpinDensity_u+d.xsf'.format(exe,vmc_path)
        out,err,rc = execute(command)

        # Assert that return code is 0
        assert(rc==0)

        # Assert that output is consistent with reference
        out_ref = '''
Norm:   tot             = 8.00000003

Non-Cumulative Value of C Species at Cutoff 0.6 is: 4.773264912162242
'''
        assert(text_eq(out,out_ref))

        # VMC cumulative
        command = '{0} -s C -r 0.6 -c {1}/vmc.avg.s000.SpinDensity_u+d.xsf'.format(exe,vmc_path)
        out,err,rc = execute(command)

        # Assert that return code is 0
        assert(rc==0)

        # Assert that output is consistent with reference
        out_ref = '''
Norm:   tot             = 8.00000003

Cumulative Value of C Species at Cutoff 0.6 is: 1.0684248641866259
'''

        # DMC extrapolated non-cumulative
        command = '{0} -s C -r 0.6 --vmc={1}/vmc.avg.s000.SpinDensity_u+d.xsf {2}/dmc.avg.s001.SpinDensity_u+d.xsf'.format(exe,vmc_path,dmc_path)
        out,err,rc = execute(command)

        # Assert that return code is 0
        assert(rc==0)

        # Assert that output is consistent with reference
        out_ref = '''
Extrapolating from VMC and DMC densities...

Norm:   tot             = 7.999999969999998

Non-Cumulative Value of C Species at Cutoff 0.6 is: 4.910093964664309
'''
        assert(text_eq(out,out_ref))

        # DMC extrapolated cumulative
        command = '{0} -s C -r 0.6 -c --vmc={1}/vmc.avg.s000.SpinDensity_u+d.xsf {2}/dmc.avg.s001.SpinDensity_u+d.xsf'.format(exe,vmc_path,dmc_path)
        out,err,rc = execute(command)

        # Assert that return code is 0
        assert(rc==0)

        # Assert that output is consistent with reference
        out_ref = '''
Extrapolating from VMC and DMC densities...

Norm:   tot             = 7.999999969999998

Cumulative Value of C Species at Cutoff 0.6 is: 1.1078267486386275
'''
        assert(text_eq(out,out_ref))


        # DMC extrapolated cumulative with error bar
        command = '{0} -s C -r 0.6 -c -n 3 --seed=0 --vmc={1}/vmc.avg.s000.SpinDensity_u+d.xsf --vmcerr={1}/vmc.avg.s000.SpinDensity_u+d+err.xsf --dmcerr={2}/dmc.avg.s001.SpinDensity_u+d+err.xsf {2}/dmc.avg.s001.SpinDensity_u+d.xsf'.format(exe,vmc_path,dmc_path)
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

Cumulative Value of C Species at Cutoff 0.6 is: 1.1078267486386275+/-0.0016066467833404942
'''
        assert(text_eq(out,out_ref))
    #end def test_radial_density
#end if
