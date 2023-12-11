
import versions
import testing
from testing import execute,text_eq



if versions.scipy_available:
    def test_fit():
        import os

        tpath = testing.setup_unit_test_output_directory('qmc_fit','test_fit')

        exe = testing.executable_path('qmc-fit')
        
        qa_files_path = testing.unit_test_file_path('qmcpack_analyzer','diamond_gamma/dmc')
        command = 'rsync -a {} {}'.format(qa_files_path,tpath)
        out,err,rc = execute(command)
        assert(rc==0)
        dmc_path = os.path.join(tpath,'dmc')
        dmc_infile = os.path.join(dmc_path,'dmc.in.xml')
        assert(os.path.exists(dmc_infile))

        command = "{} ts --noplot -e 10 -s 1 -t '0.02 0.01 0.005' -f linear {}/*scalar*".format(exe,dmc_path)

        out,err,rc = execute(command)

        out_ref = '''
            fit function  : linear
            fitted formula: (-10.5271 +/- 0.0021) + (-0.28 +/- 0.17)*t
            intercept     : -10.5271 +/- 0.0021  Ha
            '''
        
        def process_text(t):
            return t.replace('(',' ( ').replace(')',' ) ')
        #end def process_text

        out = process_text(out)
        out_ref = process_text(out_ref)

        assert(text_eq(out,out_ref,atol=1e-2,rtol=1e-2))

    #end def test_fit
#end if
