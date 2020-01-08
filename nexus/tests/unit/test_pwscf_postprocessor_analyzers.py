
import testing
from testing import value_eq,object_eq,text_eq



def test_import():
    from pwscf_postprocessors import PPAnalyzer
    from pwscf_postprocessors import DosAnalyzer
    from pwscf_postprocessors import BandsAnalyzer
    from pwscf_postprocessors import ProjwfcAnalyzer
    from pwscf_postprocessors import CpppAnalyzer
    from pwscf_postprocessors import PwexportAnalyzer
#end def test_import



def test_empty_init():
    from pwscf_postprocessors import PPAnalyzer
    from pwscf_postprocessors import DosAnalyzer
    from pwscf_postprocessors import BandsAnalyzer
    from pwscf_postprocessors import ProjwfcAnalyzer
    from pwscf_postprocessors import CpppAnalyzer
    from pwscf_postprocessors import PwexportAnalyzer

    pa = PPAnalyzer(None)
    pa = DosAnalyzer(None)
    pa = BandsAnalyzer(None)
    pa = ProjwfcAnalyzer(None)
    pa = CpppAnalyzer(None)
    pa = PwexportAnalyzer(None)
#end def test_empty_init



def test_projwfc_analyzer():
    import os
    from generic import obj
    from pwscf_postprocessors import ProjwfcAnalyzer

    tpath = testing.setup_unit_test_output_directory(
        test      = 'pwscf_postprocessor_analyzers',
        subtest   = 'test_projwfc_analyzer',
        file_sets = ['pwf.in','pwf.out'],
        )

    projwfc_in = os.path.join(tpath,'pwf.in')

    pa = ProjwfcAnalyzer(projwfc_in)

    del pa.info.path

    pa_ref = obj(
        info = obj(
            infile          = 'pwf.in',
            initialized     = True,
            outfile         = 'pwf.out',
            strict          = False,
            warn            = False,
            ),
        input = obj(
            projwfc = obj(
                lwrite_overlaps = True,
                outdir          = 'pwscf_output',
                prefix          = 'pwscf',
                ),
            ),
        )

    assert(object_eq(pa.to_obj(),pa_ref))


    pa = ProjwfcAnalyzer(projwfc_in,analyze=True)

    del pa.info.path

    pa_ref = obj(
        info = obj(
            infile      = 'pwf.in',
            initialized = True,
            outfile     = 'pwf.out',
            strict      = False,
            warn        = False,
            ),
        input = obj(
            projwfc = obj(
                lwrite_overlaps = True,
                outdir          = 'pwscf_output',
                prefix          = 'pwscf',
                ),
            ),
        lowdin = obj({
            0 : obj(
                down = obj({
                    'charge' : 1.9988,
                    'd'      : 0.0,
                    'dx2-y2' : 0.0,
                    'dxy'    : 0.0,
                    'dxz'    : 0.0,
                    'dyz'    : 0.0,
                    'dz2'    : 0.0,
                    'p'      : 0.999,
                    'px'     : 0.3318,
                    'py'     : 0.3336,
                    'pz'     : 0.3336,
                    's'      : 0.9998,
                    }),
                pol = obj({
                    'charge' : 2.0001,
                    'd'      : 0.0,
                    'dx2-y2' : 0.0,
                    'dxy'    : 0.0,
                    'dxz'    : 0.0,
                    'dyz'    : 0.0,
                    'dz2'    : 0.0,
                    'p'      : 1.9999,
                    'px'     : 0.6678,
                    'py'     : 0.666,
                    'pz'     : 0.666,
                    's'      : 0.0001,
                    }),
                tot = obj({
                    'charge' : 5.9977,
                    'd'      : 0.0,
                    'dx2-y2' : 0.0,
                    'dxy'    : 0.0,
                    'dxz'    : 0.0,
                    'dyz'    : 0.0,
                    'dz2'    : 0.0,
                    'p'      : 3.9979,
                    'px'     : 1.3314,
                    'py'     : 1.3332,
                    'pz'     : 1.3332,
                    's'      : 1.9997,
                    }),
                up = obj({
                    'charge' : 3.9989,
                    'd'      : 0.0,
                    'dx2-y2' : 0.0,
                    'dxy'    : 0.0,
                    'dxz'    : 0.0,
                    'dyz'    : 0.0,
                    'dz2'    : 0.0,
                    'p'      : 2.9989,
                    'px'     : 0.9996,
                    'py'     : 0.9996,
                    'pz'     : 0.9996,
                    's'      : 0.9999,
                    }),
                )
            }),
        states = obj(
            elem    = ['S'],
            nstates = 9,
            ),
        )

    assert(object_eq(pa.to_obj(),pa_ref))


    lowdin_file = os.path.join(tpath,'pwf.lowdin')

    pa.write_lowdin(lowdin_file)

    text = open(lowdin_file,'r').read()

    text_ref = '''
        nup+ndn = 5.9977
        nup-ndn = 2.0001
        
        tot
            0   S   6.00  s( 2.00)p( 4.00)d( 0.00)
        
        pol
            0   S   2.00  s( 0.00)p( 2.00)d( 0.00)
        
        up
            0   S   4.00  s( 1.00)p( 3.00)d( 0.00)
        
        down
            0   S   2.00  s( 1.00)p( 1.00)d( 0.00)
        '''

    def process_text(t):
        return t.replace('(',' ( ').replace(')',' ) ')
    #end def process_text

    text     = process_text(text)
    text_ref = process_text(text_ref)

    assert(text_eq(text,text_ref))


    lowdin_file = os.path.join(tpath,'pwf.lowdin_long')

    pa.write_lowdin(lowdin_file,long=True)

    text = open(lowdin_file,'r').read()

    text_ref = '''
        nup+ndn = 5.9977
        nup-ndn = 2.0001
        
        tot
            0   S   6.00  s( 2.00)px( 1.33)py( 1.33)pz( 1.33)dx2-y2( 0.00)dxy( 0.00)dxz( 0.00)dyz( 0.00)dz2( 0.00)
        
        pol
            0   S   2.00  s( 0.00)px( 0.67)py( 0.67)pz( 0.67)dx2-y2( 0.00)dxy( 0.00)dxz( 0.00)dyz( 0.00)dz2( 0.00)
        
        up
            0   S   4.00  s( 1.00)px( 1.00)py( 1.00)pz( 1.00)dx2-y2( 0.00)dxy( 0.00)dxz( 0.00)dyz( 0.00)dz2( 0.00)
        
        down
            0   S   2.00  s( 1.00)px( 0.33)py( 0.33)pz( 0.33)dx2-y2( 0.00)dxy( 0.00)dxz( 0.00)dyz( 0.00)dz2( 0.00)
        '''

    text     = process_text(text)
    text_ref = process_text(text_ref)

    assert(text_eq(text,text_ref))
#end def test_projwfc_analyzer
