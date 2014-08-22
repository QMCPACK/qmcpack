

import os
from generic import obj
from simulation import Simulation,SimulationInput,SimulationAnalyzer


class WfconvertInput(SimulationInput):
    def __init__(self,app_name='wfconvert',h5in='MISSING.h5',h5out='wfconvert.h5',spline=False,format='eshdf',factor=None):
        self.app_name = app_name
        self.h5in = h5in
        self.h5out= h5out
        self.spline = spline
        self.format = format
        self.factor = factor
    #end def __init__

#wfconvert --nospline --eshdf diamond.h5 out/diamond.pwscf.h5 >& diamond-wfconvert.out 
    def set_app_name(self,app_name):
        self.app_name = app_name
    #end def set_app_name

    def app_command(self):
        c = self.app_name+' '
        if not self.spline:
            c+= '--nospline '
        #end if
        c+='--'+self.format+' '+self.h5out+' '+self.h5in
        return c
    #end def app_command
        

    def read(self,filepath):
        None
    #end def read

    def write_contents(self):
        return self.app_command()
    #end def write_contents
#end class WfconvertInput


def generate_wfconvert_input(app_name='wfconvert',h5in='MISSING.h5',h5out='wfconvert.h5',spline=False,format='eshdf',factor=None):
    wi = WfconvertInput(
        app_name = app_name,
        h5in   = h5in,
        h5out  = h5out,
        spline = spline,
        format = format,
        factor = factor
        )
    return wi
#end def generate_wfconvert_input


class WfconvertAnalyzer(SimulationAnalyzer):
    def __init__(self,arg0):
        if isinstance(arg0,Simulation):
            sim = arg0
            self.infile = sim.infile
            self.dir    = sim.locdir
            self.h5file = os.path.join(sim.locdir,sim.input.h5out)
        else:
            self.infile = arg0
        #end if
    #end def __init__

    def analyze(self):
        if False:
            import h5py
            self.log('Fixing h5 file',n=5)
            h = h5py.File(self.h5file)
            if 'electrons' in h:
                elec = h['electrons']
                nkpoints = 0
                for name,val in elec.iteritems():
                    if name.startswith('kpoint'):
                        nkpoints+=1
                    #end for
                #end if
                nkold = elec['number_of_kpoints'][0] 
                self.log('Were',nkold,'kpoints, now',nkpoints,'kpoints',n=6)
                elec['number_of_kpoints'][0] = nkpoints
            #end for        
        #end if
    #end def analyze
#end class WfconvertAnalyzer



class Wfconvert(Simulation):
    input_type             = WfconvertInput
    analyzer_type          = WfconvertAnalyzer
    generic_identifier     = 'wfconvert'
    application            = 'wfconvert'
    application_properties = set(['serial'])
    application_results    = set(['orbitals'])

    def set_app_name(self,app_name):
        self.app_name = app_name
        self.input.set_app_name(app_name)
    #end def set_app_name

    def check_result(self,result_name,sim):
        calculating_result = False
        if result_name=='orbitals':
            calculating_result = True
        else:
            calculating_result = False
            self.error('ability to check for result '+result_name+' has not been implemented')
        #end if        
        return calculating_result
    #end def check_result

    def get_result(self,result_name,sim):
        result = obj()
        if result_name=='orbitals':
            result.h5file   = os.path.join(self.locdir,self.input.h5out)
            result.outfile  = os.path.join(self.locdir,self.outfile)
        else:
            self.error('ability to get result '+result_name+' has not been implemented')
        #end if        
        return result
    #end def get_result

    def incorporate_result(self,result_name,result,sim):
        if result_name=='orbitals':
            self.input.h5in = os.path.relpath(result.h5file,self.locdir)
            self.job.app_command = self.input.app_command()
        else:
            self.error('ability to incorporate result '+result_name+' has not been implemented')
        #end if                
    #end def incorporate_result

    def check_sim_status(self):
        outfile = os.path.join(self.locdir,self.outfile)
        errfile = os.path.join(self.locdir,self.errfile)
        fobj = open(outfile,'r')
        output = fobj.read()
        fobj.close()
        fobj = open(errfile,'r')
        errors = fobj.read()
        fobj.close()
        h5file = os.path.join(self.locdir,self.input.h5out)
        file_exists = os.path.exists(h5file)
        outfin = 'Successfully read' in errors and 'numSpins' in errors
        outfin = outfin and 'Writing laplacians' in output

        success = file_exists and outfin

        self.finished = success
    #end def check_sim_status

    def get_output_files(self):
        output_files = []
        return output_files
    #end def get_output_files

    def app_command(self):
        # app_name is passed along in post_init
        return self.input.app_command()
    #end def app_command
#end class Wfconvert









def generate_wfconvert(**kwargs):
    sim_args,inp_args = Simulation.separate_inputs(kwargs)

    if not 'input' in sim_args:
        sim_args.input = generate_wfconvert_input(**inp_args)
    #end if
    wfconvert = Wfconvert(**sim_args)

    return wfconvert
#end def generate_wfconvert
