##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  opium.py                                                          #
#    Nexus interface to the OPIUM pseudpotential generator.          #
#    The interface is not yet complete.                              #
#                                                                    #                                        
#====================================================================#


import os
from generic import obj
from developer import DevBase
from simulation import Simulation,SimulationInput,SimulationAnalyzer


# PLEASE READ THIS
#   
#   depending on what you want to do with a simulation
#     you will have to implement different functions below
#     here are a few use cases and the functions required
#     details about each function are given within the classes
#
#   use cases
#     1) standalone simulation
#         nexus drives this simulation in isolation of others
#           i.e., one performs parameter scans to drive several independent opium runs
#         in this setting, a opium simulation does not provide information to 
#           other simulations (e.g. pseudopotentials to qmcpack) 
#      
#         the input file will be read from a template file
#           and modified to obtain the desired inputs
#         one could also provide the input longhand in python 
#           in a form OpiumInput understands (this depends on your implementation)
#
#         required functions to be implemented:
#           OpiumInput: read_text, write_text
#           Opium:      app_command, check_sim_status
#
#     2) generated standalone simulation
#          as above, but with fully generated input files
#            generate functions provide a short-hand of minimal vars for input
#
#          required functions to be implemented:
#           OpiumInput: read_text, write_text
#           Opium:      app_command, check_sim_status
#           generate_opium_input
#           
#     3) simulation that provides information to subsequent chained simulations
#          as above (with or without #2)
#            other simulations can request and get information about 
#              results produced by this simulation
#            (e.g. relaxed structure data, location of orbital files, etc.)
#            this information is used by the others to populate input files
#
#         required functions to be implemented:
#           OpiumInput: read_text, write_text
#           Opium:      app_command,check_sim_status, 
#                        check_result, get_result
#
#           if required to get needed output information:
#           OpiumAnalyzer: analyze


def readval(s):
    try:
        val = int(s)
    except:
        try:
            val = float(s)
        except:
            val = s
        #end try
    #end try
    return val
#end def readval


class Section(DevBase):
    types = {int:'int',float:'float',str:'str',list:'list',obj:'obj'}
    variables = None
    def __init__(self,list_rep=None,**kwargs):
        if isinstance(list_rep,list):
            self.from_list_rep(list_rep)
        else:
            self.set(**kwargs)
        #end if
        #self.validate()
        #for name,type in self.variables.iteritems():
        #    if type in ('str','int','float'):
        #        val = None
        #    elif type=='list':
        #        val = []
        #    elif type=='obj':
        #        val = obj()
        #    #end if
        #    self[name] = val
        ##end for
    #end def __init__

    def validate(self):
        allowed = set(self.variables.keys())
        if len(allowed)>0: # some have numeric rather than named entries
            present = set(self.keys())
            missing = allowed-present
            invalid = present-allowed
            if len(missing)>0:
                self.error('the following variables are missing: '+str(sorted(missing)))
            #end if
            if len(invalid)>0:
                self.error('invalid variable names encountered\ninvalid variables: {0}\nallowed variables: {1}'.format(sorted(invalid),sorted(allowed)))
            #end if
            for name in sorted(present):
                type = self.variables[name]
                if not isinstance(self[name],type):
                    self.error('type of variable {0} is incorrect\ntype required: {1}\ncontents: {2}'.format(name,self.types[type]),str(self[name]))
                #end if
            #end for
        #end if
    #end def validate

    def read(self,text):
        list_rep = []
        lines = text.splitlines()
        for line in lines:
            tokens = line.split()
            vals = []
            for token in tokens:
                vals.append(readval(token))
            #end for
            list_rep.append(vals)
        #end for
        self.from_list_rep(list_rep)
    #end def read

    def write(self):
        s = '['+self.__class__.__name__+']\n'
        for line_list in self.list_rep():
            for val in line_list:
                s+=' '+str(val)
            #end for
            s+='\n'
        #end for
        s+='\n'
        return s
    #end def write

    def list_rep(self):
        self.not_implemented()
    #end def list_rep

    def from_list_rep(self,list_rep):
        self.not_implemented()
    #end def from_list_rep
#end class Section




class Atom(Section):
    variables = obj(symbol=str,rorbs=int,ref=obj)

    def list_rep(self):
        list_rep = [[self.symbol],[self.rorbs]]
        for index in sorted(self.ref.keys()):
            v = self.ref[index]
            list_rep.append([v.nlm,v.occ,v.eig])
        #end for
        return list_rep
    #end def list_rep

    def from_list_rep(self,list_rep):
        self.symbol = list_rep[0][0]
        self.rorbs  = list_rep[1][0]
        ref = obj()
        n = 0
        for nlm,occ,eig in list_rep[2:]:
            ref[n] = obj(nlm=nlm,occ=occ,eig=eig)
            n+=1
        #end for
        self.ref = ref
    #end def from_list_rep
#end class Atom


class Pseudo(Section):
    variables = obj(porbs=int,rcuts=list,method=str)

    def list_rep(self):
        list_rep = [[self.porbs]+list(self.rcuts),[self.method]]
        return list_rep
    #end def list_rep

    def from_list_rep(self,list_rep):
        self.porbs = list_rep[0][0]
        self.rcuts = list_rep[0][1:]
        self.method = list_rep[1][0]
    #end def from_list_rep
#end class Pseudo


class Optinfo(Section):
    variables = obj(qcuts=list,bessels=list)

    def list_rep(self):
        list_rep = zip(self.qcuts,self.bessels)
        return list_rep
    #end def list_rep

    def from_list_rep(self,list_rep):
        qc = []
        bs = []
        for q,b in list_rep:
            qc.append(q)
            bs.append(b)
        #end for
        self.qcuts   = qc
        self.bessels = bs
    #end def from_list_rep
#end class Optinfo


class XC(Section):
    variables = obj(functional=str)

    def list_rep(self):
        list_rep = [[self.functional]]
        return list_rep
    #end def list_rep

    def from_list_rep(self,list_rep):
        self.functional = list_rep[0][0]
    #end def from_list_rep
#end class XC


class Pcc(Section):
    variables = obj(radius=float,method=str)

    def list_rep(self):
        list_rep = [[self.radius],[self.method]]
        return list_rep
    #end def list_rep

    def from_list_rep(self,list_rep):
        self.radius = list_rep[0][0]
        self.method = list_rep[1][0]
    #end def from_list_rep
#end class Pcc


class Relativity(Section):
    variables = obj(rl=str)

    def list_rep(self):
        list_rep = [[self.rl]]
        return list_rep
    #end def list_rep

    def from_list_rep(self,list_rep):
        self.rl = list_rep[0][0]
    #end def from_list_rep
#end class Relativity


class Grid(Section):
    variables = obj(np=int,a=float,b=float)

    def list_rep(self):
        list_rep = [[self.np,self.a,self.b]]
        return list_rep
    #end def list_rep

    def from_list_rep(self,list_rep):
        np,a,b = list_rep[0]
        self.np = np
        self.a  = a
        self.b  = b
    #end def from_list_rep
#end class Grid


class Tol(Section):
    variables = obj(aetol=float,nltol=float)

    def list_rep(self):
        list_rep = [[self.aetol,self.nltol]]
        return list_rep
    #end def list_rep

    def from_list_rep(self,list_rep):
        ae,nl = list_rep[0]
        self.aetol = ae
        self.nltol = nl
    #end def from_list_rep
#end class Tol


class Configs(Section):
    variables = obj()

    def list_rep(self):
        list_rep = [[len(self)]]
        for index in sorted(self.keys()):
            v = self[index]
            list_rep.append([v.nlm,v.occ,v.eig])
        #end for
        return list_rep
    #end def list_rep

    def from_list_rep(self,list_rep):
        nconfigs = list_rep[0][0]
        orbs = list_rep[1:]
        orbs_per_config = len(orbs)/nconfigs
        n=1
        for nlm,occ,eig in orbs:
            if n%orbs_per_config==0:
                self.append(obj(nlm=nlm,occ=occ,eig=eig))
            #end if
            n+=1
        #end for
    #end def from_list_rep
#end class Configs


class KBdesign(Section):
    variables = obj(local=str,boxes=obj)

    def list_rep(self):
        list_rep = [[self.local]]
        if 'boxes' in self:
            list_rep.append([len(self.boxes)])
            for index in sorted(self.boxes.keys()):
                v = self.boxes[index]
                list_rep.append([v.units,v.rmin,v.rmax,v.depth])
            #end for
        #end if
        return list_rep
    #end def list_rep

    def from_list_rep(self,list_rep):
        self.local  = list_rep[0][0]
        if len(list_rep)>1:
            boxes = obj()
            for units,rmin,rmax,depth in list_rep[2:]:
                boxes.append(obj(units=units,rmin=rmin,rmax=rmax,depth=depth))
            #end for
            self.boxes = boxes
        #end if
    #end def from_list_rep
#end class KBdesign


class Loginfo(Section):
    variables = obj(config=int,radius=float,Emin=float,Emax=float)

    def list_rep(self):
        list_rep = [[self.config],[self.radius,self.Emin,self.Emax]]
        return list_rep
    #end def list_rep

    def from_list_rep(self,list_rep):
        self.config = list_rep[0][0]
        radius,Emin,Emax = list_rep[1]
        self.radius = radius
        self.Emin   = Emin
        self.Emax   = Emax
    #end def from_list_rep
#end class Loginfo



class OpiumInput(SimulationInput):

    method_map = obj(o='optimized',k='kerker',t='tm')
    section_map = obj(atom=Atom,pseudo=Pseudo,optinfo=Optinfo,xc=XC,pcc=Pcc,relativity=Relativity,grid=Grid,tol=Tol,configs=Configs,kbdesign=KBdesign,loginfo=Loginfo)

    section_order = 'Atom Pseudo Optinfo XC Pcc Relativity Grid Tol Configs KBdesign Loginfo'.split()

    def __init__(self,filepath=None,
                 atom       = None,
                 pseudo     = None,
                 optinfo    = None,
                 xc         = None,
                 pcc        = None,
                 relativity = None,
                 grid       = None,
                 tol        = None,
                 configs    = None,
                 kbdesign   = None,
                 loginfo    = None
                 ):        
        if filepath!=None:
            self.read(filepath)
        else:
            inputs = obj(atom=atom,pseudo=pseudo,optinfo=optinfo,xc=xc,pcc=pcc,relativity=relativity,grid=grid,tol=tol,configs=configs,kbdesign=kbdesign,loginfo=loginfo)
            for secname,input in inputs.iteritems():
                section_type = self.section_map[secname]
                if isinstance(input,section_type):
                    self[secname]=input
                elif isinstance(input,(dict,obj)):
                    self[secname] = section_type(**input)
                elif isinstance(input,(list,tuple)):
                    self[secname] = section_type(list(input))
                elif input!=None:
                    self.error('invalid type encountered for {0} input\nexpected types are dict, obj, list, or tuple\nvalue provided: {1}'.format(secname,input))
                #end if
            #end for
        #end if
    #end def __init__


    def read_text(self,contents,filepath=None):
        lines = contents.splitlines()
        sections = obj()
        sec=None
        secname=None
        for line in lines:
            ls = line.strip()
            if len(ls)>0 and not ls.startswith('#'):
                if ls.startswith('[') and ls.endswith(']'):
                    prevsecname = secname
                    secname = ls.strip('[]').lower()
                    if not secname in self.section_map:
                        self.error('cannot read file\n{0} is not a valid section name\nvalid options are: {1}'.format(secname,self.section_order))
                    #end if
                    if sec!=None:
                        sections[prevsecname]=sec
                    #end if
                    sec=''
                elif sec is None:
                    self.error('invalid text encountered: '+line)
                else:
                    sec+=ls+'\n'
                #end if
            #end if
        #end for
        if sec!=None and secname!=None and not secname in sections:
            sections[secname]=sec
        #end if
        for secname,sectext in sections.iteritems():
            section = self.section_map[secname]()
            section.read(sectext)
            self[secname] = section
        #end for
    #end def read_text


    def write_text(self,filepath=None):
        contents = ''
        for secname in self.section_order:
            secname = secname.lower()
            if secname in self:
                contents += self[secname].write()
            #end if
        #end for
        return contents
    #end def write_text
#end class OpiumInput



def generate_opium_input(selector,**kwargs):

    if selector=='basic':
        return generate_basic_opium_input(**kwargs)
    elif selector=='full':
        return generate_full_opium_input(**kwargs)
    else:
        OpiumInput.class_error('selection '+str(selector)+' has not been implemented for opium input generation')
    #end if
#end def generate_opium_input


def generate_basic_opium_input(
    ):
    oi = None
    return oi
#end def generate_basic_opium_input


def generate_full_opium_input(
    atom       = None,
    pseudo     = None,
    optinfo    = None,
    xc         = None,
    pcc        = None,
    relativity = None,
    grid       = None,
    tol        = None,
    configs    = None,
    kbdesign   = None,
    loginfo    = None
    ):
    
    oi = OpiumInput(
        atom       = atom      ,
        pseudo     = pseudo    ,
        optinfo    = optinfo   ,
        xc         = xc        ,
        pcc        = pcc       ,
        relativity = relativity,
        grid       = grid      ,
        tol        = tol       ,
        configs    = configs   ,
        kbdesign   = kbdesign  ,
        loginfo    = loginfo   
        )

    return oi
#end def generate_full_opium_input



class OpiumAnalyzer(SimulationAnalyzer):
    def __init__(self,arg0=None):
        # optional
        #  only necessary if you want to use results from output files
        #   to inform the inputs of subsequent simulations
        #   or if you want to have a general purpose class to scrape
        #   and process simulation data
        # below is a reasonable default implementation
        # if you don't want to implement it, just uncomment the following line
        #return

        self.path  = None
        self.input = None
        infile = None
        if isinstance(arg0,Simulation):
            sim = arg0
            infile = os.path.join(sim.locdir,sim.infile)
        else:
            infile = arg0
        #end if
        if infile!=None:
            self.path = os.path.dirname(infile)
            self.input = OpiumInput(infile)
        #end if
    #end def __init__


    def analyze(self):
        # optional
        #  only necessary if you want to use results from output files
        #   to inform the inputs of subsequent simulations
        #   or if you want to have a general purpose class to scrape
        #   and process simulation data
        # if you don't want to implement it, no action is required
        None
    #end def analyze
#end class OpiumAnalyzer



class Opium(Simulation):
    input_type         = OpiumInput
    analyzer_type      = OpiumAnalyzer
    generic_identifier = 'opium'
    application        = 'opium_exe' #replace with default name of opium executable
    application_properties = set(['serial','mpi'])
    application_results    = set(['orbitals']) #what opium produces that other simulations can use

    def check_result(self,result_name,sim):
        # optional
        #  only necessary if another simulation depends on this one
        #  e.g.
        #  other_sim.depends(opium_sim,'pseudopotential')  or similar
        # if you don't want to implement it, uncomment the line below
        #return False
        calculating_result = False
        input = self.input # a OpiumInput object
        # check the input to see if result is being calculated
        #  (e.g. result_name='pseudopotential')
        return calculating_result
    #end def check_result


    def get_result(self,result_name,sim):
        # optional
        #  only necessary if another simulation depends on this one
        #  e.g.
        #  other_sim.depends(opium_sim,'pseudopotential')  or similar
        # if you don't want to implement it, uncomment the line below
        #self.not_implemented()

        result = obj()
        input    = self.input
        #analyzer = self.load_analyzer_image()

        # package information about a result/product in the result object
        # for example, if pseudopotentials are requested, 
        # the path to the pseudopotential file might be provided:
        # result.pseudopotential_file = '/path/to/pseudopotential/file'
        return result
    #end def get_result


    def app_command(self):
        # required
        #  specify command line arguments to the executable, such as the input file
        #    e.g. command_line_args = ' '+self.infile
        command_line_args = ''
        return self.app_name + command_line_args
    #end def app_command


    def check_sim_status(self):
        # required
        #  read output/error files to check whether simulation has
        #    completed successfully
        #  one could also check whether all output files exist
        output = open(os.path.join(self.locdir,self.outfile),'r').read()
        errors = open(os.path.join(self.locdir,self.errfile),'r').read()
        
        success = False
        # check output and errors
        #  set success=True if run completed successfully

        self.finished = success and self.job.finished
    #end def check_sim_status


    def get_output_files(self):
        # optional
        #  if provided, the listed output files will be copied to the results directory
        # if you don't want to implement it, no action is required
        output_files = []
        return output_files
    #end def get_output_files
#end class Opium



def generate_opium(**kwargs):
    has_input = 'input_type' in kwargs
    if has_input:
        input_type = kwargs['input_type']
        del kwargs['input_type']
    else:
        input_type = 'basic'
    #end if    
    kw = set(kwargs.keys())
    sim_kw = kw & Simulation.allowed_inputs
    inp_kw = (kw - sim_kw)
    sim_args = dict()
    inp_args  = dict()
    for kw in sim_kw:
        sim_args[kw] = kwargs[kw]
    #end for
    for kw in inp_kw:
        inp_args[kw] = kwargs[kw]
    #end for    
    if len(inp_args)>0:
        sim_args['input'] = generate_opium_input(input_type,**inp_args)
    #end if
    opium = Opium(**sim_args)

    return opium
#end def generate_opium
