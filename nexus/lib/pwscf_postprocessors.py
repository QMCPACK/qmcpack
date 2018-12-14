##################################################################
##  (c) Copyright 2016-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  pwscf_postprocessors.py                                           #
#    Nexus interfaces for PWSCF postprocessing tools:                #
#    pp.x, dos.x, bands.x, projwfc.x, cppp.x, pw_export.x            #
#                                                                    #
#  Content summary:                                                  #
#    generate_pp                                                     #
#      User-facing function to create pp simulation objects.         #
#                                                                    #
#    generate_pp_input                                               #
#      User-facing function to create input for pp.                  #
#                                                                    #
#    PP                                                              #
#      Simulation class for pp.                                      #
#                                                                    #
#    PPInput                                                         #
#      SimulationInput class for pp.                                 #
#                                                                    #
#                                                                    #
#    generate_dos                                                    #
#      User-facing function to create dos simulation objects.        #
#                                                                    #
#    generate_dos_input                                              #
#      User-facing function to create input for dos.                 #
#                                                                    #
#    Dos                                                             #
#      Simulation class for dos.                                     #
#                                                                    #
#    DosInput                                                        #
#      SimulationInput class for dos.                                #
#                                                                    #
#                                                                    #
#    generate_bands                                                  #
#      User-facing function to create bands simulation objects.      #
#                                                                    #
#    generate_bands_input                                            #
#      User-facing function to create input for bands.               #
#                                                                    #
#    Bands                                                           #
#      Simulation class for bands.                                   #
#                                                                    #
#    BandsInput                                                      #
#      SimulationInput class for bands.                              #
#                                                                    #
#                                                                    #
#    generate_projwfc                                                #
#      User-facing function to create projwfc simulation objects.    #
#                                                                    #
#    generate_projwfc_input                                          #
#      User-facing function to create input for projwfc.             #
#                                                                    #
#    Projwfc                                                         #
#      Simulation class for projwfc.                                 #
#                                                                    #
#    ProjwfcInput                                                    #
#      SimulationInput class for projwfc.                            #
#                                                                    #
#    ProjwfcAnalyzer                                                 #
#      SimulationAnalyzer class for projwfc.                         #
#                                                                    #
#                                                                    #
#    generate_cppp                                                   #
#      User-facing function to create cppp simulation objects.       #
#                                                                    #
#    generate_cppp_input                                             #
#      User-facing function to create input for cppp.                #
#                                                                    #
#    Cppp                                                            #
#      Simulation class for cppp.                                    #
#                                                                    #
#    CpppInput                                                       #
#      SimulationInput class for cppp.                               #
#                                                                    #
#                                                                    #
#    generate_pwexport                                               #
#      User-facing function to create pwexport simulation objects.   #
#                                                                    #
#    generate_pwexport_input                                         #
#      User-facing function to create input for pwexport.            #
#                                                                    #
#    Pwexport                                                        #
#      Simulation class for pwexport.                                #
#                                                                    #
#    PwexportInput                                                   #
#      SimulationInput class for pwexport.                           #
#                                                                    #
#                                                                    #
#    NamelistInput                                                   #
#      Generic SimulationInput class to handle Fortran namelists.    #
#      Base class for all *Input classes listed above.               #
#                                                                    #
#    PostProcessSimulation                                           #
#      Generic Simulation class for all PWSCF postprocessing tools.  #
#      Base class for all *Simulation classes listed above.          #
#                                                                    #
#    generate_ppsim                                                  #
#      Generic simulation object generator function.                 #
#                                                                    #
#====================================================================#



import os   
from generic import obj
from fileio import TextFile
from simulation import Simulation,SimulationInput,SimulationAnalyzer,NullSimulationAnalyzer
from developer import DevBase,ci


booldict = {'.true.':True,'.false.':False}
def readval(val):
    if val in booldict:
        v = booldict[val]   
    else:
        try:
            v = int(val)
        except:
            try:
                v = float(val.replace('d','e').replace('D','e'))
            except:
                if ' ' in val or ',' in val:
                    v = None # fail, no array handling yet
                else:
                    v = val.strip('"').strip("'")
                #end if
            #end try
        #end try
    #end if
    return v
#end def readval


def writeval(val):
    if isinstance(val,bool):
        if val:
            sval = '.true.'
        else:
            sval = '.false.'
        #end if
    elif isinstance(val,str):
        sval =  "'"+val+"'"
    elif isinstance(val,int):
        sval = str(val)
    elif isinstance(val,float):
        sval = str(val).replace('e','d')
    else:
        sval = None # fail, no array handling yet
    #end if
    return sval
#end def writeval



class Namelist(DevBase):
    @classmethod
    def class_init(cls):
        cls.class_set_optional(
            namelist = 'unknown',
            names = [],
            )
        cls.name_set = set(cls.names)
    #end def class_init

        
    def __init__(self,text=None,**vals):
        if text!=None:
            self.read_text(text)
        #end if
        if len(vals)>0:
            self.assign_values('initialization',**vals)
        #end if
    #end def __init__


    def check_names(self,label,names):
        cls = self.__class__
        if len(cls.name_set)>0:
            invalid = set(names)-cls.name_set
            if len(invalid)>0:
                self.error('invalid names encountered in namelist during {0}\nnamelist name: {1}\ninvalid names: {2}\nvalid options are: {3}'.format(label,self.namelist,sorted(invalid),cls.names))
            #end if
        #end if
    #end def check_names


    def assign_values(self,label,**vals):
        self.check_names(label,vals.keys())
        for name,value in vals.iteritems():
            self[name] = value
        #end for
    #end def assign_values


    def read_text(self,text):
        cls = self.__class__
        if isinstance(text,str):
            lines = text.split()
        elif isinstance(text,list):
            lines = text
        else:
            self.error('read_text only accepts string or list inputs for text\nencountered invalid type for text: {0}'.format(text.__class__.__name__))
        #end if
        if len(lines)>0:
            if lines[0].strip().startswith('&'):
                lines = lines[1:]
            #end if
            if lines[-1].strip().endswith('/'):
                lines = lines[:-1]
            #end if
        #end if
        vals = obj()
        for line in lines:
            tokens = line.split(',')
            for t in tokens:
                name,value = t.split('=')
                name  = name.strip()
                value = value.strip()
                v = readval(value)
                if v!=None:
                    vals[name] = v
                else:
                    self.error('namelist read failed\nnamelist name: {0}\nvariable name: {1}\nvariable value: {2}'.format(self.namelist,name,value))
                #end if
            #end for
        #end for
        self.assign_values('read',**vals)
    #end def read_text


    def write_text(self):
        cls = self.__class__
        has_namelist = 'namelist' in self
        if has_namelist:
            namelist = self.namelist
            del self.namelist
        else:
            namelist = cls.namelist
        #end if
        self.check_names('write',self.keys())
        text = '&'+namelist+'\n'
        for name,value in self.iteritems():
            v = writeval(value)
            if v!=None:
                text += '  {0} = {1}\n'.format(name,v)
            else:
                self.error('namelist write failed\nnamelist name: {0}\nvariable name: {1}\nvariable value: {2}'.format(namelist,name,value))
            #end if
        #end for
        text += '/\n'
        if has_namelist:
            self.namelist = namelist
        #end if
        return text
    #end def write_text
#end class Namelist



class NamelistInput(SimulationInput):
    @classmethod
    def class_init(cls):
        cls.class_set_optional(
            namelists = [],
            namelist_classes = obj(),
            )
        cls.namelist_set = set(cls.namelists)
        cls.name_map     = obj()
        for namelist_name,namelist_cls in cls.namelist_classes.iteritems():
            for name in namelist_cls.names:
                cls.name_map[name] = namelist_name
            #end for
        #end for
    #end def class_init


    def __init__(self,filepath=None,**vals):
        if filepath!=None:
            self.read(filepath)
            return
        #end if
        cls = self.__class__
        if len(cls.namelists)==0:
            self.error('cannot initialize this input class as no namelists have been assigned it')
        #end if
        for name,value in vals.iteritems():
            if name in cls.name_map:
                namelist_name = cls.name_map[name]
                if namelist_name not in self:
                    namelist = cls.namelist_classes[namelist_name]()
                    self[namelist_name] = namelist
                else:
                    namelist = self[namelist_name]
                #end if
                namelist[name] = value
            else:
                self.error('encountered invalid variable name during initialization\ninvalid variable name: {0}\nthis variable does not belong to any of the following namelists: {1}'.format(name,cls.namelists))
            #end if
        #end for
    #end def __init__


    def read_text(self,text,filepath=None):
        cls      = self.__class__
        lines    = text.split('\n')
        inside   = False
        name     = None
        nl_lines = []
        for l in lines:
            ls = l.strip()
            if inside:
                nl_lines.append(l)
            #end if
            if ls.startswith('&'):
                inside=True
                name = ls[1:].lower()
            elif ls.endswith('/'):
                inside=False
                if len(cls.namelist_classes)==0:
                    self[name] = Namelist(nl_lines)
                elif name in cls.namelist_classes:
                    self[name] = cls.namelist_classes[name](nl_lines)
                else:
                    msg = 'encountered invalid namelist during read\ninvalid namelist: {0}\nvalid namelists are: {1}'.format(name,cls.namelists)
                    if filepath!=None:
                        msg += '\nfilepath: {0}'.format(filepath)
                    #end if
                    self.error(msg)
                #end if
            #end if
        #end for
    #end def read_text


    def write_text(self,filepath=None):
        cls = self.__class__
        text = ''
        for name,namelist in self.iteritems():
            if name not in cls.namelist_set:
                msg = 'encountered invalid namelist during write\ninvalid namelist: {0}\nvalid namelists are: {1}'.format(name,cls.namelists)
                if filepath!=None:
                    msg += '\nfilepath: {0}'.format(filepath)
                #end if
                self.error(msg)
            #end if
        #end for
        for name in cls.namelists:
            if name in self:
                namelist = self[name]
                text += namelist.write_text()+'\n'
            #end if
        #end for
        return text
    #end def write_text
#end class NamelistInput



class PostProcessSimulation(Simulation):
    analyzer_type = NullSimulationAnalyzer

    def check_result(self,result_name,sim):
        return False
    #end def check_result    

    def app_command(self):
        return self.app_name+'<'+self.infile
    #end def app_command

    def check_sim_status(self):
        self.finished = True
    #end def check_sim_status

    def get_output_files(self):
        return []
    #end def get_output_files
#end class PostProcessSimulation



def generate_ppsim(gen_input=None,Sim=None,**kwargs):
    sim_args,inp_args = Simulation.separate_inputs(kwargs)
    if not 'input' in sim_args:
        sim_args.input = gen_input(**inp_args)
    #end if
    sim = Sim(**sim_args)
    return sim
#end def generate_ppsim




class PPInputppNamelist(Namelist):
    namelist = 'inputpp'
    names    = ['prefix','outdir','filplot','plot_num','spin_component',
                'sample_bias','kpoint','kband','lsign','emin','emax']
#end class PPInputppNamelist

class PPPlotNamelist(Namelist):
    namelist = 'plot'
    names    = ['nfile','filepp','weight','iflag','output_format',
                'fileout','interpolation','e1','e2','e3','x0',
                'nx','ny','nz','radius']
#end class PPPlotNamelist


class PPInput(NamelistInput):
    namelists = ['inputpp','plot']
    namelist_classes = obj(
        inputpp = PPInputppNamelist,
        plot    = PPPlotNamelist,
        )
#end class PPInput


PPAnalyzer = NullSimulationAnalyzer


class PP(PostProcessSimulation):
    input_type         = PPInput
    generic_identifier = 'pp'
    application        = 'pp.x'
#end class PP


generate_pp_input = PPInput

def generate_pp(**kwargs):
    return generate_ppsim(PPInput,PP,**kwargs)
#end def generate_pp




class DosNamelist(Namelist):
    namelist = 'dos'
    names = ['prefix','outdir','ngauss','degauss',
             'Emin','Emax','DeltaE','fildos']
#end class DosNamelist


class DosInput(NamelistInput):
    namelists = ['dos']
    namelist_classes = obj(
        dos = DosNamelist,
        )
#end class DosInput


DosAnalyzer = NullSimulationAnalyzer


class Dos(PostProcessSimulation):
    input_type         = DosInput
    generic_identifier = 'dos'
    application        = 'dos.x'
#end class Dos


generate_dos_input = DosInput

def generate_dos(**kwargs):
    return generate_ppsim(DosInput,Dos,**kwargs)
#end def generate_dos




class BandsNamelist(Namelist):
    namelist = 'bands'
    names = ['prefix','outdir','filband','spin_component',
             'lsigma','lp','filp','lsym','no_overlap','plot_2d',
             'firstk','lastk']
#end class BandsNamelist


class BandsInput(NamelistInput):
    namelists = ['bands']
    namelist_classes = obj(
        bands = BandsNamelist,
        )
#end class BandsInput


BandsAnalyzer = NullSimulationAnalyzer


class Bands(PostProcessSimulation):
    input_type         = BandsInput
    generic_identifier = 'bands'
    application        = 'bands.x'
#end class Bands


generate_bands_input = BandsInput

def generate_bands(**kwargs):
    return generate_ppsim(BandsInput,Bands,**kwargs)
#end def generate_bands




class ProjwfcNamelist(Namelist):
    namelist = 'projwfc'
    names = ['ngauss','degauss','Emin','Emax','deltaE',
             'prefix','outdir','fildos','filproj',
             'lsym','pawproj','lwrite_overlaps','lbinary_data']
#end class ProjwfcNamelist


class ProjwfcInput(NamelistInput):
    namelists = ['projwfc']
    namelist_classes = obj(
        projwfc = ProjwfcNamelist,
        )
#end class ProjwfcInput


class ProjwfcAnalyzer(SimulationAnalyzer):
    def __init__(self,arg0=None,outfile=None,analyze=False,warn=False,strict=False):
        self.info = obj(
            outfile     = outfile,
            warn        = warn,
            strict      = strict,
            initialized = False,
            )

        if isinstance(arg0,Simulation):
            sim = arg0
            path    = sim.locdir
            infile  = sim.infile
            outfile = sim.outfile
            infile_path = os.path.join(path,infile)
        elif arg0!=None:
            infile_path = arg0
            path,infile = os.path.split(infile_path)
            if outfile is None:
                outfile = infile.rsplit('.',1)[0]+'.out'
            #end if
        else:
            return
        #end if

        self.input = ProjwfcInput(infile_path)

        self.info.set(
            path        = path,
            infile      = infile,
            outfile     = outfile,
            initialized = True,
            )

        if analyze:
            self.analyze()
        #end if
    #end def __init__ 


    def analyze(self):
        operations = [
            ('open_log'   ,ProjwfcAnalyzer.open_log),
            ('read_states',ProjwfcAnalyzer.read_states),
            ('read_lowdin',ProjwfcAnalyzer.read_lowdin),
            ('close_log'  ,ProjwfcAnalyzer.close_log),
            ]
        if self.info.strict:
            for name,op in operations:
                op(self)
            #end for
        else:
            failures = []
            for name,op in operations:
                try:
                    op(self)
                except:
                    failures.append(name)
                #end try
            #end for
            if len(failures)>0 and self.info.warn:
                self.warn('analysis failed, some data will not be available\noperations failed: {0}'.format(failures))
            #end if
        #end if
    #end def analyze


    def open_log(self):
        logfile = os.path.join(self.info.path,self.info.outfile)
        self.log = TextFile(logfile)
    #end def open_log

    def read_states(self):
        log = self.log
        log.seek('state #')
        nstates  = 0
        elem_ind = set()
        elem     = []
        while True:
            line = log.readline()
            tokens = line.replace('(',' ').replace(')',' ').split()
            if not (len(tokens)>0 and tokens[0]=='state'):
                break
            #end if
            ei,e = tokens[4],tokens[5] 
            if ei not in elem_ind:
                elem.append(e)
                elem_ind.add(ei)
            #end if
            nstates += 1
        #end while
        self.states = obj(nstates=nstates,elem=elem)
    #end def read_states

    def read_lowdin(self):
        log = self.log
        log.seek('Lowdin Charges')
        lowdin = obj()
        has_ud = False
        nmax = len(self.states.elem)*20
        n  = 0
        ls = ''
        cur_atom = -1
        while n<nmax and not ls.startswith('Spilling'):
            n+=1
            ls = log.readline().strip()
            if ls.startswith('Atom'):
                astr,ls = ls.split(':')
                cur_atom = int(astr.split('#')[1])-1
                if cur_atom not in lowdin:
                    lowdin[cur_atom] = obj(tot=obj(),up=obj(),down=obj())
                #end if
                lc = lowdin[cur_atom]                
            #end if
            if 'tot' in ls:
                lc_comp = lc.tot
            elif 'up' in ls:
                lc_comp = lc.up
            elif 'down' in ls:
                lc_comp = lc.down
            else:
                continue
            #end if
            tokens = ls.replace(' ','').rstrip(',').split(',')
            for t in tokens:
                name,value = t.split('=')
                if 'spin' in name:
                    has_ud = True
                    name = 'charge'
                elif 'charge' in name:
                    name = 'charge'
                #end if
                lc_comp[name] = float(value)
            #end for
        #end while
        if has_ud:
            for lc in lowdin:
                u = lc.up
                d = lc.down
                lc.pol = obj()
                for k,uv in u.iteritems():
                    dv = d[k]
                    lc.tot[k] = uv + dv
                    lc.pol[k] = uv - dv
                #end for
            #end for
        else:
            for lc in lowdin:
                del lc.up
                del lc.down
            #end for
        #end if
        self.lowdin = lowdin
    #end def read_lowdin

    def write_lowdin(self,filepath=None,sum=None,tot=None,pol=None,up=None,down=None,all=True,long=False):
        if tot is None:
            tot = all
        #end if
        if pol is None:
            pol = all
        #end if
        if up is None:
            up = all
        #end if
        if down is None:
            down = all
        #end if
        if sum is None:
            sum = all
        #end if
        sections = [('tot',tot),('pol',pol),('up',up),('down',down)]
        elem=None
        if 'states' in self:
            elem = self.states.elem
        #end if
        lowdin = self.lowdin
        text   = ''
        if sum:
            nelec = '?'
            npol  = '?'
            if len(lowdin)>0:
                if 'tot' in lowdin[0]:
                    nelec = 0
                    for lc in lowdin:
                        nelec += lc.tot.charge
                    #end for
                #end if
                if 'pol' in lowdin[0]:
                    npol = 0
                    for lc in lowdin:
                        npol += lc.pol.charge
                    #end for
                #end if
            #end if
            text += 'nup+ndn = {0}\n'.format(nelec)
            text += 'nup-ndn = {0}\n'.format(npol)
            text += '\n'
        #end if
        lvals  = 'spdfg'
        for q,qwrite in sections:
            if qwrite and len(lowdin)>0 and q in lowdin[0]:
                text+=q+'\n'
                for n in range(len(lowdin)):
                    lc = lowdin[n][q]
                    if elem is None:
                        text += '  {0:>3}  {1: 3.2f}  '.format(n,lc.charge)
                    else:
                        text += '  {0:>3}  {1:>2}  {2: 3.2f}  '.format(n,elem[n],lc.charge)
                    #end if
                    if not long:
                        for l in lvals:
                            if l in lc:
                                text += '{0}({1: 3.2f})'.format(l,lc[l])
                            #end if
                        #end for
                    else:
                        for l in lvals:
                            if l in lc:
                                lset = []
                                for k in lc.keys():
                                    if k.startswith(l) and (len(k)>1 or k=='s'):
                                        lset.append(k)
                                    #end if
                                #end for
                                for k in sorted(lset):
                                    text += '{0}({1: 3.2f})'.format(k,lc[k])
                                #end for
                            #end if
                        #end for
                    #end if
                    text+='\n'
                #end for
                text+='\n'
            #end if
        #end for
        if filepath!=None:
            open(filepath,'w').write(text)
        #end if
        return text
    #end def write_lowdin

    def close_log(self):
        if 'log' in self:
            del self.log
        #end if
    #end def close_log
        
#end class ProjwfcAnalyzer


class Projwfc(PostProcessSimulation):
    input_type         = ProjwfcInput
    analyzer_type      = ProjwfcAnalyzer
    generic_identifier = 'projwfc'
    application        = 'projwfc.x'

    def post_analyze(self,analyzer):
        # try to write lowdin output data file
        try:
            lowdin_file = self.identifier+'.lowdin'
            filepath = os.path.join(self.locdir,lowdin_file)
            analyzer.write_lowdin(filepath)
            analyzer.write_lowdin(filepath+'_long',long=True)
        except:
            None
        #end try
    #end def post_analyze
#end class Projwfc


def generate_projwfc_input(prefix='pwscf',outdir='pwscf_output',**vals):
    pp = ProjwfcInput(
        prefix = prefix,
        outdir = outdir,
        **vals
        )
    return pp
#end def generate_projwfc_input


def generate_projwfc(**kwargs):
    return generate_ppsim(generate_projwfc_input,Projwfc,**kwargs)
#end def generate_projwfc




class CpppInputppNamelist(Namelist):
    namelist = 'inputpp'
    names = ['prefix','fileout','output','outdir','lcharge',
             'lforces','ldynamics','lpdb','lrotation',
             'ns1','ns2','ns3','np1','np2','np3','nframes','ndr',
             'atomic_number','charge_density','state','lbinary']
#end class CpppInputppNamelist


class CpppInput(NamelistInput):
    namelists = ['inputpp']
    namelist_classes = obj(
        inputpp = CpppInputppNamelist,
        )
#end class CpppInput


CpppAnalyzer = NullSimulationAnalyzer


class Cppp(PostProcessSimulation):
    input_type         = CpppInput
    generic_identifier = 'cppp'
    application        = 'cppp.x'
#end class Cppp


generate_cppp_input = CpppInput

def generate_cppp(**kwargs):
    return generate_ppsim(CpppInput,Cppp,**kwargs)
#end def generate_cppp




class PwexportInputppNamelist(Namelist):
    namelist = 'inputpp'
    names = ['prefix','outdir','pseudo_dir','psfile',
             'single_file','ascii','pp_file','uspp_spsi']
#end class PwexportInputppNamelist


class PwexportInput(NamelistInput):
    namelists = ['inputpp']
    namelist_classes = obj(
        inputpp = PwexportInputppNamelist,
        )
#end class PwexportInput


PwexportAnalyzer = NullSimulationAnalyzer


class Pwexport(PostProcessSimulation):
    input_type         = PwexportInput
    generic_identifier = 'pwexport'
    application        = 'pw_export.x'
#end class Pwexport


generate_pwexport_input = PwexportInput

def generate_pwexport(**kwargs):
    return generate_ppsim(PwexportInput,Pwexport,**kwargs)
#end def generate_pwexport




namelist_classes = [
    Namelist               , NamelistInput,
    ProjwfcNamelist        , ProjwfcInput,
    PPInputppNamelist      , PPPlotNamelist, PPInput,
    DosNamelist            , DosInput,
    BandsNamelist          , BandsInput,
    CpppInputppNamelist    , CpppInput,
    PwexportInputppNamelist, PwexportInput,
    ]
for cls in namelist_classes:
    cls.class_init()
#end for




if __name__=='__main__':
    pi = generate_projwfc_input(
        prefix = 'pwscf',
        outdir = 'pwscf_output',
        )
    print pi
    print pi.write()

#end if
