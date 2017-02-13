##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  vasp_input.py                                                     #
#    Supports I/O, generation, and manipulation of VASP input files. #
#                                                                    #
#  Content summary:                                                  #
#    VaspInput                                                       #
#      SimulationInput class for VASP.                               #
#                                                                    #
#    generate_vasp_input                                             #
#      User-facing function to generate arbitrary VASP input files.  #
#                                                                    #
#    generate_poscar                                                 #
#      Function to create a Poscar object from a Structure object.   #
#                                                                    #
#    Vobj                                                            #
#      Base class for VASP input classes.                            #
#                                                                    #
#    VFile                                                           #
#      Base class for a VASP file.                                   #
#                                                                    #
#    VKeywordFile                                                    #
#      Base class for VASP input files in keyword format.            #
#      I/O handled at base class level.                              #
#      Derived classes contain the keyword spec. for each file.      #
#      See Incar and Stopcar classes.                                #
#                                                                    #
#    VFormattedFile                                                  #
#      Base class for VASP input files with strict formatting.       #
#      Derived classes handle specialized I/O for each file.         #
#      See Iconst, Kpoints, Penaltypot, Poscar, Potcar, and Exhcar.  #
#                                                                    #
#====================================================================#


import os
from numpy import array,abs,empty,ndarray,dot
from numpy.linalg import inv
from generic import obj
from periodic_table import is_element
from nexus_base import nexus_noncore
from simulation import SimulationInput
from structure import interpolate_structures,Structure
from physical_system import PhysicalSystem
from developer import DevBase,error
from debug import *
lcs = ls



# support functions for keyword files

def remove_comment(sval):
    if ' ' in sval:
        sval = sval.split(' ',1)[0]
    #end if
    return sval
#end def remove_comment

def expand_array(sval):
    sarr = []
    for v in sval.split():
        if '*' in v:
            n,vv = v.rsplit('*',1)
            sarr.extend(int(n)*[vv])
        else:
            sarr.append(v)
        #end if
    #end for
    return sarr
#end def expand_array

def read_int(sval):
    sval = remove_comment(sval)
    return int(sval)
#end def read_int

def read_real(sval):
    sval = remove_comment(sval)
    return float(sval.replace('d','e'))
#end def read_real

bool_dict = dict(true=True,false=False)
def read_bool(sval):
    sval = remove_comment(sval)
    return bool_dict[sval.lower().strip('.')]
#end def read_bool

def read_string(sval):
    return sval
#end def read_string

def read_int_array(sval):
    return array(expand_array(sval),dtype=int)
#end def read_int_array

def read_real_array(sval):
    return array(expand_array(sval),dtype=float)
#end def read_real_array

bool_array_dict = dict(T=True,F=False)
def read_bool_array(sval):
    barr = []
    for v in expand_array(sval):
        barr.append(bool_array_dict[v])
    #end for
    return array(barr,dtype=bool)
#end def read_bool_array

def write_int(v):
    return str(v)
#end def write_int

def write_real(v):
    return str(v)
#end def write_real

def write_bool(v):
    if v:
        return '.TRUE.'
    else:
        return '.FALSE.'
    #end if
#end def write_bool

def write_string(v):
    return v
#end def write_string

def equality(a,b):
    return a==b
#end def equality

def real_equality(a,b):
    return abs(a-b)<=1e-6*(abs(a)+abs(b))/2
#end def real_equality

def render_bool(v):
    if v:
        return 'T'
    else:
        return 'F'
    #end if
#end def render_bool

def write_array(arr,same=equality,render=str,max_repeat=3):
    value_counts = []
    count = 0
    value = arr[0]
    for v in arr:
        if same(v,value):
            count += 1
        else:
            value_counts.append((value,count))
            value = v
            count = 1
        #end if
    #end for
    if same(v,value):
        value_counts.append((value,count))
    else:
        value_counts.append((v,1))
    #end if
    s = ''
    for value,count in value_counts:
        if count>max_repeat:
            s += '{0}*{1} '.format(count,render(value))
        else:
            for i in xrange(count):
                s += render(value)+' '
            #end for
        #end if
    #end for
    return s
#end def write_array

def write_int_array(a):
    return write_array(a)
#end def write_int_array

def write_real_array(a):
    return write_array(a,same=real_equality)
#end def write_real_array

def write_bool_array(a):
    return write_array(a,render=render_bool)
#end def write_bool_array

assign_bool_map = {True:True,False:False,1:True,0:False}
def assign_bool(v):
    return assign_bool_map[v]
#end def assign_bool

def assign_string(v):
    if isinstance(v,str):
        return v
    else:
        raise ValueError('value must be a string')
    #end if
#end def assign_string

def assign_int_array(a):
    if isinstance(a,(tuple,list,ndarray)):
        return array(a,dtype=int)
    else:
        raise ValueError('value must be a tuple, list, or array')
    #end if
#end def assign_int_array

def assign_real_array(a):
    if isinstance(a,(tuple,list,ndarray)):
        return array(a,dtype=float)
    else:
        raise ValueError('value must be a tuple, list, or array')
    #end if
#end def assign_real_array

def assign_bool_array(a):
    if isinstance(a,(tuple,list,ndarray)):
        return array(a,dtype=bool)
    else:
        raise ValueError('value must be a tuple, list, or array')
    #end if
#end def assign_bool_array

#utility functions to convert from VASP internal objects to 
#nexus objects:
def vasp_to_nexus_elem(elem,elem_count):
  syselem=[]
  for x,count in zip(elem,elem_count):
    syselem+=[x for i in xrange(0,count)]
  return array(syselem)


read_value_functions = obj(
    ints        = read_int,
    reals       = read_real,
    bools       = read_bool,
    strings     = read_string,
    int_arrays  = read_int_array,
    real_arrays = read_real_array,
    bool_arrays = read_bool_array
    )

write_value_functions = obj(
    ints        = write_int,
    reals       = write_real,
    bools       = write_bool,
    strings     = write_string,
    int_arrays  = write_int_array,
    real_arrays = write_real_array,
    bool_arrays = write_bool_array
    )

assign_value_functions = obj(
    ints        = int,
    reals       = float,
    bools       = assign_bool,
    strings     = assign_string,
    int_arrays  = assign_int_array,
    real_arrays = assign_real_array,
    bool_arrays = assign_bool_array    
    )





class Vobj(DevBase):
    def get_path(self,filepath):
        if os.path.exists(filepath) and os.path.isdir(filepath):
            path = filepath
        else:
            path,tmp = os.path.split(filepath)
            if len(path)>0 and not os.path.exists(path):
                self.error('path {0} does not exist'.format(path))
            #end if
        #end if
        return path
    #end def get_path
#end class Vobj



class VFile(Vobj):
    def __init__(self,filepath=None):
        if filepath!=None:
            self.read(filepath)
        #end if
    #end def __init__

    def read(self,filepath):
        if not os.path.exists(filepath):
            self.error('file {0} does not exist'.format(filepath))
        #end if
        text = open(filepath,'r').read()
        self.read_text(text,filepath)
        return text
    #end def read

    def write(self,filepath=None):
        text = self.write_text(filepath)
        if filepath!=None:
            open(filepath,'w').write(text)
        #end if
        return text
    #end def write

    def read_text(self,text,filepath=''):
        self.not_implemented()
    #end def read_text

    def write_text(self,filepath=''):
        self.not_implemented()
    #end def write_text

    def remove_comment(self,line):
        cloc1 = line.find('!')
        cloc2 = line.find('#')
        has1  = cloc1!=-1
        has2  = cloc2!=-1
        if has1 or has2:
            if has1 and has2:
                cloc = min(cloc1,cloc2)
            elif has1:
                cloc = cloc1
            else:
                cloc = cloc2
            #end if
            line = line[:cloc].strip()
        #end if
        return line
    #end def remove_comment
#end class VFile



class VKeywordFile(VFile):
    kw_scalars = ['ints','reals','bools','strings']
    kw_arrays  = ['int_arrays','real_arrays','bool_arrays']
    kw_fields  = kw_scalars + kw_arrays + ['keywords','unsupported']

    keyword_classification = None

    @classmethod
    def class_init(cls):
        for kw_field in cls.kw_fields:
            if not cls.class_has(kw_field):
                cls.class_set_single(kw_field,set())
            #end if
        #end for
        #cls.check_consistency()
        cls.scalar_keywords = set()
        for scalar_field in cls.kw_scalars:
            cls.scalar_keywords |= cls.class_get(scalar_field)
        #end for
        cls.array_keywords = set()
        for array_field in cls.kw_arrays:
            cls.array_keywords |= cls.class_get(array_field)
        #end for
        cls.keywords = cls.scalar_keywords | cls.array_keywords
        cls.type = obj()
        cls.read_value   = obj()
        cls.write_value  = obj()
        cls.assign_value = obj()
        for type in cls.kw_scalars + cls.kw_arrays:
            for name in cls.class_get(type):
                cls.type[name] = type
                cls.read_value[name]   = read_value_functions[type]
                cls.write_value[name]  = write_value_functions[type]
                cls.assign_value[name] = assign_value_functions[type]
            #end for
        #end for
    #end def class_init

    @classmethod
    def check_consistency(cls):
        fail = False
        msg  = ''
        types = cls.kw_scalars+cls.kw_arrays
        untyped = cls.keywords 
        for type in types:
            untyped -= cls.class_get(type)
        #end for
        if len(untyped)>0:
            fail = True
            msg += 'variables without a type: {0}\n'.format(sorted(untyped))
        #end if
        for type in types:
            unknown = cls.class_get(type)-cls.keywords
            if len(unknown)>0:
                fail = True
                msg += 'unknown {0}: {1}\n'.format(type,sorted(unknown))
            #end if
        #end for
        if fail:
            cls.class_error(msg)
        #end if
    #end def check_consistency


    def read_text(self,text,filepath=''):
        lines = text.splitlines()
        expression = None
        continued  = False
        for line in lines:
            ls = line.strip()
            if len(ls)>0 and ls[0]!='!' and ls[0]!='#':
                ls = self.remove_comment(ls)
                this_cont = ls.endswith('\\')
                if this_cont:
                    ls = ls.rstrip('\\')
                    if continued:
                        expression += ls
                    else:
                        expression = ls
                        continued = True
                    #end if
                elif continued:
                    expression += ls
                    continued = False
                else:
                    expression = ls
                #end if
                if not continued:
                    tokens = expression.split(';')
                    for token in tokens:
                        if '=' in token:
                            name,value = token.split('=',1)
                            name  = name.lower().strip()
                            value = value.strip()
                            if name in self.keywords:
                                try:
                                    value = self.read_value[name](value)
                                    self[name] = value
                                except Exception,e:
                                    self.error('read failed for keyword {0}\nkeyword type: {1}\ninput text: {2}\nexception:\n{3}'.format(name,self.type[name],token,e))
                                #end try
                            elif name in self.unsupported:
                                self.warn('keyword {0} is not currently supported'.format(name))
                            else:
                                #ci(lcs(),gs())
                                self.error('{0} is not a keyword for the {1} file'.format(name.upper(),self.__class__.__name__.upper()))
                            #end if
                        #end if
                    #end for
                #end if
            #end if
        #end for
    #end def read_text

                                
    def write_text(self,filepath=''):
        text = ''
        maxlen=0
        for name in self.keys():
            maxlen = max(maxlen,len(name))
        #end for
        #maxlen = min(maxlen,9)
        valfmt = '{0:<'+str(maxlen)+'} = {1}\n'
        for name in sorted(self.keys()):
            value = self[name]
            try:
                svalue = self.write_value[name](value)
            except Exception,e:
                self.error('write failed for file {0} keyword {1}\nkeyword type: {2}\nvalue: {3}\nexception:\n{4}'.format(filepath,name,self.type[name],value,e))
            #end try
            text += valfmt.format(name.upper(),svalue)
        #end for
        return text
    #end def write_text


    def assign(self,**values):
        for name,value in values.iteritems():
            try:
                self[name] = self.assign_value[name](value)
            except Exception,e:
                self.error('assign failed for keyword {0}\nkeyword type: {1}\nvalue: {2}\nexception:\n{3}'.format(name,self.type[name],value,e))
            #end try
        #end for
    #end def assign
#end class VKeywordFile



class VFormattedFile(VFile):
    def read_lines(self,text,remove_empty=False):
        raw_lines = text.splitlines()
        lines = []
        for line in raw_lines:
            ls = self.remove_comment(line).strip()
            if not remove_empty or len(ls)>0:
                lines.append(ls)
            #end if
        #end for
        return lines
    #end def read_lines

    def join(self,lines,first_line,last_line):
        joined = ''
        for iline in xrange(first_line,last_line):
            joined += lines[iline]+' '
        #end for
        joined += lines[last_line]
        return joined
    #end def join

    def is_empty(self,lines,start=None,end=None):
        if start is None:
            start = 0
        #end if
        if end is None:
            end = len(lines)
        #end if
        is_empty = True
        for line in lines[start:end]: 
            is_empty &= len(line)==0
        #end for
        return is_empty
    #end def is_empty
#end class VFormattedFile


 
class Incar(VKeywordFile):

    # VTST extensions:  http://theory.cm.utexas.edu/vtsttools/index.html
    #   ichain lclimb ltangentold ldneb lnebcell jacobian timestep

    keywords = set('''
      addgrid aexx aggac aggax aldac algo amin amix amix_mag andersen_prob apaco 
      bmix bmix_mag 
      clnt cln cll clz cmbj cshift
      deper dimer_dist dipol dq
      ebreak eint ediff ediffg efield efield_pead elmin emax emin enaug encut 
      encutfock encutgw encutgwsoft enmax enmin epsilon evenonly evenonlygw
      ferdo ferwe findiff
      gga gga_compat
      hfscreen hflmaxf hills_bin hills_h hills_w 
      ialgo iband ibrion ichain icharg ichibare i_constrained_m icorelevel idipol 
      igpar images imix increm inimix iniwav ipead isif ismear ispin istart 
      isym ivdw iwavpr 
      jacobian
      kblock kgamma kpar kpuse kspacing 
      lambda langevin_gamma langevin_gamma_l lasph lasync lattice_constraints 
      lberry lblueout lcalceps lcalcpol lcharg lchimag lclimb lcorr ldau ldauj ldaul 
      ldauprint ldautype ldauu ldiag ldipol ldneb lefg lelf lepsilon lhfcalc lhyperfine 
      lkproj lmaxfock lmaxfockae lmaxfockmp2 lmaxmix lmaxmp2 lmaxpaw lmaxtau 
      lmixtau lmono lnabla lnebcell lnmr_sym_red lnoncollinear loptics lorbit lpard 
      lpead lplane lreal lrpa lscalapack lscaler0 lscalu lscsgrad lselfenergy 
      lsepb lsepk lspectral lsorbit ltangentold lthomas luse_vdw lvdw lvdw_ewald lvdwscs 
      lvhar lvtot lwave 
      magmom maxmem maxmix mbja mbjb m_constr mdalgo metagga minrot mixpre 
      nbands nbandsgw nblk nblock nbmod ncore nedos nelect nelm nelmdl nelmin 
      nfree ngx ngxf ngy ngyf ngz ngzf nkred nkredx nkredy nkredz nlspline 
      nmaxfockae nomega nomegar npaco npar nppstr nsim nsw nsubsys nupdown 
      nwrite
      oddonly oddonlygw ofield_a ofield_kappa ofield_q6_far ofield_q6_near 
      omegamax omegamin omegatl
      param1 param2 pmass pomass potim prec precfock pstress psubsys
      random_seed ropt rwigs 
      saxis scsrad shakemaxiter shaketol sigma skip_edotp smass spring step_max 
      step_size symprec system 
      tebeg teend time timestep tsubsys
      value_max value_min vdw_a1 vdw_a2 vdw_alpha vdw_cnradius vdw_c6 vdw_c6au 
      vdw_d vdw_radius vdw_r0 vdw_r0au vdw_scaling vdw_sr vdw_s6 vdw_s8 voskown 
      wc weimin 
      zab_vdw zval  
      '''.split()) # only used to check consistency of typed names below

    # some of these are mixed type arrays or other oddly formatted fields
    unsupported = set('quad_efg'.split())

    ints = set('''
      apaco 
      clnt cln cll clz 
      elmin
      findiff
      hflmaxf hills_bin
      ialgo ibrion ichain icharg ichibare i_constrained_m icorelevel idipol igpar 
      images imix inimix iniwav ipead isif ismear ispin istart isym ivdw iwavpr 
      kblock kpar 
      ldauprint ldautype lmaxfock lmaxfockae lmaxfockmp2 lmaxmix lmaxmp2 
      lmaxpaw lorbit 
      maxmem maxmix mdalgo mixpre 
      nbands nbandsgw nblk nblock nbmod ncore nedos nelm nelmdl nelmin nfree 
      ngx ngxf ngy ngyf ngz ngzf nkred nkredx nkredy nkredz nmaxfockae nomega 
      nomegar npaco npar nppstr nsim nsw nupdown nwrite 
      shakemaxiter smass spring 
      voskown
      '''.split())

    reals = set('''
      aexx aggac aggax aldac amin amix amix_mag andersen_prob
      bmix bmix_mag 
      cshift
      deper dimer_dist dq
      ebreak ediff ediffg efield emax emin enaug encut encutfock encutgw 
      encutgwsoft enmax enmin epsilon
      hfscreen hills_h hills_w
      jacobian
      kspacing 
      lambda langevin_gamma_l
      mbja mbjb minrot
      nelect 
      ofield_a ofield_kappa ofield_q6_far ofield_q6_near omegamax omegamin omegatl
      param1 param2 pmass pomass potim pstress 
      scsrad shaketol sigma step_max step_size symprec 
      tebeg teend time timestep
      vdw_a1 vdw_a2 vdw_cnradius vdw_d vdw_radius vdw_scaling vdw_sr vdw_s6 vdw_s8
      wc weimin 
      zab_vdw zval 
      '''.split())

    bools = set('''
      addgrid
      evenonly evenonlygw
      gga_compat 
      lasph lasync lberry lblueout lcalceps lcalcpol lcharg lchimag lclimb lcorr 
      ldau ldiag ldipol ldneb lefg lelf lepsilon lhfcalc lhyperfine lkproj lmaxtau 
      lmixtau lmono lnabla lnebcell lnmr_sym_red lnoncollinear loptics lpard lpead 
      lplane lrpa lscalapack lscaler0 lscalu lscsgrad lselfenergy lsepb 
      lsepk lsorbit lspectral ltangentold lthomas luse_vdw lvdw lvdw_ewald lvdwscs 
      lvhar lvtot lwave 
      kgamma 
      nlspline
      oddonly oddonlygw
      skip_edotp
      '''.split())

    strings = set('''
      algo 
      gga 
      lreal
      metagga 
      prec precfock
      system
      '''.split())

    int_arrays = set('''
      iband 
      kpuse 
      ldaul
      nsubsys 
      random_seed
      '''.split())

    real_arrays = set('''
      cmbj
      dipol 
      efield_pead eint
      ferdo ferwe 
      increm 
      langevin_gamma ldauj ldauu
      magmom m_constr
      psubsys 
      ropt rwigs 
      saxis
      tsubsys 
      value_max value_min vdw_alpha vdw_c6 vdw_c6au vdw_r0 vdw_r0au
      '''.split())

    bool_arrays = set('''
      lattice_constraints
      '''.split()) # formatted: F F T, etc


    keyword_classification = obj(
        array_dimensions = '''
        nkpts nkdim nbands nedos nions ldim lmdim nplwv irmax irdmax 
        ngx ngy ngz ngxf ngyf ngzf 
        '''.split(),
        )
#end class Incar



class Stopcar(VKeywordFile):
    keywords = set('lstop labort'.split())
    bools    = set('lstop labort'.split())
#end class Stopcar


for cls in Incar,Stopcar:
    cls.class_init()
#end for
del VKeywordFile.kw_scalars
del VKeywordFile.kw_arrays
del VKeywordFile.kw_fields



class Iconst(VFormattedFile):  # metadynamics -> 6.62.4
    None
#end class Iconst



class Kpoints(VFormattedFile):

    #  mode == explicit
    #    coord    = cartesian/reciprocal
    #    kpoints  = list of 3D kpoints
    #    kweights = list of kpoint weights
    #    tetrahedra = optional list of tetra objects (volume, degeneracy, corners)
    #
    #  mode == line
    #    coord      = cartesian/reciprocal
    #    kinsert    = number of points inserted between each pair of endpoints
    #    kendpoints = kpoint pairs forming line endpoints
    #
    #  mode == auto
    #    centering = auto/gamma/monkhorst-pack
    #    kgrid     = number of grid points for each direction (single integer)
    #    kshift    = optional shift of k-point grid
    #    
    #  mode == basis
    #    coord  = cartesian/reciprocal
    #    kbasis = 3x3 matrix of kpoint basis vectors
    #    kshift = shift of kpoint mesh

    centering_options = obj(a='auto',g='gamma',m='monkhorst-pack')

    def coord_options(self,cselect):
        if cselect=='c' or cselect=='k':
            return 'cartesian'
        else:
            return 'reciprocal'
        #end if
    #end def coord_options


    def __init__(self,filepath=None):
        self.mode = None  # explicit, line, auto, basis
        VFile.__init__(self,filepath)
    #end def __init__


    def read_text(self,text,filepath=''):
        lines = self.read_lines(text,remove_empty=True)
        if len(lines)>2:
            if not ' ' in lines[1]:
                iselect = int(lines[1])
            else: # erroneous case? (e.g. user supplies '0 0 0' instead of '0')
                iselect = int(lines[1].split()[0])
            #end if
            cselect = lines[2].lower()[0]
            if iselect==0: # auto or basis
                if cselect=='a':  # fully auto mesh
                    self.mode      = 'auto'
                    self.centering = self.centering_options[cselect]
                    self.kgrid     = int(lines[3])
                elif cselect=='g' or cselect=='m': # gamma or monkhorst mesh
                    self.mode      = 'auto'
                    self.centering = self.centering_options[cselect]
                    self.kgrid     = array(lines[3].split(),dtype=int)
                    if len(lines)>4:
                        self.kshift    = array(lines[4].split(),dtype=float)
                    else:
                        self.kshift = None
                    #end if
                else:
                    self.mode   = 'basis'  # basis generated mesh
                    self.coord  = self.coord_options(cselect)
                    self.kbasis = array(self.join(lines,3,5).split(),dtype=float)
                    self.kbasis.shape = 3,3
                    self.kshift = array(lines[6].split(),dtype=float)
                #end if
            elif cselect=='l': # line mode (band structure)
                self.mode    = 'line'
                self.kinsert = iselect
                self.coord   = self.coord_options(lines[3].lower()[0])
                endpoints = []
                for line in lines[4:]:
                    endpoints.append(line.split())
                #end for
                self.kendpoints = array(endpoints,dtype=float)
            else: # explicit kpoints
                self.mode  = 'explicit'
                self.coord = self.coord_options(cselect)
                nkpoints = iselect
                kpw = []
                for line in lines[3:3+nkpoints]:
                    kpw.append(line.split())
                #end for
                kpw = array(kpw,dtype=float)
                self.kpoints  = kpw[:,0:3]
                self.kweights = kpw[:,3].ravel()
                tetline = 3+nkpoints
                if len(lines)>tetline and lines[tetline].lower()[0]=='t':
                    self.tetrahedra = obj()
                    tokens = lines[tetline+1].split()
                    ntets      = int(tokens[0])
                    tet_volume = float(tokens[1])
                    for n in xrange(ntets):
                        tokens = lines[tetline+2+n].split()
                        self.tetrahedra.append(
                            obj(volume     = tet_volume,
                                degeneracy = int(tokens[0]),
                                corners    = array(tokens[1:],dtype=int))
                            )
                    #end for
                #end if
            #end if
        #end if
    #end def read_text


    def write_text(self,filepath=''):
        text = ''
        if self.mode=='auto':
            text+='{0} mesh\n 0\n'.format(self.centering)
            cent = self.centering.lower()
            if len(cent)>0 and cent[0] in Kpoints.centering_options:
                self.centering = Kpoints.centering_options[cent[0]]
            #end if
            if self.centering=='auto':
                text+='auto\n'
                text+=' {0:d}\n'.format(self.kgrid)
            elif self.centering=='gamma' or self.centering=='monkhorst-pack':
                text+='{0}\n'.format(self.centering)
                text+=' {0:d} {1:d} {2:d}\n'.format(*self.kgrid)
                if self.kshift!=None:
                    text+=' {0} {1} {2}\n'.format(*self.kshift)
                #end if
            else:
                self.error('invalid centering for file {0}: {1}\nvalid options are: auto, gamma, monkhorst-pack'.format(filepath,self.centering))
            #end if
        elif self.mode=='basis':
            text+='basis mesh\n 0\n'
            text+='{0}\n'.format(self.coord)
            for kb in self.kbasis:
                text+=' {0:18.14f} {1:18.14f} {2:18.14f}\n'.format(*kb)
            #end for
            text+=' {0:18.14f} {1:18.14f} {2:18.14f}\n'.format(*self.kshift)
        elif self.mode=='line':
            text+='bandstructure\n {0}\nline-mode\n'.format(self.kinsert)
            text+='{0}\n'.format(self.coord)
            npoints = len(self.kendpoints)
            for n in xrange(npoints):
                text+=' {0:18.14f} {1:18.14f} {2:18.14f}   1\n'.format(*self.kendpoints[n])
                if n!=npoints-1 and n%2==1:
                    text+='\n'
                #end if
            #end for
        elif self.mode=='explicit':
            text+='explicit kpoints\n {0}\n'.format(len(self.kpoints))
            text+='{0}\n'.format(self.coord)
            for n in xrange(len(self.kpoints)):
                kp = self.kpoints[n]
                kw = self.kweights[n]
                text+=' {0:18.14f} {1:18.14f} {2:18.14f} {3:12.8f}\n'.format(kp[0],kp[1],kp[2],kw)
            #end for
            if 'tetrahedra' in self and len(self.tetrahedra)>0:
                ntets = len(self.tetrahedra)
                tets = self.tetrahedra
                text+='tetrahedra\n'
                text+=' {0} {1}'.format(ntets,tets[0].volume)
                for n in xrange(ntets):
                    t = tets[n]
                    d = t.degeneracy
                    c = t.corners
                    text+=' {0:d} {1:d} {2:d} {3:d}\n'.format(d,*c)
                #end for
            #end if
        else:
            self.error('invalid mode: {0}\nvalid options are: auto, basis, line, explicit')
        #end if
        return text 
    #end def write_text
#end class Kpoints



class Penaltypot(VFormattedFile):  # metadynamics -> 6.62.4 (2nd one)
    None
#end class Penaltypot



class Poscar(VFormattedFile):

    bool_map = {True:'T',False:'F'}

    def __init__(self,filepath=None):
        self.description = None
        self.scale       = None
        self.axes        = None
        self.elem        = None
        self.elem_count  = None
        self.coord       = None
        self.pos         = None
        self.dynamic     = None
        self.vel_coord   = None
        self.vel         = None
        VFile.__init__(self,filepath)
    #end def __init__

    def change_specifier(self,specifier,vasp_input_class):
        axes=vasp_input_class.poscar.axes
        scale=vasp_input_class.poscar.scale
        unitcellvec=axes*scale

      #units are in angstroms.  
        pos=self.pos
        spec=self.coord  #the current specifier
        
        if spec==specifier:
            #do nothing
            pass
        else:
        
            if spec=="cartesian":
                pass
            elif spec=="direct":
                pos=dot(pos,unitcellvec)
            else:
                raise Exception("Poscar.change_specifier():  %s is not a valid coordinate specifier"%(spec))
                 
            spec=specifier  #the new specifier

            if spec=="cartesian":
                pass # already in cartesian coordinates.   
            elif spec=="direct":
                pos=dot(pos,inv(axes))
            else:
                raise Exception("Poscar.change_specifier():  %s is not a valid coordinate specifier"%(spec))

            self.coord=spec
            self.pos=pos
      
        
    def read_text(self,text,filepath=''):
        lines = self.read_lines(text,remove_empty=False)
        nlines = len(lines)
        min_lines = 8
        if nlines<min_lines:
            self.error('file {0} must have at least {1} lines\n  only {2} lines found'.format(filepath,min_lines,nlines))
        #end if
        description = lines[0]
        dim = 3
        scale = float(lines[1].strip())
        axes = empty((dim,dim))
        axes[0] = array(lines[2].split(),dtype=float)
        axes[1] = array(lines[3].split(),dtype=float)
        axes[2] = array(lines[4].split(),dtype=float)
        tokens = lines[5].split()
        if tokens[0].isdigit():
            counts = array(tokens,dtype=int)
            elem   = None
            lcur   = 6
        else:
            elem   = array(tokens,dtype=str)
            counts = array(lines[6].split(),dtype=int)
            lcur   = 7
        #end if

        if lcur<len(lines) and len(lines[lcur])>0:
            c = lines[lcur].lower()[0]
            lcur+=1
        else:
            self.error('file {0} is incomplete (missing positions)'.format(filepath))
        #end if
        selective_dynamics = c=='s'
        if selective_dynamics: # Selective dynamics
            if lcur<len(lines) and len(lines[lcur])>0:
                c = lines[lcur].lower()[0]
                lcur+=1
            else:
                self.error('file {0} is incomplete (missing positions)'.format(filepath))
            #end if
        #end if
        cartesian = c=='c' or c=='k'
        if cartesian:
            coord = 'cartesian'
        else:
            coord = 'direct'
        #end if
        npos = counts.sum()
        if lcur+npos>len(lines):
            self.error('file {0} is incomplete (missing positions)'.format(filepath))
        #end if
        spos = []
        for i in range(npos):
            spos.append(lines[lcur+i].split())
        #end for
        lcur += npos
        spos = array(spos)
        pos  = array(spos[:,0:3],dtype=float)
        if selective_dynamics:
            dynamic = array(spos[:,3:6],dtype=str)
            dynamic = dynamic=='T'
        else:
            dynamic = None
        #end if
        if lcur<len(lines) and not self.is_empty(lines,lcur):
            cline = lines[lcur].lower()
            lcur+=1
            if lcur+npos>len(lines):
                self.error('file {0} is incomplete (missing velocities)'.format(filepath))
            #end if
            cartesian = len(cline)>0 and (cline[0]=='c' or cline[0]=='k')
            if cartesian:
                vel_coord = 'cartesian'
            else:
                vel_coord = 'direct'
            #end if
            svel = []
            for i in range(npos):
                svel.append(lines[lcur+i].split())
            #end for
            lcur += npos
            vel = array(svel,dtype=float)
        else:
            vel_coord = None
            vel = None
        #end if
        self.set(
            description = description,
            scale       = scale,
            axes        = axes,
            elem        = elem,
            elem_count  = counts,
            coord       = coord,
            pos         = pos,
            dynamic     = dynamic,
            vel_coord   = vel_coord,
            vel         = vel
            )
    #end def read_text


    def write_text(self,filepath=''):
        msg = self.check_complete(exit=False)
        if msg!='':
            self.error('incomplete data to write file {0}\n{1}'.format(filepath,msg))
        #end if
        text = ''
        if self.description is None:
            text += 'System cell and coordinates\n'
        else:
            text += self.description+'\n'
        #end if
        text += ' {0}\n'.format(self.scale)
        for a in self.axes:
            text += ' {0:18.14f} {1:18.14f} {2:18.14f}\n'.format(*a)
        #end for
        if self.elem!=None:
            for e in self.elem:
                iselem,symbol = is_element(e,symbol=True)
                if not iselem:
                    self.error('{0} is not an element'.format(e))
                #end if
                text += symbol+' '
            #end for
            text += '\n'
        #end if
        for ec in self.elem_count:
            text += ' {0}'.format(ec)
        #end for
        text += '\n'
        if self.dynamic!=None:
            text += 'selective dynamics\n'
        #end if
        text += self.coord+'\n'
        if self.dynamic is None:
            for p in self.pos:
                text += ' {0:18.14f} {1:18.14f} {2:18.14f}\n'.format(*p)
            #end for
        else:
            bm = self.bool_map
            for i in xrange(len(self.pos)):
                p = self.pos[i]
                d = self.dynamic[i]
                text += ' {0:18.14f} {1:18.14f} {2:18.14f}  {3}  {4}  {5}\n'.format(p[0],p[1],p[2],bm[d[0]],bm[d[1]],bm[d[2]])
            #end for
        #end if
        if self.vel!=None:
            text += self.vel_coord+'\n'
            for v in self.vel:
                text += ' {0:18.14f} {1:18.14f} {2:18.14f}\n'.format(*v)
            #end for
        #end if
        return text
    #end def write_text


    def check_complete(self,exit=True):
        msg = ''
        if self.scale is None:
            msg += 'scale is missing\n'
        #end if
        if self.axes is None:
            msg += 'axes is missing\n'
        #end if
        if self.elem_count is None:
            msg += 'elem_count is missing\n'
        #end if
        if self.coord is None:
            msg += 'coord is missing\n'
        #end if
        if self.pos is None:
            msg += 'pos is missing\n'
        #end if
        if self.vel!=None and self.vel_coord is None:
            msg += 'vel_coord is missing\n'
        #end if
        if exit:
            self.error(msg)
        #end if
        return msg
    #end def check_complete
#end class Poscar



class NebPoscars(Vobj):
    def read(self,filepath):
        path = self.get_path(filepath)
        dirs = os.listdir(path)
        for d in dirs:
            dpath = os.path.join(path,d)
            if len(d)==2 and d.isdigit() and os.path.isdir(dpath):
                n = int(d)
                poscar = Poscar()
                poscar.read(os.path.join(dpath,'POSCAR'))
                self[n] = poscar
            #end if
        #end for
    #end def read

    def write(self,filepath):
        path = self.get_path(filepath)
        for n in xrange(len(self)):
            neb_path = os.path.join(path,str(n).zfill(2))
            if not os.path.exists(neb_path):
                os.mkdir(neb_path)
            #end if
            poscar = self[n]
            poscar.write(os.path.join(neb_path,'POSCAR'))
        #end for
    #end def write
#end class NebPoscars



class Potcar(VFormattedFile):
    def __init__(self,filepath=None,files=None):
        self.files    = files
        self.filepath = filepath
        self.pseudos  = obj()
        if not os.path.isdir(filepath):
            VFile.__init__(self,filepath)
        else:
            VFile.__init__(self)
        #end if
    #end def __init__


    def read_text(self,text,filepath=''):
        start  = 0
        end    = len(text)
        pstart = start
        pend   = end
        n      = 0
        iter   = 0
        while n<end and iter<20:
            n = text.find('End of Dataset',start,end)
            if n==-1:
                break
            #end if
            start = n
            n=text.find('\n',start,end)+1
            pend = n
            self.pseudos.append(text[pstart:pend])
            pstart = pend
            start  = pend
            iter+=1
        #end while
        if iter>=20:
            self.error('failed to read file {0}'.format(filepath))
        #end if
    #end def read_text


    def write_text(self,filepath=''):
        text = ''
        if len(self.pseudos)>0:
            for i in range(len(self.pseudos)):
                text += self.pseudos[i]
            #end for
        elif self.filepath!=None and self.files!=None:
            for file in self.files:
                text += open(os.path.join(self.filepath,file),'r').read()
            #end for
        #end if
        return text
    #end def write_text


    def pot_info(self):
        pot_info = obj()
        if len(self.pseudos)>0:
            pots = self.pseudos
        elif self.filepath!=None and self.files!=None:
            pots = obj()
            for file in self.files:
                pots.append(open(os.path.join(self.filepath,file),'r').read())
            #end for
        else:
            pots = obj()
        #end if
        for i in range(len(pots)):
            pot = pots[i]

            n1 = pot.find('\n')

            label = pot[0:n1].strip()

            n2 = pot.find('\n',n1+1)
            Zval = int(float(pot[n1:n2].strip()))

            n  = pot.find('VRHFIN')
            n1 = pot.find('=',n+1)+1
            n2 = pot.find(':',n1+1)
            element = pot[n1:n2].strip()

            pot_info.append(obj(label=label,Zval=Zval,element=element))
        #end for
        return pot_info
    #end def pot_info

    
    def label_to_potcar_name(self,label):
        func,elem = label.split()[0:2]
        tag = ''
        if '_' in elem:
            elem,tag = elem.split('_',1)
            tag = '_'+tag
        #end if
        return elem+'.'+func+tag+'.POTCAR'
    #end def label_to_potcar_name


    def load(self):
        self.pseudos.clear()
        if self.filepath!=None and self.files!=None:
            for file in self.files:
                self.pseudos.append(open(os.path.join(self.filepath,file),'r').read())
            #end for
        #end if
    #end def load
#end class Potcar



class Exhcar(VFormattedFile):
    None
#end class Exhcar



class VaspInput(SimulationInput,Vobj):

    all_inputs  = '''
      EXHCAR   ICONST  INCAR  KPOINTS  PENALTYPOT  POSCAR  POTCAR  
      STOPCAR  WAVEDER
      '''.split()
    all_outputs = '''
      CHG     CHGCAR  CONTCAR  DOSCAR   ELFCAR   EIGENVAL HILLSPOT     
      IBZKPT  LOCPOT  OSZICAR  OUTCAR   PCDAT    PRJCAR  
      PROCAR  PROOUT  REPORT   TMPCAR   WAVECAR  XDATCAR  vasprun.xml
      '''.split()# note that CHGCAR, TMPCAR, and WAVECAR sometimes contain input

    input_files = obj(
        #exhcar     = Exhcar,
        #iconst     = Iconst,
        incar      = Incar,
        kpoints    = Kpoints,
        #penaltypot = Penaltypot,
        poscar     = Poscar,
        potcar     = Potcar,
        #stopcar    = Stopcar,
        #waveder    = Waveder
        )

    keyword_files = obj(
        incar   = Incar,
        stopcar = Stopcar
        )

    vasp_save_files = all_inputs + all_outputs

    def __init__(self,filepath=None,prefix='',postfix=''):
        if filepath!=None:
            self.read(filepath,prefix,postfix)
        #end if
    #end def __init__


    def read(self,filepath,prefix='',postfix=''):
        path = self.get_path(filepath)
        for file in os.listdir(path):
            name = file.lower()
            if name in self.input_files:
                filepath = os.path.join(path,prefix+file+postfix)
                self[name] = self.input_files[name](filepath)
            #end if
        #end for
    #end def read


    def write(self,filepath,prefix='',postfix=''):
        path = self.get_path(filepath)
        for name,vfile in self.iteritems():
            filepath = os.path.join(path,prefix+name.upper()+postfix)
            vfile.write(filepath)
        #end for
    #end def write

    def incorporate_system(self,system,incorp_kpoints=True,coord='cartesian'):
        structure = system.structure

        # assign kpoints
        if len(structure.kpoints)>0 and incorp_kpoints:
            kpoints = Kpoints()
            kpoints.mode     = 'explicit'
            kpoints.coord    = 'cartesian'
            kpoints.kpoints  = structure.kpoints.copy()
            kpoints.kweights = structure.kweights.copy()
            self.kpoints = kpoints
        #end if

        # assign poscar
        species = None
        if len(structure.elem)>0:
            s = structure.copy()
            s.change_units('A')
            species,species_count = s.order_by_species()
            poscar = Poscar()
            poscar.scale      = 1.0
            poscar.axes       = s.axes
            poscar.elem       = species
            poscar.elem_count = species_count
            if coord=='cartesian':
                poscar.coord  = 'cartesian'
                poscar.pos    = s.pos
            elif coord=='direct':
                poscar.coord  = 'direct'
                poscar.pos    = s.pos_unit()
            else:
                self.error('coord must be either direct or cartesian\nyou provided: {0}'.format(coord))
            #end if
            if s.frozen!=None:
                poscar.dynamic = s.frozen==False
            #end if
            self.poscar = poscar
        #end if

        # handle charged systems
        #  warning: spin polarization is handled by the user!
        self.incar.nelect = system.particles.count_electrons()

        return species
    #end def incorporate_system
    def return_system(self,**valency):
        axes=self.poscar.axes
        scale=self.poscar.scale
        axes=scale*axes
        scale=1.0
 
        velem=self.poscar.elem
        velem_count=self.poscar.elem_count
        elem=vasp_to_nexus_elem(velem,velem_count)
 
        self.poscar.change_specifier('cartesian',self)
        pos=self.poscar.pos
        
        center=axes.sum(0)/2.0
        
        kpoints=None
        kweights=None
        kgrid=None
        kshift=None

        if self.kpoints.mode=="auto":
         
          kshift=self.kpoints.kshift
          if self.kpoints.centering=="monkhorst-pack":
            kshift=kshift+array([0.5,0.5,0.5])
          elif self.kpoints.centering=="gamma":
            pass
          #end if
        #endif   
          kgrid=self.kpoints.kgrid
         
        else:
          raise Exception("VaspInput.return_system:  Error, system generation doesn't currently\
                                                          work with manually specified k-points")
        structure = Structure(
                      axes    = axes,
                      elem    = elem,
                      scale   = scale,
                      pos     = pos,
                      center  = center,
                      kpoints = kpoints,
                      kweights = kweights,
                      kgrid   = kgrid,
                      kshift  = kshift,
                      units   = 'A',
                      rescale = False
            )
         
        structure.zero_corner()
        structure.recenter()

        ion_charge = 0
        atoms   = list(elem)
        for atom in atoms:
            if not atom in valency:
                self.error('valence charge for atom {0} has not been defined\nplease provide the valence charge as an argument to return_system()'.format(atom))
            #end if
            ion_charge += atoms.count(atom)*valency[atom]
        #end for

        ####WARNING:  Assuming that the netcharge and netspin are ZERO. 
        net_charge=0
        net_spin=0

        system = PhysicalSystem(structure,net_charge,net_spin,**valency)
 
        return system
    #end def return_system

    def set_potcar(self,pseudos,species=None):
        if species is None:
            ordered_pseudos = pseudos
        else:
            pseudo_map = obj()
            for ppname in pseudos:
                element = ppname[0:2].strip('.')
                pseudo_map[element] = ppname
            #end for
            ordered_pseudos = []
            for element in species:
                iselem,symbol = is_element(element,symbol=True)
                if not iselem:
                    self.error('{0} is not an element'.format(element))
                elif not symbol in pseudo_map:
                    self.error('pseudopotential for element {0} not found\nelements present: {1}'.format(symbol,sorted(pseudo_map.keys())))
                #end if
                ordered_pseudos.append(pseudo_map[symbol])
            #end for
        #end if
        self.potcar = Potcar(nexus_noncore.pseudo_dir,ordered_pseudos)
    #end def set_potcar


    def setup_neb(self,*structures,**interp_args):
        # check input types
        if len(structures)==1 and isinstance(structures[0],(list,tuple)):
            structures = structures[0]
        #end if
        for s in structures:
            if not isinstance(s,(Structure,PhysicalSystem)):
                self.error('arguments to setup NEB must either be structure or system objects\n  received an object of type: {0}'.format(s.__class__.__name__))
            #end if
        #end for
        interp_args['repackage'] = False

        # generate NEB image structures
        if len(structures)<2:
            self.error('must provide at least two structures to setup NEB\n  you provided: {0}'.format(len(structures)))
        elif len(structures)==2:
            incar_images = 'images' in self.incar
            kwarg_images = 'images' in interp_args
            if incar_images and kwarg_images and self.incar.images!=interp_args['images']:
                self.error('images provided in incar and setup_neb do not match\n  please ensure they match to remove ambiguity\n  incar images: {0}\n  setup_neb images: {1}'.format(self.incar.images,interp_args['images']))
            elif incar_images:
                interp_args['images'] = self.incar.images
            elif kwarg_images:
                self.incar.images = interp_args['images']
            else:
                self.error('images must be provided in INCAR to setup NEB')
            #end if
            struct1,struct2 = structures
            neb_structures = interpolate_structures(struct1,struct2,**interp_args)
        else:
            if 'images' in interp_args:
                neb_structures = interpolate_structures(structures,**interp_args)
            else:
                neb_structures = structures
            #end if
            if 'images' in self.incar and len(neb_structures)!=self.incar.images+2:
                self.error('number of structures provided to setup_neb must be consistent with number of images in INCAR\n  INCAR images: {0}\n  structures provided {1}'.format(self.incar.images,len(neb_structures)))
            #end if
            self.incar.images = len(neb_structures)-2
        #end if
        
        # create a poscar for each structure and include in input file
        neb_poscars = NebPoscars()
        for n in xrange(len(neb_structures)):
            neb_poscars[n] = generate_poscar(neb_structures[n])
        #end for
        self.poscar = neb_poscars
    #end def setup_neb


    def run_type(self):
        incar = self.incar
        # check for neb
        if 'images' in incar:
            run_type = 'neb'
        elif 'ibrion' in incar and incar.ibrion>0.5:
            run_type = 'relax'
        elif 'ibrion' in incar and incar.ibrion==0:
            run_type = 'md'
        elif 'nsw' in incar and incar.nsw>1.5:
            run_type = 'md'
        else:
            run_type = 'unknown'
        #end if
        return run_type
    #end def run_type


    def producing_structure(self):
        return self.run_type()=='relax'
    #end def producing_structure


    def performing_relax(self):
        return self.run_type()=='relax'
    #end def preforming_relax


    def performing_neb(self):
        return self.run_type()=='neb'
    #end def performing_neb
#end class VaspInput



def generate_vasp_input(**kwargs):
    if 'input_type' in kwargs:
        input_type = kwargs['input_type']
        del kwargs['input_type']
    else:
        input_type = 'general'
    #end if
    if input_type=='general' or input_type=='generic':
        vi = generate_any_vasp_input(**kwargs)
    else:
        VaspInput.class_error('input_type {0} is unrecognized\nvalid options are: general'.format(input_type))
    #end if
    return vi
#end def generate_vasp_input




generate_any_defaults = obj(
    kcenter    = None,
    kpoints    = None,
    kweights   = None,
    kbasis     = None,
    kgrid      = None,
    kshift     = (0,0,0),
    kcoord     = 'cartesian',
    kendpoints = None,
    kinsert    = 10,
    system     = None,
    pseudos    = None,
    neb        = None,
    neb_args   = obj(),
    coord      = 'cartesian'
    )

def generate_any_vasp_input(**kwargs):
    # handle 'system' name collision
    system_str = None
    if 'title' in kwargs:
        system_str = kwargs['title']
        del kwargs['title']
    #end if

    # remove keywords associated with kpoints, poscar, and any other formatted files
    vf = obj()
    for name,default in generate_any_defaults.iteritems():
        if name in kwargs:
            vf[name] = kwargs[name]
            del kwargs[name]
        else:
            vf[name] = default
        #end if
    #end for
    gen_kpoints = not 'kspacing' in kwargs

    # create an empty input file
    vi = VaspInput()

    # assign values to incar and any other keyword files
    keywords = set(kwargs.keys())
    for name,keyword_file in VaspInput.keyword_files.iteritems():
        keys = keywords & keyword_file.keywords
        if len(keys)>0:
            kw = obj()
            kw.move_from(kwargs,keys)
            vfile = keyword_file()
            vfile.assign(**kw)
            vi[name] = vfile
        #end if
    #end for

    # check for leftover keywords
    if len(kwargs)>0:
        VaspInput.class_error('unrecognized keywords: {0}'.format(sorted(kwargs.keys())),'generate_vasp_input')
    #end if

    # incorporate system information
    species = None
    if vf.system!=None:
        species = vi.incorporate_system(vf.system,gen_kpoints,vf.coord)
    #end if

    # set potcar
    if vf.pseudos!=None:
        vi.set_potcar(vf.pseudos,species)
    #end if

    # add kpoints information (override anything provided by system)
    if gen_kpoints and (vf.kpoints!=None or vf.kweights!=None or vf.kbasis!=None or vf.kgrid!=None or vf.kcenter!=None or vf.kendpoints!=None):
        if 'kpoints' in vi:
            kp = vi.kpoints
            kp.clear()
        else:
            kp = Kpoints()
            vi.kpoints = kp
        #end if
        if vf.kpoints!=None:
            kp.mode     = 'explicit'
            kp.kpoints  = vf.kpoints
            kp.kweights = vf.kweights
            kp.coord    = vf.kcoord
        elif vf.kgrid!=None:
            kp.mode      = 'auto'
            kp.centering = vf.kcenter
            if vf.kgrid!=None:
                kp.kgrid = vf.kgrid
            #end if
            if vf.kshift!=None:
                kp.kshift = vf.kshift
            #end if
        elif vf.kendpoints!=None:
            kp.mode       = 'line'
            kp.coord      = vf.kcoord
            kp.kinsert    = vf.kinsert
            kp.kendpoints = vf.kendpoints
        else:
            VaspInput.class_error('could not set kpoints from user inputs','generate_vasp_input')
        #end if
    #end if

    # create many poscars if doing nudged elastic band
    if vf.neb!=None:
        vi.setup_neb(*vf.neb,**vf.neb_args)
    #end if

    # handle 'system' name collision
    if system_str!=None:
        vi.incar.system = system_str
    #end if

    return vi
#end def generate_any_vasp_input




def generate_poscar(structure,coord='cartesian'):
    s = structure.copy()
    s.change_units('A')
    species,species_count = s.order_by_species()
    poscar = Poscar()
    poscar.scale      = 1.0
    poscar.axes       = s.axes
    poscar.elem       = species
    poscar.elem_count = species_count
    if coord=='cartesian':
        poscar.coord  = 'cartesian'
        poscar.pos    = s.pos
    elif coord=='direct':
        poscar.coord  = 'direct'
        poscar.pos    = s.pos_unit()
    else:
        error('coord must be either direct or cartesian\nyou provided: {0}'.format(coord),'generate_poscar')
    #end if
    if s.frozen!=None:
        poscar.dynamic = s.frozen==False
    #end if
    return poscar
#end def generate_poscar
