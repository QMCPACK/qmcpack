
import os
from numpy import array,abs,empty,ndarray
from generic import obj
from simulation import SimulationInput
from developer import DevBase
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




class VFile(DevBase):
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
            if not kw_field in cls.__dict__:
                cls.__dict__[kw_field] = set()
            #end if
        #end for
        #cls.check_consistency()
        cls.scalar_keywords = set()
        for scalar_field in cls.kw_scalars:
            cls.scalar_keywords |= cls.__dict__[scalar_field]
        #end for
        cls.array_keywords = set()
        for array_field in cls.kw_arrays:
            cls.array_keywords |= cls.__dict__[array_field]
        #end for
        cls.keywords = cls.scalar_keywords | cls.array_keywords
        cls.type = obj()
        cls.read_value   = obj()
        cls.write_value  = obj()
        cls.assign_value = obj()
        for type in cls.kw_scalars + cls.kw_arrays:
            for name in cls.__dict__[type]:
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
            untyped -= cls.__dict__[type]
        #end for
        if len(untyped)>0:
            fail = True
            msg += 'variables without a type: {0}\n'.format(sorted(untyped))
        #end if
        for type in types:
            unknown = cls.__dict__[type]-cls.keywords
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


# Dimension of arrays:
#   k-points           NKPTS =     27   k-points in BZ     NKDIM =     27   number of bands    NBANDS=    378
#   number of dos      NEDOS =    301   number of ions     NIONS =     40
#   non local maximal  LDIM  =      6   non local SUM 2l+1 LMDIM =     18
#   total plane-waves  NPLWV = 201600
#   max r-space proj   IRMAX =      1   max aug-charges    IRDMAX= 149001
#   dimension x,y,z NGX =    60 NGY =   60 NGZ =   56
#   dimension x,y,z NGXF=   120 NGYF=  120 NGZF=  112
#   support grid    NGXF=   240 NGYF=  240 NGZF=  224
#   ions per type =               8   2   6  24
# NGX,Y,Z   is equivalent  to a cutoff of  12.66, 12.66, 12.06 a.u.
# NGXF,Y,Z  is equivalent  to a cutoff of  25.33, 25.33, 24.12 a.u.
#
#
# I would recommend the setting:
#   dimension x,y,z NGX =    57 NGY =   57 NGZ =   56
# SYSTEM =  BFO                                     
# POSCAR =   O    Fe   Bi                           
#
# Startparameter for this run:
#   NWRITE =      2    write-flag & timer
#   PREC   = accura    normal or accurate (medium, high low for compatibility)
#   ISTART =      0    job   : 0-new  1-cont  2-samecut
#   ICHARG =      2    charge: 1-file 2-atom 10-const
#   ISPIN  =      1    spin polarized calculation?
#   LNONCOLLINEAR =      T non collinear calculations
#   LSORBIT =      T    spin-orbit coupling
#   INIWAV =      1    electr: 0-lowe 1-rand  2-diag
#   LASPH  =      F    aspherical Exc in radial PAW
#   METAGGA=      F    non-selfconsistent MetaGGA calc.
#
# Electronic Relaxation 1
#   ENCUT  =  500.0 eV  36.75 Ry    6.06 a.u.  14.36 14.36 14.07*2*pi/ulx,y,z
#   ENINI  =  500.0     initial cutoff
#   ENAUG  =  605.4 eV  augmentation charge cutoff
#   NELM   =    300;   NELMIN=  2; NELMDL= -5     # of ELM steps 
#   EDIFF  = 0.1E-07   stopping-criterion for ELM
#   LREAL  =      F    real-space projection
#   NLSPLINE    = F    spline interpolate recip. space projectors
#   LCOMPAT=      F    compatible to vasp.4.4
#   GGA_COMPAT  = T    GGA compatible to vasp.4.4-vasp.4.6
#   LMAXPAW     = -100 max onsite density
#   LMAXMIX     =    4 max onsite mixed and CHGCAR
#   VOSKOWN=      0    Vosko Wilk Nusair interpolation
#   ROPT   =    0.00000   0.00000   0.00000   0.00000
# Ionic relaxation
#   EDIFFG = 0.1E-06   stopping-criterion for IOM
#   NSW    =      0    number of steps for IOM
#   NBLOCK =      1;   KBLOCK =      1    inner block; outer block 
#   IBRION =     -1    ionic relax: 0-MD 1-quasi-New 2-CG
#   NFREE  =      0    steps in history (QN), initial steepest desc. (CG)
#   ISIF   =      2    stress and relaxation
#   IWAVPR =     10    prediction:  0-non 1-charg 2-wave 3-comb
#   ISYM   =      0    0-nonsym 1-usesym 2-fastsym
#   LCORR  =      T    Harris-Foulkes like correction to forces
#
#   POTIM  = 0.5000    time-step for ionic-motion
#   TEIN   =    0.0    initial temperature
#   TEBEG  =    0.0;   TEEND  =   0.0 temperature during run
#   SMASS  =  -3.00    Nose mass-parameter (am)
#   estimated Nose-frequenzy (Omega)   =  0.10E-29 period in steps =****** mass=  -0.142E-26a.u.
#   SCALEE = 1.0000    scale energy and forces
#   NPACO  =    256;   APACO  = 16.0  distance and # of slots for P.C.
#   PSTRESS=    0.0 pullay stress
#
#  Mass of Ions in am
#   POMASS = 208.98 55.85 26.98 16.00
#  Ionic Valenz
#   ZVAL   =  15.00 14.00  3.00  6.00
#  Atomic Wigner-Seitz radii
#   RWIGS  =   1.97  1.97  1.97  1.97
#  virtual crystal weights 
#   VCA    =   1.00  1.00  1.00  1.00
#   NELECT =     310.0000    total number of electrons
#   NUPDOWN=      -1.0000    fix difference up-down
#
# DOS related values:
#   EMIN   =  10.00;   EMAX   =-10.00  energy-range for DOS
#   EFERMI =   0.00
#   ISMEAR =    -5;   SIGMA  =   0.05  broadening in eV -4-tet -1-fermi 0-gaus
#
# Electronic relaxation 2 (details)
#   IALGO  =     38    algorithm
#   LDIAG  =      T    sub-space diagonalisation (order eigenvalues)
#   LSUBROT=      F    optimize rotation matrix (better conditioning)
#   TURBO    =      0    0=normal 1=particle mesh
#   IRESTART =      0    0=no restart 2=restart with 2 vectors
#   NREBOOT  =      0    no. of reboots
#   NMIN     =      0    reboot dimension
#   EREF     =   0.00    reference energy to select bands
#   IMIX   =      4    mixing-type and parameters
#     AMIX     =   0.40;   BMIX     =  1.00
#     AMIX_MAG =   1.60;   BMIX_MAG =  1.00
#     AMIN     =   0.10
#     WC   =   100.;   INIMIX=   1;  MIXPRE=   1;  MAXMIX= -45
#
# Intra band minimization:
#   WEIMIN = 0.0000     energy-eigenvalue tresh-hold
#   EBREAK =  0.66E-11  absolut break condition
#   DEPER  =   0.30     relativ break condition  
#
#   TIME   =   0.40     timestep for ELM
#
#  volume/ion in A,a.u.               =      11.97        80.80
#  Fermi-wavevector in a.u.,A,eV,Ry     =   1.416120  2.676079 27.285078  2.005397
#  Thomas-Fermi vector in A             =   2.537488
# 
# Write flags
#   LWAVE  =      F    write WAVECAR
#   LCHARG =      T    write CHGCAR
#   LVTOT  =      F    write LOCPOT, total local potential
#   LVHAR  =      F    write LOCPOT, Hartree potential only
#   LELF   =      F    write electronic localiz. function (ELF)
#   LORBIT =     11    0 simple, 1 ext, 2 COOP (PROOUT)
#
#
# Dipole corrections
#   LMONO  =      F    monopole corrections only (constant potential shift)
#   LDIPOL =      F    correct potential (dipole corrections)
#   IDIPOL =      0    1-x, 2-y, 3-z, 4-all directions 
#   EPSILON=  1.0000000 bulk dielectric constant
#
# LDA+U is selected, type is set to LDAUTYPE =  1
#   angular momentum for each species LDAUL =    -1    2   -1   -1
#   U (eV)           for each species LDAUU =   0.0  2.0  0.0  0.0
#   J (eV)           for each species LDAUJ =   0.0  0.0  0.0  0.0
# Exchange correlation treatment:
#   GGA     =    --    GGA type
#   LEXCH   =     2    internal setting for exchange type
#   VOSKOWN=      0    Vosko Wilk Nusair interpolation
#   LHFCALC =     F    Hartree Fock is set to
#   LHFONE  =     F    Hartree Fock one center treatment
#   AEXX    =    0.0000 exact exchange contribution
#
# Linear response parameters
#   LEPSILON=     F    determine dielectric tensor
#   LRPA    =     F    only Hartree local field effects (RPA)
#   LNABLA  =     F    use nabla operator in PAW spheres
#   LVEL    =     F    velocity operator in full k-point grid
#   LINTERFAST=   F  fast interpolation
#   KINTER  =     0    interpolate to denser k-point grid
#   CSHIFT  =0.1000    complex shift for real part using Kramers Kronig
#   OMEGAMAX=  -1.0    maximum frequency
#   DEG_THRESHOLD= 0.2000000E-02 threshold for treating states as degnerate
#   RTIME   =    0.100 relaxation time in fs
#
# Orbital magnetization related:
#   ORBITALMAG=     F  switch on orbital magnetization
#   LCHIMAG   =     F  perturbation theory with respect to B field
#   DQ        =  0.001000  dq finite difference perturbation B field





 
class Incar(VKeywordFile):

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
      ialgo iband ibrion icharg ichibare i_constrained_m icorelevel idipol 
      igpar images imix increm inimix iniwav ipead isif ismear ispin istart 
      isym ivdw iwavpr 
      kblock kgamma kpar kpuse kspacing 
      lambda langevin_gamma langevin_gamma_l lasph lasync lattice_constraints 
      lberry lblueout lcalceps lcalcpol lcharg lchimag lcorr ldau ldauj ldaul 
      ldauprint ldautype ldauu ldiag ldipol lefg lelf lepsilon lhfcalc lhyperfine 
      lkproj lmaxfock lmaxfockae lmaxfockmp2 lmaxmix lmaxmp2 lmaxpaw lmaxtau 
      lmixtau lmono lnabla lnmr_sym_red lnoncollinear loptics lorbit lpard 
      lpead lplane lreal lrpa lscalapack lscaler0 lscalu lscsgrad lselfenergy 
      lsepb lsepk lspectral lsorbit lthomas luse_vdw lvdw lvdw_ewald lvdwscs 
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
      tebeg teend time tsubsys
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
      ialgo ibrion icharg ichibare i_constrained_m icorelevel idipol igpar 
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
      kspacing 
      lambda langevin_gamma_l
      mbja mbjb minrot
      nelect 
      ofield_a ofield_kappa ofield_q6_far ofield_q6_near omegamax omegamin omegatl
      param1 param2 pmass pomass potim pstress 
      scsrad shaketol sigma step_max step_size symprec 
      tebeg teend time 
      vdw_a1 vdw_a2 vdw_cnradius vdw_d vdw_radius vdw_scaling vdw_sr vdw_s6 vdw_s8
      wc weimin 
      zab_vdw zval 
      '''.split())

    bools = set('''
      addgrid
      evenonly evenonlygw
      gga_compat 
      lasph lasync lberry lblueout lcalceps lcalcpol lcharg lchimag lcorr 
      ldau ldiag ldipol lefg lelf lepsilon lhfcalc lhyperfine lkproj lmaxtau 
      lmixtau lmono lnabla lnmr_sym_red lnoncollinear loptics lpard lpead 
      lplane lrpa lscalapack lscaler0 lscalu lscsgrad lselfenergy lsepb 
      lsepk lsorbit lspectral lthomas luse_vdw lvdw lvdw_ewald lvdwscs lvhar 
      lvtot lwave 
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
    #    coord    = cartesian/reciprocal
    #    ninsert   = number of points inserted between each set of endpoints
    #    endpoints = kpoint pairs forming line endpoints
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
                self.ninsert = iselect
                self.coord   = self.coord_options(lines[3].lower()[0])
                endpoints = []
                for line in lines[4:]:
                    endpoints.append(line.split())
                #end for
                self.endpoints = array(endpoints,dtype=float)
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
            text+='kpoints along lines\n {0}\nline-mode\n'.format(self.ninsert)
            text+='{0}\n'.format(self.coord)
            npoints = len(self.endpoints)
            for n in xrange(npoints):
                text+=' {0:18.14f} {1:18.14f} {2:18.14f}\n'.format(*self.endpoints[n])
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
                text += e+' '
            #end for
            text += '\n'
        #end if
        for ec in self.elem_count:
            text += ' {0}'.format(ec)
        #end for
        text += '\n'
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
            n2 = pot.find('\n',n1+1)
            Zval = int(float(pot[n1:n2].strip()))

            n  = pot.find('VRHFIN')
            n1 = pot.find('=',n+1)+1
            n2 = pot.find(':',n1+1)
            element = pot[n1:n2].strip()

            pot_info.append(obj(Zval=Zval,element=element))
        #end for
        return pot_info
    #end def pot_info


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




class VaspInput(SimulationInput):

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
        path,tmp = os.path.split(filepath)
        if len(path)>0 and not os.path.exists(path):
            self.error('path {0} does not exist'.format(path))
        #end if
        for file in os.listdir(path):
            name = file.lower()
            if name in self.input_files:
                filepath = os.path.join(path,prefix+file+postfix)
                self[name] = self.input_files(filepath)
            #end if
        #end for
    #end def read


    def write(self,filepath,prefix='',postfix=''):
        path,tmp = os.path.split(filepath)
        if len(path)>0 and not os.path.exists(path):
            self.error('path {0} does not exist'.format(path))
        #end if
        for name,vfile in self.iteritems():
            filepath = os.path.join(path,prefix+name.upper()+postfix)
            vfile.write(filepath)
        #end for
    #end def write


    def incorporate_system(self,system,incorp_kpoints=True):
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
            species,species_count = s.order_by_species()
            poscar = Poscar()
            poscar.scale      = 1.0
            poscar.axes       = s.axes
            poscar.elem       = species
            poscar.elem_count = species_count
            poscar.coord      = 'cartesian'
            poscar.pos        = s.pos
            if 'frozen' in structure:
                poscar.dynamic = s.frozen==False
            #end if
            self.poscar = poscar
        #end if

        # handle charged systems
        #  warning: spin polarization is handled by the user!
        self.incar.nelect = system.particles.count_electrons()

        return species
    #end def incorporate_system


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
                if not element in pseudo_map:
                    self.error('pseudopotential for element {0} not found\nelements present: {1}\n'.format(element,sorted(pseudo_map.keys())))
                #end if
                ordered_pseudos.append(pseudo_map[element])
            #end for
        #end if
        self.potcar = Potcar(VaspInput.pseudo_dir,ordered_pseudos)
    #end def set_potcar
#end class VaspInput



def generate_vasp_input(**kwargs):
    if 'input_type' in kwargs:
        input_type = kwargs['input_type']
        del kwargs['input_type']
    else:
        input_type = 'general'
    #end if
    if input_type=='general':
        vi = generate_any_vasp_input(**kwargs)
    else:
        VaspInput.class_error('input_type {0} is unrecognized\nvalid options are: general'.format(input_type))
    #end if
    return vi
#end def generate_vasp_input




generate_any_defaults = obj(
    kcenter  = None,
    kpoints  = None,
    kweights = None,
    kbasis   = None,
    kgrid    = None,
    kshift   = (0,0,0),
    kcoord   = 'cartesian',
    system   = None,
    pseudos  = None
    )

def generate_any_vasp_input(**kwargs):
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

    gen_kpoints = not 'kspacing' in vf

    # incorporate system information
    species = None
    if vf.system!=None:
        species = vi.incorporate_system(vf.system,gen_kpoints)
    #end if

    # set potcar
    if vf.pseudos!=None:
        vi.set_potcar(vf.pseudos,species)
    #end if

    # add kpoints information (override anything provided by system)
    if gen_kpoints and (vf.kpoints!=None or vf.kweights!=None or vf.kbasis!=None or vf.kgrid!=None or vf.kcenter!=None):
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
        else:
            VaspInput.class_error('could not set kpoints from user inputs','generate_vasp_input')
        #end if
    #end if

    return vi
#end def generate_any_vasp_input
