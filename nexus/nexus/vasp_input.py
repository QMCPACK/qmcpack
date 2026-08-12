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
from copy import deepcopy
from functools import partial
from numbers import Integral,Real
from types import MappingProxyType
import numpy as np
from .periodic_table import Elements
from .nexus_base import nexus_noncore
from .simulation import SimulationInput
from .structure import interpolate_structures, Structure
from .physical_system import PhysicalSystem
from .developer import DevBase, obj, error
from .utilities import path_string
from . import numpy_extensions as npe

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
            n = int(n)
            if n<0:
                raise ValueError('array repeat count must be non-negative')
            #end if
            sarr.extend(n*[vv])
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
    return float(sval.lower().replace('d','e'))
#end def read_real


bool_dict = dict(true=True,false=False,t=True,f=False)
def read_bool(sval):
    sval = remove_comment(sval)
    return bool_dict[sval.lower().strip('.')]
#end def read_bool


def read_string(sval):
    return sval
#end def read_string


def read_int_array(sval):
    return np.array(expand_array(sval),dtype=int)
#end def read_int_array


def read_real_array(sval):
    values = [v.lower().replace('d','e') for v in expand_array(sval)]
    return np.array(values,dtype=float)
#end def read_real_array


bool_array_dict = dict(T=True,F=False,TRUE=True,FALSE=False)
def read_bool_array(sval):
    barr = []
    for v in expand_array(sval):
        barr.append(bool_array_dict[v.upper().strip('.')])
    #end for
    return np.array(barr,dtype=bool)
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
    quote = (
        '\n' in v or ';' in v or '#' in v or '!' in v or
        v.endswith('\\') or v != v.strip()
        )
    if quote:
        return '"'+v+'"'
    else:
        return v
    #end if
#end def write_string


def equality(a,b):
    return a==b
#end def equality


def render_bool(v):
    if v:
        return 'T'
    else:
        return 'F'
    #end if
#end def render_bool


def write_array(arr,same=equality,render=str,max_repeat=3):
    if len(arr)==0:
        return ''
    #end if
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
            for i in range(count):  # noqa: B007
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
    return write_array(a)
#end def write_real_array


def write_bool_array(a):
    return write_array(a,render=render_bool)
#end def write_bool_array


assign_bool_map = {True:True,False:False,1:True,0:False}
def assign_bool(v):
    return assign_bool_map[v]
#end def assign_bool


integer_tolerance = 1e-8
max_exact_float_integer = 2**53
def assign_int(v):
    if isinstance(v,(bool,np.bool_)):
        raise TypeError('value must be an integer, not a boolean')
    elif isinstance(v,(Integral,np.integer)):
        return int(v)
    elif isinstance(v,(Real,np.floating)):
        value = float(v)
        if not np.isfinite(value):
            raise ValueError('value must be finite')
        elif abs(value)>max_exact_float_integer:
            raise ValueError(
                'floating-point integer values must not exceed {0}'
                .format(max_exact_float_integer)
                )
        #end if
        nearest = round(value)
        if not np.isclose(value,nearest,rtol=0.0,atol=integer_tolerance):
            raise ValueError(
                'real value {0} is not within {1} of an integer'
                .format(v,integer_tolerance)
                )
        #end if
        return int(nearest)
    else:
        raise TypeError('value must be an integer or real number')
    #end if
#end def assign_int


def assign_string(v):
    if isinstance(v,str):
        return v
    else:
        raise ValueError('value must be a string')
    #end if
#end def assign_string


def assign_int_array(a):
    if isinstance(a,(tuple,list,np.ndarray)):
        array = np.asarray(a,dtype=object)
        if array.ndim!=1:
            raise ValueError('value must be a one-dimensional array')
        #end if
        values = []
        for index,value in enumerate(array):
            try:
                values.append(assign_int(value))
            except (TypeError,ValueError) as e:
                raise type(e)(
                    'invalid integer array element at index {0}: {1}'
                    .format(index,e)
                    ) from e
            #end try
        #end for
        try:
            return np.array(values,dtype=int)
        except (OverflowError,TypeError,ValueError) as e:
            raise ValueError(
                'integer array values cannot be represented: {0}'.format(e)
                ) from e
        #end try
    else:
        raise ValueError('value must be a tuple, list, or array')
    #end if
#end def assign_int_array


def assign_real_array(a):
    if isinstance(a,(tuple,list,np.ndarray)):
        return np.array(a,dtype=float)
    else:
        raise ValueError('value must be a tuple, list, or array')
    #end if
#end def assign_real_array


def assign_bool_array(a):
    if isinstance(a,(tuple,list,np.ndarray)):
        return np.array(a,dtype=bool)
    else:
        raise ValueError('value must be a tuple, list, or array')
    #end if
#end def assign_bool_array


#utility functions to convert from VASP internal objects to 
#nexus objects:
def vasp_to_nexus_elem(elem,elem_count):
    syselem=[]
    for x,count in zip(elem,elem_count):
        syselem+=[x for i in range(0,count)]
    #end for
    return np.array(syselem)
#end def vasp_to_nexus_elem


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
    ints        = assign_int,
    reals       = float,
    bools       = assign_bool,
    strings     = assign_string,
    int_arrays  = assign_int_array,
    real_arrays = assign_real_array,
    bool_arrays = assign_bool_array    
    )


def mixed_type_matches(value,value_type):
    if value_type=='strings':
        return isinstance(value,str)
    elif value_type=='bools':
        return isinstance(value,(bool,np.bool_,Integral,np.integer))
    elif value_type in ('ints','reals'):
        return (
            isinstance(value,(Real,np.integer,np.floating)) and
            not isinstance(value,(bool,np.bool_))
            )
    elif value_type in ('int_arrays','real_arrays','bool_arrays'):
        return isinstance(value,(tuple,list,np.ndarray))
    else:
        raise ValueError('unknown keyword value type: {0}'.format(value_type))
    #end if
#end def mixed_type_matches


def assign_mixed(types,value):
    errors = []
    for value_type in types:
        if mixed_type_matches(value,value_type):
            try:
                return assign_value_functions[value_type](value)
            except Exception as e:
                errors.append('{0}: {1}'.format(value_type,e))
            #end try
        #end if
    #end for
    message = 'value does not match permitted types: {0}'.format(types)
    if len(errors)>0:
        message += '\n'+'\n'.join(errors)
    #end if
    raise ValueError(message)
#end def assign_mixed


def read_mixed(types,sval):
    scalar_types = {'ints','reals','bools'}
    array_types = {'int_arrays','real_arrays','bool_arrays'}
    try:
        nvalues = len(expand_array(sval))
    except Exception:
        nvalues = len(sval.split())
    #end try
    if nvalues>1:
        ordered_types = (array_types,)
    else:
        ordered_types = scalar_types,array_types
    #end if
    errors = []
    for type_group in ordered_types:
        for value_type in types:
            if value_type in type_group:
                try:
                    return read_value_functions[value_type](sval)
                except Exception as e:
                    errors.append('{0}: {1}'.format(value_type,e))
                #end try
            #end if
        #end for
    #end for
    if nvalues==1 and 'strings' in types:
        return read_string(sval)
    #end if
    raise ValueError(
        'value does not match permitted types: {0}\n{1}'
        .format(types,'\n'.join(errors))
        )
#end def read_mixed


def write_mixed(types,value):
    converted = assign_mixed(types,value)
    for value_type in types:
        if mixed_type_matches(converted,value_type):
            return write_value_functions[value_type](converted)
        #end if
    #end for
    raise ValueError('value does not match permitted types: {0}'.format(types))
#end def write_mixed





class Vobj(DevBase):
    def get_path(self,filepath):
        filepath = path_string(filepath)
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
        if filepath is not None:
            filepath = path_string(filepath)
            self.read(filepath)
        #end if
    #end def __init__


    def read(self,filepath):
        if not os.path.exists(filepath):
            self.error('file {0} does not exist'.format(filepath))
        #end if
        with open(filepath, "r") as f:
            text = f.read()
        self.read_text(text,filepath)
        return text
    #end def read


    def write(self,filepath=None):
        text = self.write_text(filepath)
        if filepath is not None:
            with open(filepath, "w") as f:
                f.write(text)
        #end if
        return text
    #end def write


    def read_text(self,text,filepath=''):
        raise NotImplementedError
    #end def read_text


    def write_text(self,filepath=''):
        raise NotImplementedError
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

    def preprocess_multiline_strings(self,text):
        mvals = obj()
        if '"' in text:
            text_in = text
            text = ''
            plocs = []
            i = 0
            istart = 0
            n = 0
            nqmax = 101
            while i!=-1 and n<nqmax:
                i = text_in.find('"',istart)
                plocs.append(i)
                istart = i+1
                n+=1
            #end while
            if len(plocs)>0:
                plocs.pop()
            #end if
            if n>=nqmax:
                self.error('max number of multi-line strings exceeded.\nOver {} quotation marks found in file.'.format(nqmax-1))
            #end if
            if len(plocs)%2!=0:
                self.error('quotation marks for multi-line strings are not paired')
            #end if
            nlabel = 0
            istart = 0
            for n,i in enumerate(plocs):
                if n%2==0:
                    text += text_in[istart:i]
                else:
                    q = text_in[istart:i]
                    label = '--multiline{}--'.format(str(nlabel).zfill(3))
                    text += label+'\n'
                    mvals[label] = q
                    nlabel += 1
                #end if
                istart = i+1
            #end for
            text += text_in[istart:]
        #end if
        return text,mvals
    #end def preprocess_multiline_strings
#end class VFile



class VKeywordFile(VFile):
    kw_scalars = ('ints','reals','bools','strings')
    kw_arrays  = ('int_arrays','real_arrays','bool_arrays')
    kw_fields  = kw_scalars + kw_arrays + ('keywords','unsupported')

    keyword_classification = None
    mixed_types = MappingProxyType({})

    @classmethod
    def class_init(cls):
        cls.kw_scalars = VKeywordFile.kw_scalars
        cls.kw_arrays  = VKeywordFile.kw_arrays
        cls.kw_fields  = VKeywordFile.kw_fields
        for kw_field in cls.kw_fields:
            if not hasattr(cls,kw_field):
                setattr(cls,kw_field,set())
        #end for
        cls.scalar_keywords = set()
        for scalar_field in cls.kw_scalars:
            cls.scalar_keywords |= getattr(cls,scalar_field)
        #end for
        cls.array_keywords = set()
        for array_field in cls.kw_arrays:
            cls.array_keywords |= getattr(cls,array_field)
        #end for
        cls.keywords = cls.scalar_keywords | cls.array_keywords
        cls.type = obj()
        cls.read_value   = obj()
        cls.write_value  = obj()
        cls.assign_value = obj()
        for type in cls.kw_scalars + cls.kw_arrays:
            for name in getattr(cls,type):
                cls.type[name] = type
                cls.read_value[name]   = read_value_functions[type]
                cls.write_value[name]  = write_value_functions[type]
                cls.assign_value[name] = assign_value_functions[type]
            #end for
        #end for
        for name,types in cls.mixed_types.items():
            if name not in cls.keywords:
                raise ValueError(
                    'mixed-type keyword is not classified: {0}'.format(name)
                    )
            #end if
            classified_type = cls.type[name]
            if classified_type not in types:
                raise ValueError(
                    'classified type {0} is not permitted for keyword {1}'
                    .format(classified_type,name)
                    )
            #end if
            for value_type in types:
                if value_type not in read_value_functions:
                    raise ValueError(
                        'invalid type {0} for keyword {1}'
                        .format(value_type,name)
                        )
                #end if
            #end for
            cls.type[name] = types
            cls.read_value[name] = partial(read_mixed,types)
            cls.write_value[name] = partial(write_mixed,types)
            cls.assign_value[name] = partial(assign_mixed,types)
        #end for
    #end def class_init


    def read_text(self,text,filepath=''):
        text,multiline_values = self.preprocess_multiline_strings(text)
        lines = text.splitlines()
        expression = None
        continued  = False
        for line in lines:
            ls = line.strip()
            if len(ls)>0 and ls[0]!='!' and ls[0]!='#':
                ls = self.remove_comment(ls)
                this_cont = ls.endswith('\\')
                if this_cont:
                    ls = ls.rstrip('\\')+' '
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
                                if value.startswith('--multiline'):
                                    value = multiline_values[value]
                                #end if
                                try:
                                    value = self.read_value[name](value)
                                    self[name] = value
                                except Exception as e:
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
        if continued:
            self.error('incomplete line continuation at end of file')
        #end if
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
            except Exception as e:
                self.error('write failed for file {0} keyword {1}\nkeyword type: {2}\nvalue: {3}\nexception:\n{4}'.format(filepath,name,self.type[name],value,e))
            #end try
            text += valfmt.format(name.upper(),svalue)
        #end for
        return text
    #end def write_text


    def assign(self,**values):
        for name,value in values.items():
            if name not in self.keywords:
                self.error(
                    '{0} is not a keyword for the {1} file'
                    .format(name.upper(),self.__class__.__name__.upper())
                    )
            #end if
            try:
                self[name] = self.assign_value[name](value)
            except Exception as e:
                self.error('assign failed for keyword {0}\nkeyword type: {1}\nvalue: {2}\nexception:\n{3}'.format(name,self.type[name],value,e))
            #end try
        #end for
    #end def assign
#end class VKeywordFile



class VFormattedFile(VFile):
    def read_lines(self,text,*,remove_empty=False):
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
        for iline in range(first_line,last_line):
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

    # VASP wiki with incar keys/tags
    #   https://www.vasp.at/wiki/index.php/Category:INCAR_tag

    # VTST extensions:  https://henkelmangroup.github.io/vtsttools/

    # Block-style INCAR constructs cannot be represented by the flat
    # keyword-file reader/writer.
    unsupported = frozenset({
        'image_1', 'kernel_truncation/factor',
        'kernel_truncation/idimensionality', 'kernel_truncation/ipad',
        'kernel_truncation/isurface', 'kernel_truncation/lcoarsen',
        'kernel_truncation/ltruncate', 'plugins/force_and_stress',
        'plugins/local_potential', 'plugins/machine_learning',
        'plugins/ml_mode', 'plugins/ml_outblock', 'plugins/ml_output_mode',
        'plugins/neighbor_cutoff', 'plugins/occupancies', 'plugins/structure',
        })

    # Pages still present in Nexus history but absent or malformed on the wiki.
    broken_docs = frozenset({
        'dmft_basis',
        })

    ints = frozenset({
        'antires', 'ch_nedos', 'cll', 'cln', 'clnt', 'drotmax', 'efermi_nedos',
        'elmin', 'elph_fermi_nedos', 'elph_ismear', 'elph_lr', 'elph_nbands',
        'elph_selfen_band_start_kp', 'elph_selfen_band_stop_kp',
        'elph_selfen_nw', 'elph_transport_driver', 'elph_transport_nedos_plot',
        'elph_wf_comm_opt', 'exxoep', 'findiff', 'fmp_direction', 'fmp_period',
        'fmp_snumber', 'fmp_swapnum', 'fnmin', 'fockcorr', 'hflmax', 'hflmaxf',
        'hills_bin', 'i_constrained_m', 'ialgo', 'iall_in_one', 'ibrion',
        'ibse', 'ichain', 'icharg', 'ichibare', 'icorelevel', 'idipol',
        'iepsilon', 'ifc_asr', 'ifc_lr', 'igpar', 'ilbfgsmem', 'images',
        'imix', 'inimix', 'iniwav', 'inmrprint', 'iopt', 'ipead',
        'irc_direction', 'irc_stop', 'isearch', 'isif', 'ismear', 'ispin',
        'istart', 'isym', 'ivdw', 'ivdw_nl', 'iwavpr', 'kblock', 'kpar',
        'kpoints_opt_mode', 'kpoints_opt_nkbatch', 'ldauprint', 'ldautype',
        'libmbd_n_omega_grid', 'lmaxfock', 'lmaxfockae', 'lmaxfockmp2',
        'lmaxmix', 'lmaxmp2', 'lmaxpaw', 'lmaxtau', 'lorbit', 'maxmem',
        'maxmix', 'mdalgo', 'mixpre',
        'ml_calgo', 'ml_desc_type', 'ml_estblock', 'ml_ff_icouple_mb',
        'ml_ff_ireg_mb', 'ml_ff_istart', 'ml_ff_lmax2_mb', 'ml_ff_mrb1_mb',
        'ml_ff_mrb2_mb', 'ml_ff_natom_coupled_mb', 'ml_iafilt2',
        'ml_ialgo_linreg', 'ml_icriteria', 'ml_ireg', 'ml_iscale_toten',
        'ml_istart', 'ml_iweight', 'ml_lmax2', 'ml_mb', 'ml_mb_min',
        'ml_mconf', 'ml_mconf_new', 'ml_mhis', 'ml_mopot_nm', 'ml_mrb1',
        'ml_mrb2', 'ml_natom_coupled', 'ml_ncshmem', 'ml_nhyp', 'ml_nmdint',
        'ml_nrank_sparsdes', 'ml_outblock', 'ml_output_mode', 'ml_srpot_n0',
        'naturalo', 'nbands', 'nbands_wave', 'nbandsexact', 'nbandsgw',
        'nbandso', 'nbandsv', 'nblk', 'nblock', 'nblock_fock', 'nbmod',
        'nbseblocko', 'nbseblockv', 'nbseeig', 'ncore', 'ncore_in_image1',
        'ncshmem', 'nedos', 'nelm', 'nelmall', 'nelmdl', 'nelmgw', 'nelmin',
        'nfree', 'ngx', 'ngxf', 'ngy', 'ngyf', 'ngz', 'ngzf', 'nhc_nchains',
        'nhc_nrespa', 'nhc_ns', 'nkred', 'nkredx', 'nkredy', 'nkredz',
        'nmaxfockae', 'nomega', 'nomega_dump', 'nomegapar', 'nomegar', 'npaco',
        'npar', 'nppstr', 'nrmm', 'nsim', 'nstorb', 'nsw', 'ntaupar',
        'ntemper', 'num_wann', 'nwrite', 'phon_dos', 'phon_nedos',
        'phon_nstruct', 'phon_ntlist', 'phon_nwrite', 'plevel', 'proutine',
        'pyamff_maxepoch', 'pyamff_maxiter', 'shakemaxiter', 'snl',
        'transport_nedos', 'voskown', 'wrt_nmrcur',
        })

    reals = frozenset({
        'aexx', 'aggac', 'aggax', 'aldac', 'aldax', 'alpha_vdw', 'amggac',
        'amggax', 'amin', 'amix', 'amix_mag', 'andersen_prob', 'apaco', 'bexx',
        'bmix', 'bmix_mag', 'bparam', 'ch_amplification', 'ch_sigma', 'clz',
        'cmbja', 'cmbjb', 'cmbje', 'cparam', 'cshift', 'ddr', 'deg_threshold',
        'deper', 'dfnmax', 'dfnmin', 'dimer_dist', 'dq', 'ebreak', 'ediff',
        'ediffg', 'efield', 'elph_kspacing', 'elph_selfen_band_start',
        'elph_selfen_band_stop', 'elph_selfen_broad_tol', 'elph_selfen_wrange',
        'elph_transport_dfermi_tol', 'elph_transport_emax',
        'elph_transport_emax_plot', 'elph_transport_emin',
        'elph_transport_emin_plot', 'elph_wf_cache_mb', 'emax', 'emin',
        'enaug', 'encut', 'encutfock', 'encutgw', 'encutgwsoft', 'encutlr',
        'enini', 'enmax', 'enmin', 'epsilon', 'estop', 'ewald_cutoff',
        'falpha', 'fdstep', 'ftimedec', 'ftimeinc', 'ftimemax', 'gamma_vdw',
        'hfalpha', 'hfrcut', 'hfscreen', 'hills_h', 'hills_w', 'hitoler',
        'invcurv', 'irc_delta0', 'irc_maxstep', 'irc_minstep', 'irc_vnorm0',
        'jacobian', 'kspacing', 'kspacing_opt', 'lambda', 'lanczosthr',
        'langevin_gamma_l', 'libmbd_k_grid_shift', 'libmbd_mbd_a',
        'libmbd_mbd_beta', 'libmbd_ts_d', 'libmbd_ts_sr', 'libxc1_pn',
        'libxc2_pn', 'maxdis', 'maxmove', 'mbja', 'mbjb', 'minrot',
        'ml_afilt2', 'ml_cdoub', 'ml_csig', 'ml_cslope', 'ml_ctifor', 'ml_cx',
        'ml_emppot_rcut', 'ml_eps_low', 'ml_eps_reg', 'ml_ff_rcouple_mb',
        'ml_ff_rcut1_mb', 'ml_ff_rcut2_mb', 'ml_ff_sion1_mb', 'ml_ff_sion2_mb',
        'ml_ff_w1_mb', 'ml_ff_w2_mb', 'ml_mopot_dm', 'ml_mopot_rkm',
        'ml_mopot_rm', 'ml_rcouple', 'ml_rcut1', 'ml_rcut2',
        'ml_rdes_sparsdes', 'ml_sclc_ctifor', 'ml_sigv0', 'ml_sigw0',
        'ml_sion1', 'ml_sion2', 'ml_srpot_b0', 'ml_srpot_s0', 'ml_w1',
        'ml_wtifor', 'ml_wtoten', 'ml_wtsif', 'msdgw_f', 'nelect', 'nupdown',
        'ofield_a', 'ofield_kappa', 'ofield_q6_far', 'ofield_q6_near',
        'omegamax', 'omegamin', 'omegatl', 'param1', 'param2', 'phon_g_cutoff',
        'phon_sigma', 'pmass', 'potim', 'pstress', 'pthreshold', 'pyamff_etol',
        'pyamff_fcoeff', 'pyamff_ftol', 'pyamff_maxmove', 'pyamff_swftol',
        'pyamff_tol', 'rsmbj', 'scalee', 'scissor', 'scsrad', 'sdalpha', 'sdr',
        'shaketol', 'shaketolsoft', 'sigma', 'sltol', 'smass', 'smbj',
        'spring', 'step_max', 'step_size', 'symprec', 'tebeg', 'teend',
        'tilambda', 'time', 'timestep', 'transport_relaxation_time',
        'vacpotflat', 'vcaimages', 'vcutoff', 'vdw_a1', 'vdw_a2', 'vdw_beta',
        'vdw_cnradius', 'vdw_d', 'vdw_radius', 'vdw_s6', 'vdw_s8', 'vdw_s9',
        'vdw_scaling', 'vdw_sr', 'vdw_sr8', 'wc', 'weimin', 'xcm_pn',
        'zab_vdw', 'zval',
        })

    bools = frozenset({
        'addgrid', 'ch_lspec', 'elph_ignore_imag_phonons', 'elph_pot_generate',
        'elph_prepare', 'elph_run', 'elph_selfen_dw', 'elph_selfen_fan',
        'elph_selfen_g_skip', 'elph_selfen_gaps', 'elph_selfen_imag_skip',
        'elph_selfen_static', 'elph_transport', 'elph_useblas',
        'elph_userecip', 'elph_wf_cache_prefill', 'elph_wf_redistribute',
        'elph_write_hdf5vel', 'elph_write_textvel', 'evenonly', 'evenonlygw',
        'gga_compat', 'interactive', 'kgamma', 'kpoints_opt', 'ladder',
        'laechg', 'lall_in_one', 'lasph', 'lasync', 'lautoscale', 'lberry',
        'lblueout', 'lbone', 'lcalceps', 'lcalcpol', 'lcharg', 'lchargh5',
        'lchimag', 'lclimb', 'lcorr', 'ldau', 'ldiag', 'ldipol',
        'ldisentangle', 'ldisentangled', 'ldmatrix', 'ldneb', 'ldownsample',
        'lefg', 'lelf', 'lepsilon', 'lfermigw', 'lfinite_temperature',
        'lfockace', 'lfockaedft', 'lfockstd', 'lfxc', 'lglobal', 'lh5',
        'lhartree', 'lhfcalc', 'lhyperfine', 'lintpol_kpath', 'lkpoints_opt',
        'lkpoints_wan', 'lkproj', 'llineopt', 'llraug', 'lmixtau', 'lmodelhf',
        'lmono', 'lmp2lt', 'lnabla', 'lnebcell', 'lnicsall', 'lnlrpa',
        'lnmr_sym_red', 'lnmrcar', 'lnmrleg', 'lnmrshield', 'lnoaugxc',
        'lnoncollinear', 'loptics', 'lorbitalreal', 'lorbmom', 'lpard',
        'lpardh5', 'lpead', 'lpead_sym_red', 'lphon_dispersion', 'lphon_polar',
        'lphon_read_force_constants', 'lplane', 'lposnics', 'lreal_compat',
        'lrhfcalc', 'lrpa', 'lrpaforce', 'lscaaware', 'lscalapack', 'lscaler0',
        'lscalu', 'lscdm', 'lsck', 'lscrpa', 'lscsgrad', 'lselfenergy',
        'lsepb', 'lsepk', 'lsfbxc', 'lsingles', 'lsmp2lt', 'lsorbit',
        'lsoshift', 'lspectral', 'lspectralgw', 'lspin_vdw', 'lspiral',
        'lsubrot', 'lsynch5', 'ltangentold', 'ltau', 'ltboundlibxc', 'ltemper',
        'lthomas', 'ltriplet', 'ltssurf', 'ltwo_centre', 'luse_vdw',
        'lusenccl', 'lvacpotav', 'lvdw', 'lvdw_ewald', 'lvdwexpansion',
        'lvdwscs', 'lvgvappl', 'lvgvcalc', 'lvhar', 'lvtot', 'lwannier90',
        'lwannier90_auto_window', 'lwannier90_run', 'lwave', 'lwaveh5',
        'lweighted', 'lwrite_mmn_amn', 'lwrite_spn', 'lwrite_unk',
        'lwrite_wannier_xsf', 'lwrite_wanproj', 'lwrt_augmented_density',
        'lzeroz', 'lzora', 'ml_ff_lcouple_mb', 'ml_ff_lheat_mb', 'ml_ff_lmlff',
        'ml_ff_lnorm1_mb', 'ml_ff_lnorm2_mb', 'ml_ff_lsic_mb',
        'ml_ff_lsupervec_mb', 'ml_lafilt2', 'ml_lbasis_discard', 'ml_lcouple',
        'ml_leatom', 'ml_lemppot', 'ml_lerr', 'ml_lfast', 'ml_lheat', 'ml_lib',
        'ml_lmlff', 'ml_lsparsdes', 'ml_luse_names', 'nlspline', 'nucind',
        'oddonly', 'oddonlygw', 'pflat', 'phon_lbose', 'phon_lmc',
        'skip_edotp', 'velocity',
        })

    strings = frozenset({
        'algo', 'bandgap', 'bseprec', 'checkpoint_fd', 'cutoff_type',
        'dftd4_model', 'dftd4_xc', 'efermi', 'elph_decompose', 'elph_driver',
        'elph_mode', 'elph_scattering_approx', 'fftwmakeplan', 'gga',
        'libmbd_method', 'libmbd_parallel_mode', 'libmbd_vdw_params_kind',
        'libmbd_xc', 'libxc1', 'libxc2', 'locproj', 'lreal', 'metagga',
        'ml_grace_model', 'ml_mode', 'ml_type', 'nthreads_hi', 'nthreads_lo',
        'nthreads_mu', 'prec', 'precfock', 'pyamff_conv', 'pyamff_model',
        'pyamff_opt', 'random_generator', 'sdftd3_damping', 'sdftd3_xc',
        'stop_on', 'system', 'wannier90_win', 'wrt_density', 'wrt_potential',
        'xc',
        })

    int_arrays = frozenset({
        'elph_nbands_sum', 'iband', 'k_multiply', 'kpoint_bse', 'kpuse',
        'ldaul', 'libmbd_k_grid', 'ml_icouple', 'ml_mopot_ijm', 'ncrpa_bands',
        'nsubsys', 'ntarget_states', 'random_seed',
        })

    real_arrays = frozenset({
        'bext', 'bseelectron', 'bsehole', 'cmbj', 'cutoff_mu', 'cutoff_sigma',
        'dipol', 'efield_pead', 'efor', 'eint', 'elph_pot_fft_mesh',
        'elph_pot_lattice', 'elph_selfen_carrier_den',
        'elph_selfen_carrier_den_range', 'elph_selfen_carrier_per_cell',
        'elph_selfen_delta', 'elph_selfen_energy_window', 'elph_selfen_ikpt',
        'elph_selfen_kpts', 'elph_selfen_mu', 'elph_selfen_mu_range',
        'elph_selfen_temps', 'elph_selfen_temps_range', 'fbias_a', 'fbias_d',
        'fbias_r0', 'ferdo', 'ferwe', 'increm', 'langevin_gamma', 'ldauj',
        'ldauu', 'libmbd_alpha', 'libmbd_c6au', 'libmbd_r0au', 'm_constr',
        'magmom', 'ml_eatom_ref', 'ml_ff_eatom', 'ngyromag',
        'phon_born_charges', 'phon_dielectric', 'phon_tlist', 'pomass',
        'psubsys', 'qmaxfockae', 'qspiral', 'quad_efg', 'ropt', 'rwigs',
        'saxis', 'smearings', 'spring_k', 'spring_r0', 'spring_v0', 'tsubsys',
        'value_max', 'value_min', 'vca', 'vdw_alpha', 'vdw_c6', 'vdw_c6au',
        'vdw_r0', 'vdw_r0au', 'xc_c',
        })

    bool_arrays = frozenset({
        'fmp_active', 'lattice_constraints', 'lvdw_onecell',
        })

    mixed_types = MappingProxyType({
        'bext':   ('reals','real_arrays'),
        'efermi': ('strings','reals'),
        'libxc1': ('strings','ints'),
        'libxc2': ('strings','ints'),
        'lreal':  ('strings','bools'),
        })

    keyword_classification = obj(
        array_dimensions = (
            'nbands', 'ngzf', 'nplwv', 'ngxf', 'lmdim', 'ngy', 'ngz',
            'irdmax', 'irmax',
            'nions', 'nkpts', 'nkdim', 'ldim', 'ngyf', 'nedos', 'ngx'
            ),
        )

    # Retained for backwards compatibility.
    deprecated = frozenset({
        'elmin', 'enmax', 'enmin', 'hflmaxf', 'ichain', 'jacobian', 'lclimb',
        'ldneb', 'lmaxfockmp2', 'lmaxmp2', 'lnebcell', 'ltangentold', 'lvdw',
        'lvdwscs', 'mbja', 'mbjb', 'skip_edotp', 'timestep', 'vdw_scaling',
        })

#end class Incar



class Stopcar(VKeywordFile):
    keywords = frozenset({'lstop', 'labort'})
    bools    = frozenset({'lstop', 'labort'})
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
            if ' ' not in lines[1]:
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
                    self.kgrid     = np.array(lines[3].split(),dtype=int)
                    if len(lines)>4:
                        self.kshift    = np.array(lines[4].split(),dtype=float)
                    else:
                        self.kshift = None
                    #end if
                else:
                    self.mode   = 'basis'  # basis generated mesh
                    self.coord  = self.coord_options(cselect)
                    self.kbasis = np.array(self.join(lines,3,5).split(),dtype=float)
                    npe.reshape_inplace(self.kbasis, (3, 3))
                    self.kshift = np.array(lines[6].split(),dtype=float)
                #end if
            elif cselect=='l': # line mode (band structure)
                self.mode    = 'line'
                self.kinsert = iselect
                self.coord   = self.coord_options(lines[3].lower()[0])
                endpoints = []
                for line in lines[4:]:
                    endpoints.append(line.split())
                #end for
                self.kendpoints = np.array(endpoints,dtype=float)
            else: # explicit kpoints
                self.mode  = 'explicit'
                self.coord = self.coord_options(cselect)
                nkpoints = iselect
                kpw = []
                for line in lines[3:3+nkpoints]:
                    kpw.append(line.split())
                #end for
                kpw = np.array(kpw,dtype=float)
                self.kpoints  = kpw[:,0:3]
                self.kweights = kpw[:,3].ravel()
                tetline = 3+nkpoints
                if len(lines)>tetline and lines[tetline].lower()[0]=='t':
                    self.tetrahedra = obj()
                    tokens = lines[tetline+1].split()
                    ntets      = int(tokens[0])
                    tet_volume = float(tokens[1])
                    for n in range(ntets):
                        tokens = lines[tetline+2+n].split()
                        self.tetrahedra[len(self.tetrahedra)] = obj(
                            volume     = tet_volume,
                            degeneracy = int(tokens[0]),
                            corners    = np.array(tokens[1:],dtype=int),
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
                if self.kshift is not None:
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
            for n in range(npoints):
                text+=' {0:18.14f} {1:18.14f} {2:18.14f}   1\n'.format(*self.kendpoints[n])
                if n!=npoints-1 and n%2==1:
                    text+='\n'
                #end if
            #end for
        elif self.mode=='explicit':
            text+='explicit kpoints\n {0}\n'.format(len(self.kpoints))
            text+='{0}\n'.format(self.coord)
            for n in range(len(self.kpoints)):
                kp = self.kpoints[n]
                kw = self.kweights[n]
                text+=' {0:18.14f} {1:18.14f} {2:18.14f} {3:12.8f}\n'.format(kp[0],kp[1],kp[2],kw)
            #end for
            if 'tetrahedra' in self and len(self.tetrahedra)>0:
                ntets = len(self.tetrahedra)
                tets = self.tetrahedra
                text+='tetrahedra\n'
                text+=' {0} {1}\n'.format(ntets,tets[0].volume)
                for n in range(ntets):
                    t = tets[n]
                    d = t.degeneracy
                    c = t.corners
                    text+=' {0:d} {1:d} {2:d} {3:d} {4:d}\n'.format(d,*c)
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

    bool_map = MappingProxyType({True:'T',False:'F'})

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

        pos=self.pos
        spec=self.coord  #the current specifier
        
        if spec==specifier:
            return
        #end if
        if spec=="cartesian":
            pass
        elif spec=="direct":
            pos=np.dot(pos,axes)
        else:
            self.error("Poscar.change_specifier():  %s is not a valid coordinate specifier"%(spec))
        #end if
        spec=specifier  #the new specifier

        if spec=="cartesian":
            pass # already in cartesian coordinates.
        elif spec=="direct":
            pos=np.dot(pos,np.linalg.inv(axes))
        else:
            self.error("Poscar.change_specifier():  %s is not a valid coordinate specifier"%(spec))
        #end if

        self.coord=spec
        self.pos=pos
    #end def change_specifier


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
        axes = np.empty((dim,dim))
        axes[0] = np.array(lines[2].split(),dtype=float)
        axes[1] = np.array(lines[3].split(),dtype=float)
        axes[2] = np.array(lines[4].split(),dtype=float)
        tokens = lines[5].split()
        if tokens[0].isdigit():
            counts = np.array(tokens,dtype=int)
            elem   = None
            lcur   = 6
        else:
            elem   = np.array(tokens,dtype=str)
            counts = np.array(lines[6].split(),dtype=int)
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
        spos = np.array(spos)
        pos  = np.array(spos[:,0:3],dtype=float)
        if selective_dynamics:
            dynamic = np.array(spos[:,3:6],dtype=str)
            dynamic = np.char.upper(dynamic)=='T'
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
            vel = np.array(svel,dtype=float)
        else:
            vel_coord = None
            vel = None
        #end if
        self.update(
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
        if self.elem is not None:
            for e in self.elem:
                iselem, element = Elements.is_element(e, return_element=True)
                if not iselem:
                    self.error('{0} is not an element'.format(e))
                #end if
                text += element.symbol+' '
            #end for
            text += '\n'
        #end if
        for ec in self.elem_count:
            text += ' {0}'.format(ec)
        #end for
        text += '\n'
        if self.dynamic is not None:
            text += 'selective dynamics\n'
        #end if
        text += self.coord+'\n'
        if self.dynamic is None:
            for p in self.pos:
                text += ' {0:18.14f} {1:18.14f} {2:18.14f}\n'.format(*p)
            #end for
        else:
            bm = self.bool_map
            for i in range(len(self.pos)):
                p = self.pos[i]
                d = self.dynamic[i]
                text += ' {0:18.14f} {1:18.14f} {2:18.14f}  {3}  {4}  {5}\n'.format(p[0],p[1],p[2],bm[d[0]],bm[d[1]],bm[d[2]])
            #end for
        #end if
        if self.vel is not None:
            text += self.vel_coord+'\n'
            for v in self.vel:
                text += ' {0:18.14f} {1:18.14f} {2:18.14f}\n'.format(*v)
            #end for
        #end if
        return text
    #end def write_text


    def check_complete(self,*,exit=True):
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
        if self.vel is not None and self.vel_coord is None:
            msg += 'vel_coord is missing\n'
        #end if
        if exit and len(msg)>0:
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
        for n in range(len(self)):
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
        if filepath is not None and not os.path.isdir(filepath):
            VFile.__init__(self,filepath)
        else:
            VFile.__init__(self)
        #end if
    #end def __init__


    def read_text(self,text,filepath=''):
        marker = 'End of Dataset'
        pstart = 0
        end = len(text)
        while pstart<end:
            n = text.find(marker,pstart,end)
            if n==-1:
                break
            #end if
            line_end = text.find('\n',n+len(marker),end)
            if line_end==-1:
                pend = end
            else:
                pend = line_end+1
            #end if
            self.pseudos[len(self.pseudos)] = text[pstart:pend]
            pstart = pend
        #end while
    #end def read_text


    def write_text(self,filepath=''):
        text = ''
        if len(self.pseudos)>0:
            for i in range(len(self.pseudos)):
                text += self.pseudos[i]
            #end for
        elif self.filepath is not None and self.files is not None:
            for file in self.files:
                with open(os.path.join(self.filepath, file), "r") as f:
                    text += f.read()
            #end for
        #end if
        return text
    #end def write_text


    def pot_info(self):
        pot_info = obj()
        if len(self.pseudos)>0:
            pots = self.pseudos
        elif self.filepath is not None and self.files is not None:
            pots = obj()
            for file in self.files:
                with open(os.path.join(self.filepath,file),'r') as f:
                    pots[len(pots)] = f.read()
                #end with
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

            pot_info[len(pot_info)] = obj(label=label,Zval=Zval,element=element)
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
        if self.filepath is not None and self.files is not None:
            for file in self.files:
                with open(os.path.join(self.filepath,file),'r') as f:
                    self.pseudos[len(self.pseudos)] = f.read()
                #end with
            #end for
        #end if
    #end def load
#end class Potcar



class Exhcar(VFormattedFile):
    None
#end class Exhcar



class VaspInput(SimulationInput,Vobj):

    all_inputs  = (
        'PENALTYPOT', 'KPOINTS', 'WAVEDER', 'ICONST',
        'POSCAR', 'STOPCAR', 'EXHCAR', 'POTCAR', 'INCAR',
        )
    all_outputs = (
        'IBZKPT', 'LOCPOT', 'PRJCAR', 'PCDAT', 'ELFCAR', 'TMPCAR', 'DOSCAR', 'CHGCAR',
        'CONTCAR', 'HILLSPOT', 'PROCAR', 'EIGENVAL', 'WAVECAR', 'REPORT', 'CHG',
        'OUTCAR', 'PROOUT', 'XDATCAR', 'vasprun.xml', 'OSZICAR',
        )# note that CHGCAR, TMPCAR, and WAVECAR sometimes contain input

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
        if filepath is not None:
            self.read(filepath,prefix,postfix)
        #end if
    #end def __init__


    def read(self,filepath,prefix='',postfix=''):
        filepath = path_string(filepath)
        path = self.get_path(filepath)
        for file in os.listdir(path):
            name = str(file)
            if len(prefix)>0:
                if name.startswith(prefix):
                    name = name.split(prefix,1)[1]
                else:
                    continue
                #end if
            #end if
            if len(postfix)>0:
                if name.endswith(postfix):
                    name = name.rsplit(postfix,1)[0]
                else:
                    continue
                #end if
            #end if
            name = name.lower()
            if name in self.input_files:
                filepath = os.path.join(path,file)
                self[name] = self.input_files[name](filepath)
            #end if
        #end for
    #end def read


    def write(self,filepath,prefix='',postfix=''):
        path = self.get_path(filepath)
        for name,vfile in self.items():
            filepath = os.path.join(path,prefix+name.upper()+postfix)
            vfile.write(filepath)
        #end for
    #end def write


    def incorporate_system(self,system,coord='cartesian',*,incorp_kpoints=True,set_nelect=True):
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
            s = deepcopy(structure)
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
            if s.frozen is not None:
                poscar.dynamic = s.frozen==False
            #end if
            self.poscar = poscar
        #end if

        # handle charged systems
        if set_nelect or system.net_charge!=0:
            #  warning: spin polarization is handled by the user!
            if 'incar' not in self:
                self.incar = Incar()
            #end if
            self.incar.nelect = system.n_elec
        #end if

        return species
    #end def incorporate_system


    def return_system(self,*,structure_only=False,**valency):
        raw_axes = self.poscar.axes
        input_scale = self.poscar.scale
        if input_scale<0:
            scale_factor = (abs(input_scale)/abs(np.linalg.det(raw_axes)))**(1.0/3.0)
        else:
            scale_factor = input_scale
        #end if
        axes  = scale_factor*raw_axes
        scale = 1.0
 
        velem       = self.poscar.elem
        velem_count = self.poscar.elem_count
        elem        = vasp_to_nexus_elem(velem,velem_count)
 
        if self.poscar.coord=='direct':
            pos = np.dot(self.poscar.pos,axes)
        elif self.poscar.coord=='cartesian':
            pos = scale_factor*self.poscar.pos
        else:
            self.error('invalid POSCAR coordinate specifier: {0}'.format(self.poscar.coord))
        #end if
        
        center=axes.sum(0)/2.0
        
        kpoints  = None
        kweights = None
        kgrid    = None
        kshift   = None

        if structure_only:
            kpoints_file = None
        elif 'kpoints' in self:
            kpoints_file = self.kpoints
        else:
            self.error('KPOINTS is required to generate a physical system')
        #end if

        if kpoints_file is not None and kpoints_file.mode=="auto":
            kshift=kpoints_file.kshift
            if kshift is None:
                kshift = np.zeros(3,dtype=float)
            #end if
            if kpoints_file.centering=="monkhorst-pack":
                kshift=kshift+np.array([0.5,0.5,0.5])
            elif kpoints_file.centering=="gamma":
                pass
            #end if
            kgrid=kpoints_file.kgrid
        elif kpoints_file is not None:
            self.error('system generation does not currently work with manually specified k-points')
        #end if
        structure = Structure(
            axes     = axes,
            elem     = elem,
            scale    = scale,
            pos      = pos,
            center   = center,
            kpoints  = kpoints,
            kweights = kweights,
            kgrid    = kgrid,
            kshift   = kshift,
            units    = 'A',
            rescale  = False,
            )
         
        structure.zero_corner()
        structure.recenter()

        if structure_only:
            return structure
        #end if

        ion_charge = 0
        atoms      = list(elem)
        for atom in atoms:
            if atom not in valency:
                self.error('valence charge for atom {0} has not been defined\nplease provide the valence charge as an argument to return_system()'.format(atom))
            #end if
            ion_charge += atoms.count(atom)*valency[atom]
        #end for

        ####WARNING:  Assuming that the netcharge and netspin are ZERO. 
        net_charge = 0
        net_spin   = 0

        system = PhysicalSystem(
            structure  = structure,
            net_charge = net_charge,
            net_spin   = net_spin,
            **valency
            )
 
        return system
    #end def return_system


    def set_potcar(self,pseudos,species=None):
        if species is None:
            ordered_pseudos = pseudos
        else:
            pseudo_map = obj()
            for ppname in pseudos:
                element = ppname[0:2].strip('._')
                pseudo_map[element] = ppname
            #end for
            ordered_pseudos = []
            for element in species:
                iselem, elem = Elements.is_element(element, return_element=True)
                symbol = elem.symbol
                if not iselem:
                    self.error('{0} is not an element'.format(element))
                elif symbol not in pseudo_map:
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
        for n in range(len(neb_structures)):
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
        error('input_type {0} is unrecognized\nvalid options are: general'.format(input_type))
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
    coord      = 'cartesian',
    set_nelect = True,
    )

def generate_any_vasp_input(**kwargs):
    # handle 'system' name collision
    system_str = kwargs.pop('title',None)

    # remove keywords associated with kpoints, poscar, and any other formatted files
    vf = obj()
    for name,default in generate_any_defaults.items():
        if name in kwargs:
            vf[name] = kwargs[name]
            del kwargs[name]
        else:
            vf[name] = default
        #end if
    #end for
    gen_kpoints = 'kspacing' not in kwargs

    # create an empty input file
    vi = VaspInput()

    # assign values to incar and any other keyword files
    keywords = set(kwargs.keys())
    for name,keyword_file in VaspInput.keyword_files.items():
        keys = keywords & keyword_file.keywords
        if len(keys)>0:
            kw = obj()
            for k in keys:
                if k in kwargs:
                    kw[k] = kwargs[k]
                    del kwargs[k]
            vfile = keyword_file()
            vfile.assign(**kw)
            vi[name] = vfile
        #end if
    #end for

    # check for leftover keywords
    if len(kwargs)>0:
        error('unrecognized keywords: {0}'.format(sorted(kwargs.keys())),'generate_vasp_input')
    #end if

    # incorporate system information
    species = None
    if system_str is not None and 'incar' not in vi:
        vi.incar = Incar()
    #end if
    if vf.system is not None:
        species = vi.incorporate_system(
            system         = vf.system,
            incorp_kpoints = gen_kpoints,
            coord          = vf.coord,
            set_nelect     = vf.set_nelect,
            )
    #end if

    # set potcar
    if vf.pseudos is not None:
        vi.set_potcar(vf.pseudos,species)
    #end if

    # add kpoints information (override anything provided by system)
    if gen_kpoints and (vf.kpoints is not None or vf.kweights is not None or vf.kbasis is not None or vf.kgrid is not None or vf.kcenter is not None or vf.kendpoints is not None):
        if 'kpoints' in vi:
            kp = vi.kpoints
            kp.clear()
        else:
            kp = Kpoints()
            vi.kpoints = kp
        #end if
        if vf.kpoints is not None:
            kp.mode     = 'explicit'
            kp.kpoints  = vf.kpoints
            kp.kweights = vf.kweights
            kp.coord    = vf.kcoord
        elif vf.kgrid is not None:
            kp.mode      = 'auto'
            kp.centering = vf.kcenter
            if vf.kgrid is not None:
                kp.kgrid = vf.kgrid
            #end if
            if vf.kshift is not None:
                kp.kshift = vf.kshift
            #end if
        elif vf.kendpoints is not None:
            kp.mode       = 'line'
            kp.coord      = vf.kcoord
            kp.kinsert    = vf.kinsert
            kp.kendpoints = vf.kendpoints
        else:
            error('could not set kpoints from user inputs','generate_vasp_input')
        #end if
    #end if

    # create many poscars if doing nudged elastic band
    if vf.neb is not None:
        vi.setup_neb(*vf.neb,**vf.neb_args)
    #end if

    # handle 'system' name collision
    if system_str is not None:
        vi.incar.system = system_str
    #end if

    return vi
#end def generate_any_vasp_input




def generate_poscar(structure,*,coord='cartesian'):
    s = deepcopy(structure)
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
    if s.frozen is not None:
        poscar.dynamic = s.frozen==False
    #end if
    return poscar
#end def generate_poscar
