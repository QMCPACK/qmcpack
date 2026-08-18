##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################

"""Support I/O, generation, and manipulation of VASP input files.

``VaspInput`` collects the files for a VASP calculation, while
``generate_vasp_input`` and ``generate_poscar`` provide user-facing input
generation. Keyword files share parsing and writing through ``VKeywordFile``;
strictly formatted files derive from ``VFormattedFile``.
"""


import os
import re
from collections.abc import Mapping
from copy import deepcopy
from functools import partial
from types import MappingProxyType
import numpy as np
from .periodic_table import Elements
from .nexus_base import nexus_noncore
from .simulation import SimulationInput
from .structure import interpolate_structures, Structure
from .physical_system import PhysicalSystem
from .pseudopotential import PseudoSet
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
                msg = 'array repeat count must be non-negative'
                raise ValueError(msg)
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
        not set(v).isdisjoint(set('\n;#!'))
        or v.endswith('\\')
        or v != v.strip()
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
        msg = 'value must be an integer, not a boolean'
        raise TypeError(msg)
    elif isinstance(v,(int,np.integer)):
        return int(v)
    elif isinstance(v,(float,np.floating)):
        value = float(v)
        if not np.isfinite(value):
            msg = 'value must be finite'
            raise ValueError(msg)
        elif abs(value)>max_exact_float_integer:
            msg = 'floating-point integer values must not exceed {}'.format(
                max_exact_float_integer
                )
            raise ValueError(msg)
        #end if
        nearest = round(value)
        if not np.isclose(value,nearest,rtol=0.0,atol=integer_tolerance):
            msg = 'real value {} is not within {} of an integer'.format(
                v,integer_tolerance
                )
            raise ValueError(msg)
        #end if
        return nearest
    else:
        msg = 'value must be an integer or real number'
        raise TypeError(msg)
    #end if
#end def assign_int


def assign_string(v):
    if isinstance(v,str):
        return v
    else:
        msg = 'value must be a string'
        raise ValueError(msg)
    #end if
#end def assign_string


def assign_int_array(a):
    if isinstance(a,(tuple,list,np.ndarray)):
        array = np.asarray(a,dtype=object)
        if array.ndim!=1:
            msg = 'value must be a one-dimensional array'
            raise ValueError(msg)
        #end if
        values = []
        for index,value in enumerate(array):
            try:
                values.append(assign_int(value))
            except (TypeError,ValueError) as e:
                msg = 'invalid integer array element at index {}: {}'.format(
                    index,e
                    )
                raise type(e)(msg) from None
            #end try
        #end for
        return np.array(values,dtype=int)
    else:
        msg = 'value must be a tuple, list, or array'
        raise ValueError(msg)
    #end if
#end def assign_int_array


def assign_real_array(a):
    if isinstance(a,(tuple,list,np.ndarray)):
        return np.array(a,dtype=float)
    else:
        msg = 'value must be a tuple, list, or array'
        raise ValueError(msg)
    #end if
#end def assign_real_array


def assign_bool_array(a):
    if isinstance(a,(tuple,list,np.ndarray)):
        return np.array(a,dtype=bool)
    else:
        msg = 'value must be a tuple, list, or array'
        raise ValueError(msg)
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


block_type_names = MappingProxyType({
    'int':        'ints',
    'real':       'reals',
    'bool':       'bools',
    'string':     'strings',
    'int_array':  'int_arrays',
    'real_array': 'real_arrays',
    'bool_array': 'bool_arrays',
    })


def mixed_type_matches(value,value_type):
    if value_type=='strings':
        return isinstance(value,str)
    elif value_type=='bools':
        return isinstance(value,(bool,np.bool_,int,np.integer))
    elif value_type in ('ints','reals'):
        return (
            isinstance(value,(int,float,np.integer,np.floating))
            and not isinstance(value,(bool,np.bool_))
            )
    elif value_type in ('int_arrays','real_arrays','bool_arrays'):
        return isinstance(value,(tuple,list,np.ndarray))
    else:
        msg = (
            'unknown value for keyword value_type: {}, must be one of {}'
            .format(value_type,list(block_type_names.values()))
            )
        raise ValueError(msg)
    #end if
#end def mixed_type_matches


def assign_mixed(value,*,types):
    errors = []
    for value_type in types:
        if mixed_type_matches(value,value_type):
            try:
                return assign_value_functions[value_type](value)
            except Exception as e:  # noqa: BLE001
                errors.append(
                    '{0} for value {1}: {2}'.format(value_type,value,e)
                    )
            #end try
        #end if
    #end for
    message = 'value does not match permitted types: {0}'.format(types)
    if len(errors)>0:
        message += '\n'+'\n'.join(errors)
    #end if
    raise ValueError(message)
#end def assign_mixed


def read_mixed(sval,*,types):
    scalar_types = {'ints','reals','bools'}
    array_types = {'int_arrays','real_arrays','bool_arrays'}
    try:
        nvalues = len(expand_array(sval))
    except Exception:  # noqa: BLE001
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
                except Exception as e:  # noqa: BLE001
                    errors.append('{0}: {1}'.format(value_type,e))
                #end try
            #end if
        #end for
    #end for
    if nvalues==1 and 'strings' in types:
        return read_string(sval)
    #end if
    msg = 'value does not match permitted types: {}\n{}'.format(
        types,'\n'.join(errors)
        )
    raise ValueError(msg)
#end def read_mixed


def write_mixed(value,*,types):
    converted = assign_mixed(value,types=types)
    for value_type in types:
        if mixed_type_matches(converted,value_type):
            return write_value_functions[value_type](converted)
        #end if
    #end for
    msg = 'value does not match permitted types: {}'.format(types)
    raise ValueError(msg)
#end def write_mixed





class Vobj(DevBase):
    """Base class for VASP input objects."""

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
    """Base class for VASP input files."""

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
                self.error(
                    'max number of multi-line strings exceeded.\n'
                    'Over {} quotation marks found in file.'
                    .format(nqmax-1)
                    )
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
    """Base class for keyword-based VASP files."""

    kw_scalars = ('ints','reals','bools','strings')
    kw_arrays  = ('int_arrays','real_arrays','bool_arrays')
    kw_fields  = kw_scalars + kw_arrays + ('keywords','unsupported')

    keyword_classification = None
    mixed_types = MappingProxyType({})
    block_constructs = MappingProxyType({})

    @classmethod
    def class_init(cls):
        cls.kw_scalars = VKeywordFile.kw_scalars
        cls.kw_arrays  = VKeywordFile.kw_arrays
        cls.kw_fields  = VKeywordFile.kw_fields
        for kw_field in cls.kw_fields:
            if not hasattr(cls,kw_field):
                setattr(cls,kw_field,set())
            #end if
            setattr(cls,kw_field,frozenset(getattr(cls,kw_field)))
        #end for
        cls.scalar_keywords = set()
        for scalar_field in cls.kw_scalars:
            cls.scalar_keywords |= getattr(cls,scalar_field)
        #end for
        cls.scalar_keywords = frozenset(cls.scalar_keywords)
        cls.array_keywords = set()
        for array_field in cls.kw_arrays:
            cls.array_keywords |= getattr(cls,array_field)
        #end for
        cls.array_keywords = frozenset(cls.array_keywords)
        cls.keywords = frozenset(
            cls.scalar_keywords | cls.array_keywords
            )
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
                msg = 'mixed-type keyword is not classified: {}'.format(name)
                raise ValueError(msg)
            #end if
            classified_type = cls.type[name]
            if classified_type not in types:
                msg = (
                    'classified type {} is not permitted for keyword {}'
                    .format(classified_type,name)
                    )
                raise ValueError(msg)
            #end if
            for value_type in types:
                if value_type not in read_value_functions:
                    msg = 'invalid type {} for keyword {}'.format(
                        value_type,name
                        )
                    raise ValueError(msg)
                #end if
            #end for
            cls.type[name] = types
            cls.read_value[name] = partial(read_mixed,types=types)
            cls.write_value[name] = partial(write_mixed,types=types)
            cls.assign_value[name] = partial(assign_mixed,types=types)
        #end for
        for block_name,schema in cls.block_constructs.items():
            if block_name=='image_*' and schema=='keywords':
                continue
            elif not isinstance(schema,MappingProxyType):
                msg = (
                    'schema for block construct {} must be read-only'
                    .format(block_name)
                    )
                raise TypeError(msg)
            #end if
            for field,value_type in schema.items():
                if value_type not in block_type_names:
                    msg = 'invalid type {} for block field {}/{}'.format(
                        value_type,block_name,field
                        )
                    raise ValueError(msg)
                #end if
            #end for
        #end for
    #end def class_init


    @classmethod
    def block_schema(cls,name):
        if name in cls.block_constructs:
            return cls.block_constructs[name]
        elif re.fullmatch(r'image_[1-9][0-9]*',name):
            return cls.block_constructs.get('image_*',None)
        else:
            return None
        #end if
    #end def block_schema


    @classmethod
    def is_block_name(cls,name):
        return cls.block_schema(name) is not None
    #end def is_block_name


    def extract_block_constructs(self,text):
        block_names = [
            re.escape(name) for name in self.block_constructs
            if name!='image_*'
            ]
        if 'image_*' in self.block_constructs:
            block_names.append(r'image_[1-9][0-9]*')
        #end if
        if len(block_names)==0:
            return text,obj()
        #end if
        pattern = re.compile(
            r'\b('+'|'.join(block_names)+r')\s*(?:=\s*)?\{',
            re.IGNORECASE,
            )
        blocks = obj()
        output = ''
        start = 0
        while True:
            match = pattern.search(text,start)
            if match is None:
                output += text[start:]
                break
            #end if
            output += text[start:match.start()]
            depth = 1
            end = match.end()
            while end<len(text) and depth>0:
                if text[end]=='{':
                    depth += 1
                elif text[end]=='}':
                    depth -= 1
                #end if
                end += 1
            #end while
            name = match.group(1).lower()
            if depth!=0:
                self.error('block construct {0} is not closed'.format(name))
            #end if
            label = '--block{0}--'.format(str(len(blocks)).zfill(3))
            blocks[label] = text[match.end():end-1]
            output += name+' = '+label
            start = end
        #end while
        return output,blocks
    #end def extract_block_constructs


    def block_field_type(self,block_name,field,schema):
        if schema=='keywords':
            if field not in self.keywords:
                msg = '{} is not an INCAR keyword'.format(field.upper())
                raise ValueError(msg)
            #end if
            return self.type[field]
        elif field not in schema:
            msg = '{} is not a field of block construct {}'.format(
                field.upper(),block_name.upper()
                )
            raise ValueError(msg)
        else:
            return block_type_names[schema[field]]
        #end if
    #end def block_field_type


    def read_block_construct(self,name,text,multiline_values):
        schema = self.block_schema(name)
        value = obj()
        expression = ''
        for line in text.splitlines():
            line = self.remove_comment(line).strip()
            if len(line)==0:
                continue
            #end if
            expression += line
            if expression.endswith('\\'):
                expression = expression[:-1]+' '
                continue
            #end if
            for token in expression.split(';'):
                if len(token.strip())==0:
                    continue
                elif '=' not in token:
                    msg = 'missing assignment in block {}: {}'.format(
                        name.upper(),token.strip()
                        )
                    raise ValueError(msg)
                #end if
                field,sval = token.split('=',1)
                field = field.lower().strip()
                sval = sval.strip()
                if sval.startswith('--multiline'):
                    sval = multiline_values[sval]
                #end if
                value_type = self.block_field_type(name,field,schema)
                if isinstance(value_type,tuple):
                    value[field] = read_mixed(sval,types=value_type)
                else:
                    value[field] = read_value_functions[value_type](sval)
                #end if
            #end for
            expression = ''
        #end for
        if len(expression)>0:
            msg = 'incomplete line continuation in block {}'.format(name)
            raise ValueError(msg)
        #end if
        return value
    #end def read_block_construct


    def assign_block_construct(self,name,value):
        if not isinstance(value,Mapping):
            msg = (
                'block construct value should be a mapping, but is {}'
                .format(type(value).__name__)
                )
            raise TypeError(msg)
        #end if
        schema = self.block_schema(name)
        assigned = obj()
        for field,field_value in value.items():
            field = field.lower()
            value_type = self.block_field_type(name,field,schema)
            if isinstance(value_type,tuple):
                assigned[field] = assign_mixed(
                    field_value,types=value_type
                    )
            else:
                assigned[field] = assign_value_functions[value_type](
                    field_value
                    )
            #end if
        #end for
        return assigned
    #end def assign_block_construct


    def write_block_construct(self,name,value):
        value = self.assign_block_construct(name,value)
        schema = self.block_schema(name)
        fields = value.keys()
        maxlen = max(map(len,fields),default=0)
        text = name.upper()+' = {\n'
        for field in sorted(fields):
            value_type = self.block_field_type(name,field,schema)
            if isinstance(value_type,tuple):
                sval = write_mixed(value[field],types=value_type)
            else:
                sval = write_value_functions[value_type](value[field])
            #end if
            text += '  {0:<{fmt}} = {1}\n'.format(
                field.upper(),sval,fmt=maxlen
                )
        #end for
        return text+'}\n'
    #end def write_block_construct


    def read_text(self,text,filepath=''):
        text,multiline_values = self.preprocess_multiline_strings(text)
        text = '\n'.join(
            self.remove_comment(line) for line in text.splitlines()
            )
        text,blocks = self.extract_block_constructs(text)
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
                            if value.startswith('--multiline'):
                                value = multiline_values[value]
                            #end if
                            if self.is_block_name(name):
                                block_value = self.read_block_construct(
                                    name,blocks[value],multiline_values
                                    )
                                if name not in self:
                                    self[name] = obj()
                                #end if
                                for field,field_value in block_value.items():
                                    self[name][field] = field_value
                                #end for
                                continue
                            elif '/' in name:
                                block_name,field = name.split('/',1)
                                schema = self.block_schema(block_name)
                                if schema is not None:
                                    value_type = self.block_field_type(
                                        block_name,field,schema
                                        )
                                    if isinstance(value_type,tuple):
                                        value = read_mixed(
                                            value,types=value_type
                                            )
                                    else:
                                        value = read_value_functions[
                                            value_type
                                            ](value)
                                    #end if
                                    if block_name not in self:
                                        self[block_name] = obj()
                                    #end if
                                    self[block_name][field] = value
                                    continue
                                #end if
                            #end if
                            if name in self.keywords:
                                try:
                                    value = self.read_value[name](value)
                                    self[name] = value
                                except Exception as e:  # noqa: BLE001
                                    self.error(
                                        'read failed for keyword {}\n'
                                        'keyword type: {}\n'
                                        'input text: {}\n'
                                        'exception:\n'
                                        '{}'.format(
                                            name,self.type[name],token,e
                                            )
                                        )
                                #end try
                            elif name in self.unsupported:
                                self.warn(
                                    'keyword {0} is not currently supported'
                                    .format(name)
                                    )
                            else:
                                #ci(lcs(),gs())
                                self.error(
                                    '{0} is not a keyword for the {1} file'
                                    .format(
                                        name.upper(),
                                        self.__class__.__name__.upper(),
                                        )
                                    )
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
            if self.is_block_name(name):
                text += self.write_block_construct(name,value)
                continue
            #end if
            try:
                svalue = self.write_value[name](value)
            except Exception as e:  # noqa: BLE001
                self.error(
                    'write failed for file {} keyword {}\n'
                    'keyword type: {}\n'
                    'value: {}\n'
                    'exception:\n'
                    '{}'.format(filepath,name,self.type[name],value,e)
                    )
            #end try
            text += valfmt.format(name.upper(),svalue)
        #end for
        return text
    #end def write_text


    def assign(self,**values):
        for name,value in values.items():
            if self.is_block_name(name):
                self[name] = self.assign_block_construct(name,value)
                continue
            elif name not in self.keywords:
                self.error(
                    '{0} is not a keyword for the {1} file'
                    .format(name.upper(),self.__class__.__name__.upper())
                    )
            #end if
            try:
                self[name] = self.assign_value[name](value)
            except Exception as e:  # noqa: BLE001
                self.error(
                    'assign failed for keyword {}\n'
                    'keyword type: {}\n'
                    'value: {}\n'
                    'exception:\n'
                    '{}'.format(name,self.type[name],value,e)
                    )
            #end try
        #end for
    #end def assign
#end class VKeywordFile



class VFormattedFile(VFile):
    """Base class for structured VASP files."""

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
    """Represent a VASP INCAR file."""

    # VASP wiki with incar keys/tags
    #   https://www.vasp.at/wiki/index.php/Category:INCAR_tag

    # VTST extensions:  https://henkelmangroup.github.io/vtsttools/

    unsupported = frozenset()

    # Pages still present in Nexus history but absent or malformed on the wiki.
    broken_docs = frozenset({
        'dmft_basis',
        })

    ints = frozenset({
        'antires', 'ch_nedos', 'cll', 'cln', 'clnt', 'drotmax', 'efermi_nedos',
        'elmin', 'elph_fermi_nedos', 'elph_ismear', 'elph_lr', 'elph_nbands',
        'elph_selfen_band_start_kp', 'elph_selfen_band_stop_kp',
        'elph_selfen_nw', 'elph_transport_driver', 'elph_transport_nedos',
        'elph_transport_nedos_plot', 'elph_wf_comm_opt', 'exxoep', 'findiff',
        'fmp_direction', 'fmp_period', 'fmp_snumber', 'fmp_swapnum', 'fnmin',
        'fockcorr', 'hflmax', 'hflmaxf', 'hills_bin', 'i_constrained_m',
        'ialgo', 'iall_in_one', 'ibrion',
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
        'cmbja', 'cmbjb', 'cmbje', 'cparam', 'cshift', 'csvr_period', 'ddr',
        'deg_threshold', 'deper', 'dfnmax', 'dfnmin', 'dimer_dist', 'dq',
        'ebreak', 'ediff', 'ediffg', 'efield', 'elph_kspacing',
        'elph_selfen_band_start',
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
        'elph_prepare', 'elph_rotateprojectors', 'elph_run', 'elph_selfen_dw',
        'elph_selfen_fan', 'elph_selfen_g_skip', 'elph_selfen_gaps',
        'elph_selfen_imag_skip', 'elph_selfen_static', 'elph_transport',
        'elph_useblas', 'elph_userecip', 'elph_wf_cache_prefill',
        'elph_wf_redistribute', 'elph_write_hdf5vel', 'elph_write_textvel',
        'evenonly', 'evenonlygw',
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

    block_constructs = MappingProxyType({
        'kernel_truncation': MappingProxyType({
            'factor':          'real',
            'idimensionality': 'int',
            'ipad':            'int',
            'isurface':        'int',
            'lcoarsen':        'bool',
            'ltruncate':       'bool',
            }),
        'plugins': MappingProxyType({
            'force_and_stress': 'bool',
            'local_potential':  'bool',
            'machine_learning': 'bool',
            'ml_mode':          'string',
            'ml_outblock':      'int',
            'ml_output_mode':   'int',
            'neighbor_cutoff':  'real',
            'occupancies':      'bool',
            'structure':        'bool',
            }),
        'image_*': 'keywords',
        })

    # Retained for backwards compatibility.
    deprecated = frozenset({
        'elmin', 'enmax', 'enmin', 'hflmaxf', 'ichain', 'jacobian', 'lclimb',
        'ldneb', 'lmaxfockmp2', 'lmaxmp2', 'lnebcell', 'ltangentold', 'lvdw',
        'lvdwscs', 'mbja', 'mbjb', 'skip_edotp', 'timestep', 'vdw_scaling',
        })

#end class Incar



class Stopcar(VKeywordFile):
    """Represent a VASP STOPCAR file."""

    keywords = frozenset({'lstop', 'labort'})
    bools    = frozenset({'lstop', 'labort'})
#end class Stopcar


for cls in Incar,Stopcar:
    cls.class_init()
#end for



class Iconst(VFormattedFile):  # metadynamics -> 6.62.4
    """Represent geometric constraints in an ICONST file."""

    def __init__(self,filepath=None):
        self.coordinates = obj()
        VFile.__init__(self,filepath)
    #end def __init__


    def read_text(self,text,filepath=''):
        self.coordinates.clear()
        lines = self.read_lines(text,remove_empty=True)
        for line in lines:
            tokens = line.split()
            if len(tokens)<2:
                self.error('invalid ICONST line: {0}'.format(line))
            #end if
            values = []
            for token in tokens[1:-1]:
                try:
                    value = int(token)
                except ValueError:
                    value = read_real(token)
                #end try
                values.append(value)
            #end for
            self.coordinates[len(self.coordinates)] = obj(
                flag   = tokens[0],
                items  = tuple(values),
                status = int(tokens[-1]),
                )
        #end for
    #end def read_text


    def write_text(self,filepath=''):
        text = ''
        for index in sorted(self.coordinates.keys()):
            coordinate = self.coordinates[index]
            values = [coordinate.flag]
            values.extend(map(str,coordinate['items']))
            values.append(str(assign_int(coordinate.status)))
            text += ' '.join(values)+'\n'
        #end for
        return text
    #end def write_text
#end class Iconst



class Kpoints(VFormattedFile):
    """Represent Brillouin-zone sampling in a KPOINTS file."""

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
        lines_in = self.read_lines(text,remove_empty=False)
        if len(lines_in)>0:
            lines = [lines_in[0]] + [
                line for line in lines_in[1:] if len(line)>0
                ]
        else:
            lines = []
        #end if
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
                    self.kgrid     = read_real(lines[3])
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
                labels = []
                for line in lines[4:]:
                    tokens = line.split()
                    endpoints.append(tokens[:3])
                    labels.append(' '.join(tokens[3:]))
                #end for
                self.kendpoints = np.array(endpoints,dtype=float)
                if any(len(label)>0 for label in labels):
                    self.labels = np.array(labels,dtype=str)
                #end if
            else: # explicit kpoints
                self.mode  = 'explicit'
                self.coord = self.coord_options(cselect)
                nkpoints = iselect
                kpw = []
                labels = []
                for line in lines[3:3+nkpoints]:
                    tokens = line.split()
                    kpw.append(tokens[:4])
                    labels.append(' '.join(tokens[4:]))
                #end for
                kpw = np.array(kpw,dtype=float)
                self.kpoints  = kpw[:,0:3]
                self.kweights = kpw[:,3].ravel()
                if any(len(label)>0 for label in labels):
                    self.labels = np.array(labels,dtype=object)
                #end if
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
                text+=' {0}\n'.format(self.kgrid)
            elif self.centering=='gamma' or self.centering=='monkhorst-pack':
                text+='{0}\n'.format(self.centering)
                text+=' {0} {1} {2}\n'.format(*self.kgrid)
                if self.kshift is not None:
                    text+=' {0} {1} {2}\n'.format(*self.kshift)
                #end if
            else:
                self.error(
                    'invalid centering for file {0}: {1}\n'
                    'valid options are: auto, gamma, monkhorst-pack'
                    .format(filepath,self.centering)
                    )
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
                text += ' {0:18.14f} {1:18.14f} {2:18.14f}'.format(
                    *self.kendpoints[n]
                    )
                if 'labels' in self and len(self.labels[n])>0:
                    text += '  {0}'.format(self.labels[n])
                #end if
                text += '\n'
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
                text += (
                    ' {0:18.14f} {1:18.14f} {2:18.14f} {3:12.8f}'
                    .format(kp[0],kp[1],kp[2],kw)
                    )
                if 'labels' in self and len(self.labels[n])>0:
                    text += '  {0}'.format(self.labels[n])
                #end if
                text += '\n'
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
                    text+=' {0} {1} {2} {3} {4}\n'.format(d,*c)
                #end for
            #end if
        else:
            self.error(
                'invalid mode: {0}\n'
                'valid options are: auto, basis, line, explicit'
                )
        #end if
        return text
    #end def write_text
#end class Kpoints


class KpointsOpt(Kpoints):
    """Represent a VASP KPOINTS_OPT file."""
#end class KpointsOpt


class KpointsElph(Kpoints):
    """Represent a VASP KPOINTS_ELPH file."""
#end class KpointsElph


class KpointsWan(Kpoints):
    """Represent a VASP KPOINTS_WAN file."""
#end class KpointsWan


class Qpoints(Kpoints):
    """Represent phonon sampling in a QPOINTS file."""
#end class Qpoints



class Penaltypot(VFormattedFile):  # metadynamics -> 6.62.4 (2nd one)
    """Represent bias potentials in a PENALTYPOT file."""

    def __init__(self,filepath=None):
        self.hills = np.empty((0,0),dtype=float)
        VFile.__init__(self,filepath)
    #end def __init__


    def read_text(self,text,filepath=''):
        lines = self.read_lines(text,remove_empty=True)
        rows = [line.split() for line in lines]
        if len(rows)==0:
            self.hills = np.empty((0,0),dtype=float)
        elif len({len(row) for row in rows})!=1:
            self.error('PENALTYPOT rows must have equal lengths')
        else:
            self.hills = np.array(rows,dtype=float)
        #end if
    #end def read_text


    def write_text(self,filepath=''):
        text = ''
        hills = np.asarray(self.hills,dtype=float)
        if hills.ndim!=2:
            self.error('PENALTYPOT hills must be a two-dimensional array')
        #end if
        for hill in hills:
            text += ' '.join(map(str,hill))+'\n'
        #end for
        return text
    #end def write_text
#end class Penaltypot


class Irccar(VFormattedFile):
    """Represent a discretized path in an IRCCAR file."""

    def __init__(self,filepath=None):
        self.points = np.empty((0,0),dtype=float)
        VFile.__init__(self,filepath)
    #end def __init__


    def read_text(self,text,filepath=''):
        lines = self.read_lines(text,remove_empty=True)
        if len(lines)==0:
            self.error('IRCCAR is empty')
        #end if
        npoints = int(lines[0])
        if len(lines[1:])!=npoints:
            self.error(
                'IRCCAR declares {0} points but contains {1}'
                .format(npoints,len(lines[1:]))
                )
        #end if
        rows = [line.split() for line in lines[1:]]
        if len(rows)>0 and len({len(row) for row in rows})!=1:
            self.error('IRCCAR point rows must have equal lengths')
        #end if
        self.points = np.array(rows,dtype=float)
    #end def read_text


    def write_text(self,filepath=''):
        points = np.asarray(self.points,dtype=float)
        if points.ndim!=2:
            self.error('IRCCAR points must be a two-dimensional array')
        #end if
        text = str(len(points))+'\n'
        for point in points:
            text += ' '.join(map(str,point))+'\n'
        #end for
        return text
    #end def write_text
#end class Irccar


# Retained for backwards compatibility.
Ircar = Irccar


class VRawFile(VFormattedFile):
    """Preserve a VASP text file without interpreting it."""

    def __init__(self,filepath=None):
        self.text = ''
        VFile.__init__(self,filepath)
    #end def __init__


    def read_text(self,text,filepath=''):
        self.text = text
    #end def read_text


    def write_text(self,filepath=''):
        return self.text
    #end def write_text
#end class VRawFile


class Hessemat(VRawFile):
    """Represent a VASP HESSEMAT file."""
#end class Hessemat


class Gamma(VRawFile):
    """Represent a VASP GAMMA file."""
#end class Gamma


class Wanproj(VRawFile):
    """Represent a VASP WANPROJ file."""
#end class Wanproj


class MlAb(VRawFile):
    """Represent a VASP ML_AB file."""
#end class MlAb


class MlFf(VRawFile):
    """Represent a VASP ML_FF file."""
#end class MlFf


class Dynmatfull(VRawFile):
    """Represent a VASP DYNMATFULL file."""
#end class Dynmatfull


class Chgcar(VRawFile):
    """Represent a VASP CHGCAR file."""
#end class Chgcar


class Taucar(VRawFile):
    """Represent a VASP TAUCAR file."""
#end class Taucar



class Poscar(VFormattedFile):
    """Represent a VASP POSCAR file."""

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
            self.error(
                'Poscar.change_specifier():  {0} is not a valid coordinate '
                'specifier'.format(spec)
                )
        #end if
        spec=specifier  #the new specifier

        if spec=="cartesian":
            pass # already in cartesian coordinates.
        elif spec=="direct":
            pos=np.dot(pos,np.linalg.inv(axes))
        else:
            self.error(
                'Poscar.change_specifier():  {0} is not a valid coordinate '
                'specifier'.format(spec)
                )
        #end if

        self.coord=spec
        self.pos=pos
    #end def change_specifier


    def read_text(self,text,filepath=''):
        lines = self.read_lines(text,remove_empty=False)
        nlines = len(lines)
        min_lines = 8
        if nlines<min_lines:
            self.error(
                'file {0} must have at least {1} lines\n'
                '  only {2} lines found'
                .format(filepath,min_lines,nlines)
                )
        #end if
        description = text.split('\n',1)[0].strip()
        dim = 3
        scale_values = np.array(lines[1].split(),dtype=float)
        if len(scale_values)==1:
            scale = float(scale_values[0])
        elif len(scale_values)==3:
            scale = scale_values
        else:
            self.error(
                'file {0} must contain one or three scaling factors'
                .format(filepath)
                )
        #end if
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
        pos = np.array([tokens[:3] for tokens in spos],dtype=float)
        if selective_dynamics:
            dynamic = np.array([tokens[3:6] for tokens in spos],dtype=str)
            dynamic = np.char.upper(dynamic)=='T'
            label_start = 6
        else:
            dynamic = None
            label_start = 3
        #end if
        labels = []
        for tokens in spos:
            labels.append(' '.join(tokens[label_start:]))
        #end for
        if not any(len(label)>0 for label in labels):
            labels = None
        else:
            labels = np.array(labels,dtype=str)
        #end if

        lattice_vel_init = None
        lattice_vel = None
        lattice_vectors = None
        if lcur<len(lines) and lines[lcur].lower().startswith('l'):
            if lcur+8>len(lines):
                self.error(
                    'file {0} is incomplete (missing lattice velocities)'
                    .format(filepath)
                    )
            #end if
            lcur += 1
            lattice_vel_init = int(lines[lcur].split()[0])
            lcur += 1
            lattice_vel = np.array(
                [lines[lcur+i].split()[:3] for i in range(3)],
                dtype=float,
                )
            lcur += 3
            lattice_vectors = np.array(
                [lines[lcur+i].split()[:3] for i in range(3)],
                dtype=float,
                )
            lcur += 3
        #end if
        vector_header = None
        vectors = None
        if lcur<len(lines) and not self.is_empty(lines,lcur):
            header = lines[lcur]
            cline = header.lower()
            lcur+=1
            if lcur+npos>len(lines):
                self.error(
                    'file {0} is incomplete (missing post-position vectors)'
                    .format(filepath)
                    )
            #end if
            is_velocity = (
                len(cline)==0 or cline[0] in ('c','k','d')
                )
            if is_velocity:
                cartesian = len(cline)==0 or cline[0] in ('c','k')
                vel_coord = 'cartesian' if cartesian else 'direct'
                vector_header = header
            else:
                vel_coord = None
                vector_header = header
            #end if
            svectors = []
            for i in range(npos):
                svectors.append(lines[lcur+i].split()[:3])
            #end for
            lcur += npos
            vectors = np.array(svectors,dtype=float)
            vel = vectors if is_velocity else None
        else:
            vel_coord = None
            vel = None
        #end if
        if lcur<len(lines) and not self.is_empty(lines,lcur):
            md_extra = '\n'.join(lines[lcur:])+'\n'
        else:
            md_extra = None
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
            vel         = vel,
            )
        if vectors is not None:
            self.vector_header = vector_header
            self.vectors = vectors
        #end if
        if labels is not None:
            self.labels = labels
        #end if
        if lattice_vel is not None:
            self.lattice_vel_init = lattice_vel_init
            self.lattice_vel = lattice_vel
            self.lattice_vectors = lattice_vectors
        #end if
        if md_extra is not None:
            self.md_extra = md_extra
        #end if
    #end def read_text


    def write_text(self,filepath=''):
        msg = self.check_complete(exit=False)
        if msg!='':
            self.error(
                'incomplete data to write file {0}\n'
                '{1}'.format(filepath,msg)
                )
        #end if
        text = ''
        if self.description is None:
            text += 'System cell and coordinates\n'
        else:
            text += self.description+'\n'
        #end if
        if np.isscalar(self.scale):
            text += ' {0}\n'.format(self.scale)
        else:
            text += ' {0} {1} {2}\n'.format(*self.scale)
        #end if
        for a in self.axes:
            text += ' {0:18.14f} {1:18.14f} {2:18.14f}\n'.format(*a)
        #end for
        if self.elem is not None:
            for e in self.elem:
                text += str(e)+' '
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
            for i,p in enumerate(self.pos):
                text += ' {0:18.14f} {1:18.14f} {2:18.14f}'.format(*p)
                if 'labels' in self and len(self.labels[i])>0:
                    text += '  {0}'.format(self.labels[i])
                #end if
                text += '\n'
            #end for
        else:
            bm = self.bool_map
            for i in range(len(self.pos)):
                p = self.pos[i]
                d = self.dynamic[i]
                text += (
                    ' {0:18.14f} {1:18.14f} {2:18.14f}'
                    '  {3}  {4}  {5}'
                    .format(p[0],p[1],p[2],bm[d[0]],bm[d[1]],bm[d[2]])
                    )
                if 'labels' in self and len(self.labels[i])>0:
                    text += '  {0}'.format(self.labels[i])
                #end if
                text += '\n'
            #end for
        #end if
        if 'lattice_vel' in self:
            text += 'Lattice velocities and vectors\n'
            text += ' {0}\n'.format(self.lattice_vel_init)
            for vector in self.lattice_vel:
                text += ' {0:18.14f} {1:18.14f} {2:18.14f}\n'.format(
                    *vector
                    )
            #end for
            for vector in self.lattice_vectors:
                text += ' {0:18.14f} {1:18.14f} {2:18.14f}\n'.format(
                    *vector
                    )
            #end for
        #end if
        if self.vel is not None:
            vectors = self.vel
            vector_header = (
                '' if self.vel_coord=='cartesian' else self.vel_coord
                )
        else:
            vectors = self.vectors if 'vectors' in self else None
            vector_header = (
                self.vector_header if 'vector_header' in self else None
                )
        #end if
        if vectors is not None:
            if vector_header is None:
                self.error(
                    'post-position vector header is missing for file {0}'
                    .format(filepath)
                    )
            #end if
            text += vector_header+'\n'
            for v in vectors:
                text += ' {0:18.14f} {1:18.14f} {2:18.14f}\n'.format(*v)
            #end for
        #end if
        if 'md_extra' in self:
            text += self.md_extra
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
        if 'vectors' in self and self.vectors is not None:
            if 'vector_header' not in self or self.vector_header is None:
                msg += 'vector_header is missing\n'
            elif np.shape(self.vectors)!=(int(np.sum(self.elem_count)),3):
                msg += 'vectors must have shape (number of atoms, 3)\n'
            #end if
        #end if
        if exit and len(msg)>0:
            self.error(msg)
        #end if
        return msg
    #end def check_complete
#end class Poscar



class NebPoscars(Vobj):
    """Manage POSCAR files for a chain of NEB images."""

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
    """Represent concatenated datasets in a POTCAR file."""

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
            pend = end if line_end==-1 else line_end+1
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



class VaspInput(SimulationInput,Vobj):
    """Collect the input files for a VASP calculation."""

    all_inputs  = (
        'CHGCAR', 'DYNMATFULL', 'GAMMA', 'HESSEMAT', 'ICONST', 'INCAR',
        'IRCCAR', 'KPOINTS', 'KPOINTS_ELPH', 'KPOINTS_OPT', 'KPOINTS_WAN',
        'ML_AB', 'ML_FF', 'PENALTYPOT', 'POSCAR', 'POTCAR', 'QPOINTS',
        'STOPCAR', 'TAUCAR', 'WANPROJ', 'WAVEDER',
        )
    all_outputs = (
        'IBZKPT', 'LOCPOT', 'PRJCAR', 'PCDAT', 'ELFCAR', 'TMPCAR', 'DOSCAR', 'CHGCAR',
        'CONTCAR', 'HILLSPOT', 'PROCAR', 'EIGENVAL', 'WAVECAR', 'REPORT', 'CHG',
        'OUTCAR', 'PROOUT', 'XDATCAR', 'vasprun.xml', 'OSZICAR',
        )# note that CHGCAR, TMPCAR, and WAVECAR sometimes contain input

    input_files = obj(
        chgcar       = Chgcar,
        dynmatfull   = Dynmatfull,
        gamma        = Gamma,
        hessemat     = Hessemat,
        iconst       = Iconst,
        incar        = Incar,
        irccar       = Irccar,
        kpoints      = Kpoints,
        kpoints_elph = KpointsElph,
        kpoints_opt  = KpointsOpt,
        kpoints_wan  = KpointsWan,
        ml_ab        = MlAb,
        ml_ff        = MlFf,
        penaltypot   = Penaltypot,
        poscar       = Poscar,
        potcar       = Potcar,
        qpoints      = Qpoints,
        stopcar      = Stopcar,
        taucar       = Taucar,
        wanproj      = Wanproj,
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
        image_dirs = [
            name for name in os.listdir(path)
            if len(name)==2
            and name.isdigit()
            and os.path.isfile(os.path.join(path,name,'POSCAR'))
            ]
        if len(image_dirs)>0:
            self.poscar = NebPoscars()
            self.poscar.read(path)
        #end if
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
            kpoints.coord    = 'reciprocal'
            kpoints.kpoints  = np.dot(
                structure.kpoints,np.linalg.inv(structure.kaxes)
                )
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
            if s.vel is not None:
                poscar.vel_coord = 'cartesian'
                poscar.vel = s.vel.copy()
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
        if 'poscar' not in self:
            self.error('POSCAR is required to generate a physical system')
        elif isinstance(self.poscar,NebPoscars):
            self.error(
                'return_system cannot select a structure from an NEB path; '
                'select an image POSCAR first'
                )
        #end if
        raw_axes = self.poscar.axes
        input_scale = self.poscar.scale
        if np.isscalar(input_scale) and input_scale<0:
            scale_factor = (abs(input_scale)/abs(np.linalg.det(raw_axes)))**(1.0/3.0)
        else:
            scale_factor = input_scale
        #end if
        axes  = scale_factor*raw_axes
        scale = 1.0

        pot_info = None
        velem = self.poscar.elem
        if velem is None and 'potcar' in self:
            pot_info = self.potcar.pot_info()
            if len(pot_info)>0:
                velem = np.array(
                    [pot_info[n].element for n in range(len(pot_info))]
                    )
            #end if
        #end if
        if velem is None:
            self.error(
                'POSCAR does not contain species names and they could not be '
                'recovered from POTCAR'
                )
        #end if
        velem_count = self.poscar.elem_count
        if len(velem)!=len(velem_count):
            self.error(
                'the number of POSCAR species does not match the number of '
                'species counts'
                )
        #end if
        elem        = vasp_to_nexus_elem(velem,velem_count)

        if self.poscar.coord=='direct':
            pos = np.dot(self.poscar.pos,axes)
        elif self.poscar.coord=='cartesian':
            pos = scale_factor*self.poscar.pos
        else:
            self.error('invalid POSCAR coordinate specifier: {0}'.format(self.poscar.coord))
        #end if
        vel = None
        if self.poscar.vel is not None:
            if self.poscar.vel_coord=='direct':
                vel = np.dot(self.poscar.vel,axes)
            elif self.poscar.vel_coord=='cartesian':
                vel = self.poscar.vel.copy()
            else:
                self.error(
                    'invalid POSCAR velocity coordinate specifier: {0}'
                    .format(self.poscar.vel_coord)
                    )
            #end if
        #end if

        center=axes.sum(0)/2.0

        kgrid    = None
        kshift   = None

        if structure_only:
            kpoints_file = None
        elif 'kpoints' in self:
            kpoints_file = self.kpoints
        elif 'incar' in self and 'kspacing' in self.incar:
            kpoints_file = None
        else:
            kpoints_file = None
        #end if

        if kpoints_file is not None and kpoints_file.mode=="auto":
            if 'kshift' in kpoints_file and kpoints_file.kshift is not None:
                kshift = kpoints_file.kshift
            else:
                kshift = np.zeros(3,dtype=float)
            #end if
            centering = kpoints_file.centering.lower()
            if centering.startswith('m'):
                kshift=kshift+np.array([0.5,0.5,0.5])
            elif centering.startswith('g'):
                pass
            #end if
            if centering.startswith('a'):
                reciprocal_axes = 2*np.pi*np.linalg.inv(axes).T
                kgrid = tuple(
                    max(1,int(
                        kpoints_file.kgrid*np.linalg.norm(kaxis)/(2*np.pi)+0.5
                        ))
                    for kaxis in reciprocal_axes
                    )
            else:
                kgrid=kpoints_file.kgrid
            #end if
        elif kpoints_file is not None and kpoints_file.mode not in (
            'explicit',
            ):
            self.error(
                'system generation does not currently work with KPOINTS '
                'mode: {0}'.format(kpoints_file.mode)
                )
        #end if
        structure = Structure(
            axes     = axes,
            elem     = elem,
            scale    = scale,
            pos      = pos,
            vel      = vel,
            center   = center,
            kgrid    = kgrid,
            kshift   = kshift,
            units    = 'A',
            rescale  = False,
            )

        if kpoints_file is not None and kpoints_file.mode=='explicit':
            if kpoints_file.coord=='reciprocal':
                kpoints = np.dot(kpoints_file.kpoints,structure.kaxes)
            elif kpoints_file.coord=='cartesian':
                kpoints = 2*np.pi*kpoints_file.kpoints/scale_factor
            else:
                self.error(
                    'invalid KPOINTS coordinate specifier: {0}'
                    .format(kpoints_file.coord)
                    )
            #end if
            structure.add_kpoints(
                kpoints,kpoints_file.kweights,recenter=False
                )
        elif not structure_only and kpoints_file is None:
            if 'incar' in self and 'kspacing' in self.incar:
                kspacing = self.incar.kspacing
            else:
                kspacing = 0.5
            #end if
            if 'incar' in self and 'kgamma' in self.incar:
                kgamma = self.incar.kgamma
            else:
                kgamma = True
            #end if
            kshift = (0,0,0) if kgamma else (0.5,0.5,0.5)
            structure.add_kmesh(
                kspacing=kspacing,kshift=kshift
                )
        #end if

        structure.zero_corner()
        structure.recenter()

        if structure_only:
            return structure
        #end if

        system_valency = dict(valency)
        if pot_info is None and 'potcar' in self:
            pot_info = self.potcar.pot_info()
        #end if
        if pot_info is not None and len(pot_info)==len(velem):
            for n,species in enumerate(velem):
                info = pot_info[n]
                if species not in system_valency:
                    system_valency[species] = info.Zval
                #end if
            #end for
        #end if

        for atom in set(elem):
            if atom not in system_valency:
                self.error(
                    'valence charge for atom {0} has not been defined\n'
                    'please provide the valence charge as an argument to '
                    'return_system()'.format(atom)
                    )
            #end if
        #end for

        ion_charge = sum(system_valency[atom] for atom in elem)
        if 'incar' in self and 'nelect' in self.incar:
            net_charge = ion_charge-self.incar.nelect
        else:
            net_charge = 0
        #end if
        if 'incar' in self and 'nupdown' in self.incar:
            net_spin = self.incar.nupdown
        else:
            net_spin = 0
        #end if

        system = PhysicalSystem(
            structure  = structure,
            net_charge = net_charge,
            net_spin   = net_spin,
            **system_valency
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
                    self.error(
                        'pseudopotential for element {0} not found\n'
                        'elements present: {1}'
                        .format(symbol,sorted(pseudo_map.keys()))
                        )
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
                self.error(
                    'arguments to setup NEB must either be structure or '
                    'system objects\n'
                    '  received an object of type: {0}'
                    .format(s.__class__.__name__)
                    )
            #end if
        #end for
        interp_args['repackage'] = False

        # generate NEB image structures
        if len(structures)<2:
            self.error(
                'must provide at least two structures to setup NEB\n'
                '  you provided: {0}'.format(len(structures))
                )
        elif len(structures)==2:
            incar_images = 'images' in self.incar
            kwarg_images = 'images' in interp_args
            if incar_images and kwarg_images and self.incar.images!=interp_args['images']:
                self.error(
                    'images provided in incar and setup_neb do not match\n'
                    '  please ensure they match to remove ambiguity\n'
                    '  incar images: {0}\n'
                    '  setup_neb images: {1}'
                    .format(self.incar.images,interp_args['images'])
                    )
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
                self.error(
                    'number of structures provided to setup_neb must be '
                    'consistent with number of images in INCAR\n'
                    '  INCAR images: {0}\n'
                    '  structures provided {1}'
                    .format(self.incar.images,len(neb_structures))
                    )
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
        if 'incar' not in self:
            return 'unknown'
        #end if
        incar = self.incar
        ibrion = incar.ibrion if 'ibrion' in incar else -1
        nsw = incar.nsw if 'nsw' in incar else 0
        # check for neb
        if 'images' in incar and incar.images>0:
            run_type = 'neb'
        elif ibrion in {5,6,7,8}:
            run_type = 'phonon'
        elif nsw<=0 or ibrion==-1:
            run_type = 'static'
        elif ibrion==0:
            run_type = 'md'
        elif ibrion>0:
            run_type = 'relax'
        else:
            run_type = 'unknown'
        #end if
        return run_type
    #end def run_type


    def validate(self,*,runnable=True,exit=True):
        """Check VASP input completeness and cross-file consistency."""
        messages = []
        required = ('incar','poscar','potcar') if runnable else ()
        for name in required:
            if name not in self:
                messages.append('{0} is missing'.format(name.upper()))
            #end if
        #end for

        poscars = []
        if 'poscar' in self:
            if isinstance(self.poscar,NebPoscars):
                indices = sorted(self.poscar.keys())
                if indices!=list(range(len(indices))):
                    messages.append(
                        'NEB image directories are not consecutively numbered from 00'
                        )
                #end if
                poscars = [self.poscar[n] for n in indices]
                if len(poscars)<2:
                    messages.append('NEB input requires at least two POSCARs')
                #end if
                if 'incar' not in self or 'images' not in self.incar:
                    messages.append('NEB POSCARs require IMAGES in INCAR')
                elif self.incar.images!=len(poscars)-2:
                    messages.append(
                        'INCAR IMAGES does not match the number of NEB POSCARs'
                        )
                #end if
            else:
                poscars = [self.poscar]
                if (
                    'incar' in self
                    and 'images' in self.incar
                    and self.incar.images>0
                    ):
                    messages.append(
                        'INCAR IMAGES requires POSCARs in NEB image directories'
                        )
                #end if
            #end if
        #end if

        for n,poscar in enumerate(poscars):
            prefix = 'POSCAR' if len(poscars)==1 else 'NEB POSCAR {0:02d}'.format(n)
            complete = poscar.check_complete(exit=False).strip()
            if len(complete)>0:
                messages.append('{0}: {1}'.format(
                    prefix,complete.replace('\n','; ')
                    ))
                continue
            #end if
            natoms = int(np.sum(poscar.elem_count))
            if len(poscar.pos)!=natoms:
                messages.append(
                    '{0}: atom counts do not match the number of positions'
                    .format(prefix)
                    )
            #end if
            if poscar.elem is not None and len(poscar.elem)!=len(poscar.elem_count):
                messages.append(
                    '{0}: species names do not match species counts'
                    .format(prefix)
                    )
            #end if
            if (
                poscar.dynamic is not None
                and np.shape(poscar.dynamic)!=(natoms,3)
                ):
                messages.append(
                    '{0}: selective-dynamics flags must have shape ({1}, 3)'
                    .format(prefix,natoms)
                    )
            #end if
        #end for

        if len(poscars)>1:
            reference = poscars[0]
            for n,poscar in enumerate(poscars[1:],1):
                if not np.array_equal(poscar.elem_count,reference.elem_count):
                    messages.append(
                        'NEB POSCAR {0:02d} has different species counts'
                        .format(n)
                        )
                elif (
                    reference.elem is not None
                    and poscar.elem is not None
                    and not np.array_equal(poscar.elem,reference.elem)
                    ):
                    messages.append(
                        'NEB POSCAR {0:02d} has different species names'
                        .format(n)
                        )
                #end if
                if (
                    reference.scale is not None
                    and poscar.scale is not None
                    and reference.axes is not None
                    and poscar.axes is not None
                    ):
                    reference_axes = np.asarray(reference.scale)*reference.axes
                    image_axes = np.asarray(poscar.scale)*poscar.axes
                    if not np.allclose(image_axes,reference_axes):
                        messages.append(
                            'NEB POSCAR {0:02d} has different lattice vectors'
                            .format(n)
                            )
                    #end if
                #end if
            #end for
        #end if

        if 'kpoints' in self:
            kpoints = self.kpoints
            if kpoints.mode not in ('auto','basis','line','explicit'):
                messages.append('KPOINTS mode is missing or invalid')
            elif kpoints.mode=='explicit':
                if 'kpoints' not in kpoints or 'kweights' not in kpoints:
                    messages.append(
                        'explicit KPOINTS require points and weights'
                        )
                else:
                    nkpoints = len(kpoints.kpoints)
                    if np.shape(kpoints.kpoints)!=(nkpoints,3):
                        messages.append(
                            'KPOINTS points must have shape (N, 3)'
                            )
                    #end if
                    if len(kpoints.kweights)!=nkpoints:
                        messages.append(
                            'KPOINTS weights do not match the number of points'
                            )
                    #end if
                #end if
                if (
                    'coord' not in kpoints
                    or kpoints.coord not in ('cartesian','reciprocal')
                    ):
                    messages.append(
                        'explicit KPOINTS coordinates are missing or invalid'
                        )
                #end if
            elif kpoints.mode=='auto':
                if 'centering' not in kpoints:
                    messages.append('automatic KPOINTS centering is missing')
                else:
                    centering = str(kpoints.centering).lower()
                #end if
                if 'centering' in kpoints and (
                    len(centering)==0 or centering[0] not in ('a','g','m')
                    ):
                    messages.append('automatic KPOINTS centering is invalid')
                elif 'centering' in kpoints and centering.startswith('a'):
                    if (
                        'kgrid' not in kpoints
                        or not np.isscalar(kpoints.kgrid)
                        or kpoints.kgrid<=0
                        ):
                        messages.append(
                            'fully automatic KPOINTS length must be positive'
                            )
                    #end if
                elif 'centering' in kpoints and (
                    'kgrid' not in kpoints or np.shape(kpoints.kgrid)!=(3,)
                    ):
                    messages.append('automatic KPOINTS grid must have length 3')
                #end if
            elif kpoints.mode=='basis':
                if 'kbasis' not in kpoints or np.shape(kpoints.kbasis)!=(3,3):
                    messages.append('KPOINTS basis must have shape (3, 3)')
                #end if
            elif (
                kpoints.mode=='line'
                and (
                    'kendpoints' not in kpoints
                    or np.ndim(kpoints.kendpoints)!=2
                    or np.shape(kpoints.kendpoints)[1:]!=(3,)
                    )
                ):
                messages.append('KPOINTS points must have shape (N, 3)')
            #end if
        #end if

        if (
            'incar' in self
            and 'kspacing' in self.incar
            and self.incar.kspacing<=0
            ):
            messages.append('INCAR KSPACING must be positive')
        #end if

        if 'potcar' in self and len(poscars)>0:
            try:
                pot_info = self.potcar.pot_info()
            except Exception as exception:  # noqa: BLE001
                messages.append('POTCAR metadata could not be read: {0}'.format(
                    exception
                    ))
            else:
                if len(pot_info)>0:
                    nspecies = len(poscars[0].elem_count)
                    if len(pot_info)!=nspecies:
                        messages.append(
                            'POTCAR datasets do not match POSCAR species counts'
                            )
                    elif poscars[0].elem is not None:
                        pot_elem = [
                            pot_info[n].element for n in range(len(pot_info))
                            ]
                        pos_elem = []
                        for species in poscars[0].elem:
                            iselem,element = Elements.is_element(
                                species,return_element=True
                                )
                            pos_elem.append(
                                element.symbol if iselem else species
                                )
                        #end for
                        if pot_elem!=pos_elem:
                            messages.append(
                                'POTCAR species do not match POSCAR species'
                                )
                        #end if
                    #end if
                elif runnable:
                    messages.append('POTCAR contains no datasets')
                #end if
            #end try
        #end if

        if len(messages)>0 and exit:
            self.error(
                'VASP input is invalid:\n  '+'\n  '.join(messages)
                )
        #end if
        return len(messages)==0
    #end def validate


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
        error(
            'input_type {0} is unrecognized\n'
            'valid options are: general'.format(input_type)
            )
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
    if vf.pseudos is not None:
        vf.pseudos = PseudoSet.pseudo_remap('vasp',vf.pseudos,vf.system)
    #end if
    gen_kpoints = 'kspacing' not in kwargs

    # create an empty input file
    vi = VaspInput()

    # assign values to incar and any other keyword files
    keywords = set(kwargs.keys())
    for name,keyword_file in VaspInput.keyword_files.items():
        keys = keywords & keyword_file.keywords
        keys |= {key for key in keywords if keyword_file.is_block_name(key)}
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
        elif vf.kbasis is not None:
            kp.mode   = 'basis'
            kp.kbasis = vf.kbasis
            kp.coord  = vf.kcoord
            kp.kshift = vf.kshift
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
    if isinstance(structure,PhysicalSystem):
        structure = structure.structure
    #end if
    if not isinstance(structure,Structure):
        error(
            'structure must be a Structure or PhysicalSystem\n'
            'you provided: {}'.format(type(structure).__name__),
            'generate_poscar'
            )
    #end if
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
        error(
            'coord must be either direct or cartesian\n'
            'you provided: {0}'.format(coord),
            'generate_poscar',
            )
    #end if
    if s.frozen is not None:
        poscar.dynamic = s.frozen==False
    #end if
    if s.vel is not None:
        poscar.vel_coord = 'cartesian'
        poscar.vel = s.vel.copy()
    #end if
    return poscar
#end def generate_poscar
