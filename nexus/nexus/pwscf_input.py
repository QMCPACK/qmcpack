##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  pwscf_input.py                                                    #
#    Supports I/O for PWSCF input files.                             #
#                                                                    #
#  Content summary:                                                  #
#    PwscfInput                                                      #
#      SimulationInput class for PWSCF.                              #
#      Can read/write most PWSCF input files.                        #
#                                                                    #
#    generate_pwscf_input                                            #
#      User-facing function to create arbitrary input files.         #
#                                                                    #
#    PwscfInputBase                                                  #
#      Base class for all other PWSCF input classes.                 #
#      Contains a listing of all keyword variable names and types.   #
#                                                                    #
#    Element                                                         #
#      Base class for two different types of PWSCF input sections:   #
#      namelists (Section class) and 'cards' (Card class)            #
#                                                                    #
#    Section                                                         #
#      Class representing keyword (namelist) input sections.         #
#      Base class for other keyword input section classes.           #
#      See control, system, electrons, ions, cell, phonon, and ee.   #
#                                                                    #
#    Card                                                            #
#      Class representing formatted (card) input sections.           #
#      Base class for other formatted section classes.               #
#      See atomic_species, atomic_positions, k_points,               #
#        cell_parameters, climbing_images, constraints,              #
#        collective_vars, occupations and hubbard.                   #
#                                                                    #
#    QEXML                                                           #
#      Class to represent an xml element.                            #
#                                                                    #
#    readval                                                         #
#      Function converts an attribute value string to numeric form.  #
#                                                                    #
#====================================================================#


import inspect
import os
import sys
from copy import deepcopy
from types import MappingProxyType
from typing import ClassVar, Literal, TypeAlias

import numpy as np
from numpy import pi
from numpy.linalg import inv

from . import numpy_extensions as npe
from .developer import DevBase, error, log, obj, warn
from .periodic_table import Elements
from .physical_system import PhysicalSystem
from .pseudoset import pp_elem_label, PseudoSet
from .pwscf_input_defs import (
    CellDefinitions,
    ControlDefinitions,
    ElectronsDefinitions,
    FcpDefinitions,
    IonsDefinitions,
    PwscfInputType,
    RismDefinitions,
    SystemDefinitions,
)
from .simulation import SimulationInput
from .structure import Structure, kmesh
from .unit_converter import convert

"""Union of the namelist definition enums."""
NamelistType: TypeAlias = (
    type[CellDefinitions]
    | type[ControlDefinitions]
    | type[ElectronsDefinitions]
    | type[FcpDefinitions]
    | type[IonsDefinitions]
    | type[RismDefinitions]
    | type[SystemDefinitions]
)

"""Tuple of all of the namelist definition enums."""
NAMELIST_DEFINITIONS = (
    ControlDefinitions,
    SystemDefinitions,
    ElectronsDefinitions,
    IonsDefinitions,
    CellDefinitions,
    FcpDefinitions,
    RismDefinitions,
    )

def read_str(sv):
    return sv.strip('"').strip("'")
#end def read_str


def read_int(sv):
    return int(sv)
#end def read_int


def read_float(sv):
    return float(sv.replace('d','e').replace('D','e'))
#end def read_float


bconv = {'.true.':True,'.false.':False}
def read_bool(sv):
    return bconv[sv.lower()]
#end def read_bool


def write_str(val):
    return "'"+val+"'"
#end def write_str


def write_int(val):
    return str(val)
#end def write_int


def write_float(val):
    return str(val)
    #return '{0:20.18e}'.format(val)
#end def write_float


def write_bool(val):
    if val:
        return '.true.'
    else:
        return '.false.'
    #end if
#end def write_bool


readval={str:read_str,int:read_int,float:read_float,bool:read_bool}
writeval={str:write_str,int:write_int,float:write_float,bool:write_bool}


def write_scalar(var,val):
    if isinstance(val,str):
        vtype = str
    elif isinstance(val,float):
        vtype = float
    elif var in PwscfInputBase.bools:
        vtype = bool
    elif isinstance(val,int):
        vtype = int
    else:
        error('cannot write pwscf input file\nattempted to write variable with unknown scalar type\nvariable: {0}\ndata type: {1}'.format(var,val.__class__.__name__))
    #end if
    return writeval[vtype](val)
#end def write_scalar


def noconv(v):
    return v
#end def noconv


def array_from_lines(lines):
    s=''
    for l in lines:
        s+=l+' '
    #end for
    a = np.fromstring(s,sep=' ')
    nelem = len(a)
    nlines = len(lines)
    dim = nelem//nlines
    npe.reshape_inplace(a, (nlines, dim))
    return a
#end def array_from_lines


pwscf_precision = '16.8f'
pwscf_array_format = '{0:'+pwscf_precision+'}' 
def array_to_string(a,pad='   ',format=pwscf_array_format,converter=noconv,rowsep='\n'):
    s=''
    if len(a.shape)==1:
        s+=pad
        for v in a:
            s += format.format(converter(v))+' '
        #end for
        s+=rowsep
    else:
        for r in a:
            s+=pad
            for v in r:
                s += format.format(converter(v))+' '
            #end for
            s+=rowsep
        #end for
    #end if
    return s
#end def array_to_string


def _get_var_types() -> (
    tuple[
        frozenset[str], # ints
        frozenset[str], # floats
        frozenset[str], # strs
        frozenset[str], # bools
        frozenset[str], # real_arrays
        frozenset[str], # species_arrays
        frozenset[str], # multidimensional_arrays
        ]
    ):
    """Get all variable types from the namelist definitions.

    Returns
    -------
    ints : frozenset of str
        All integer-type variables. (scalars and arrays)
    floats : frozenset of str
        All float-type variables. (scalars and arrays)
    strs : frozenset of str
        All string-type variables. (scalars and arrays)
    bools : frozenset of str
        All boolean-type variables. (scalars and arrays)
    real_arrays : frozenset of str
        All array variables. (all datatypes)
    species_arrays : frozenset of str
        All array variables whose size depends on species. (all datatypes)
    multidimensional_arrays : frozenset of str
        All multidimensional array variables. (all datatypes, fixed-size and species arrays)
    """
    ints = set()
    floats = set()
    strs = set()
    bools = set()
    real_arrays = set()
    species_arrays = set()
    multidimensional_arrays = set()

    for nmlist in NAMELIST_DEFINITIONS:
        for var in nmlist:
            if var.shape is not None:
                real_arrays.add(var.name)

                if isinstance(var.shape[1], str):
                    species_arrays.add(var.name)
                elif (
                    isinstance(var.shape[1], tuple)
                    and (
                        isinstance(var.shape[1][0], str)
                        or isinstance(var.shape[1][1], str)
                        or isinstance(var.shape[1][-1], str)
                        )
                    ):
                    species_arrays.add(var.name)

                if isinstance(var.shape[0], tuple) and len(var.shape[0]) >= 2:
                    multidimensional_arrays.add(var.name)
            #end if var.shape is not None

            if var.datatype is int:
                ints.add(var.name)
            elif var.datatype is float:
                floats.add(var.name)
            elif var.datatype is str:
                strs.add(var.name)
            elif var.datatype is bool:
                bools.add(var.name)
            else:
                msg = f"Variable {var.name} has an invalid datatype `{var.datatype.__name__}`!"
                raise ValueError(msg)

    return (
        frozenset(ints),
        frozenset(floats),
        frozenset(strs),
        frozenset(bools),
        frozenset(real_arrays),
        frozenset(species_arrays),
        frozenset(multidimensional_arrays),
        )
#end def get_var_types


class PwscfInputBase(DevBase):
    (ints, floats, strs, bools, real_arrays,
     species_arrays, multidimensional_arrays) = _get_var_types()

    species_array_indices = obj(
        hubbard_j=1,
        starting_ns_eigenvalue=2,
        )

    all_variables = ints | floats | strs | bools

    section_aliases = MappingProxyType(dict(
        celldm1='celldm(1)',
        celldm2='celldm(2)',
        celldm3='celldm(3)',
        celldm4='celldm(4)',
        celldm5='celldm(5)',
        celldm6='celldm(6)'
        ))

    var_types = dict()  # noqa: RUF012
    for v in ints:
        var_types[v]=int
    #end for
    for v in floats:
        var_types[v]=float
    #end for
    for v in strs:
        var_types[v]=str
    #end for
    for v in bools:
        var_types[v]=bool
    #end for
    var_types: MappingProxyType[str, type] = MappingProxyType(var_types)
#end class PwscfInputBase




class Element(PwscfInputBase):
    name = None
    def add(self,**variables):
        self.update(**variables)
    #end def add

    def read(self,lines):
        raise NotImplementedError
    #end def read

    def write(self,parent):
        raise NotImplementedError
    #end def write

    def post_process_read(self,parent):
        None
    #end def post_process_read
#end class Element


class Section(Element):
    defs: ClassVar[NamelistType]
    @classmethod
    def class_init(cls):
        cls.variables = frozenset(cls.defs.__members__.keys())
        cls.case_map  = {name: val.input_name for name, val in cls.defs.__members__.items()}
    #end def class_init

    def assign(self,**variables):
        self.update(**variables)
    #end def assign


    def read(self,lines):
        for l in lines:
            # exclude comments
            cloc = l.find('!')
            if cloc!=-1:
                l = l[:cloc]
            #end if
            # parse tokens accounting for significant whitespace
            #  replace commas outside array brackets with spaces
            if '(' not in l:
                lout = l
            else:
                lout = ''
                inparen=False
                for c in l:
                    if c=='(':
                        inparen=True
                        lout+=c
                    elif c==')':
                        inparen=False
                        lout+=c
                    elif inparen and c==',':
                        lout+='|' # new delimiter
                    else:
                        lout += c
                    #end if
                #end for
            #end if
            tokens = lout.split(',')
            for t in tokens:
                if len(t)>0:
                    tsplt = t.split('=')
                    if len(tsplt)!=2:
                        self.error('attempted to read misformatted line\nmisformatted line: {0}\ntokens: {1}'.format(l,tsplt))
                    #end if
                    var,val = tsplt
                    var     = var.strip().lower()
                    val     = val.strip()
                    if '(' not in var:
                        varname = var
                        index   = None
                    else:
                        var = var.replace('|',',')
                        varname = var.split('(')[0]
                        index   = var.replace(varname,'').strip('()')
                        if ',' not in index:
                            index = int(index)
                        else:
                            index = tuple(np.array(index.split(','),dtype=int))
                        #end if
                    #end if
                    if varname not in self.variables:
                        self.error('pwscf input section {0} does not have a variable named "{1}", please check your input\nif correct, please add a new variable ({1}) to the {0} PwscfInput class'.format(self.__class__.__name__,varname),trace=False)
                    #end if
                    if varname not in self.var_types:
                        self.error('a type has not been specified for variable "{0}"\nplease add it to PwscfInputBase'.format(varname),trace=False)
                    #end if
                    vtype = self.var_types[varname]
                    val = readval[vtype](val)
                    if varname not in self.real_arrays:
                        self[var] = val
                    else:
                        if index is None and '(' not in var:
                            index = 1
                        #end if
                        if varname not in self:
                            self[varname] = obj()
                        #end if
                        self[varname][index] = val
                    #end if
                #end for
            #end for
        #end for
    #end def read


    def write(self,parent):
        atom_index = obj()
        if 'atomic_species' in parent and 'atoms' in parent.atomic_species:
            for i,a in enumerate(parent.atomic_species.atoms):
                atom_index[a] = i+1
            #end for
        #end if
        cls = self.__class__
        c='&'+self.name.upper()+'\n'
        vars = list(self.keys())
        vars.sort()
        for var in vars:
            val = self[var]
            vname = var
            if var in cls.case_map:
                vname = cls.case_map[var]
            #end if
            if var not in self.real_arrays:
                # write scalar values
                sval = write_scalar(var,val)
                c+='   '+'{0:<15} = {1}\n'.format(vname,sval)
            else:
                # write array values
                allow_spec = var in self.species_arrays
                index_map = obj()
                for index in val.keys():
                    index_map[index] = index
                #end for
                index_map_inv = index_map
                if var in self.species_arrays:
                    for index in val.keys():
                        if isinstance(index,str):
                            if index in atom_index:
                                index_map[index] = atom_index[index]
                            else:
                                self.error('cannot write pwscf input\ninvalid array species index encountered\nspecies index provided is not in the set of species present\nspecies present: {0}\nspecies used as index: {1}\narray variable: {2}'.format(sorted(atom_index.keys()),index),var)
                            #end if
                        elif isinstance(index,tuple):
                            if var not in self.species_array_indices:
                                self.error('cannot write pwscf input\ninvalid multidimensional array species index encountered\narray variable "{0}" does not support multidimensional species indices\nindex received: {1}'.format(var,index))
                            #end if
                            indloc = self.species_array_indices[var]
                            atom = index[indloc]
                            if not isinstance(atom,str):
                                continue
                            #end if
                            if atom not in atom_index:
                                self.error('cannot write pwscf input\ninvalid array species index encountered\nspecies index provided is not in the set of species present\nspecies present: {0}\nspecies used as index: {1}\nfull index provided: {2}\narray variable: {3}'.format(sorted(atom_index.keys()),atom,index),var)
                            #end if
                            indlist = list(index)
                            indlist[indloc] = atom_index[atom]
                            index_map[index] = tuple(indlist)
                        #end if
                    #end for
                    index_map_inv = obj({v:k for k,v in index_map.items()})
                #end if

                for index in sorted(index_map_inv.keys()):
                    value = val[index_map_inv[index]]
                    if isinstance(index,int):
                        sind = '({0})'.format(index)
                    elif isinstance(index,tuple):
                        if not allow_spec:
                            None
                        #end if
                        sind = str(index).replace(' ','')
                    else:
                        self.error('cannot write pwscf input\ninvalid array index encountered\nmust be an integer or tuple of integers\nindex received: {0}\narray variable: {1}'.format(str(index)),var)
                    #end if
                    svar = vname+sind
                    sval = write_scalar(var,value)
                    c+='   '+'{0:<15} = {1}\n'.format(svar,sval)
                #end for
            #end if
        #end for
        c+='/'+'\n\n'
        return c
    #end def write
#end class Section


class Card(Element):
    def __init__(self):
        self.specifier = ''
    #end def __init__

    def get_specifier(self,line):
        tokens = line.split()
        if len(tokens)>1:
            self.specifier = tokens[1].strip('{}()').lower()
        #end if
    #end def get_specifier

    def read(self,lines):
        self.get_specifier(lines[0])
        self.read_text(lines[1:])
    #end def read

    def write(self,parent):
        c = self.name.upper()+' '+self.specifier+'\n'
        c += self.write_text()+'\n'
        return c
    #end def write

    def read_text(self,lines):
        raise NotImplementedError
    #end def read_text

    def write_text(self):
        raise NotImplementedError
    #end def write_text

    def change_specifier(self,new_specifier):
        raise NotImplementedError
    #end def change_specifier

    def change_option(self,*args,**kwargs):
        self.change_specifier(*args,**kwargs)
    #end def change_option
#end class Card


class control(Section):
    name = 'control'
    defs = ControlDefinitions
#end class control



class system(Section):
    name = 'system'
    defs = SystemDefinitions
    atomic_variables = obj(
        hubbard_u = 'Hubbard_U',
        start_mag = 'starting_magnetization',
        hubbard_j = 'Hubbard_J',
        angle1    = 'angle1',
        angle2    = 'angle2',
        )

    # specialized read for partial array variables (hubbard U, starting mag, etc)
    def post_process_read(self,parent):
        if 'atomic_species' in parent:
            keys = self.keys()
            for alias,name in self.atomic_variables.items():
                has_var = False
                avals = obj()
                akeys = []
                for key in keys:
                    if key.startswith(name):
                        akeys.append(key)
                        if '(' not in key:
                            index=1
                        else:
                            index = int(key.replace(name,'').strip('()'))
                        #end if
                        avals[index] = self[key]
                    #end if
                #end for
                if has_var:
                    for key in akeys:
                        del self[key]
                    #end for
                    atoms = parent.atomic_species.atoms
                    value = obj()
                    for i in range(len(atoms)):
                        index = i+1
                        if index in avals:
                            value[atoms[i]] = avals[i+1]
                        #end if
                    #end for
                    self[alias] = value
                #end if
            #end for
        #end if
    #end def post_process_read


    # specialized write for odd handling of hubbard U
    def write_old(self,parent):
        cls = self.__class__
        c='&'+self.name.upper()+'\n'
        vars = list(self.keys())
        vars.sort()
        for var in vars:
            val = self[var]
            if var in self.atomic_variables:
                if val is None: # jtk mark: patch fix for generate_pwscf_input
                    continue    #    i.e. hubbard_u = None
                #end if
                if 'atomic_species' in parent:
                    atoms = parent.atomic_species.atoms
                    avar = self.atomic_variables[var]
                    for i in range(len(atoms)):
                        index = i+1
                        vname = '{0}({1})'.format(avar,index)
                        atom = atoms[i]
                        if atom in val:
                            sval = writeval[float](val[atom])
                            c+='   '+'{0:<15} = {1}\n'.format(vname,sval)
                        #end if
                    #end for
                else:
                    self.error('cannot write {0}, atomic_species is not present'.format(var))
                #end if
            else:
                #vtype = type(val)
                #sval = writeval[vtype](val)

                vtype = None
                if isinstance(val,str):
                    vtype = str
                elif isinstance(val,float):
                    vtype = float
                #elif isinstance(val,bool):
                #    vtype = bool
                elif var in self.bools:
                    vtype = bool
                elif isinstance(val,int):
                    vtype = int
                else:
                    self.error('Type "{0}" is not known as a value of variable "{1}".\nThis may reflect a need for added developer attention to support this type.  Please contact a developer.'.format(vtype.__class__.__name__,var))
                #end if
                sval = writeval[vtype](val)

                if var in self.section_aliases.keys():
                    vname = self.section_aliases[var]
                else:
                    vname = var
                #end if
                if vname in cls.case_map:
                    vname = cls.case_map[vname]
                #end if
                #c+='   '+vname+' = '+sval+'\n'
                c+='   '+'{0:<15} = {1}\n'.format(vname,sval)
            #end if
        #end for
        c+='/'+'\n\n'
        return c
    #end def write
#end class system


class electrons(Section):
    name = 'electrons'
    defs = ElectronsDefinitions
#end class electrons


class ions(Section):
    name = 'ions'
    defs = IonsDefinitions
#end class ions


class cell(Section):
    name = 'cell'
    defs = CellDefinitions
#end class cell


class fcp(Section):
    name = "fcp"
    defs = FcpDefinitions
#end class fcp


class rism(Section):
    name = "rism"
    defs = RismDefinitions
#end class rism


section_classes = (
    control, system, electrons, ions, cell, fcp, rism
    )
for sec in section_classes:
    sec.class_init()


def check_section_classes(*,exit=True):
    sections = section_classes
    all_variables = PwscfInputBase.all_variables
    global_missing = set(all_variables)
    local_missing = obj()
    locs_missing = False
    secs = obj()
    for section in sections:
        variables = section.variables
        global_missing -= variables
        loc_missing = variables - all_variables
        local_missing[section.name] = loc_missing
        locs_missing |= len(loc_missing)>0
        secs[section.name] = section
    #end for
    if len(global_missing)>0 or locs_missing:
        msg = 'PwscfInput: variable information is not consistent for section classes\n'
        if len(global_missing)>0:
            msg+='  some typed variables have not been assigned to a section:\n    {0}\n'.format(sorted(global_missing))
        #end if
        if locs_missing:
            for name in sorted(local_missing.keys()):
                lmiss = local_missing[name]
                if len(lmiss)>0:
                    vmiss = []
                    for vname in secs[name].variables:
                        if vname in lmiss:
                            vmiss.append(vname)
                        #end if
                    #end for
                    msg+='  some variables in section {0} have not been assigned a type\n    missing variable counts: {1} {2}\n    missing variables: {3}\n'.format(name,len(lmiss),len(vmiss),vmiss)
                #end if
            #end for
        #end if
        error(msg)
    else:
        log('pwscf input checks passed')
    #end if
    if exit:
        sys.exit()
    #end if
#end def check_section_classes
#check_section_classes()



class atomic_species(Card):
    name = 'atomic_species'

    def read_text(self,lines):
        atoms = []
        masses   = obj()
        pseudopotentials = obj()
        for l in lines:
            tokens = l.split()
            atom = tokens[0]
            atoms.append(tokens[0])
            masses[atom] = read_float(tokens[1])
            pseudopotentials[atom] = tokens[2]
        #end for
        self.add(atoms=atoms,masses=masses,pseudopotentials=pseudopotentials)
    #end def read_text

    def write_text(self):
        c = ''
        for at in self.atoms:
            c += '   '+'{0:2}'.format(at)+' '+str(self.masses[at])+' '+self.pseudopotentials[at]+'\n'
        #end for
        return c
    #end def write_text
#end class atomic_species



class atomic_positions(Card):
    name = 'atomic_positions'

    def read_text(self,lines):
        npos = len(lines)
        dim = 3
        atoms = []
        positions = np.empty((npos,dim))
        relax_directions = np.empty((npos,dim),dtype=int)
        i=0
        has_relax_directions = False
        for l in lines:
            tokens = l.split()
            atoms.append(tokens[0])
            positions[i,:] = np.array(tokens[1:4],dtype=np.float64)
            has_relax_directions = len(tokens)==7
            if has_relax_directions:
                relax_directions[i,:] = np.array(tokens[4:7],dtype=int)
            #end if
            i+=1
        #end for
        self.add(atoms=atoms,positions=positions)
        if has_relax_directions:
            self.add(relax_directions=relax_directions)
        #end if
    #end def read_text


    def write_text(self):
        c = ''
        has_relax_directions = 'relax_directions' in self
        if has_relax_directions:
            rowsep = ''
        else:
            rowsep = '\n'
        #end if
        for i in range(len(self.atoms)):
            c +='   '+'{0:2}'.format(self.atoms[i])+' '
            c += array_to_string(self.positions[i],pad='',rowsep=rowsep)
            if has_relax_directions:
                c += array_to_string(self.relax_directions[i],pad='',format='{0}')
            #end if
        #end for
        return c
    #end def write_text


    def change_specifier(self,new_specifier,pwi):
        scale = pwi.get_common_vars('scale')

        pos = self.positions

        spec = self.specifier
        if spec=='alat' or spec=='':
            pos *= scale
        elif spec=='bohr':
            None
        elif spec=='angstrom':
            pos *= convert(1.,'A','B')
        elif spec=='crystal':
            axes = pwi.get_common_vars('axes')
            pos = np.dot(pos,axes)
        else:
            self.error('old specifier for atomic_positions is invalid\n  old specifier: '+spec+'\n  valid options: alat, bohr, angstrom, crystal')
        #end if

        spec = new_specifier
        if spec=='alat' or spec=='':
            pos /= scale
        elif spec=='bohr':
            None
        elif spec=='angstrom':
            pos /= convert(1.,'A','B')
        elif spec=='crystal':
            axes = pwi.get_common_vars('axes')
            pos = np.dot(pos,inv(axes))
        else:
            self.error('new specifier for atomic_positions is invalid\n  new specifier: '+spec+'\n  valid options: alat, bohr, angstrom, crystal')
        #end if
            
        self.positions = pos
        self.specifier = new_specifier
    #end def change_specifier

#end class atomic_positions



class atomic_forces(Card):
    name = 'atomic_forces'

    def read_text(self,lines):
        npos = len(lines)
        dim = 3
        atoms = []
        forces = np.empty((npos,dim))
        i=0
        for l in lines:
            tokens = l.split()
            atoms.append(tokens[0])
            forces[i,:] = np.array(tokens[1:4],dtype=np.float64)
            i+=1
        #end for
        self.add(atoms=atoms,forces=forces)
    #end def read_text


    def write_text(self):
        c = ''
        rowsep = '\n'
        for i in range(len(self.atoms)):
            c +='   '+'{0:2}'.format(self.atoms[i])+' '
            c += array_to_string(self.forces[i],pad='',rowsep=rowsep)
        #end for
        return c
    #end def write_text
#end class atomic_forces



class k_points(Card):
    name = 'k_points'

    def read_text(self,lines):
        if self.specifier in {'tpiba','crystal','tpiba_b','crystal_b',''}:
            self.nkpoints = int(lines[0])
            a = array_from_lines(lines[1:])
            self.kpoints = a[:,0:3]
            self.weights = a[:,3]
            npe.reshape_inplace(self.weights, (self.nkpoints,))
        elif self.specifier == 'automatic':
            a = np.fromstring(lines[0],sep=' ')
            self.grid  = a[0:3]
            self.shift = a[3:]
        elif self.specifier == 'gamma':
            None
        else:
            self.error('k_points specifier '+self.specifier+' is unrecognized')
        #end if
    #end def read_text


    def write_text(self):
        c = ''        
        if self.specifier in {'tpiba','crystal','tpiba_b','crystal_b',''}:
            self.nkpoints = len(self.kpoints)
            c+='   '+str(self.nkpoints)+'\n'
            a = np.empty((self.nkpoints,4))
            a[:,0:3] = self.kpoints
            a[:,3]   = self.weights[:]
            c+=array_to_string(a)
        elif self.specifier == 'automatic':
            c+='   '
            c+=array_to_string(np.array(self.grid),pad='',format='{0}',converter=int,rowsep='')
            c+=array_to_string(np.array(self.shift),pad=' ',format='{0}',converter=int)
        elif self.specifier == 'gamma':
            None
        else:
            self.error('k_points specifier '+self.specifier+' is unrecognized')
        #end if
        return c
    #end def write_text


    def change_specifier(self,new_specifier,pwi):
        scale,kaxes = pwi.get_common_vars('scale','kaxes')

        spec = self.specifier
        if spec=='tpiba' or spec=='':
            kpoints = self.kpoints*(2*pi)/scale
        elif spec=='gamma':
            kpoints = np.array([[0,0,0]])
        elif spec=='crystal':
            kpoints = np.dot(self.kpoints,kaxes)
        elif spec=='automatic':
            grid  = np.array(self.grid,dtype=int)
            shift = .5*np.array(self.shift)
            kpoints = kmesh(kaxes,grid,shift)
        elif spec=='tpiba_b' or spec=='crystal_b':
            self.error('specifiers tpiba_b and crystal_b have not yet been implemented in change_specifier')
        else:
            self.error('old specifier for k_points is invalid\n  old specifier: '+spec+'\n  valid options: tpiba, gamma, crystal, automatic, tpiba_b, crystal_b')
        #end if

        spec = new_specifier
        if spec=='tpiba' or spec=='':
            kpoints = kpoints/((2*pi)/scale)
        elif spec=='gamma':
            kpoints = np.array([[0,0,0]])
        elif spec=='crystal':
            kpoints = np.dot(kpoints,inv(kaxes))
        elif spec=='automatic':
            if self.specifier!='automatic':
                self.error('cannot map arbitrary kpoints into a Monkhorst-Pack mesh')
            #end if
        elif spec=='tpiba_b' or spec=='crystal_b':
            self.error('specifiers tpiba_b and crystal_b have not yet been implemented in change_specifier')
        else:
            self.error('new specifier for k_points is invalid\n  new specifier: '+spec+'\n  valid options: tpiba, gamma, crystal, automatic, tpiba_b, crystal_b')
        #end if
            
        self.kpoints   = kpoints
        self.specifier = new_specifier
    #end def change_specifier
#end class k_points




class cell_parameters(Card):
    name = 'cell_parameters'

    def read_text(self,lines):
        self.vectors = array_from_lines(lines)
    #end def read_text

    def write_text(self):
        return array_to_string(self.vectors)
    #end def write_text

    def change_specifier(self,new_specifier,pwi):
        scale = pwi.get_common_vars('scale')

        vec = self.vectors

        spec = self.specifier
        if spec=='alat' or spec=='':
            vec *= scale
        elif spec=='bohr':
            None
        elif spec=='angstrom':
            vec *= convert(1.,'A','B')
        else:
            self.error('old specifier for cell_parameters is invalid\nold specifier: '+spec+'\nvalid options: alat, bohr, angstrom')
        #end if

        spec = new_specifier
        if spec=='alat' or spec=='':
            vec /= scale
        elif spec=='bohr':
            None
        elif spec=='angstrom':
            vec /= convert(1.,'A','B')
        else:
            self.error('new specifier for cell_parameters is invalid\nnew specifier: '+spec+'\nvalid options: alat, bohr, angstrom')
        #end if
            
        self.vectors   = vec
        self.specifier = new_specifier
    #end def change_specifier
#end class cell_parameters




class climbing_images(Card):
    name = 'climbing_images'

    def read_text(self,lines):
        self.images = array_from_lines(lines)
    #end def read_text

    def write_text(self):
        c='   '
        for n in self.images:
            c+=str(int(n))+' '
        #end for
        return c
    #end def write_text
#end class climbing_images




class constraints(Card):
    name = 'constraints'

    def read_text(self,lines):
        tokens = lines[0].split()
        self.ncontraints = int(tokens[0])
        if len(tokens)>1:
            self.tolerance = float(tokens[1])
        #end if
        self.constraints = obj()
        for i in range(len(lines)-1):
            tokens = lines[i+1].split()
            cons = obj()
            cons.type = tokens[0]
            cons.parameters = np.array(tokens[1:],dtype=np.float64)
            self.constraints[i] = cons
        #end for
    #end def read_text

    def write_text(self):
        c = '   '+str(self.nconstraints)
        if 'tolerance' in self:
            c+=' '+str(self.tolerance)
        #end if
        for cons in self.constraints:
            c+='   '+cons.type+' '+array_to_string(cons.parameters,pad='')
        #end for
        return c
    #end def write_text
#end class constraints




class collective_vars(Card):
    name = 'collective_vars'

    def read_text(self,lines):
        tokens = lines[0].split()
        self.ncontraints = int(tokens[0])
        if len(tokens)>1:
            self.tolerance = float(tokens[1])
        #end if
        self.collective_vars = obj()
        for i in range(len(lines)-1):
            tokens = lines[i+1].split()
            collv = obj()
            collv.type = tokens[0]
            collv.parameters = np.array(tokens[1:],dtype=np.float64)
            self.collective_vars[i] = collv
        #end for
    #end def read_text

    def write_text(self):
        c= '   '+str(self.ncollective_vars)
        if 'tolerance' in self:
            c+=' '+str(self.tolerance)
        #end if
        for collv in self.collective_vars:
            c+='   '+collv.type+' '+array_to_string(collv.parameters,pad='')
        #end for        
        return c
    #end def write_text
#end class collective_vars




class occupations(Card):
    name = 'occupations'

    def read_text(self,lines):
        self.occupations = array_from_lines(lines)
    #end def read_text
 
    def write_text(self):
        return array_to_string(self.occupations)
    #end def write_text
#end class occupations


class hubbard(Card):
    name = 'hubbard'
    available_specifiers = ('atomic', 'ortho-atomic', 'norm-atomic', 'wf', 'pseudo')
    default_specifier = 'atomic'
    system = None
    def read_text(self, lines):        
        contents = ''
        self.hubbard = {}
        for line in lines:
            line = line.strip().split()
            if line:
                intrxn = line[0]
                if len(line) == 3:
                    specie = line[1]
                    val = float(line[2])
                    self.hubbard[intrxn] = {specie:val}
                elif len(line) == 6:
                    specie1 = line[1]
                    specie2 = line[2]
                    ind1 = int(line[3])
                    ind2 = int(line[4])
                    val = float(line[5])
                    if intrxn not in self.hubbard.keys():
                        self.hubbard[intrxn] = {(specie1, specie2):[{'indices':(ind1, ind2), 'value':val}]}
                    else:
                        if (specie1, specie2) not in self.hubbard[intrxn].keys():
                            self.hubbard[intrxn][(specie1, specie2)]=[{'indices':(ind1, ind2), 'value':val}]
                        else:
                            self.hubbard[intrxn][(specie1, specie2)].append({'indices':(ind1, ind2), 'value':val})
                        #end if 
                    #end if
            #end for
        #end for
    #end def read_text

    def write_text(self):
        manifold_dict = {} 
        contents = ''
        for param, interaction in self.hubbard.items():
            valid_format = True
            assert(param in {'U', 'J', 'V'})
            assert(isinstance(interaction, dict))
            for label_manifold, value in interaction.items():
                if isinstance(label_manifold, str):
                    # Ex: {'U':{'C-2p': 1.0}}
                    assert(isinstance(value, (int, float)))
                    contents += f"{param} {label_manifold} {value} \n"
                elif isinstance(label_manifold, tuple):
                    assert(len(label_manifold) == 2)
                    assert(all([isinstance(_, str) for _ in label_manifold]))
                    if isinstance(value, (int, float)):
                        # Ex: {'V' : {('C-2p', 'C-2p'): 1e-8}}
                        atom1, manifold1 = label_manifold[0].split('-')
                        atom2, manifold2 = label_manifold[1].split('-')
                        if atom1 not in manifold_dict.keys():
                            manifold_dict[atom1] = [manifold1]
                        elif manifold1 not in manifold_dict[atom1]:
                            manifold_dict[atom1].append(manifold1)
                        #end if
                        if atom2 not in manifold_dict.keys():
                            manifold_dict[atom2] = manifold2
                        elif manifold2 not in manifold_dict[atom2]:
                            manifold_dict[atom2].append(manifold2)
                        #end if
                        elem = self.system.structure.elem
                        index1 = np.where(elem == atom1)[0] + 1
                        index2 = np.where(elem == atom2)[0] + 1
                        combs = []
                        for in1 in index1:
                            for in2 in index2:
                                combs.append([in1, in2])
                            #end for
                        #end for
                        for ind1, ind2 in combs:
                            contents += f"{param} {label_manifold[0]} {label_manifold[1]} {ind1} {ind2} {value}\n"
                        #end for
                    elif isinstance(value, list):
                        for val in value:
                            if 'indices' in val.keys() and 'value' in val.keys():
                                # Ex: {'V' : {('C-2p', 'C-2p'): [{'indices':(1,2), 'value':1e-8}]}
                                ind1 = val['indices'][0]
                                ind2 = val['indices'][1]
                                val = val['value']
                                contents += f"{param} {label_manifold[0]} {label_manifold[1]} {ind1} {ind2} {val}\n"
                            elif 'radius' in val.keys() and 'value' in val.keys():
                                # Ex: {'V' : {('C-2p', 'C-2p'): [{'radius':4.5, 'value':1e-8}]} radius is in Bohr
                                atom1, manifold1 = label_manifold[0].split('-')
                                atom2, manifold2 = label_manifold[1].split('-')
                                elem = self.system.structure.elem
                                index1 = np.where(elem == atom1)[0]
                                index2 = np.where(elem == atom2)[0]
                                nn = self.system.structure.nearest_neighbors(rmax = val['radius'])
                                combs = []
                                for in1 in index1:
                                    for in2 in index2:
                                        if in2 in nn[in1]:
                                            combs.append([in1+1, in2+1])
                                        #end if
                                    #end for
                                #end for
                                # import pdb
                                # pdb.set_trace()
                                for ind1, ind2 in combs:
                                    contents += f"{param} {label_manifold[0]} {label_manifold[1]} {ind1} {ind2} {val['value']}\n"
                                #end for
                            else:
                                valid_format = False
                            #end if 
                        #end for
                    else:
                        valid_format = False
                    #end if 
                else:
                    valid_format = False
                #end for
                if not valid_format:
                    self.error('Hubbard card unknown input format')
            #end for
        #end for
        for key, value in manifold_dict.items():
            if len(value) > 2:
                self.error('Element "{}" has more than 2 Hubbard manifolds "{}". Up to 3 manifolds are allowed in QE 7.1, but in that case \
                2nd and 3rd manifolds must be defined as one effective manifold, e.g. "U Mn-3d 5.0" and "U Mn-3p-3s 3.0"'.format(key, value))
            #end if
        #end for
        contents += '\n'
        return contents
    #end def write_text
#end class hubbard







class PwscfInput(SimulationInput):

    sections = ('control','system','electrons','ions','cell','fcp','rism')
    cards    = ('atomic_species','atomic_positions','atomic_forces',
                'k_points','cell_parameters','climbing_images','constraints',
                'collective_vars','occupations', 'hubbard')

    section_types = obj(
        control   = control,
        system    = system,
        electrons = electrons,
        ions      = ions,
        cell      = cell,
        fcp       = fcp,
        rism      = rism,
        )
    card_types = obj(
        atomic_species   = atomic_species,
        atomic_positions = atomic_positions,
        atomic_forces    = atomic_forces,
        k_points         = k_points,
        cell_parameters  = cell_parameters,
        climbing_images  = climbing_images,
        constraints      = constraints,
        collective_vars  = collective_vars,
        occupations      = occupations,
        hubbard          = hubbard,
        )

    element_types = obj(**section_types)
    element_types.update(**card_types)

    required_elements = ('control','system','electrons','atomic_species','atomic_positions','k_points')
    def __init__(self,*elements):
        elements = list(elements)
        if len(elements)==1 and os.path.exists(elements[0]):
            self.read(elements[0])
        elif len(elements)==1 and ('.' in elements[0] or '/' in elements[0]):
            self.error('input file '+elements[0]+' does not exist')
        else:
            for element in self.required_elements:
                if element not in elements:
                    elements.append(element)
                #end if
            #end for
            for element in elements:
                if element in self.element_types:
                    self[element] = self.element_types[element]()
                else:
                    self.error('  Error: '+element+' is not a pwscf element\n  valid options are '+str(list(self.element_types.keys())))
                #end if
            #end for
        #end if
    #end def __init__


    def read_text(self,contents,filepath=None):
        lines = contents.splitlines()
        in_element = False
        elem_type = None
        for line in lines:
            l = line.strip()
            if len(l)>0 and l[0]!='!':
                tokens = l.split()
                if l.startswith('&'):
                    l=l.strip('/').strip(" ") # allow all default section such as &ions /
                    if l[1:].lower() in self.sections:
                        prev_type = elem_type
                        in_element = True
                        elem_name = l[1:].lower()
                        elem_type = 'section'
                        c=[]
                    else:
                        self.error('encountered unrecognized input section during read\n{0} is not a recognized pwscf section\nfile read failed'.format(l[1:]))
                    #end if
                elif tokens[0].lower() in self.cards and '=' not in l:
                    if elem_type == 'card':
                        if elem_name not in self:
                            self[elem_name] = self.element_types[elem_name]()
                        #end if
                        self[elem_name].read(c)
                    #end if
                    in_element = True
                    elem_name = tokens[0].lower()
                    elem_type = 'card'
                    c=[l]
                elif l=='/':
                    in_element = False
                    if elem_name not in self:
                        self[elem_name] = self.element_types[elem_name]()
                    #end if
                    self[elem_name].read(c)
                elif in_element:
                    c.append(l)
                else:
                    self.error('invalid line encountered during read\ninvalid line: {0}\nfile read failed'.format(l))
                #end if
            #end if
        #end for
        if elem_type == 'card':
            if elem_name not in self:
                self[elem_name] = self.element_types[elem_name]()
            #end if
            self[elem_name].read(c)
        #end if
        #post-process hubbard u and related variables
        #for element in self:
        #    element.post_process_read(self)
        ##end for
    #end def read_text


    def write_text(self,filepath=None):
        contents = ''
        for s in self.sections:
            if s in self:
                contents += self[s].write(self)
            #end if
        #end for
        contents+='\n'
        for c in self.cards:
            if c in self:
                contents += self[c].write(self)
            #end if
        #end for
        contents+='\n'
        return contents
    #end def write_text


    def get_common_vars(self,*vars):
        scale = 1.0
        axes  = None
        kaxes = None

        if 'celldm(1)' in self.system:
            scale = self.system['celldm(1)']
        #end if
        if 'cell_parameters' in self:
            axes = scale*np.array(self.cell_parameters.vectors)
            kaxes = 2*pi*inv(axes).T
        #end if

        vals = []
        loc = locals()
        msg = ""
        for var in vars:
            if var in loc:
                val = loc[var]
                if val is None:
                    msg += 'requested variable '+var+' was not found\n'
                #end if
            else:
                msg += 'requested variable '+var+' is not computed by get_common_vars\n'
                val = None
            #end if
            vals.append(val)
        #end for
        if len(msg) > 0:
            self.error(f'could not get requested variables:\n{msg}')
        #end if
        return vals
    #end def get_common_vars

    def incorporate_hubbard(self, hubbard_result):
        hub_obj = hubbard()
        hubbard_result = hubbard_result.split('\n')
        hub_obj.specifier = hubbard_result[1].split()[-1]
        hub_obj.read_text(hubbard_result[2:])
        self.hubbard = hub_obj
    #end def incorporate_hubbard

    def incorporate_system(self,system,elem_order=None):
        system.check_folded_system()
        system.change_units('B')
        s  = system.structure
        nc = system.net_charge
        ns = system.net_spin

        nup = system.n_up
        ndn = system.n_down

        self.system.ibrav        = 0
#        self.system['celldm(1)'] = 1.0e0
        nions = system.n_ions
        nspecies = system.n_species
        self.system.nat          = nions
        self.system.ntyp         = nspecies
        self.system.tot_charge   = nc

        if 'cell_parameters' not in self:
            self.cell_parameters = self.element_types['cell_parameters']()
        #end if
        self.cell_parameters.specifier = 'bohr' 
        self.cell_parameters.vectors   = s.axes.copy()

        self.k_points.clear()
        nkpoints = len(s.kpoints)
        if nkpoints>0:
            if s.at_Gpoint():
                self.k_points.specifier = 'gamma'
            elif s.at_Lpoint():
                self.k_points.update(
                    specifier = 'automatic',
                    grid  = (1,1,1),
                    shift = (1,1,1)
                    )
            else:
                kpoints = s.kpoints/(2*pi)
                self.k_points.specifier = 'tpiba'
                self.k_points.nkpoints  = nkpoints
                self.k_points.kpoints   = kpoints
                self.k_points.weights   = s.kweights.copy()
                self.k_points.change_specifier('crystal',self) #added to make debugging easier
            #end if
        #end if

        if 'masses' not in self.atomic_species:
            self.atomic_species.masses = obj()
        #end if
        for name in system.ion_labels:
            is_elem, element = Elements.is_element(name, return_element=True)
            self.atomic_species.masses[name] = element.atomic_weight
        #end for
        if elem_order is None:
            self.atomic_species.atoms = list(sorted(system.ion_labels))
        else:
            if set(elem_order)!=set(system.ion_labels):
                self.error('elem_order is missing some atomic species\natomic species present: {0}\nelem_order: {1}'.format(sorted(system.ion_labels),elem_order))
            elif len(elem_order)!=system.n_ions:
                self.error('elem_order has repeated elements\nelem_order: {0}'.format(elem_order))
            #end if
            self.atomic_species.atoms = list(elem_order)
        #end if

        # set pseudopotentials for renamed atoms (e.g. Cu3 is same as Cu)
        pp = self.atomic_species.pseudopotentials
        for atom in self.atomic_species.atoms:
            if atom not in pp:
                iselem,element = Elements.is_element(atom,return_element=True)
                symbol = element.symbol
                if iselem and symbol in pp:
                    pp[atom] = str(pp[symbol])
                #end if
            #end if
        #end for

        self.atomic_positions.specifier = 'bohr'
        self.atomic_positions.positions = s.pos.copy()
        self.atomic_positions.atoms     = list(s.elem)
        if s.frozen is not None:
            frozen = s.frozen
            if 'relax_directions' in self.atomic_positions:
                relax_directions = self.atomic_positions.relax_directions
            else:
                relax_directions = np.ones(s.pos.shape,dtype=int)
            #end if
            for i in range(len(s.pos)):
                relax_directions[i,0] = int(not frozen[i,0] and relax_directions[i,0])
                relax_directions[i,1] = int(not frozen[i,1] and relax_directions[i,1])
                relax_directions[i,2] = int(not frozen[i,2] and relax_directions[i,2])
            #end for
            self.atomic_positions.relax_directions = relax_directions
        #end if                    
    #end def incorporate_system


    def incorporate_system_old(self,system,spin_polarized=None):
        system.check_folded_system()
        system.change_units('B')
        s  = system.structure
        nc = system.net_charge
        ns = system.net_spin

        nup = system.n_up
        ndn = system.n_down

        self.system.ibrav        = 0
#        self.system['celldm(1)'] = 1.0e0
        nions = system.n_ions
        nspecies = system.n_species
        self.system.nat          = nions
        self.system.ntyp         = nspecies
        #self.system.nelec        = nup+ndn
        self.system.tot_charge   = nc
        mag = nup-ndn
        if (mag!=0 and spin_polarized!=False) or spin_polarized==True:
            self.system.nspin = 2
            self.system.tot_magnetization = mag
        #end if

        if 'cell_parameters' not in self:
            self.cell_parameters = self.element_types['cell_parameters']()
        #end if
        self.cell_parameters.specifier = 'bohr'
        self.cell_parameters.vectors   = s.axes.copy()

        self.k_points.clear()
        nkpoints = len(s.kpoints)
        if nkpoints>0:
            if s.at_Gpoint():
                self.k_points.specifier = 'gamma'
            elif s.at_Lpoint():
                self.k_points.update(
                    specifier = 'automatic',
                    grid  = (1,1,1),
                    shift = (1,1,1)
                    )
            else:
                kpoints = s.kpoints/(2*pi)
                self.k_points.specifier = 'tpiba'
                self.k_points.nkpoints  = nkpoints
                self.k_points.kpoints   = kpoints
                self.k_points.weights   = s.kweights.copy()
                self.k_points.change_specifier('crystal',self) #added to make debugging easier

            #end if
        #end if

        masses = obj()
        for name in system.ion_labels:
            is_elem, element = Elements.is_element(name, return_element=True)
            masses[name] = element.atomic_weight
        #end for
        self.atomic_species.atoms  = list(sorted(system.ion_labels))
        self.atomic_species.masses = masses
        # set pseudopotentials for renamed atoms (e.g. Cu3 is same as Cu)
        pp = self.atomic_species.pseudopotentials
        for atom in self.atomic_species.atoms:
            if atom not in pp:
                iselem,element = Elements.is_element(atom,return_element=True)
                symbol = element.symbol
                if iselem and symbol in pp:
                    pp[atom] = str(pp[symbol])
                #end if
            #end if
        #end for

        self.atomic_positions.specifier = 'bohr'
        self.atomic_positions.positions = s.pos.copy()
        self.atomic_positions.atoms     = list(s.elem)
        if s.frozen is not None:
            frozen = s.frozen
            if 'relax_directions' in self.atomic_positions:
                relax_directions = self.atomic_positions.relax_directions
            else:
                relax_directions = np.ones(s.pos.shape,dtype=int)
            #end if
            for i in range(len(s.pos)):
                relax_directions[i,0] = int(not frozen[i,0] and relax_directions[i,0])
                relax_directions[i,1] = int(not frozen[i,1] and relax_directions[i,1])
                relax_directions[i,2] = int(not frozen[i,2] and relax_directions[i,2])
            #end for
            self.atomic_positions.relax_directions = relax_directions
        #end if                    
    #end def incorporate_system_old

        
    # test needed
    def return_system(self,*,structure_only=False,**valency):
        ibrav = self.system.ibrav
        if ibrav!=0:
            self.error('ability to handle non-zero ibrav not yet implemented')
        #end if

        scale,axes,kaxes = self.get_common_vars('scale','axes','kaxes')

        elem = list(self.atomic_positions.atoms)
        ap = deepcopy(self.atomic_positions)
        ap.change_specifier('bohr',self)
        pos = ap.positions

        kp = deepcopy(self.k_points)
        kp.change_specifier('tpiba',self)
        kpoints = kp.kpoints*(2*pi)/scale

        center = axes.sum(0)/2
        structure = Structure(
            axes    = axes,
            elem    = elem,
            scale   = scale,
            pos     = pos,
            center  = center,
            kpoints = kpoints,
            units   = 'B',
            rescale = False
            )
        structure.zero_corner()
        structure.recenter()

        if structure_only:
            return structure
        #end if
  
        ion_charge = 0
        atoms   = list(self.atomic_positions.atoms)
        for atom in self.atomic_species.atoms:
            if atom not in valency:
                self.error('valence charge for atom {0} has not been defined\nplease provide the valence charge as an argument to return_system()'.format(atom))
            #end if
            ion_charge += atoms.count(atom)*valency[atom]
        #end for

        if 'nelup' in self.system:
            nup = self.system.nelup
            ndn = self.system.neldw
            net_charge = ion_charge - nup - ndn
            net_spin   = nup - ndn
        elif 'tot_magnetization' in self.system:
            net_spin = self.system.tot_magnetization
            if 'nelec' in self.system:
                net_charge = ion_charge - self.system.nelec
            else:
                net_charge = 0
            #end if
        else:
            net_spin = 0
            if 'nelec' in self.system:
                net_charge = ion_charge - self.system.nelec
            else:
                net_charge = 0
            #end if
        #end if

        system = PhysicalSystem(structure,net_charge,net_spin,**valency)

        return system
    #end def return_system


    def standardize_types(self):
        for s in self.values():
            if isinstance(s,Section):
                array_keys = []
                for k in s.keys():
                    if '(' in k:
                        array_keys.append(k)
                    #end if
                #end for
                for k in array_keys:
                    name  = k.split('(',1)[0]
                    index = k.replace(name,'').strip('()')
                    if ',' not in index:
                        index = int(index)
                    else:
                        index = tuple(np.array(index.split(','),dtype=int))
                    #end if
                    if name not in s:
                        s[name] = obj()
                    #end if
                    s[name][index] = s[k]
                    del s[k]
                #end for
                for k in PwscfInputBase.real_arrays:
                    if k in s and not isinstance(s[k],obj):
                        s[k] = obj(s[k])
                    #end if
                #end for
            #end if
        #end for
    #end def standardize_types

#end class PwscfInput



def generate_pwscf_input(selector,**kwargs):
    if selector=='generic':
        return generate_any_pwscf_input(**kwargs)
    if selector=='scf':
        return generate_scf_input(**kwargs)
    elif selector=='nscf':
        return generate_nscf_input(**kwargs)
    elif selector=='relax':
        return generate_relax_input(**kwargs)
    elif selector=='vc-relax':
        return generate_vcrelax_input(**kwargs)
    else:
        error('selection '+str(selector)+' has not been implemented for pwscf input generation')
    #end if
#end def generate_pwscf_input


generate_any_defaults = obj(
    standard = obj(
        prefix     = 'pwscf',
        outdir     = 'pwscf_output',
        pseudo_dir = './',
        kgrid      = None,
        kshift     = (0,0,0),
        use_folded = True,
        ),
    oldscf   = obj(
        calculation      = 'scf',
        prefix           = 'pwscf',
        outdir           = 'pwscf_output',
        pseudo_dir       = './',
        ecutwfc          = 200.,
        occupations      = 'smearing',
        smearing         = 'fermi-dirac',
        degauss          = 0.0001,
        nosym            = False,
        wf_collect       = True,
        restart_mode     = 'from_scratch',
        tstress          = True,
        tprnfor          = True,
        disk_io          = 'low',
        verbosity        = 'high',
        ibrav            = 0,
        conv_thr         = 1e-10,
        electron_maxstep = 1000,
        mixing_mode      = 'plain',
        mixing_beta      = .7,
        diagonalization  = 'david',
        kgrid            = None,
        kshift           = None,
        use_folded       = True,
        ),
    oldrelax = obj(
        calculation       = 'relax',
        prefix            = 'pwscf',
        outdir            = 'pwscf_output',
        pseudo_dir        = './',
        ecutwfc           = 50.,
        occupations       = 'smearing',
        smearing          = 'fermi-dirac',
        degauss           = 0.0001,
        nosym             = True,
        wf_collect        = False,
        restart_mode      = 'from_scratch',
        disk_io           = 'low',
        verbosity         = 'high',
        ibrav             = 0,
        conv_thr          = 1e-6,
        electron_maxstep  = 1000,
        mixing_mode       = 'plain',
        mixing_beta       = .7,
        diagonalization   = 'david',
        ion_dynamics      = 'bfgs',
        upscale           = 100,
        pot_extrapolation = 'second_order',
        wfc_extrapolation = 'second_order',
        kgrid             = None,
        kshift            = None,
        use_folded        = False,
        ),
    )

def generate_any_pwscf_input(**kwargs):
    #move values into a more convenient representation
    #kwargs = obj(**kwargs)
    # enforce lowercase internally, but remain case insensitive to user input
    kw_in = kwargs
    kwargs = obj()
    for name,value in kw_in.items():
        kwargs[name.lower()] = value
    #end for

    # setup for k-point symmetry run
    #   ecutwfc is set to 1 so that pwscf will crash after initialization
    #   symmetrized k-points will still be written to log output
    ksymm_run = kwargs.pop('ksymm_run',False)
    if ksymm_run:
        kwargs.ecutwfc    = 1
        kwargs.nosym      = False
        kwargs.use_folded = False
    #end if

    #assign default values
    if 'defaults' in kwargs:
        defaults = kwargs.defaults
        del kwargs.defaults
        if isinstance(defaults,str):
            if defaults in generate_any_defaults:
                defaults = generate_any_defaults[defaults]
            else:
                error('invalid default set requested: {0}\n  valid options are {1}'.format(defaults,sorted(generate_any_defaults.keys())))
            #end if
        #end if
    else:
        defaults = generate_any_defaults.standard
    #end if
    for name,default in defaults.items():
        if name not in kwargs:
            deftype = type(default)
            if inspect.isclass(default) or inspect.isfunction(default):
                kwargs[name] = default()
            else:
                kwargs[name] = default
            #end if
        #end if
    #end for

    if ksymm_run and 'calculation' in kwargs and kwargs.calculation!='scf':
        error('input parameter "calculation" must be set to "scf" when ksymm_run is requested')
    #end if

    #copy certain keywords
    tot_magnetization = kwargs.get('tot_magnetization',None)
    nspin             = kwargs.get('nspin',None)
    nbnd              = kwargs.get('nbnd',None)
    hubbard_u         = kwargs.get('hubbard_u',None)
    # Pre 7.2 Hubbard tags
    hub_keys_pre72 = 'hubbard_u hubbard_j0 hubbard_j U_projection_type'.lower().split()
    has_pre72_keys = any(([_ in kwargs.keys() for _ in hub_keys_pre72]))
    # QE >=7.2 Hubbard tags
    hub_keys_v72 = 'hubbard hubbard_proj'.lower().split()
    has_v72_keys = any(([_ in kwargs.keys() for _ in hub_keys_v72]))
    if has_pre72_keys + has_v72_keys > 1:
        error('Please use {} for QE version <7.2 and {} for QE version >=7.2'.format(hub_keys_pre72, hub_keys_v72))
    #end if     
    occ               = kwargs.get('occupations',None)
    
    #make an empty input file
    pw = PwscfInput()

    #pull out section keywords
    section_keywords = obj()
    for section_name,section_type in PwscfInput.section_types.items():
        keys = set(kwargs.keys()) & section_type.variables
        if len(keys)>0:
            kw = obj()
            for k in keys:
                if k in kwargs:
                    kw[k] = kwargs.pop(k)
            section = section_type()
            section.assign(**kw)
            pw[section_name] = section
        #end if
    #end for

    #process other keywords
    use_folded       = kwargs.pop('use_folded')
    kgrid            = kwargs.pop('kgrid')
    kshift           = kwargs.pop('kshift')
    system           = kwargs.pop('system',None)
    pseudos          = kwargs.pop('pseudos',[])
    elem_order       = kwargs.pop('elem_order',None)
    mass             = kwargs.pop('mass',None)
    elem             = kwargs.pop('elem',None)
    pos              = kwargs.pop('pos',None)
    totmag_sys       = kwargs.pop('totmag_sys',False)
    start_mag        = kwargs.pop('start_mag',None)
    bandfac          = kwargs.pop('bandfac',None)
    nogamma          = kwargs.pop('nogamma',False)
    positions_option = kwargs.pop('pos_specifier',None)
    if positions_option is None:
        positions_option = kwargs.pop('positions_option',None)
    #end if
    if positions_option is None:
        positions_option = kwargs.pop('atomic_positions_option',None)
    #end if
    kpoints_option = kwargs.pop('kpoints_option',None)
    if kpoints_option is None:
        kpoints_option = kwargs.pop('k_points_option',None)
    #end if
    cell_option = kwargs.pop('cell_option',None)
    if cell_option is None:
        cell_option = kwargs.pop('cell_parameters_option',None)
    #end if
    hubbard_input  = kwargs.pop('hubbard', None)
    hubbard_option = kwargs.pop('hubbard_proj',None)
    
    #  pseudopotentials
    pseudopotentials = obj()
    atom_species = []
    if system is not None:
        pseudos = PseudoSet.pseudo_remap('pwscf',pseudos,system)
    for ppname in pseudos:
        #element = ppname[0:2].strip('.')
        label,element = pp_elem_label(ppname,guard=True)
        atom_species.append(element)
        pseudopotentials[element] = ppname
    #end for
    pw.atomic_species.update(
        atoms            = list(sorted(atom_species)),
        pseudopotentials = pseudopotentials,
        )

    #  physical system information
    if system is None:
        if elem is None:
            error('system must be provided','generate_pwscf_input')
        else:
            if mass is None:
                error('"mass" must be provided when "elem" is given','generate_pwscf_input')
            #end if
            if pos is None:
                error('"pos" must be provided when "elem" is given','generate_pwscf_input')
            #end if
            if positions_option is None:
                error('"atomic_positions_option" must be provided when "elem" is given','generate_pwscf_input')
            #end if

            # fill in atomic_species
            species = set(elem)
            if elem_order is not None:
                if set(elem_order)!=species:
                    error('elem_order is missing some atomic species\natomic species present: {0}\nelem_order: {1}'.format(sorted(species),elem_order),'generate_pwscf_input')
                elif len(elem_order)!=len(species):
                    error('elem_order has repeated elements\nelem_order: {0}'.format(elem_order),'generate_pwscf_input')
                #end if
                pw.atomic_species.atoms = list(elem_order)
            else:
                pw.atomic_species.atoms = list(sorted(species))
            #end if
            pw.atomic_species.masses = obj(mass)
            pp = pw.atomic_species.pseudopotentials
            for atom in pw.atomic_species.atoms:
                if atom not in pp:
                    iselem, element = Elements.is_element(atom, return_element=True)
                    if iselem and element.symbol in pp:
                        pp[atom] = str(pp[element.symbol])
                    #end if
                #end if
            #end for

            # fill in atomic_positions
            pw.atomic_positions.specifier = positions_option
            pw.atomic_positions.positions = np.array(pos,dtype=float)
            pw.atomic_positions.atoms     = list(elem)
        #end if
    else:
        system.check_folded_system()
        system.change_units('B')
        s = system.structure
        #setting the 'lattice' (cell axes) requires some delicate care
        #  qmcpack will fail if this is even 1e-10 off of what is in 
        #  the wavefunction hdf5 file from pwscf
        if s.folded_structure is not None:
            fs = s.folded_structure
            axes = np.array(array_to_string(fs.axes).split(),dtype=float)
            npe.reshape_inplace(axes, fs.axes.shape)
            axes = np.dot(s.tmatrix,axes)
            if abs(axes-s.axes).sum()>1e-5:
                error('supercell axes do not match tiled version of folded cell axes\nyou may have changed one set of axes (super/folded) and not the other\nfolded cell axes:\n'+str(fs.axes)+'\nsupercell axes:\n'+str(s.axes)+'\nfolded axes tiled:\n'+str(axes),'generate_pwscf_input')
            #end if
        else:
            axes = np.array(array_to_string(s.axes).split(),dtype=float)
            npe.reshape_inplace(axes, s.axes.shape)
        #end if
        s.adjust_axes(axes)
        if use_folded:
            system = system.get_smallest()
        #end if
        pw.incorporate_system(system,elem_order)
    #end if

    #  tot_magnetization from system
    if totmag_sys and 'tot_magnetization' not in pw.system:
        tot_magnetization = system.net_spin
        pw.system.tot_magnetization = tot_magnetization
    #end if

    # set the number of spins (totmag=0 still needs nspin=2)
    if nspin is not None:
        pw.system.nspin = nspin
    elif start_mag is None and tot_magnetization is None:
        pw.system.nspin = 1
    else:
        pw.system.nspin = 2
    #end if

    # set nbnd using bandfac, if provided
    if nbnd is None and bandfac is not None:
        nocc = max(system.n_up, system.n_down)
        pw.system.nbnd = int(np.ceil(nocc*bandfac))
    #end if

    #  Hubbard U
    if hubbard_u is not None:
        if not isinstance(hubbard_u,(dict,obj)):
            error('input hubbard_u must be of type dict or obj','generate_pwscf_input')
        #end if
        pw.system.hubbard_u = deepcopy(hubbard_u)
        pw.system.lda_plus_u = True
    #end if

    #  starting magnetization
    if start_mag is not None:
        if not isinstance(start_mag,(dict,obj)):
            error('input start_mag must be of type dict or obj','generate_pwscf_input')
        #end if
        pw.system.starting_magnetization = deepcopy(start_mag)
    #end if

    #  kpoints
    if kshift is None:
        zero_shift = False
    else:
        zero_shift = tuple(kshift)==(0,0,0)
    #end if

    single_point = kgrid is not None and tuple(kgrid)==(1,1,1)
    no_info      = kgrid is None and system is None
    at_gamma  = zero_shift and (single_point or no_info)
    sys_gamma = 'specifier' in pw.k_points and pw.k_points.specifier=='gamma'
    auto      = kgrid is not None and kshift is not None
    shifted   = not zero_shift and kshift is not None

    if at_gamma and not nogamma:
        pw.k_points.clear()
        pw.k_points.specifier = 'gamma'
    elif auto:
        pw.k_points.clear()
        pw.k_points.update(
            specifier = 'automatic',
            grid      = kgrid,
            shift     = kshift
            )
    elif sys_gamma and not nogamma:
        pw.k_points.clear()
        pw.k_points.specifier = 'gamma'
    elif (at_gamma or sys_gamma) and nogamma:
        pw.k_points.clear()
        pw.k_points.update(
            specifier = 'automatic',
            grid      = (1,1,1),
            shift     = (0,0,0)
            )
    #elif shifted:
    #    pw.k_points.clear()
    #    pw.k_points.set(
    #        specifier = 'automatic',
    #        grid      = (1,1,1),
    #        shift     = kshift
    #        )
    #end if

    # occupations card
    if isinstance(occ,(list,np.ndarray)):
        pw.system.occupations = 'from_input'
        occ_card = occupations()
        occ_card.occupations = np.array(occ,dtype=float)
        pw.occupations = occ_card
    #end if
    
    # hubbard card
    if hubbard_input is not None:
        hubbard_card = hubbard()
        hubbard_card.system  = system
        hubbard_card.hubbard = hubbard_input
        pw.hubbard = hubbard_card
        if hubbard_option is None:
            hubbard_option = hubbard_card.default_specifier
        else:
            if hubbard_option not in hubbard_card.available_specifiers:
                error('HUBBARD card specifier "{}" is not valid. Available specifiers: {}'.format(hubbard_option, hubbard_card.available_specifiers))
            #end if
        #end if
        pw.hubbard.specifier = hubbard_option
    #end if 

    # adjust card options, if requested
    options = obj(
        atomic_positions = positions_option,
        k_points         = kpoints_option,
        cell_parameters  = cell_option,
        )
    for card_name,option in options.items():
        if option is not None:
            if card_name not in pw:
                error('Card option provided for card "{}" but card is not present\noption provided: {}'.format(card_name,option))
            #end if
            pw[card_name].change_option(option,pw)
        #end if
    #end for

    # check for misformatted kpoints
    if len(pw.k_points)==0:
        error('k_points section has not been filled in\nplease provide k-point information in either of\n  1) the kgrid input argument\n  2) in the PhysicalSystem object (system input argument)','generate_pwscf_input')
    #end if

    # check for leftover keywords
    if len(kwargs)>0:
        error('unrecognized keywords: {0}\nthese keywords are not known to belong to any namelist for PWSCF'.format(sorted(kwargs.keys())),'generate_pwscf_input')
    #end if  
    
    return pw
#end def generate_any_pwscf_input



def generate_scf_input(*,
                       prefix       = 'pwscf',
                       outdir       = 'pwscf_output',
                       input_dft    = None,
                       exx_fraction = None,
                       exxdiv_treatment = None,
                       ecut         = 200.,
                       ecutrho      = None,
                       ecutfock     = None,
                       occupations  = 'smearing',
                       smearing     = 'fermi-dirac',
                       degauss      = 0.0001,
                       nosym        = False,
                       spin_polarized = None,
                       assume_isolated = None,
                       wf_collect   = True,
                       hubbard_u    = None,
                       start_mag    = None,
                       restart_mode = 'from_scratch',
                       tstress      = True,
                       tprnfor      = True,
                       disk_io      = 'low',
                       verbosity    = 'high',
                       ibrav        = 0,
                       conv_thr     = 1e-10,
                       electron_maxstep = 1000,
                       mixing_mode  = 'plain',
                       mixing_beta  = .7,
                       diagonalization = 'david',
                       kgrid        = None,
                       kshift       = None,
                       pseudos      = None,
                       system       = None,
                       use_folded   = True,
                       group_atoms  = False,
                       la2F         = None,
                       nbnd         = None,
                       lspinorb     = False,
                       noncolin     = False,
                       ):
    if pseudos is None:
        pseudos = []
    #end if
    if system is not None:
        pseudos = PseudoSet.pseudo_remap('pwscf',pseudos,system)
    #end if
    pseudopotentials = obj()
    atoms = []
    for ppname in pseudos:
        element = ppname[0:2].strip('.')
        atoms.append(element)
        pseudopotentials[element] = ppname
    #end for

    if ecutrho is None:
        ecutrho = 4*ecut
    #end if

    pw = PwscfInput()
    pw.control.update(
        calculation  = 'scf',
        prefix       = prefix,
        restart_mode = restart_mode,
        tstress      = tstress,
        tprnfor      = tprnfor,
        pseudo_dir   = './',
        outdir       = outdir,
        disk_io      = disk_io,
        verbosity    = verbosity,
        wf_collect   = wf_collect
        )
    pw.system.update(
        ibrav       = ibrav,
        ecutwfc     = ecut,
        ecutrho     = ecutrho,
        nosym       = nosym,
        )
    if la2F is not None:
        pw.system.la2F = la2F
    #end if
    if nbnd is not None:
        pw.system.nbnd = nbnd
    # end if
    if assume_isolated is not None:
        pw.system.assume_isolated = assume_isolated
    #end if
    if occupations is not None:
        if occupations=='smearing':
            pw.system.update(
                occupations = occupations,
                smearing    = smearing,
                degauss     = degauss,
                )
        else:
            pw.system.occupations = occupations
        #end if
    #end if
    pw.electrons.update(
        electron_maxstep = electron_maxstep,
        conv_thr    = conv_thr,
        mixing_mode = mixing_mode,
        mixing_beta = mixing_beta,
        diagonalization = diagonalization,
        )
    pw.atomic_species.update(
        atoms            = atoms,
        pseudopotentials = pseudopotentials
        )

    if noncolin or lspinorb:
        pw.system.update(
            noncolin = noncolin or lspinorb,
            lspinorb = lspinorb
            )
    #end if
    if input_dft is not None:
        pw.system.input_dft = input_dft
    #end if
    if exx_fraction is not None:
        pw.system.exx_fraction = exx_fraction
    #end if
    if exxdiv_treatment is not None:
        pw.system.exxdiv_treatment = exxdiv_treatment
    #end if
    if ecutfock is not None:
        pw.system.ecutfock = ecutfock
    #end if

    system.check_folded_system()
    system.change_units('B')
    s = system.structure
    #setting the 'lattice' (cell axes) requires some delicate care
    #  qmcpack will fail if this is even 1e-10 off of what is in 
    #  the wavefunction hdf5 file from pwscf
    if s.folded_structure is not None:
        fs = s.folded_structure
        axes = np.array(array_to_string(fs.axes).split(),dtype=float)
        npe.reshape_inplace(axes, fs.axes.shape)
        axes = np.dot(s.tmatrix,axes)
        if abs(axes-s.axes).sum()>1e-5:
            error('supercell axes do not match tiled version of folded cell axes\n  you may have changed one set of axes (super/folded) and not the other\n  folded cell axes:\n'+str(fs.axes)+'\n  supercell axes:\n'+str(s.axes)+'\n  folded axes tiled:\n'+str(axes))
        #end if
    else:
        axes = np.array(array_to_string(s.axes).split(),dtype=float)
        npe.reshape_inplace(axes, s.axes.shape)
    #end if
    s.adjust_axes(axes)

    if use_folded:
        system = system.get_smallest()
    #end if
        
    if start_mag is not None:
        spin_polarized=True
    #end if

    if system is not None:
        pw.incorporate_system_old(system,spin_polarized=spin_polarized)
    #end if

    if hubbard_u is not None:
        if not isinstance(hubbard_u,(dict,obj)):
            error('input hubbard_u must be of type dict or obj')
        #end if
        pw.system.hubbard_u = deepcopy(hubbard_u)
        pw.system.lda_plus_u = True
    #end if
    if start_mag is not None:
        if not isinstance(start_mag,(dict,obj)):
            error('input start_mag must be of type dict or obj')
        #end if
        pw.system.starting_magnetization = deepcopy(start_mag)
        #if 'tot_magnetization' in pw.system:
        #    del pw.system.tot_magnetization
        ##end if
        if 'multiplicity' in pw.system:
            del pw.system.multiplicity
        #end if
    #end if

    if kshift is None:
        kshift = (1,1,1)
    #end if

    if system is not None:
        structure = system.structure
        if group_atoms:
            warn('requested grouping by atomic species, but pwscf does not group atoms anymore!')
        #end if
        #if group_atoms:  # disabled, hopefully not needed for qmcpack
        #    structure.group_atoms()
        ##end if
        if structure.at_Gpoint():
            pw.k_points.clear()
            pw.k_points.update(
                specifier = 'automatic',
                grid  = (1,1,1),
                shift = (0,0,0)
                )
        elif structure.at_Lpoint():
            pw.k_points.clear()
            pw.k_points.update(
                specifier = 'automatic',
                grid  = (1,1,1),
                shift = (1,1,1)
                )
        #end if
    #end if

    if kgrid is not None:
        pw.k_points.clear()
        pw.k_points.update(
            specifier = 'automatic',
            grid     = kgrid,
            shift    = kshift
            )
    elif system is None:
        pw.k_points.clear()
        pw.k_points.update(
            specifier = 'automatic',
            grid     = (1,1,1),
            shift    = kshift
            )
    #end if


    return pw
#end def generate_scf_input




def generate_nscf_input(**kwargs):
    pw = generate_scf_input(**kwargs)
    pw.control.update(
        calculation = 'nscf'
        )
    return pw
#end def generate_nscf_input




def generate_relax_input(*,
                         prefix       = 'pwscf',
                         outdir       = 'pwscf_output',
                         input_dft    = None,
                         exx_fraction = None,
                         ecut         = 50.,
                         ecutrho      = None,
                         ecutfock     = None,
                         conv_thr     = 1e-6,
                         mixing_mode  = 'plain',
                         mixing_beta  = .7,
                         diagonalization = 'david',
                         occupations  = 'smearing',
                         smearing     = 'fermi-dirac',
                         degauss      = 0.0001,
                         nosym        = True,
                         spin_polarized = None,
                         assume_isolated = None,
                         upscale      = 100,
                         pot_extrapolation = 'second_order',
                         wfc_extrapolation = 'second_order',
                         hubbard_u    = None,
                         start_mag    = None,
                         restart_mode = 'from_scratch',
                         kgrid        = None,
                         kshift       = None,
                         pseudos      = None,
                         system       = None,
                         use_folded   = False,
                         group_atoms  = False,
                         forc_conv_thr= None,
                         disk_io      = 'low',
                         wf_collect   = False,
                         verbosity    = 'high',
                         ):
    if pseudos is None:
        pseudos = []
    #end if
    if system is not None:
        pseudos = PseudoSet.pseudo_remap('pwscf',pseudos,system)
    #end if
    
    pseudopotentials = obj()
    atoms = []
    for ppname in pseudos:
        element = ppname[0:2].strip('.')
        atoms.append(element)
        pseudopotentials[element] = ppname
    #end for

    if ecutrho is None:
        ecutrho = 4*ecut
    #end if

    pw = PwscfInput('ions')
    pw.control.update(
        calculation  = 'relax',
        prefix       = prefix,
        restart_mode = 'from_scratch',
        #tstress      = True,
        #tprnfor      = True,
        pseudo_dir   = './',
        outdir       = outdir,
        disk_io      = disk_io,
        verbosity    = verbosity,
        wf_collect   = wf_collect
        )
    pw.system.update(
        ibrav       = 0,
        ecutwfc     = ecut,
        ecutrho     = ecutrho,
        nosym       = nosym
        )
    if assume_isolated is not None:
        pw.system.assume_isolated = assume_isolated
    #end if
    if occupations is not None:
        if occupations=='smearing':
            pw.system.update(
                occupations = occupations,
                smearing    = smearing,
                degauss     = degauss,
                )
        #end if
    #end if
    pw.electrons.update(
        electron_maxstep = 1000,
        conv_thr         = conv_thr,
        mixing_beta      = mixing_beta,
        mixing_mode      = mixing_mode,
        diagonalization  = diagonalization
        )
    pw.atomic_species.update(
        atoms            = atoms,
        pseudopotentials = pseudopotentials
        )
    pw.ions.update(
        ion_dynamics      = 'bfgs',
        upscale           = upscale,
        pot_extrapolation = pot_extrapolation,
        wfc_extrapolation = wfc_extrapolation
        )

    if input_dft is not None:
        pw.system.input_dft = input_dft
    #end if
    if exx_fraction is not None:
        pw.system.exx_fraction = exx_fraction
    #end if
    if ecutfock is not None:
        pw.system.ecutfock = ecutfock
    #end if
    if system is not None:
        if use_folded:
            system = system.get_smallest()
        #end if
        pw.incorporate_system_old(system,spin_polarized=spin_polarized)
    #end if

    if hubbard_u is not None:
        if not isinstance(hubbard_u,(dict,obj)):
            error('input hubbard_u must be of type dict or obj')
        #end if
        pw.system.hubbard_u = deepcopy(hubbard_u)
        pw.system.lda_plus_u = True
    #end if
    if start_mag is not None:
        if not isinstance(start_mag,(dict,obj)):
            error('input start_mag must be of type dict or obj')
        #end if
        pw.system.starting_magnetization = deepcopy(start_mag)
        #if 'tot_magnetization' in pw.system:
        #    del pw.system.tot_magnetization
        ##end if
        if 'multiplicity' in pw.system:
            del pw.system.multiplicity
        #end if
    #end if

    if kshift is None:
        kshift = (1,1,1)
    #end if


    if system is not None:
        structure = system.structure
        if group_atoms:
            warn('requested grouping by atomic species, but pwscf does not group atoms anymore!','generate_relax_input')
        #end if
        #if group_atoms: #don't group atoms for pwscf, any downstream consequences?
        #    structure.group_atoms()
        ##end if
        if structure.at_Gpoint():
            pw.k_points.clear()
            pw.k_points.update(
                specifier = 'automatic',
                grid  = (1,1,1),
                shift = (0,0,0)
                )
        elif structure.at_Lpoint():
            pw.k_points.clear()
            pw.k_points.update(
                specifier = 'automatic',
                grid  = (1,1,1),
                shift = (1,1,1)
                )
        #end if
    #end if

    if kgrid is not None:
        pw.k_points.clear()
        pw.k_points.update(
            specifier = 'automatic',
            grid     = kgrid,
            shift    = kshift
            )
    elif system is None:
        pw.k_points.clear()
        pw.k_points.update(
            specifier = 'automatic',
            grid     = (1,1,1),
            shift    = kshift
            )
    #end if

    if forc_conv_thr is not None:
        pw.control.forc_conv_thr = forc_conv_thr
    # end if

    return pw
#end def generate_relax_input


def generate_vcrelax_input(
    press          = None, # None = use pw.x default
    cell_factor    = None, 
    cell_dofree    = None,
    forc_conv_thr  = None,
    ion_dynamics   = None,
    press_conv_thr = None,
    **kwargs):

    pw = generate_scf_input(**kwargs)
    pw.control.update(
        calculation = 'vc-relax'
        )
    pw['ions'] = pw.element_types['ions']()
    pw['cell'] = pw.element_types['cell'](
        press       = press,
        )

    # expand this section if you need more control over the input
    if forc_conv_thr is not None:
        pw.control.forc_conv_thr = forc_conv_thr
    # end if
    if cell_factor is not None:
        pw.cell.update(cell_factor=cell_factor)
    # end if
    if ion_dynamics is not None:
        pw.ions.update(ion_dynamics=ion_dynamics)
    # end if
    if press_conv_thr is not None:
        pw.cell.update(press_conv_thr=press_conv_thr)
    # end if
    if cell_dofree is not None:
        pw.cell.update(cell_dofree=cell_dofree)
    # end if

    return pw
# end def


#def generate_nscf_input(prefix='pwscf',outdir='pwscf_output',ecut=200.,kpoints=None,weights=None,pseudos=None,system=None):
#    if pseudos is None:
#        pseudos = []
#    #end if
#    
#    pseudopotentials = obj()
#    atoms = []
#    for ppname in pseudos:
#        element = ppname[0:2]
#        atoms.append(element)
#        pseudopotentials[element] = ppname
#    #end for
#
#    pw = PwscfInput()
#    pw.control.update(
#        calculation  = 'nscf',
#        prefix       = prefix,
#        restart_mode = 'from_scratch',
#        tstress      = True,
#        tprnfor      = True,
#        pseudo_dir   = './',
#        outdir       = outdir,
#        disk_io      = 'low',
#        wf_collect   = True
#        )
#    pw.system.update(
#        ibrav       = 0,
#        degauss     = 0.001,
#        smearing    = 'mp',
#        occupations = 'smearing',
#        ecutwfc     = ecut,
#        ecutrho     = 4*ecut
#        )
#    pw.electrons.update(
#        conv_thr    = 1.e-10,
#        mixing_beta = 0.7
#        )
#    pw.atomic_species.update(
#        atoms            = atoms,
#        pseudopotentials = pseudopotentials
#        )
#
#    if system is not None:
#        pw.incorporate_system(system)
#    #end if
#
#    overwrite_kpoints = kpoints is not None or system==None
#    if kpoints==None:
#        kpoints = np.array([[0.,0,0]])
#    else:
#        kpoints = np.array(kpoints)
#    #end if
#    if weights==None:
#        weights = np.ones((len(kpoints),),dtype=float)
#    #end if
#    if overwrite_kpoints:
#        pw.k_points.clear()
#        pw.k_points.update(
#            specifier = 'tpiba',
#            kpoints   = kpoints,
#            weights   = weights
#            )
#    #end if
#
#    return pw
##end def generate_nscf_input
