##################################################################
##  (c) Copyright 2018-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  quantum_package_input.py                                          #
#    Supports read/write of PySCF input.                             #
#                                                                    #
#  Content summary:                                                  #
#    PyscfInput                                                      #
#      SimulationInput class for PySCF.                              #
#                                                                    #
#    generate_pyscf_input                                            #
#      User-facing function to create template based input.          #
#====================================================================#


from copy import deepcopy
import numpy as np
from .developer import obj
from .simulation import SimulationInputTemplateDev


def render_string(s_in,n):
    indent = n*' '
    s = ''
    lines = s_in.splitlines()
    s += 3*"'"+'\n'
    for line in lines:
        s += indent+line.strip()+'\n'
    #end for
    s += indent+3*"'"
    return s
#end def render_string


def render_array(a,n):
    indent = n*' '
    s = ''
    if len(a.shape)==1:
        s += 'array('+repr(list(a))+')'
    elif len(a.shape)==2:
        s += 'array([\n'
        for a1 in a:
            s += indent+repr(list(a1))+',\n'
        #end for
        s = s[:-2]+'])'
    else:
        msg = 'render_array only supports up to 2D arrays'
        raise NotImplementedError(msg)
    #end if
    return s
#end def render_array



class PyscfInput(SimulationInputTemplateDev):

    basic_types = bool,int,float,str,tuple,list,dict

    allowed_types = basic_types+(np.ndarray,)

    mole_order = (
        'dump_input',
        'parse_arg',
        'verbose',
        'output',
        'max_memory',
        'atom',
        'basis',
        'unit',
        'nucmod',
        'ecp',
        'charge',
        'spin',
        'symmetry',
        'symmetry_subgroup',
        'cart',
        'nelec',
        'nelectron',
        'multiplicity',
        'ms',
        )

    cell_order = [  # noqa: RUF012
        'dump_input',
        'parse_arg',
        'a',
        'mesh',
        'ke_cutoff',
        'precision',
        'nimgs',
        'ew_eta',
        'ew_cut',
        'pseudo',
        'basis',
        'h',
        'dimension',
        'rcut',
        'ecp',
        'low_dim_ft_type',
        'unit',
        'atom',
        'gs',
        'h',
        'drop_exponent',
        'nimgs',
        ]

    for k in mole_order:
        if k not in cell_order:
            cell_order.append(k)
        #end if
    #end for
    cell_order: tuple[str] = tuple(cell_order)
    mole_allowed = frozenset(mole_order)
    cell_allowed = frozenset(cell_order)
    

    def __init__(self,
                 template    = None,     # path to template input file
                 prefix      = None,     # $prefix var for file prefixes
                 custom      = None,     # obj w/ $ prefixed vars in template
                 system      = None,     # physical system object
                 units       = None,     # input units desired
                 *,
                 use_folded  = True,     # use folded system/primitive cell
                 mole        = None,     # obj w/ Mole variables
                 cell        = None,     # obj w/ Cell variables
                 sys_var     = None,     # local var name for Mole/Cell
                 mole_var    = 'mol',    # local var name for Mole in written input
                 cell_var    = 'cell',   # local var name for Cell in written input
                 save_qmc    = False,    # convert to QMCPACK format
                 checkpoint  = False,    # set $chkfile variable
                 mf_var      = 'mf',     # local var name for mf, used for convert
                 kpts_var    = 'kpts',   # local var name for kpts, used for convert
                 filepath    = None,     # alias for template
                 text        = None,     # full text of (and alternate to) template 
                 calculation = None,     # obj w/ Calculation variables
                 chkfile     = None,     # obj w/ Calculation variables
                 twist_num   = None,     # Twist index
                 python_exe  = 'python3', # Python executable
                 ):
        if filepath is None and template is not None:
            filepath = template
        elif calculation is not None:
            self.calculation = calculation
            text='''#! /usr/bin/env $python_exe

import numpy as np
from numpy import array
$pyscfimport


$system

$calculation
'''
        #end if
        SimulationInputTemplateDev.__init__(self,filepath,text)

        self.prefix     = prefix
        self.save_qmc   = save_qmc
        self.checkpoint = checkpoint
        self.addendum   = None     # used for save2qmcpack

        if custom is not None:
            self.assign(**custom)
        #end if

        if sys_var is not None:
            mole_var = sys_var
            cell_var = sys_var
        #end if

        sys_name    = None
        sys_inputs  = obj()
        sys_allowed = None
        sys_kpoints = None
        is_mole = mole is not None
        is_cell = cell is not None
        if is_mole and is_cell:
            msg = (
                'both mole and cell provided\n'
                'please provide only one of them\n'
                'system cannot be both molecule and periodic solid'
                )
            raise ValueError(msg)
        elif is_mole and not isinstance(mole,(dict,obj)):
            msg = (
                'mole input must be a dict or obj\n'
                'you provided input of type: {0}'.format(
                    mole.__class__.__name__
                    )
                )
            raise TypeError(msg)
        elif is_cell and not isinstance(cell,(dict,obj)):
            msg = (
                'cell input must be a dict or obj\n'
                'you provided input of type: {0}'.format(
                    cell.__class__.__name__
                    )
                )
            raise TypeError(msg)
        #end if
        extra = ''
        if filepath is not None:
            extra = '\ntemplate located at: {0}'.format(filepath)
        #end if
        if system is not None and 'system' not in self.keywords:
            msg = 'system input is provided, but $system is not present in template input'+extra
            raise ValueError(msg)
        #end if
        folded_structure = None
        tiled_structure  = None
        if system is not None and 'system' not in self.values:
            if 'system' not in self.keywords:
                msg = (
                    'cannot incorporate "system" input\n'
                    '$system is not present in template input'+extra
                    )
                raise ValueError(msg)
            #end if
            system = deepcopy(system) # make a local copy
            if use_folded and system.has_folded():
                tiled_structure = deepcopy(system.structure)
                system = system.folded_system
            #end if
            s = system.structure
            if units is not None:
                s.change_units(units)
            #end if
            is_solid    = s.has_axes()
            is_molecule = not is_solid
            if is_solid and is_mole:
                msg = 'mole input provided, but provided system is not a molecule (cell axes are present)'
                raise ValueError(msg)
            elif is_molecule and is_cell:
                msg = 'cell input provided, but provided system is a molecule (cell axes are not present)'
                raise ValueError(msg)
            #end if
            is_mole |= is_molecule
            is_cell |= is_solid
            sys_inputs.atom   = s.pos_to_str(units=s.units, with_elem=True)
            sys_inputs.unit   = s.units
            sys_inputs.charge = system.net_charge
            sys_inputs.spin   = system.net_spin
            if is_solid:
                folded_structure     = deepcopy(s)
                sys_inputs.dimension = len(s.axes)
                sys_inputs.a         = s.write_axes()
                if len(s.kpoints)>0:
                    skp = deepcopy(s)
                    skp.change_units('B')
                    sys_kpoints = skp.kpoints.copy()
                #end if
            #end if
        #end if

        if is_mole:
            sys_name    = 'mole'
            sys_var     = mole_var
            sys_allowed = PyscfInput.mole_allowed
            sys_order   = PyscfInput.mole_order
            if mole is not None:
                sys_inputs.update(**mole)
            #end if
        elif is_cell:
            sys_name    = 'cell'
            sys_var     = cell_var
            sys_allowed = PyscfInput.cell_allowed
            sys_order   = PyscfInput.cell_order
            if cell is not None:
                sys_inputs.update(**cell)
            #end if
        else:
            None # no action needed if not molecule or periodic solid
        #end if

        if calculation is not None and 'calculation' not in self.values:

            calc = deepcopy(calculation) # make a local copy
            d = dict(
                method       = 'RKS',
                df_fitting   = True,
                xc           = 'pbe',
                tol          = '1e-10',
                df_method    = 'GDF',
                exxdiv       = 'ewald',
                u_idx        = None,
                max_cycle    = None,
                level_shift  = None,
                chkfile      = None,
                u_val        = None,
                C_ao_lo      = 'minao',
                )
            for k,v in d.items():
                if k not in calc:
                    calc[k] = v
            if calc.u_val is not None:
                calc.u_val = np.array(calc.u_val)
            #end if
            if calc.u_idx is not None:
                calc.u_idx = np.array(calc.u_idx)
            #end if

            c = '\n### generated calculation text ###\n'
            if sys_name is not None:
                df_str = '.density_fit()' if calc.df_fitting else ''
                if sys_name == 'mole':
                    if calc.u_idx is None:
                        c += 'mf = scf.{}({}){}\n'.format(calc.method,sys_var,df_str)
                    else:
                        c += 'mf = dft.{}({},U_idx={},U_val={},C_ao_lo=\'{}\'){}\n'.format(calc.method,sys_var,render_array(calc.u_idx,1),render_array(calc.u_val,1),calc.C_ao_lo,df_str)    
                    #end if
                elif sys_name == 'cell':
                    c += 'mydf          = df.{}({})\n'.format(calc.df_method,sys_var,'kpts')
                    c += 'mydf.auxbasis = \'weigend\'\n'
                    c += 'dfpath = \'df_ints.h5\'\n'
                    c += 'mydf._cderi_to_save = dfpath\n'
                    c += 'mydf.build()\n\n'
                    if calc.u_idx is None:
                        c += 'mf = scf.{}({},{}){}\n'.format(calc.method,sys_var,'kpts',df_str)
                    else:
                        c += 'mf = dft.{}({},{},U_idx={},U_val={},C_ao_lo=\'{}\'){}\n'.format(calc.method,sys_var,'kpts',render_array(calc.u_idx,1),render_array(calc.u_val,1),calc.C_ao_lo,df_str)    
                    c += 'mf.exxdiv      = \'{}\'\n'.format(calc.exxdiv)
                #end if
                if calc.max_cycle is not None: 
                    c += 'mf.max_cycle={}\n'.format(calc.max_cycle)
                #end if
                if calc.level_shift is not None: 
                    c += 'mf.level_shift={}\n'.format(calc.level_shift)
                #end if
                if calc.chkfile is not None: 
                    c += 'mf.chkfile=\'{}\'\n'.format(calc.chkfile)
                #end if
            #end if
            if 'KS' in calc.method:
                c += 'mf.xc          = \'{}\'\n'.format(calc.xc)
            #end if
            c += 'mf.tol         = \'{}\'\n'.format(calc.tol)
            if calc.df_fitting and not is_mole:
                c += 'mf.with_df     = mydf\n'
            #end if
            c += 'e_scf = mf.kernel()\n'

            c += '### end generated calculation text ###\n\n'
            self.assign(calculation=c)
            cimp = '\n### generated pyscfimport text ###\n'
            if is_mole:
                cimp += 'from pyscf import df, scf, dft\n'
            elif is_cell:
                cimp += 'from pyscf.pbc import df, scf\n'
            #endif
            cimp += '### end generated pyscfimport text ###\n\n'
            self.assign(pyscfimport=cimp)
            self.assign(python_exe=python_exe)
        #end if

        if sys_name is not None:
            invalid = set(sys_inputs.keys())-sys_allowed
            if len(invalid)>0:
                msg = (
                    'invalid {0} inputs\n'
                    'invalid inputs: {1}\n'
                    'valid options are: {2}'.format(
                        sys_name, sorted(invalid), sorted(sys_allowed)
                        )
                    )
                raise ValueError(msg)
            #end if
            klen = 0
            has_array = False
            for k,v in sys_inputs.items():
                if not isinstance(v,PyscfInput.allowed_types):
                    tlist = ''
                    for t in PyscfInput.allowed_types:
                        tlist += t.__name__+','
                    #end for
                    tlist = tlist[:-1]
                    msg = (
                        '{0} input "{1}" has an invalid type\n'
                        'invalid type: {2}\n'
                        'allowed types are: {3}'.format(
                            sys_name, k, v.__class__.__name__, tlist
                            )
                        )
                    raise TypeError(msg)
                #end if
                klen = max(klen,len(k))
                has_array |= isinstance(v,np.ndarray)
            #end for
            has_array |= isinstance(sys_kpoints,np.ndarray)
            c = '\n### generated system text ###\n'
            if has_array:
                c += 'from numpy import array\n'
            #end if
            if is_mole:
                c += 'from pyscf import gto as gto_loc\n'
                c += '{0} = gto_loc.Mole()\n'.format(sys_var)
            elif is_cell:
                c += 'from pyscf.pbc import gto as gto_loc\n'
                c += '{0} = gto_loc.Cell()\n'.format(sys_var)
            #end if
            fmt = sys_var+'.{0:<'+str(klen)+'} = {1}\n'
            nalign = klen+len(sys_var)+4
            for k in sys_order:
                if k in sys_inputs:
                    v = sys_inputs[k]
                    if isinstance(v,str) and '\n' in v:
                        vs = render_string(v,nalign)
                    elif isinstance(v,np.ndarray):
                        if len(v.shape)>2:
                            msg = (
                                'cannot write system input variable {0}\n'
                                '{0} is an array with more than two dimensions\n'
                                'only two dimensions are currently supported for writing\n'
                                'array contents: {1}'.format(k,v)
                                )
                            raise ValueError(msg)
                        #end if
                        vs = render_array(v,nalign)
                    elif isinstance(v,PyscfInput.basic_types):
                        vs = repr(v)
                    else:
                        vs = None # should not be possible
                    #end if
                    c += fmt.format(k,vs)
                #end if
            #end for
            c += '{0}.build()\n'.format(sys_var)
            if sys_kpoints is not None:
                c += '{0} = {1}\n'.format(kpts_var,render_array(sys_kpoints,4))
            #end if
            c += '### end generated system text ###\n\n'
            self.assign(system=c)
        #end if

        if prefix is not None:
            self.allow_no_assign('prefix')
            self.assign(prefix=prefix)
        #end if

        if save_qmc:
            if prefix is None:
                msg = (
                    'cannot generate save2qmcpack text\n'
                    'please provide input variable "prefix"\n'
                    '(used to set "title" in save2qmcpack)'
                    )
                raise ValueError(msg)
            elif sys_var is None:
                msg = (
                    'cannot generate save2qmcpack text\n'
                    'please provide input variable "sys_var"\n'
                    '(used to set "cell" in save2qmcpack) '
                    )
                raise ValueError(msg)
            #end if
            if folded_structure is not None:
                folded_structure.change_units('B') # always use Bohr units for k-points
                nkpoints = len(folded_structure.kpoints)
                if sys_kpoints is None and nkpoints!=0 or sys_kpoints is not None and len(sys_kpoints)!=nkpoints:
                    msg = 'inconsistency in written and saved k-points'
                    raise ValueError(msg)
                #end if
                if tiled_structure is None:
                    tiling  = np.array([1,1,1],dtype=int)
                    ntwists = nkpoints
                else:
                    if nkpoints==0:
                        msg = 'k-points must be present for save2qmcpack to write supercell wavefunction'
                        raise ValueError(msg)
                    #end if
                    tiled_structure.change_units('B') # always use Bohr units for twists/k-points
                    tiled_kpoints = tiled_structure.kpoints.copy()
                    ntwists = len(tiled_kpoints)
                    if ntwists==0:
                        msg = 'supercell k-points must be present for save2qmcpack to write supercell wavefunction'
                        raise ValueError(msg)
                    #end if
                    kmap = tiled_structure.kmap()
                    if len(kmap)!=ntwists:
                        msg = 'inconsistency between supercell twist mapping and supercell twist count'
                        raise ValueError(msg)
                    #end if
                    assert len(kmap)==ntwists
                    kmap_array = []
                    folded_indices = []
                    for itwist in sorted(kmap.keys()):
                        ifolded = list(kmap[itwist])
                        kmap_array.append(ifolded)
                        folded_indices.extend(ifolded)
                    #end for
                    kmap_array = np.array(kmap_array,dtype=int)
                    if len(set(folded_indices))!=nkpoints:
                        msg = 'inconsistency in mapping between supercell twists and folded cell k-points'
                        raise ValueError(msg)
                    #end if
                    if not tiled_structure.has_tmatrix():
                        msg = 'tiling matrix is missing'
                        raise ValueError(msg)
                    #end if
                    tmatrix  = tiled_structure.tmatrix.copy()
                    tiling   = np.diag(tmatrix)
                    diagonal = np.abs(tmatrix-np.diag(tiling)).sum()==0
                    if not diagonal:
                        msg = 'non-diagonal tilings are not yet supported by save2qmcpack'
                        raise NotImplementedError(msg)
                    #end if
                #end if
            #end if
            s = '### generated conversion text ###\n'
            s += 'from PyscfToQmcpack import savetoqmcpack\n'
            self.kpoints       = None
            self.tiled_kpoints = None
            self.tiled_kmap    = None
            if sys_kpoints is None:
                s += "savetoqmcpack({0},{1},'{2}')\n".format(sys_var,mf_var,prefix)
            elif tiled_structure is None:
                s += 'tiling = [{},{},{}]\n'.format(*tiling)
                s += 'for n,kp in enumerate({}):\n'.format(kpts_var)
                s += "    savetoqmcpack({0},{1},'{2}.twistnum_{{}}'.format(str(n).zfill(3)),kmesh=tiling,kpts=[kp],sp_twist=kp)\n".format(sys_var,mf_var,prefix)
                s += "#end for\n"
                self.kpoints = sys_kpoints.copy()
            else:
                s += 'tiling = [{},{},{}]\n'.format(*tiling)
                s += "sp_kpoints = {}\n".format(render_array(tiled_kpoints,4))
                s += "sp_kmap    = {}\n".format(render_array(kmap_array,4))
                s += 'for n,kp in enumerate(sp_kpoints):\n'
                s += "    savetoqmcpack({0},{1},'{2}.twistnum_{{}}'.format(str(n).zfill(3)),kmesh=tiling,kpts={3}[sp_kmap[n]],sp_twist=kp,kmap=sp_kmap[n])\n".format(sys_var,mf_var,prefix,kpts_var)
                s += "#end for\n"
                self.kpoints       = sys_kpoints.copy()
                self.tiled_kpoints = tiled_kpoints.copy()
                self.tiled_kmap    = kmap_array.copy()
            #end if
            s += '### end generated conversion text ###\n'
            self.addendum = '\n'+s+'\n'
        #end if

        if checkpoint and 'chkfile' not in self.values:
            if prefix is None:
                msg = (
                    'cannot set $chkpoint variable\n'
                    'please provide input variable "prefix"'
                    )
                raise ValueError(msg)
            #end if
            chkfile = '{}.chk'.format(prefix)
            self.chkfile = chkfile
            self.assign(chkfile="'"+chkfile+"'")
        #end if

    #end def __init__


    def write_text(self,filepath=None):
        text = SimulationInputTemplateDev.write_text(self,filepath)
        if self.addendum is not None:
            text += self.addendum
        #end if
        return text
    #end def write_text
#end class PyscfInput



def generate_pyscf_input(*args,**kwargs):
    return PyscfInput(*args,**kwargs)
#end def generate_pyscf_input
