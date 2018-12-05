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


from numpy import ndarray
from generic import obj
from developer import error
from simulation import SimulationInputTemplateDev



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
        error('render_array only supports up to 2D arrays')
    #end if
    return s
#end def render_array


class PyscfInput(SimulationInputTemplateDev):

    basic_types = bool,int,float,str,tuple,list,dict

    allowed_types = basic_types+(ndarray,)

    mole_order = '''
        dump_input
        parse_arg
        verbose
        output
        max_memory
        atom
        basis
        unit
        nucmod
        ecp
        charge
        spin
        symmetry
        symmetry_subgroup
        cart
        nelec
        nelectron
        multiplicity
        ms
        '''.split()

    cell_order = '''
        dump_input
        parse_arg
        a
        mesh
        ke_cutoff
        precision
        nimgs
        ew_eta
        ew_cut
        pseudo
        basis
        h
        dimension
        rcut
        ecp
        low_dim_ft_type
        unit
        atom
        gs
        h
        drop_exponent
        nimgs
        '''.split()
    for k in mole_order:
        if k not in cell_order:
            cell_order.append(k)
        #end if
    #end for

    mole_allowed = set(mole_order)
    cell_allowed = set(cell_order)
    

    def __init__(self,
                 filepath = None,
                 text     = None,
                 custom   = None,
                 system   = None,
                 mole     = None,
                 cell     = None,
                 mole_var = 'mol',
                 cell_var = 'cell',
                 ):
        SimulationInputTemplateDev.__init__(self,filepath,text)

        if custom is not None:
            self.assign(**custom)
        #end if

        sys_name    = None
        sys_inputs  = obj()
        sys_allowed = None
        sys_kpoints = None
        is_mole = mole is not None
        is_cell = cell is not None
        if is_mole and is_cell:
            self.error('both mole and cell provided\nplease provide only one of them\nsystem cannot be both molecule and periodic solid')
        elif is_mole and not isinstance(mole,(dict,obj)):
            self.error('mole input must be a dict or obj\nyou provided input of type: {0}'.format(mole.__class__.__name__))
        elif is_cell and not isinstance(cell,(dict,obj)):
            self.error('cell input must be a dict or obj\nyou provided input of type: {0}'.format(cell.__class__.__name__))
        #end if
        extra = ''
        if filepath is not None:
            extra = '\ntemplate located at: {0}'.format(filepath)
        #end if
        if system is not None and 'system' not in self.keywords:
            self.error('system input is provided, but $system is not present in template input'+extra)
        #end if
        if system is not None and 'system' not in self.values:
            if 'system' not in self.keywords:
                self.error('cannot incorporate "system" input\n$system is not present in template input'+extra)
            #end if
            if system.has_folded():
                system = system.folded_system
            #end if
            s = system.structure
            is_solid    = s.has_axes()
            is_molecule = not is_solid
            if is_solid and is_mole:
                self.error('mole input provided, but provided system is not a molecule (cell axes are present)')
            elif is_molecule and is_cell:
                self.error('cell input provided, but provided system is a molecule (cell axes are not present)')
            #end if
            is_mole |= is_molecule
            is_cell |= is_solid
            sys_inputs.atom   = s.write_xyz(header=False,units=s.units)
            sys_inputs.unit   = s.units
            sys_inputs.charge = system.net_charge
            sys_inputs.spin   = system.net_spin
            if is_solid:
                sys_inputs.dimension = len(s.axes)
                sys_inputs.a         = s.write_axes()
                if len(s.kpoints)>0:
                    sys_kpoints = s.kpoints.copy()
                #end if
            #end if
        #end if
        if is_mole:
            sys_name    = 'mole'
            sys_var     = mole_var
            sys_allowed = PyscfInput.mole_allowed
            sys_order   = PyscfInput.mole_order
            if mole is not None:
                sys_inputs.set(**mole)
            #end if
        elif is_cell:
            sys_name    = 'cell'
            sys_var     = cell_var
            sys_allowed = PyscfInput.cell_allowed
            sys_order   = PyscfInput.cell_order
            if cell is not None:
                sys_inputs.set(**cell)
            #end if
        else:
            None # no action needed if not molecule or periodic solid
        #end if

        if sys_name is not None:
            invalid = set(sys_inputs.keys())-sys_allowed
            if len(invalid)>0:
                self.error('invalid {0} inputs\ninvalid inputs: {1}\nvalid options are: {2}'.format(sys_name,sorted(invalid),sorted(sys_allowed)))
            #end if
            klen = 0
            has_array = False
            for k,v in sys_inputs.iteritems():
                if not isinstance(v,PyscfInput.allowed_types):
                    tlist = ''
                    for t in PyscfInput.allowed_types:
                        tlist += t.__name__+','
                    #end for
                    tlist = tlist[:-1]
                    self.error('{0} input "{1}" has an invalid type\ninvalid type: {2}\nallowed types are: {3}'.format(sys_name,k,v.__class__.__name__,tlist))
                #end if
                klen = max(klen,len(k))
                has_array |= isinstance(v,ndarray)
            #end for
            c = '\n### generated system text ###\n'
            if has_array:
                c += 'from numpy import array\n'
            #end if
            if is_mole:
                c += 'from pyscf import gto_loc\n'
                c += '{0} = gto_loc.Mole()\n'.format(sys_var)
            elif is_cell:
                c += 'from pyscf.pbc import gto_loc\n'
                c += '{0} = gto_loc.Cell()\n'.format(sys_var)
            #end if
            fmt = sys_var+'.{0:<'+str(klen)+'} = {1}\n'
            nalign = klen+len(sys_var)+4
            for k in sys_order:
                if k in sys_inputs:
                    v = sys_inputs[k]
                    if isinstance(v,str) and '\n' in v:
                        vs = render_string(v,nalign)
                    elif isinstance(v,ndarray):
                        if len(v.shape)>2:
                            self.error('cannot write system input variable {0}\n{0} is an array with more than two dimensions\nonly two dimensions are currently supported for writing\narray contents: {1}'.format(k,v))
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
            c += '### end generated system text ###\n\n'
            self.assign(system=c)
        #end if

    #end def __init__
#end class PyscfInput


def generate_pyscf_input(*args,**kwargs):
    return PyscfInput(*args,**kwargs)
#end def generate_pyscf_input
