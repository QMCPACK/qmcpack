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
                 template   = None,   # path to template input file
                 prefix     = None,   # $prefix var for file prefixes
                 custom     = None,   # obj w/ $ prefixed vars in template
                 system     = None,   # physical system object
                 units      = None,   # input units desired
                 use_folded = True,   # use folded system/primitive cell
                 mole       = None,   # obj w/ Mole variables
                 cell       = None,   # obj w/ Cell variables
                 sys_var    = None,   # local var name for Mole/Cell
                 mole_var   = 'mol',  # local var name for Mole in written input
                 cell_var   = 'cell', # local var name for Cell in written input
                 save_qmc   = False,  # convert to QMCPACK format
                 checkpoint = False,  # set $chkfile variable
                 mf_var     = 'mf',   # local var name for mf, used for convert
                 kpts_var   = 'kpts', # local var name for kpts, used for convert
                 filepath   = None,   # alias for template
                 text       = None,   # full text of (and alternate to) template 
                 ):
        if filepath is None and template is not None:
            filepath = template
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
            system = system.copy() # make a local copy
            if use_folded and system.has_folded():
                system = system.folded_system
            #end if
            s = system.structure
            if units is not None:
                s.change_units(units)
            #end if
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
                    skp = s.copy()
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
            for k,v in sys_inputs.items():
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
            has_array |= isinstance(sys_kpoints,ndarray)
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
                self.error('cannot generate save2qmcpack text\nplease provide input variable "prefix"\n(used to set "title" in save2qmcpack)')
            elif sys_var is None:
                self.error('cannot generate save2qmcpack text\nplease provide input variable "sys_var"\n(used to set "cell" in save2qmcpack) ')
            #end if
            s = '### generated conversion text ###\n'
            s += 'from PyscfToQmcpack import savetoqmcpack\n'
            if sys_kpoints is None:
                s += "savetoqmcpack({0},{1},'{2}')\n".format(sys_var,mf_var,prefix)
            else:
                s += "savetoqmcpack({0},{1},'{2}',{3})\n".format(sys_var,mf_var,prefix,kpts_var)
            #end if
            s += '### end generated conversion text ###\n'
            self.addendum = '\n'+s+'\n'
        #end if

        if checkpoint:
            if prefix is None:
                self.error('cannot set $chkpoint variable\nplease provide input variable "prefix"')
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
