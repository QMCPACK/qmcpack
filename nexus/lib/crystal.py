##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


import os
from generic import obj
from periodic_table import pt as ptable
from pseudopotential import GaussianPP
from nexus_base import nexus_noncore
from simulation import Simulation
from developer import error


molecule_symm_text = '''
1    1       C1     
2    1b      Ci     
3    2x      C2x    
4    2y      C2y    
5    2z      C2z    
6    mx      Csx    
7    my      Csy    
8    mz      Csz    
9    2/mx    C2hx   
10   2/my    C2hy   
11   2/mz    C2hz   
12   222     D2     
13   2mm     C2vx   
14   m2m     C2vy   
15   mm2     C2vz   
16   mmm     D2h    
17   4       C4     
18   4b      S4     
19   4/m     C4h    
20   422     D4     
21   4mm     C4v    
22   4b2m    D2d    
23   4bm2    D2d    
24   4/mmm   D4h    
25   3       C3     
26   3b      C3i    
27   321     D3     
28   312     D3     
29   3m1     C3v    
30   31m     C3v    
31   3bm1    D3d    
32   3b1m    D3d    
33   6       C6     
34   6b      C3h    
35   6/m     C6h    
36   622     D6     
37   6mm     C6v    
38   6bm2    D3h    
39   6b2m    D3h    
40   6/mmm   D6h    
41   23      T      
42   m3b     Th     
43   432     O      
44   4b3m    Td     
45   m3bm    Oh     
46   235     I      
47   m3b5b   Ih     
'''

symm_map = obj(
    molecule = obj(),
    polymer  = obj(),
    slab     = obj(),
    crystal  = obj(),
    )

smap = symm_map.molecule
for line in molecule_symm_text.splitlines():
    ls = line.strip()
    if len(ls)>0:
        index,hm,schoen = line.split()
        index = int(index)
        smap[hm]    = index
        smap[schoen] = index
        # default to x direction if exists
        if hm.endswith('x'):
            smap[hm[:-1]] = index
        #end if
        if schoen.endswith('x'):
            smap[schoen[:-1]] = index
        #end if
    #end if
#end for







def write_geometry(title,bcond,system,symmetry=1,pseudo=True):
    t = '{0}\n{1}\n'.format(title,bcond.upper())
    bcond = bcond.lower()
    if not bcond in symm_map:
        error('unknown boundary conditions: {0}'.format(bcond))
    #end if
    smap = symm_map[bcond]
    if isinstance(symmetry,int):
        t += str(symmetry)+'\n'
    elif symmetry in symm_map:
        t += str(smap[symmetry])+'\n'
    else:
        error('symmetry {0} is unknown'.format(symmetry),'write_geometry')
    #end if
    s = system.structure.copy()
    s.change_units('A')
    t += '{0}\n'.format(len(s.elem))
    if pseudo:
        an_base = 200 # for 'conventional' atomic numbers
    else:
        an_base = 0
    #end if
    for n in xrange(len(s.elem)):
        e = s.elem[n]
        p = s.pos[n]
        conv_atomic_number = ptable[e].atomic_number+an_base
        t += '{0} {1: 12.8f} {2: 12.8f} {3: 12.8f}\n'.format(conv_atomic_number,p[0],p[1],p[2])
    #end for
    t+='END\n'
    return t
#end def write_geometry


def write_basis(pseudos,occupations,formats):
    s = ''
    if len(pseudos)!=len(occupations):
        error('must provide one set of occupations for each pseudopotential','write_basis')
    #end if
    if len(pseudos)!=len(formats):
        error('must specify file format for each pseudopotential','write_basis')
    #end if
    for n in range(len(pseudos)):
        pp = GaussianPP(pseudos[n],format=formats[n])
        s += pp.write_text(format='crystal',occ=occupations[n])
    #end for
    s += '99 0\n'
    s += 'END\n'
    return s
#end def write_basis


def write_hamiltonian(theory          = None,
                      system          = None,
                      spinlock_cycles = 100,
                      levshift        = None,
                      levshift_lock   = 1,
                      maxcycle        = None,
                      ):
    s = ''
    theory = theory.lower().strip()
    dft_functionals = set(['pbe0'])
    known_theories = set(['rhf','uhf']) | dft_functionals
    if theory in known_theories:
        if theory not in dft_functionals:
            s+=theory.upper()+'\n'
        else:
            s+='DFT\n'+theory.upper()+'\nSPIN\nEND\n'
        #end if
    else:
        error('unknown theory: '+theory,'write_hamiltonian')
    #end if
    if system!=None:
        s+='SPINLOCK\n{0} {1}\n'.format(system.net_spin,spinlock_cycles)
    #end if
    if levshift!=None:
        s+='LEVSHIFT\n{0} {1}\n'.format(levshift,levshift_lock)
    #end if
    if maxcycle!=None:
        s+='MAXCYCLE\n{0}\n'.format(maxcycle)
    #end if
    s += 'END\n'
    return s
#end def write_hamiltonian


class CrystalSim: None
class PropertiesSim: None


from simulation import GenericSimulation
class Crystal_1(GenericSimulation,CrystalSim):
    generic_identifier = 'crystal'
    application = 'crystal'
    infile_extension = '.d12'
    application_results = set(['orbitals'])
    
    def check_result(self,result_name,sim):
        return result_name=='orbitals'
    #end def check_result

    def get_result(self,result_name,sim):
        result = obj() # its up to the other application how to handle crystal
        if result_name!='orbitals':
            self.error('ability to get result '+result_name+' has not been implemented')
        #end if
        return result
    #end def get_result

    def incorporate_result(self,result_name,result,sim):
        self.error('ability to incorporate result '+result_name+' has not been implemented')
    #end def incorporate_result

    def app_command(self):
        return self.app_name+'<'+self.infile
    #end def app_command
#end class Crystal_1



from simulation import generate_simulation,input_template_dev
def gen_crystal_molecule_sim(
    identifier  = 'crystal',
    path        = None,
    job         = None,
    title       = 'atom',
    system      = None,
    pseudos     = None,
    occupations = None,
    formats     = None,
    theory      = None,
    levshift    = None,
    maxcycle    = None,
    ):
    required = obj(path=path,system=system,pseudos=pseudos,occupations=occupations,formats=formats,theory=theory,job=job)
    for k,v in required.iteritems():
        if v is None:
            error(k+' is a required input','gen_crystal_molecule_sim')
        #end if
    #end for
    pseudos_in = pseudos
    pseudos = []
    for pp in pseudos_in:
        pseudos.append(os.path.join(nexus_noncore.pseudo_dir,pp))
    #end for

    crys_input = input_template_dev()
    crys_input.read_text(
        write_geometry(
            title = title,
            bcond = 'molecule',
            system = system
            )+
        write_basis(
            pseudos     = pseudos,
            occupations = occupations,
            formats     = formats,
            )+
        write_hamiltonian(
            theory   = theory,
            system   = system,
            levshift = levshift,
            maxcycle = maxcycle
            )
        )

    #crys = generate_simulation(
    #    identifier = identifier,
    #    path       = path,
    #    job        = job(cores=1,serial=True,app_command='crystal<{0}.d12>&{0}.out&'.format(identifier)),
    #    input      = crys_input,
    #    infile     = '{0}.d12'.format(identifier),
    #    outfile    = '{0}.out'.format(identifier),
    #    )


    crys = Crystal_1(
        identifier = identifier,
        path       = path,
        job        = job,
        input      = crys_input,
        )


    return crys
#end def gen_crystal_molecule_sim






class Properties_1(GenericSimulation,PropertiesSim): # same implementation, not 'is a'
    generic_identifier = 'properties'
    application = 'properties'
    infile_extension = '.d3'
    application_results = set(['orbitals'])
    
    def check_result(self,result_name,sim):
        return result_name=='orbitals'
    #end def check_result

    def get_result(self,result_name,sim):
        result = obj() # its up to the other application how to handle crystal
        if result_name!='orbitals':
            self.error('ability to get result '+result_name+' has not been implemented')
        #end if
        return result
    #end def get_result

    def incorporate_result(self,result_name,result,sim):
        if result_name!='orbitals' or not isinstance(sim,Crystal_1):
            self.error('ability to incorporate result '+result_name+' has not been implemented')
        #end if
    #end def incorporate_result

    def app_command(self):
        return self.app_name+'<'+self.infile
    #end def app_command
#end class Properties_1


molecule_text = '''
NEWK
1 0
CRYAPI_OUT
END
'''
periodic_text = '''
NEWK
0 0
1 0
CRYAPI_OUT
END
'''

def gen_properties(**kwargs):

    if 'systype' not in kwargs:
        error('systype is a required input','gen_properties')
    #end if
    systype = kwargs['systype']
    del kwargs['systype']

    if systype=='molecule_qmc':
        text = molecule_text
    elif systype=='periodic_qmc':
        text = periodic_text
    else:
        error('invalid systype encountered\nsystype provided: {0}\nvalid options are: molecule_qmc, periodic_qmc'.format(systype))
    #end if

    sim_args,inp_args = Simulation.separate_inputs(kwargs)
    if len(inp_args)>0:
        error('invalid arguments encountered\ninvalid keywords: {0}'.format(sorted(inp_args.keys())),'gen_properties')
    #end if
    if not 'input' in sim_args:
        sim_args.input = input_template_dev(text=text.strip())
    #end if
    prop = Properties_1(**sim_args)
    return prop
#end def gen_properties
