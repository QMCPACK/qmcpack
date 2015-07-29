##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


from generic import obj
from periodic_table import pt as ptable
from pseudopotential import GaussianPP


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


def write_hamiltonian(theory):
    s = ''
    theory = theory.lower()
    if theory=='uhf':
        s+=theory.upper()+'\n'
    else:
        error('unknown theory: '+theory)
    #end if
    s += 'END\n'
    return s
#end def write_hamiltonian



