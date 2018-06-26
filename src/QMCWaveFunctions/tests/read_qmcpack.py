
# Read parts of the QMCPACK XML input file

import xml.etree.ElementTree as ET
from collections import namedtuple, defaultdict
import math
import gaussian_orbitals
import numpy as np


# Get the total number of basis functions after expanding the angular parts
def get_total_basis_functions(basis_set):
    ijk_list = gaussian_orbitals.get_ijk()
    angular_funcs = defaultdict(int)
    for i,j,k,name in ijk_list:
        L = i + j + k
        angular_funcs[L] += 1

    total = 0
    for cg in basis_set:
        total += angular_funcs[cg.orbtype]
    return total

def read_basis_groups(atomic_basis_set):
    basis_groups =  atomic_basis_set.findall('basisGroup')
    #print basis_groups
    basis_set = []
    for basis_group in basis_groups:
        if basis_group.attrib['type'] != 'Gaussian':
            print 'Expecting Gaussian type basisGroup'
        #print basis_group.attrib['n']
        n = int(basis_group.attrib['n'])
        #print basis_group.attrib['l']
        ang_mom_l = int(basis_group.attrib['l'])
        #print basis_group.attrib['type']
        zeta_list = []
        coef_list = []
        radfuncs = basis_group.findall('radfunc')
        for radfunc in radfuncs:
            zeta = float(radfunc.attrib['exponent'])
            contraction_coef =  float(radfunc.attrib['contraction'])
            zeta_list.append(zeta)
            coef_list.append(contraction_coef)

        cg = gaussian_orbitals.CG_basis(ang_mom_l, len(zeta_list), zeta_list, coef_list)
        basis_set.append(cg)

    return basis_set


#  Read the atomic basis set and MO coefficients

def parse_qmc_wf(fname):
    tree = ET.parse(fname)

    atomic_basis_sets = tree.findall('.//atomicBasisSet')
    basis_sets = dict()
    for atomic_basis_set in atomic_basis_sets:
      element = atomic_basis_set.attrib['elementType']
      basis_set = read_basis_groups(atomic_basis_set)
      basis_sets[element] = basis_set


    basis_size = 0
    for element, basis_set in basis_sets.iteritems():
      basis_size += get_total_basis_functions(basis_set)
    #print 'total basis size',basis_size

    #  Just use the first one for now - assume up and down MO's are the same
    MO_coeff_node = tree.find('.//determinant/coefficient')
    MO_matrix = None
    if MO_coeff_node is None:
        print 'Molecular orbital coefficients not found'
    else:
        #print 'MO_coeff = ',MO_coeff_node
        MO_size = int(MO_coeff_node.attrib['size'])
        #print 'MO coeff size = ',MO_size

        MO_text = MO_coeff_node.text
        MO_text_entries = MO_text.split()
        MO_values = [float(a) for a in MO_text_entries]

        #MO_matrix = np.array(MO_values).reshape( (basis_size, MO_size) )
        MO_matrix = np.array(MO_values).reshape( (MO_size, basis_size) )
        #print 'MO_values = ',MO_values

    return basis_sets, MO_matrix


# Read the ion positions and types

def parse_structure(node):
    particleset = node.find("particleset[@name='ion0']")
    npos = int(particleset.attrib['size'])
    pos_node = particleset.find("attrib[@name='position']")
    pos_values = [float(a) for a in pos_node.text.split()]
    pos = np.array(pos_values).reshape( (3, npos) )

    elements_node = particleset.find("attrib[@name='ionid']")
    elements = elements_node.text.split()
    return pos, elements


def read_structure_file(fname):
  tree = ET.parse(fname)
  return parse_structure(tree)



if __name__ == '__main__':
    # For He
    basis_set, MO_matrix = parse_qmc_wf('he_sto3g.wfj.xml')


    # For HCN - need ion positions as well
    basis_sets, MO_matrix = parse_qmc_wf('hcn.wfnoj.xml')

    pos_list, elements = read_structure_file("hcn.structure.xml")

    gtos = gaussian_orbitals.GTO_centers(pos_list, elements, basis_sets)
    atomic_orbs =  gtos.eval_v(1.0, 0.0, 0.0)
    #print np.dot(MO_matrix, atomic_orbs)

