#!/bin/env python
from xml.dom.minidom import *
import sys
from check_qmcsystem import *

class simulation:
    system = qmcsystem()


def check_project (proj, sim):
    print 'checking project section'
    return
    
def check_random (rand, sim):
    print 'checking random section'
    return

#def check_qmcsystem(system):
#    print 'checking qmcsystem section'
#    return

def check_wavefunction(wf, sim):
    print 'checking wavefunction section'
    return

def check_hamiltonian(ham, sim):
    print 'checking hamiltonian section'
    return

def check_loop(loop, sim):
    print 'checking loop section'
    return

def check_qmc(qmc, sim):
    print 'checking qmc section'
    return


def check_simulation(sim_elem):
    mysim = simulation()

    sim_dict = {'project'      : check_project     , \
                'random'       : check_random      , \
                'qmcsystem'    : check_qmcsystem   , \
                'wavefunction' : check_wavefunction, \
                'hamiltonian'  : check_hamiltonian , \
                'loop'         : check_loop        , \
                'qmc'          : check_qmc         }

    for elem in sim_elem.childNodes:
        if elem.nodeType == elem.ELEMENT_NODE:
            n = elem.localName
            if (n in sim_dict):
                sim_dict[elem.localName](elem, mysim)
            else:
                print "Error:  unrecognized subelement \"" + n + "\" in "\
                    "element \"simulation\""


##########################
# Main execution routine #
##########################

if (len(sys.argv) < 2):
    print "Usage:  check_input.py  myinput.xml"
    sys.exit()

doc = parse(sys.argv[1])

for elem in doc.childNodes:
    if elem.nodeType == elem.ELEMENT_NODE:
        if (elem.localName != 'simulation'):
            print "Error:  Top-level element should be 'simulation'"
        else:
            check_simulation (elem)
            
