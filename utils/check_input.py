#!/bin/env python
from xml.dom.minidom import *
import sys

def check_project (proj):
    return
    
def check_random (rand):
    print 'checking random section'
    return

def check_qmcsystem(system):
    print 'checking qmcsystem section'
    return

def check_wavefunction(wf):
    print 'checking wavefunction section'
    return

def check_hamiltonian(ham):
    print 'checking hamiltonaina section'
    return

def check_loop(loop):
    print 'checking loop section'
    return

def check_qmc(qmc):
    print 'checking qmc section'
    return


def check_simulation(sim):
    sim_dict = {'project'      : check_project     , \
                'random'       : check_random      , \
                'qmcsystem'    : check_qmcsystem   , \
                'wavefunction' : check_wavefunction, \
                'hamiltonian'  : check_hamiltonian , \
                'loop'         : check_loop        , \
                'qmc'          : check_qmc         }

    for elem in sim.childNodes:
        if elem.nodeType == elem.ELEMENT_NODE:
            n = elem.localName
            if (n in sim_dict):
                sim_dict[elem.localName](elem)
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
            
