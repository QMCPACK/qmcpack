#!/bin/env python
from xml.dom.minidom import *
import sys
from pylab import *
from numpy import *


class channel:
    def __init__(self):
        self.channel_dict = { 's':0, 'p':1, 'd':2, 'f':3, 'g':4, 'h':5 }

    def parse_grid(self, elem):
        self.ri = float(elem.attributes['ri'].value)
        self.rf   = float(elem.attributes['rf'].value)
        self.npts = int(elem.attributes['npts'].value)
        delta = (self.rf-self.ri)/(self.npts-1)
        self.grid = arange(self.ri, self.rf+delta, delta)
        #print self.grid

    def parse_data(self, elem):
        for e in elem.childNodes:
            if (e.nodeType == elem.TEXT_NODE):
                data = e.nodeValue.split()
                list = []
                for d in data:
                    list.append(float(d))
                self.V = array(list)
                # print self.V.shape
        
    def parse_radfunc(self, elem):
        for e in elem.childNodes:
            if (e.nodeType == elem.ELEMENT_NODE):
                if (e.localName == 'grid'):
                    self.parse_grid(e)
                elif (e.localName == 'data'):
                    self.parse_data(e)

    def parse(self, elem):
        lstr = elem.attributes['l'].value
        #print 'l = ' + lstr
        self.l = self.channel_dict[lstr]
        self.cutoff = float(elem.attributes['cutoff'].value)
        #print ('l=%d  cutoff=%1.5f' % (self.l, self.cutoff))
        for e in elem.childNodes:
            if (e.nodeType == elem.ELEMENT_NODE):
                if (e.localName == 'radfunc'):
                    self.parse_radfunc(e)

class pseudo:
    def __init__(self):
        self.channels = []
        print 'pseudo constructor'

    def parse_semilocal(self, semi_elem):
        for elem in semi_elem.childNodes:
           if elem.nodeType == elem.ELEMENT_NODE:
                n = elem.localName
                if(n == 'vps'):
                    ch = channel()
                    ch.parse(elem)
                    self.channels.append(ch)
        
    def plot (self):
        spd = ['s', 'p', 'd', 'f', 'g', 'h', 'i']
        leglist = []
        r = []
        for ch in self.channels:
            Vofr = ch.V / ch.grid
            Vofr[0] = Vofr[1]
            plot (ch.grid, Vofr, lw=1.5)
            r = ch.grid
            leglist.append(spd[ch.l])
        a = axis()
        plot (r, -1.0/r * self.Z)
        axis(a)
        leg = legend(leglist, loc='lower right')
        setp(leg.get_texts(), 'FontSize', 16)
        xlabel (r'$r$ (bohr)', size=16)
        ylabel (r'$V_l(r)$ (Hartree)', size=16)
        show()

    def parse_grid (self, grid_elem):
        print 'parsing grid'

    def parse_header(self, header_elem):
        self.Z = int(header_elem.attributes['zval'].value)

    def parse_wavefunctions(self, header_elem):
        print 'parsing wavefunctions'

    def parse(self, pseudo_elem):
        pseudo_dict = {'semilocal' : self.parse_semilocal, \
                       'grid'      : self.parse_grid     , \
                       'header'    : self.parse_header   , \
                       'pseudowave-functions' : self.parse_wavefunctions}
        for elem in pseudo_elem.childNodes:
            if elem.nodeType == elem.ELEMENT_NODE:
                n = elem.localName
                if (n in pseudo_dict):
                    pseudo_dict[elem.localName](elem)
                else:
                    print "Error:  unrecognized subelement \"" + n + "\" in "\
                        "element \"simulation\""



##########################
# Main execution routine #
##########################
if (len(sys.argv) < 2):
    print "Usage:  plot_fsatom.py  my_pp.xml"
    sys.exit()

doc = parse(sys.argv[1])

for elem in doc.childNodes:
    if elem.nodeType == elem.ELEMENT_NODE:
        if (elem.localName != 'pseudo'):
            print "Error:  Top-level element should be 'pseudo'"
        else:
            pp = pseudo()
            pp.parse(elem)
            pp.plot()
