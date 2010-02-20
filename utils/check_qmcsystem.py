

class qmcsystem:
    lattice = [[0.0, 0.0, 0.0],\
               [0.0, 0.0, 0.0],\
               [0.0, 0.0, 0.0]]\

    def get_parameter (self, par_elem):
        name = par_elem.attributes["name"]
        print "name = " + name.value
        val = ""
        for elem in par_elem.childNodes:
            print elem
            if elem.nodeType == elem.TEXT_NODE:
                print elem
                val = elem.text
        return (name,val)

    def check_cell(self, cell_elem):
        for elem in cell_elem.childNodes:
            if elem.nodeType == elem.ELEMENT_NODE:
                n = elem.localName
                print n
                if n == 'parameter':
                    self.get_parameter(elem)
        return

    def check_particleset (self, pset_elem):
        print "checking particleset"
        return

    mydict = {"simulationcell" : check_cell,\
              "particleset"     : check_particleset }

    def check (self, sys_elem, sim):
        print 'checking qmcsystem element'
        print 'lattice = ' + repr (sim.system.lattice)
        for elem in sys_elem.childNodes:
            if elem.nodeType == elem.ELEMENT_NODE:
                n = elem.localName
                print n
                if (n in self.mydict):
                    self.mydict[n](self, elem)




def check_qmcsystem(sys_elem, sim):
    sim.system.check (sys_elem, sim)
    return
