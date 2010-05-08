import pygtk
import gtk
import IO
import numpy

class Orbitals(gtk.Frame):
    def h5_chosen_callback(self, fileDialog, response):
        if (response == gtk.RESPONSE_ACCEPT):
            filename = self.FileButton.get_filename()
            # Check to see if the file is a valid PP file
            okay = self.read_h5_file(filename)

    def read_eshdf (self, io):
        # Read primitive lattice
        io.OpenSection('supercell')
        self.prim_vecs = io.ReadVar('primitive_vectors')
        a = numpy.max(numpy.abs(self.prim_vecs))
        io.CloseSection()
        self.Geometry.LatticeFrame.set_lattice(self.prim_vecs)
        self.Geometry.LatticeFrame.ArbRadio.set_active(True)

        # Read atom species
        io.OpenSection('atoms')
        num_species = io.ReadVar ('number_of_species')
        oldtypes = self.Geometry.Types.GetElementTypes()
#        for t in oldtypes:
#            self.Geometry.Types.Remove
        TypeList = []
        for isp in range(0,num_species):
            io.OpenSection('species')
            Z = io.ReadVar('atomic_number')
            Zion = io.ReadVar('valence_charge')
                
            symbol = self.Geometry.Types.Elements.ElementList[Z-1][1]
            TypeList.append(symbol)
            row = self.Geometry.Types.AddRow(None)
            row.set_elem (symbol, Z)
            if (Zion != Z):
                row.combo.set_active(1)
            io.CloseSection()


        # Read atom positions
        N = io.ReadVar('number_of_atoms')
        self.Geometry.AtomPos.set_num_atoms(N)
        pos = io.ReadVar('reduced_positions')
        self.Geometry.AtomPos.set_atom_positions(pos)
        for symbol in TypeList:
            self.Geometry.AtomPos.AddTypeCallback(None, symbol)


        io.CloseSection()

    def read_h5_file(self, filename):
        io = IO.IOSectionClass()
        if (not io.OpenFile(filename)):
            return False
        format = io.ReadVar('format')
        if (format == 'ES-HDF'):
            return self.read_eshdf (io)
        
        return False

    def __init__(self, geometry):
        self.Geometry = geometry
        gtk.Frame.__init__(self, "Orbitals")

        # Setup orbital HDF5 file chooser
        filter = gtk.FileFilter()
        filter.add_pattern("*.h5")
        filter.set_name ("XML files")
        buttons = (gtk.STOCK_CANCEL,gtk.RESPONSE_CANCEL,gtk.STOCK_OPEN,gtk.RESPONSE_ACCEPT)
        self.FileDialog = gtk.FileChooserDialog("Select orbital file", buttons=buttons)
        self.FileDialog.set_action(gtk.FILE_CHOOSER_ACTION_OPEN)
        self.FileDialog.connect("response", self.h5_chosen_callback)
        self.FileButton = gtk.FileChooserButton(self.FileDialog)
        self.FileButton.add_filter(filter)
        self.FileButton.set_sensitive(True)
        self.FileButton.set_action (gtk.FILE_CHOOSER_ACTION_OPEN)
        filebox = gtk.HBox(True)
        filebox.pack_start(self.FileButton)
        self.add(filebox)

class Jastrows(gtk.Frame):
    def __init__(self):
        gtk.Frame.__init__(self, "Jastrow correlation functions")
        
        

class Wavefunction(gtk.VBox):
    def __init__(self, geometry):
        gtk.VBox.__init__(self)
        self.OrbitalsFrame = Orbitals(geometry)
        self.pack_start (self.OrbitalsFrame, False, False, 4)
        self.Geometry = geometry

        self.JastrowsFrame = Jastrows()
        self.pack_start (self.JastrowsFrame, False, False, 4)


        # self.Widgets.append (self.FileButton)
