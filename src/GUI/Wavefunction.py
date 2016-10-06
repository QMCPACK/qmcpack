#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign   
#//
#// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign 
#//////////////////////////////////////////////////////////////////////////////////////


import pygtk
import gtk
import IO
import numpy

class TilingMatrix(gtk.Frame):
    def __init__(self):
        gtk.Frame.__init__(self, 'Orbital tiling')
        self.matrix = numpy.array([[1,0,0],[0,1,0],[0,0,1]])
        self.set_label('Orbital tiling')
        self.TileTable = gtk.Table(3,3)
        self.TileOrbitals = gtk.CheckButton('Tile orbitals')
        self.TileOrbitals.set_active(False)
        self.TileOrbitals.connect("toggled", self.tile_toggled)

        self.TileButtons = []
        for i in range(0,3):
            TileList = []
            for j in range(0,3):
                tile = gtk.SpinButton\
                    (gtk.Adjustment(0.0, -100.0, 100.0, 1.0, 2.0))
                if (i == j):
                    tile.set_value(self.matrix[i,j])
                tile.set_digits(0)
                tile.set_width_chars(2)
                tile.connect('value_changed', self.matrix_changed)
                self.TileTable.attach(tile, i, i+1, j, j+1)
                TileList.append(tile)
            self.TileButtons.append(TileList)
        vbox = gtk.VBox()
        self.UnitLabel = gtk.Label()
        self.UnitLabel.set_text('Unit cells:  1')
        vbox.pack_start(self.TileOrbitals)
        vbox.pack_start(self.TileTable)
        vbox.pack_start(self.UnitLabel)
        self.TileTable.set_sensitive(False)
        self.add(vbox)

    def matrix_changed(self, button):
        units = self.get_units()
        self.UnitLabel.set_text('Unit cells:  %d' %(units))

    def get_units(self):
       mat = self.get_matrix()
       units = numpy.abs(numpy.linalg.det(mat))
       return units

    def get_matrix(self):
        mat = []
        for i in range(0,3):
            row = []
            for j in range(0,3):
                row.append(int(self.TileButtons[i][j].get_value()))
            mat.append(row)
        return numpy.array(mat)

    def set_matrix(self, mat):
        for i in range(0,3):
            for j in range(0,3):
                TileButtons[i,j].set_value(mat[i,j])

    def tile_toggled(self, button):
        self.TileTable.set_sensitive(button.get_active())
    



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
        buttons = (gtk.STOCK_CANCEL,gtk.RESPONSE_CANCEL,\
                   gtk.STOCK_OPEN,gtk.RESPONSE_ACCEPT)
        self.FileDialog = gtk.FileChooserDialog \
            ("Select orbital file", buttons=buttons)
        self.FileDialog.set_action(gtk.FILE_CHOOSER_ACTION_OPEN)
        self.FileDialog.connect("response", self.h5_chosen_callback)
        self.FileButton = gtk.FileChooserButton(self.FileDialog)
        self.FileButton.add_filter(filter)
        self.FileButton.set_sensitive(True)
        self.FileButton.set_action (gtk.FILE_CHOOSER_ACTION_OPEN)
        filebox = gtk.HBox(True)
        vbox = gtk.VBox(True)
        self.TileFrame = TilingMatrix()
        filebox.pack_start(self.FileButton, True, False)
        filebox.pack_start(self.TileFrame , True, False)

        self.add(filebox)

    def tile_matrix_changed(self, button):
        print

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
