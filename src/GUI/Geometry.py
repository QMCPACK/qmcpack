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
pygtk.require('2.0')
import gtk
import gobject
import pango
from xml.dom import getDOMImplementation
import sys
from xml.dom.ext.reader import Sax2
from xml.dom.ext import PrettyPrint
from xml import xpath


######################################################################
# class Lattice:                                                     #
#   This class allows the user to select among the standard crystal  #
# types and set their parameters.  Alternatively, the user may set   #
# the simulation lattice vectors to arbitrary values.                #
######################################################################
class Lattice(gtk.Frame):
    def __init__(self):
        gtk.Frame.__init__(self, "Lattice")
        self.set_label_align (0.5, 0.5);
        TypeBox = gtk.VBox()
        HBox = gtk.HBox()
        TypeFrame  = gtk.Frame ("Type")
        ParamFrame = gtk.Frame ("Parameters")
        self.CubicRadio = gtk.RadioButton (None,            "Simple Cubic")
        self.TetraRadio = gtk.RadioButton (self.CubicRadio, "Tetragonal")
        self.FCCRadio   = gtk.RadioButton (self.CubicRadio, "FCC")
        self.BCCRadio   = gtk.RadioButton (self.CubicRadio, "BCC")
        self.OrthoRadio = gtk.RadioButton (self.CubicRadio, "Orthorhombic")
        self.ArbRadio   = gtk.RadioButton (self.CubicRadio, "Arbitrary")
        self.CubicRadio.connect ("toggled", self.RadioCallback, "Cubic")
        self.FCCRadio.connect   ("toggled", self.RadioCallback, "FCC")
        self.BCCRadio.connect   ("toggled", self.RadioCallback, "BCC")
        self.OrthoRadio.connect ("toggled", self.RadioCallback, "Ortho")
        self.ArbRadio.connect   ("toggled", self.RadioCallback, "Arb")
        TypeBox.pack_start (self.CubicRadio)
        TypeBox.pack_start (self.FCCRadio)
        TypeBox.pack_start (self.BCCRadio)
        TypeBox.pack_start (self.OrthoRadio)
        TypeBox.pack_start (self.ArbRadio)
        TypeFrame.add (TypeBox)

        # Setup the lattice vector table
        VectorTable=gtk.Table(3, 4, False)
        a0label = gtk.Label(); a0label.set_markup("a<small><sub>0</sub></small>")
        a1label = gtk.Label(); a1label.set_markup("a<small><sub>1</sub></small>")
        a2label = gtk.Label(); a2label.set_markup("a<small><sub>2</sub></small>")
        VectorTable.attach(a0label, 0, 1, 0, 1, gtk.SHRINK, gtk.SHRINK, 5);
        VectorTable.attach(a1label, 0, 1, 1, 2, gtk.SHRINK, gtk.SHRINK, 5);
        VectorTable.attach(a2label, 0, 1, 2, 3, gtk.SHRINK, gtk.SHRINK, 5);
        self.SpinButtons = []
        for row in range(0,3):
            spinlist = []
            for col in range(1,4):
                spin = gtk.SpinButton(\
                    gtk.Adjustment(0.0, -1.0e5, 1.00e5, 0.01, 0.1))
                spinlist.append(spin)
                spin.set_digits(5)
                spin.set_width_chars(8)
                VectorTable.attach (spin, col, col+1, row, row+1)
            self.SpinButtons.append(spinlist)
        VectorFrame = gtk.Frame ("Lattice Vectors");
        VectorFrame.add(VectorTable)

        # Setup the parameters
        ParamBox = gtk.VBox()
        # Cubic parameters
        ParamTable = gtk.Table(6,4)
        self.aButton = gtk.SpinButton\
                       (gtk.Adjustment(1.0, 0.0, 1e5, 0.01, 0.1))
        self.aButton.set_digits(5); self.aButton.set_width_chars(8);
        self.bButton = gtk.SpinButton(gtk.Adjustment(1.0, 0.0, 1e5, 0.01, 0.1));
        self.bButton.set_digits(5); self.bButton.set_width_chars(8);
        self.cButton = gtk.SpinButton(gtk.Adjustment(1.0, 0.0, 1e5, 0.01, 0.1));
        self.cButton.set_digits(5); self.cButton.set_width_chars(8);
        self.alphaButton = gtk.SpinButton(gtk.Adjustment(90.0, 0.0, 180.0, 0.01, 0.1));
        self.alphaButton.set_digits(3)
        self.betaButton = gtk.SpinButton(gtk.Adjustment(90.0, 0.0, 180.0, 0.01, 0.1));
        self.betaButton.set_digits(3)
        self.gammaButton = gtk.SpinButton(gtk.Adjustment(90.0, 0.0, 180.0, 0.01, 0.1));
        self.gammaButton.set_digits(3)
        
        # Set parameters changing callback
        self.aButton.connect    ("value_changed", self.parameters_callback)
        self.bButton.connect    ("value_changed", self.parameters_callback)
        self.cButton.connect    ("value_changed", self.parameters_callback)        
        self.alphaButton.connect("value_changed", self.parameters_callback)
        self.betaButton.connect ("value_changed", self.parameters_callback)
        self.gammaButton.connect("value_changed", self.parameters_callback)        
                             
        ParamTable.attach(gtk.Label("a"), 0, 1, 0, 1)
        ParamTable.attach(gtk.Label("b"), 0, 1, 1, 2)
        ParamTable.attach(gtk.Label("c"), 0, 1, 2, 3)
        alphaLabel = gtk.Label();
        alphaLabel.set_markup("<span foreground=\"blue\" font_family=\"Standard Symbols L\">&#945;</span>")
        betaLabel = gtk.Label();
        betaLabel.set_markup("<span foreground=\"blue\" font_family=\"Standard Symbols L\">&#946;</span>")
        gammaLabel = gtk.Label();
        gammaLabel.set_markup("<span foreground=\"blue\" font_family=\"Standard Symbols L\">&#947;</span>")
        ParamTable.attach(alphaLabel,  2, 3, 0, 1, gtk.SHRINK, gtk.SHRINK, 5, 2)
        ParamTable.attach(betaLabel ,  2, 3, 1, 2, gtk.SHRINK, gtk.SHRINK, 5, 2)
        ParamTable.attach(gammaLabel,  2, 3, 2, 3, gtk.SHRINK, gtk.SHRINK, 5, 2)
        ParamTable.attach(self.aButton,     1, 2, 0, 1, gtk.SHRINK, gtk.SHRINK, 5, 2)
        ParamTable.attach(self.bButton,     1, 2, 1, 2, gtk.SHRINK, gtk.SHRINK, 5, 2)
        ParamTable.attach(self.cButton,     1, 2, 2, 3, gtk.SHRINK, gtk.SHRINK, 5, 2)
        ParamTable.attach(self.alphaButton, 3, 4, 0, 1, gtk.SHRINK, gtk.SHRINK, 5, 2)
        ParamTable.attach(self.betaButton,  3, 4, 1, 2, gtk.SHRINK, gtk.SHRINK, 5, 2)
        ParamTable.attach(self.gammaButton, 3, 4, 2, 3, gtk.SHRINK, gtk.SHRINK, 5, 2)
        ParamBox.pack_start(ParamTable)
        ParamFrame.add (ParamBox)

        # Pack the main Lattice box
        HBox.pack_start (TypeFrame,   False, False, 3)
        HBox.pack_start (ParamFrame,  False, False, 3)
        HBox.pack_start (VectorFrame, False, False, 3)
        self.set_cubic()
        self.lattice_sensitive(False)
        self.params_sensitive([True, False, False, False, False, False])
        self.add (HBox)

    # Takes a 2D list and sets the lattice spin buttons to those values
    def set_lattice (self, lattice):
        for i in range(0,3):
            for j in range(0,3):
                self.SpinButtons[i][j].set_value (lattice[i][j])

    def get_lattice (self):
        lattice = []
        for i in range(0,3):
            row = []
            for j in range (0,3):
                row.append(self.SpinButtons[i][j].get_value())
            lattice.append(row)
        return lattice

    # The following routines set the lattice spin buttons to the appropriate values
    # based on the parameter values.
    def set_cubic(self):
        a = self.aButton.get_value()
        self.set_lattice ([ [a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]])
        self.alphaButton.set_value(90.0)
        self.betaButton.set_value (90.0)
        self.gammaButton.set_value(90.0)

    def set_FCC(self):
        a = self.aButton.get_value()
        h = 0.5*a
        self.set_lattice([ [0.0, h, h], [h, 0.0, h], [h, h, 0.0]])

    def set_BCC(self):
        a = self.aButton.get_value()
        h = 0.5*a
        self.set_lattice([[-h, h, h],[h, -h, h], [h, h, -h]])

    def set_ortho(self):
        a = self.aButton.get_value()
        b = self.bButton.get_value()
        c = self.cButton.get_value()
        self.set_lattice([[a, 0.0, 0.0,], [0.0, b, 0.0], [0.0, 0.0, c]])

    # Sets the spin buttons insensitive or sensitive.  They are only sensitive
    # when the Arbitray lattice radio button is selected.
    def lattice_sensitive(self, sensitive):
        for i in range(0,3):
            for j in range(0,3):
                self.SpinButtons[i][j].set_sensitive (sensitive)
                
    def params_sensitive (self, sensitive):
        self.aButton.set_sensitive (sensitive[0])
        self.bButton.set_sensitive (sensitive[1])
        self.cButton.set_sensitive (sensitive[2])
        self.alphaButton.set_sensitive (sensitive[3])
        self.betaButton.set_sensitive (sensitive[4])
        self.gammaButton.set_sensitive (sensitive[5])

    def parameters_callback (self, button):
        if (self.CubicRadio.get_active()):
            self.set_cubic()
        elif (self.FCCRadio.get_active()):
            self.set_FCC()
        elif (self.FCCRadio.get_active()):
            self.set_FCC()
        elif (self.BCCRadio.get_active()):
            self.set_BCC()
        elif (self.OrthoRadio.get_active()):
            self.set_ortho()
        
    def RadioCallback(self, widget, name):
        if (widget.get_active()):
            if (name == "Cubic" or name=="BCC" or name=="FCC"):
                self.lattice_sensitive(False)
                self.params_sensitive(\
                    [True, False, False, False, False, False])
                if (name=="Cubic"):
                    self.set_cubic()
                if (name=="FCC"):
                    self.set_FCC()
                if (name=="BCC"):
                    self.set_BCC()
            if (name == "Ortho"):
                self.lattice_sensitive(False)
                self.params_sensitive(\
                    [True, True, True, False, False, False])
                self.set_ortho()
            if (name == "Arb"):
                self.lattice_sensitive(True)
                self.params_sensitive(\
                    [False, False, False, False, False, False])


######################################################################
# class PeriodicTable:                                               #
######################################################################
class PeriodicTable(gtk.Dialog):
    def __init__(self, parent):
        gtk.Dialog.__init__(self, "Please select an element", parent, \
                            gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,\
                            (gtk.STOCK_OK, gtk.RESPONSE_ACCEPT))
        self.ElementList = []
        ElementTable = gtk.Table(10, 18, True)
        # List row consists of Z, name, row, column
        self.ElementList.append ([1,  "H"  , 1,  1])
        self.ElementList.append ([2,  "He" , 1, 18])
        self.ElementList.append ([3,  "Li" , 2,  1])
        self.ElementList.append ([4,  "Be" , 2,  2])
        self.ElementList.append ([5,  "B"  , 2, 13])
        self.ElementList.append ([6,  "C"  , 2, 14])
        self.ElementList.append ([7,  "N"  , 2, 15])
        self.ElementList.append ([8,  "O"  , 2, 16])
        self.ElementList.append ([9,  "F"  , 2, 17])
        self.ElementList.append ([10, "Ne" , 2, 18])
        self.ElementList.append ([11, "Na" , 3,  1])
        self.ElementList.append ([12, "Mg" , 3,  2])
        self.ElementList.append ([13, "Al" , 3, 13])
        self.ElementList.append ([14, "Si" , 3, 14])
        self.ElementList.append ([15, "P"  , 3, 15])
        self.ElementList.append ([16, "S"  , 3, 16])
        self.ElementList.append ([17, "Cl" , 3, 17])
        self.ElementList.append ([18, "Ar" , 3, 18])
        self.ElementList.append ([19, "K"  , 4,  1])
        self.ElementList.append ([20, "Ca" , 4,  2])
        self.ElementList.append ([21, "Sc" , 4,  3])
        self.ElementList.append ([22, "Ti" , 4,  4])
        self.ElementList.append ([23, "V"  , 4,  5])
        self.ElementList.append ([24, "Cr" , 4,  6])
        self.ElementList.append ([25, "Mn" , 4,  7])
        self.ElementList.append ([26, "Fe" , 4,  8])
        self.ElementList.append ([27, "Co" , 4,  9])
        self.ElementList.append ([28, "Ni" , 4, 10])
        self.ElementList.append ([29, "Cu" , 4, 11])
        self.ElementList.append ([30, "Zn" , 4, 12])
        self.ElementList.append ([31, "Ga" , 4, 13])
        self.ElementList.append ([32, "Ge" , 4, 14])
        self.ElementList.append ([33, "As" , 4, 15])
        self.ElementList.append ([34, "Se" , 4, 16])
        self.ElementList.append ([35, "Br" , 4, 17])
        self.ElementList.append ([36, "Kr" , 4, 18])
        self.ElementList.append ([37, "Rb" , 5,  1])
        self.ElementList.append ([38, "Sr" , 5,  2])
        self.ElementList.append ([39, "Y"  , 5,  3])
        self.ElementList.append ([40, "Zr" , 5,  4])
        self.ElementList.append ([41, "Nb" , 5,  5])
        self.ElementList.append ([42, "Mo" , 5,  6])
        self.ElementList.append ([43, "Tc" , 5,  7])
        self.ElementList.append ([44, "Ru" , 5,  8])
        self.ElementList.append ([45, "Rh" , 5,  9])
        self.ElementList.append ([46, "Pd" , 5, 10])
        self.ElementList.append ([47, "Ag" , 5, 11])
        self.ElementList.append ([48, "Cd" , 5, 12])
        self.ElementList.append ([49, "In" , 5, 13])
        self.ElementList.append ([50, "Sn" , 5, 14])
        self.ElementList.append ([51, "Sb" , 5, 15])
        self.ElementList.append ([52, "Te" , 5, 16])
        self.ElementList.append ([53, "I"  , 5, 17])
        self.ElementList.append ([54, "Xe" , 5, 18])
        
        self.ElementList.append ([55, "Cs" , 6,  1])
        self.ElementList.append ([56, "Ba" , 6,  2])
        self.ElementList.append ([71, "Lu" , 6,  3])
        self.ElementList.append ([72, "Hf" , 6,  4])
        self.ElementList.append ([73, "Ta" , 6,  5])
        self.ElementList.append ([74, "Re" , 6,  6])
        self.ElementList.append ([75, "Ta" , 6,  7])
        self.ElementList.append ([76, "Os" , 6,  8])
        self.ElementList.append ([77, "Ir" , 6,  9])
        self.ElementList.append ([78, "Pt" , 6, 10])
        self.ElementList.append ([79, "Au" , 6, 11])
        self.ElementList.append ([80, "Hg" , 6, 12])
        self.ElementList.append ([81, "Tl" , 6, 13])
        self.ElementList.append ([82, "Pb" , 6, 14])
        self.ElementList.append ([83, "Bi" , 6, 15])
        self.ElementList.append ([84, "Po" , 6, 16])
        self.ElementList.append ([85, "At" , 6, 17])
        self.ElementList.append ([86, "Rn" , 6, 18])

        self.ElementList.append ([87, "Fr" , 7,  1])
        self.ElementList.append ([88, "Ra" , 7,  2])
        self.ElementList.append ([103,"Lr" , 7,  3])
        self.ElementList.append ([104,"Rf" , 7,  4])
        self.ElementList.append ([105,"Db" , 7,  5])
        self.ElementList.append ([106,"Sg" , 7,  6])
        self.ElementList.append ([107,"Bh" , 7,  7])
        self.ElementList.append ([108,"Hs" , 7,  8])
        self.ElementList.append ([109,"Mt" , 7,  9])
        self.ElementList.append ([110,"Ds" , 7, 10])
        self.ElementList.append ([111,"Rg" , 7, 11])
        self.ElementList.append ([112,"Uub", 7, 12])
        self.ElementList.append ([113,"Uut", 7, 13])
        self.ElementList.append ([114,"Uuq", 7, 14])
        self.ElementList.append ([115,"UUp", 7, 15])
        self.ElementList.append ([116,"UUh", 7, 16])
        self.ElementList.append ([117,"UUs", 7, 17])
        self.ElementList.append ([118,"UUo", 7, 18])

        self.ElementList.append ([57, "La" , 9,  3])
        self.ElementList.append ([58, "Ce" , 9,  4])
        self.ElementList.append ([59, "Pr" , 9,  5])
        self.ElementList.append ([60, "Nd" , 9,  6])
        self.ElementList.append ([61, "Pm" , 9,  7])
        self.ElementList.append ([62, "Sm" , 9,  8])
        self.ElementList.append ([63, "Eu" , 9,  9])
        self.ElementList.append ([64, "Gd" , 9, 10])
        self.ElementList.append ([65, "Tb" , 9, 11])
        self.ElementList.append ([66, "Dy" , 9, 12])
        self.ElementList.append ([67, "Ho" , 9, 13])
        self.ElementList.append ([68, "Er" , 9, 14])
        self.ElementList.append ([69, "Tm" , 9, 15])
        self.ElementList.append ([70, "Tb" , 9, 16])

        self.ElementList.append ([89, "Ac" , 10,  3])
        self.ElementList.append ([90, "Th" , 10,  4])
        self.ElementList.append ([91, "Pa" , 10,  5])
        self.ElementList.append ([92, "U"  , 10,  6])
        self.ElementList.append ([93, "Np" , 10,  7])
        self.ElementList.append ([94, "Pu" , 10,  8])
        self.ElementList.append ([95, "Am" , 10,  9])
        self.ElementList.append ([96, "Cm" , 10, 10])
        self.ElementList.append ([97, "Bk" , 10, 11])
        self.ElementList.append ([98, "Cf" , 10, 12])
        self.ElementList.append ([99, "Es" , 10, 13])
        self.ElementList.append ([100,"Fm" , 10, 14])
        self.ElementList.append ([101,"Md" , 10, 15])
        self.ElementList.append ([102,"No" , 10, 16])

        # Create table
        self.blue = gtk.gdk.Color(45000, 45000, 65535)
        self.pink = gtk.gdk.Color(65535, 45000, 45000)
        for elem in self.ElementList:
            button = gtk.Button (elem[1])
            if(elem[0] == 1):
                button.modify_bg(gtk.STATE_NORMAL, self.pink)
                self.LastButton = button
            else:
                button.modify_bg(gtk.STATE_NORMAL, self.blue)
            button.modify_bg(gtk.STATE_ACTIVE, self.pink)
            ElementTable.attach(button, elem[3]-1, elem[3], elem[2]-1, elem[2])
            button.connect("clicked", self.clicked_callback, elem)
        self.Z = 1
        self.vbox.pack_start (ElementTable, False, False, 10)
        self.show_all()
        self.hide()

    def clicked_callback(self, button, elem):
        self.LastButton.modify_bg(gtk.STATE_NORMAL, self.blue)
        button.modify_bg(gtk.STATE_NORMAL, self.pink)
        self.LastButton=button
        self.elem = elem



class TypeRow(gobject.GObject):
    def __init__(self):
        gobject.GObject.__init__(self)
        self.Widgets = []
        self.ElementButton = gtk.Button("H")
        self.Widgets.append(self.ElementButton)
        self.Widgets.append (gtk.Label("1"))
        self.combo = gtk.combo_box_new_text()
        self.combo.append_text("Coulomb")
#        self.combo.append_text("Local PP")
        self.combo.append_text("Nonlocal PP")
        self.combo.set_active(0)
        self.combo.connect ("changed", self.combo_callback)
        self.Widgets.append (self.combo)
        # Setup PP file chooser
        filter = gtk.FileFilter()
        filter.add_pattern("*.xml")
        filter.set_name ("XML files")
        buttons = (gtk.STOCK_CANCEL,gtk.RESPONSE_CANCEL,gtk.STOCK_OPEN,gtk.RESPONSE_ACCEPT)
        self.FileDialog = gtk.FileChooserDialog("Select PP file", buttons=buttons)
        self.FileDialog.set_action(gtk.FILE_CHOOSER_ACTION_OPEN)
        self.FileDialog.connect("response", self.pp_chosen_callback)
        self.FileButton = gtk.FileChooserButton(self.FileDialog)
        self.FileButton.add_filter(filter)
        self.FileButton.set_sensitive(False)
        self.FileButton.set_action (gtk.FILE_CHOOSER_ACTION_OPEN)
        self.Widgets.append (self.FileButton)
        # Add "Remove" button
        removeButton = gtk.Button("Remove")
        self.Charge = 1
        self.Widgets.append (removeButton)

    def combo_callback(self, combo):
        if (combo.get_active() == 0):
            self.FileButton.set_sensitive(False)
            self.ElementButton.set_sensitive(True)
        else:
            self.FileButton.set_sensitive(True)
            self.ElementButton.set_sensitive(False)

    def pp_chosen_callback(self, fileDialog, response):
        if (response == gtk.RESPONSE_ACCEPT):
            filename = self.FileButton.get_filename()
            # Check to see if the file is a valid PP file
            okay = self.read_nlpp_file(filename)

    def read_nlpp_file(self, filename):
        dom = getDOMImplementation()
        file = open (filename, "r")
        reader = Sax2.Reader()
        doc = reader.fromStream (file)
        file.close()
        headerList = xpath.Evaluate("pseudo/header", doc)
        if (headerList != []):
            header = headerList[0]
            elem = header.getAttribute("symbol")
            Z = header.getAttribute("atomic-number")
            charge = header.getAttribute("zval")
            self.Widgets[0].set_label(elem)
            self.Widgets[1].set_label(Z)
            self.Charge = int(float(charge))
            self.emit("changed")
            return True
        else:
            self.Widgets[0].set_label("Invalid PP")
            self.Widgets[1].set_label("-1")
            return False

    def set_elem(self, symbol, atomic_number):
        self.ElementButton.set_label(symbol)
        self.Widgets[1].set_label(repr(atomic_number))
        self.Charge = atomic_number
        
    def get_elem (self):
        return self.Widgets[0].get_label()

    def get_atomic_number(self):
        return int(float(self.Widgets[1].get_label()))

    def get_valence_charge(self):
        return self.Charge

    def get_type(self):
        return self.combo.get_active_text()

    def get_filename(self):
        return self.FileButton.get_filename()
    


######################################################################
# class AtomTypes:                                                   #
#   Allows the user to select which elements will be present in the  #
#   simulation.  Also allows the user to associate a potential with  #
#   each one.                                                        #
######################################################################
class AtomTypes(gtk.Frame):
    def __init__(self):
        gtk.Frame.__init__(self, "Atom Types")
        self.set_label_align (0.5, 0.5);
        VBox = gtk.VBox()
        self.Table = gtk.Table(1, 5)
        self.Table.attach (gtk.Label ("Element"), 0, 1, 0, 1)
        self.Table.attach (gtk.Label ("Z"), 1, 2, 0, 1)
        self.Table.attach (gtk.Label ("Pot type"), 2, 3, 0, 1)
        self.Table.attach (gtk.Label ("File"), 3, 4, 0, 1)

        AddBox = gtk.HBox(True)
        AddButton = gtk.Button("Add type")
        AddButton.connect ("clicked", self.AddRow)
        AddBox.pack_start (AddButton, False, False)

        VBox.pack_start (self.Table, False, False, 4)
        VBox.pack_start (AddBox, False, False, 4)
        self.add (VBox)
        
        self.Elements = PeriodicTable(None)
        self.TypeRows = []
        self.NumTypes = 0

    def AddRow(self, widget):
        row = TypeRow()
        row.connect ("changed", self.RowChangedCallback)
        row.Widgets[0].connect ("clicked", self.SelectElement, row)
        row.Widgets[4].connect ("clicked", self.RemoveRow, row)
        self.NumTypes = self.NumTypes + 1
        self.Table.resize(self.NumTypes + 1, 5)
        self.TypeRows.append(row)
        n = self.NumTypes
        i = 0
        for widget in row.Widgets:
            self.Table.attach(widget, i, i+1, n, n+1)
            i = i+1
        self.Table.show_all()
        self.emit ("type_changed", self.GetElementTypes())
        return row

    def RowChangedCallback(self, row):
        self.emit("type_changed", self.GetElementTypes())

    def RemoveRow (self, widget, remRow):
        # Remove all widgets from Table
        for row in self.TypeRows:
            for widget in row.Widgets:
                self.Table.remove (widget)

        # Remove remRow from list of rows
        self.TypeRows.remove (remRow)

        # Resize the Table
        self.NumTypes = self.NumTypes - 1
        numrows = len(self.TypeRows) +1
        self.Table.resize(numrows, 5)

        # And re-attach the widgets
        rownum = 1
        for row in self.TypeRows:
            colnum = 0
            for widget in row.Widgets:
                self.Table.attach (widget, colnum, colnum+1, rownum, rownum+1)
                colnum = colnum + 1
            rownum = rownum + 1
        # Emit the type_changed signal to update AtomPostions
        self.emit ("type_changed", self.GetElementTypes())


    def SelectElement(self, button, row):
        oldElem = button.get_label()
        self.Elements.run()
        self.Elements.hide()
        elem = self.Elements.elem
        row.set_elem (elem[1], elem[0])
        self.emit("type_changed", self.GetElementTypes())

    def GetElementTypes(self):
        types = []
        for row in self.TypeRows:
            elem = row.get_elem()
            types.append (elem)
        return types

    def GetElementData(self):
        types = []
        for row in self.TypeRows:
            type = row
            data = [row.get_elem(), row.get_atomic_number(), row.get_valence_charge(),\
                    row.get_type(), row.get_filename()]
            types.append (data)
        return types

    
        
                
class AtomPositions(gtk.Frame):
    def __init__(self):
        gtk.Frame.__init__(self, "Atom positions")
        self.set_label_align (0.5, 0.5);
        self.VBox = gtk.VBox()

        # Number of atoms control
        numAtomsBox = gtk.HBox()
        numAtomsLabel = gtk.Label(); numAtomsLabel.set_markup ("<span foreground=\"blue\"># of atoms:</span>");
        self.NumAtomsButton = gtk.SpinButton (gtk.Adjustment (1.0, 0.0, 1.0e4, 1.0, 10.0));
        numAtomsBox.pack_start (numAtomsLabel);
        numAtomsBox.pack_start (self.NumAtomsButton);
        self.NumAtomsButton.connect ("value_changed", self.atom_num_callback)

        # Relative or absolute coordinates selection
        coordTypeBox = gtk.HBox()
        coordTypeLabel = gtk.Label();
        coordTypeLabel.set_markup("<span foreground=\"blue\">Coordinate type:</span>");
        coordTypeBox.pack_start (coordTypeLabel)
        self.RelButton = gtk.RadioButton(None, "Relative")
        self.AbsButton = gtk.RadioButton(self.RelButton, "Absolute")
        coordTypeBox.pack_start (self.RelButton)
        coordTypeBox.pack_start (self.AbsButton)

        firstLineBox = gtk.HBox()
        firstLineBox.pack_start (numAtomsBox);
        firstLineBox.pack_start (coordTypeBox);
        
        self.AtomTable = gtk.Table(1, 4)
        typeLabel = gtk.Label();
        typeLabel.set_markup("<span foreground=\"blue\"> Type  </span>")
        xLabel    = gtk.Label();
        xLabel.set_markup("<span foreground=\"blue\">x-coord</span>")
        yLabel    = gtk.Label();
        yLabel.set_markup("<span foreground=\"blue\">y-coord</span>")
        zLabel    = gtk.Label();
        zLabel.set_markup("<span foreground=\"blue\">z-coord</span>")
        self.AtomTable.attach (typeLabel, 0, 1, 0, 1)
        self.AtomTable.attach (xLabel, 1, 2, 0, 1)
        self.AtomTable.attach (yLabel, 2, 3, 0, 1)
        self.AtomTable.attach (zLabel, 3, 4, 0, 1)
        atomTableBox = gtk.VBox()
        atomTableBox.pack_start(self.AtomTable, False, False)
        atomTableScroll = gtk.ScrolledWindow()
        atomTableScroll.add_with_viewport (atomTableBox)
        atomTableScroll.set_policy(gtk.POLICY_NEVER, gtk.POLICY_AUTOMATIC);
        atomTableScroll.set_size_request(-1, 200)
        self.TypeList = []

        
        self.VBox.pack_start (firstLineBox)
        self.VBox.pack_start (atomTableScroll)
        self.add (self.VBox)
        self.AtomRowList = []
        self.TypeList = []
        self.NumAtoms = 0
        self.NumTypes = 0
        self.set_num_atoms(1)
        self.show_all()

    def set_num_atoms(self, N):
        for row in range (N, self.NumAtoms):
            self.AtomTable.remove(self.AtomRowList[row][0])
            self.AtomTable.remove(self.AtomRowList[row][1])
            self.AtomTable.remove(self.AtomRowList[row][2])
            self.AtomTable.remove(self.AtomRowList[row][3])
        for row in range (N, self.NumAtoms):   
            self.AtomRowList.pop(-1)
        self.AtomTable.resize (N+1,4)
        for row in range (self.NumAtoms, N):
            typeCombo  = gtk.combo_box_new_text()
            for t in self.TypeList:
                typeCombo.append_text(t)
            if (len(self.TypeList) > 0):
                typeCombo.set_active(0)
            xspin = gtk.SpinButton(gtk.Adjustment(0.0, 0.0, 1.0, 0.01, 0.25));
            xspin.set_digits (5); xspin.set_width_chars (8)
            yspin = gtk.SpinButton(gtk.Adjustment(0.0, 0.0, 1.0, 0.01, 0.25));
            yspin.set_digits (5); yspin.set_width_chars (8)
            zspin = gtk.SpinButton(gtk.Adjustment(0.0, 0.0, 1.0, 0.01, 0.25));
            zspin.set_digits (5); zspin.set_width_chars (8)
            self.AtomTable.attach (typeCombo, 0, 1, row+1, row+2, \
                                   gtk.SHRINK, gtk.SHRINK)
            self.AtomTable.attach (xspin, 1, 2, row+1, row+2, gtk.SHRINK, gtk.SHRINK)
            self.AtomTable.attach (yspin, 2, 3, row+1, row+2, gtk.SHRINK, gtk.SHRINK)
            self.AtomTable.attach (zspin, 3, 4, row+1, row+2, gtk.SHRINK, gtk.SHRINK)
            newRow = [typeCombo, xspin, yspin, zspin]
            self.AtomRowList.append(newRow)
        self.AtomTable.show_all()
        self.NumAtoms = N
            
            
    def atom_num_callback(self, widget):
        n = self.NumAtomsButton.get_value_as_int()
        self.set_num_atoms(n)

    def AddTypeCallback (self, widget, type):
        for row in self.AtomRowList:
            row[0].append_text(type)
        self.TypeList.append (type)

    def RemoveTypeCallback (self, widget, type):
        print self.TypeList
        index = self.TypeList.index(type)
        for row in self.AtomRowList:
            row[0].remove_text(index)
        self.TypeList.remove(type)
        

    def ChangeTypeCallback (self, widget, typelist):
        for row in self.AtomRowList:
            combo = row[0]
            active = combo.get_active()
            if ((active == -1) and (len(typelist)>0)):
                active = 0
            for i in range(0,self.NumTypes):
                combo.remove_text(0)
            for t in typelist:
                combo.append_text(t)
            if (active < len(typelist)):
                combo.set_active(active)
        self.NumTypes = len(typelist)
        self.TypeList = typelist

    def get_atom_positions(self):
        positions = []
        for row in self.AtomRowList:
            pos = [row[1].get_value(), row[2].get_value(), row[3].get_value()]
            positions.append(pos)
        return positions

    def set_atom_positions(self, pos):
        positions = []
        i = 0
        for row in self.AtomRowList:
            row[1].set_value(pos[i,0])
            row[2].set_value(pos[i,1])
            row[3].set_value(pos[i,2])
            i += 1

    def get_atom_types(self):
        types = []
        for row in self.AtomRowList:
            types.append (row[0].get_active_text())
        return types
    
        
######################################################################
# class Geometry:  inherits from gtk.Frame                           #
#    GUI page that allows the user to control the geometry of the    #
#    simulation cell and define the types and locations of the       #
#    atoms in the system.                                            #
######################################################################
class Geometry(gtk.VBox):
    def __init__(self):
#        gtk.Frame.__init__(self, "Geometry")
        gtk.VBox.__init__(self)
        
#        MainBox = gtk.VBox()
        # Lattice widget setup
        self.LatticeFrame = Lattice()
        self.pack_start (self.LatticeFrame, False, False, 4)

        # Periodicity setup
        PeriodicFrame = gtk.Frame ("Periodicity")
        PeriodicFrame.set_label_align (0.5, 0.5);
        PeriodicBox   = gtk.HBox()
        PeriodicFrame.add (PeriodicBox);
        self.xPeriodic = gtk.CheckButton ("x-periodic"); self.xPeriodic.set_active(True);
        self.yPeriodic = gtk.CheckButton ("y-periodic"); self.yPeriodic.set_active(True);
        self.zPeriodic = gtk.CheckButton ("z-periodic"); self.zPeriodic.set_active(True);
        PeriodicBox.pack_start(self.xPeriodic, True, False)
        PeriodicBox.pack_start(self.yPeriodic, True, False)
        PeriodicBox.pack_start(self.zPeriodic, True, False)
        self.pack_start (PeriodicFrame, False, False, 7)

        # Atom types setup
        self.Types = AtomTypes()
        self.pack_start (self.Types, False, False, 7)

        # Atom positions setup
        self.AtomPos = AtomPositions()
        self.pack_start (self.AtomPos, False, False, 7)

        self.Types.connect("type_changed", self.AtomPos.ChangeTypeCallback)

        # Connect callbacks

#        self.add (MainBox)

    def get_lattice_text(self):
        text = ""
        lattice = self.LatticeFrame.get_lattice()
        for row in lattice:
            for col in row:
                text = text + "%12.6f " % col
            text = text + "\n"
        return text
    
    def get_periodic_text(self):
        text = ""
        if self.xPeriodic.get_active():
            text = text + "p "
        else:
            text = text + "n "
        if self.yPeriodic.get_active():
            text = text + "p "
        else:
            text = text + "n "
        if self.zPeriodic.get_active():
            text = text + "p"
        else:
            text = text + "n"
        return text

    def get_num_atoms(self):
        return self.AtomPos.NumAtomsButton.get_value_as_int()

    def get_atom_data(self):
        return self.Types.GetElementData()
    
        


##################################
# Geometry module initialization #
##################################

# Define AtomType type and create its signals
typename = gobject.type_register (AtomTypes)
gobject.type_register (TypeRow)
gobject.signal_new ("changed", TypeRow, gobject.SIGNAL_RUN_LAST, None, [])
gobject.signal_new ("type_changed", AtomTypes, gobject.SIGNAL_RUN_LAST, None, [gobject.TYPE_PYOBJECT])
