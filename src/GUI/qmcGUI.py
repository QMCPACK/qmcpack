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


#!/bin/env python

import pygtk
pygtk.require('2.0')
import gtk
import pango
from Geometry import *
from Wavefunction import *
from Run import *
from xml.dom import getDOMImplementation
from xml.dom.ext import PrettyPrint
from os.path import basename
import os.path

viewerOkay = True
try:
    import gtk.gtkgl
    from OpenGL.GL import *
    from OpenGL.GLU import *
except ImportError:
    viewerOkay = False

if (viewerOkay):
    from StructViewer import *

class MainManager (gtk.UIManager):
    def __init__(self):
        gtk.UIManager.__init__(self)

    def Setup(self):
        FileGroup = gtk.ActionGroup("File")
        FileAction = gtk.Action("File", "_File", "File operations", None)
        self.SaveAction   = gtk.Action ("Save", "_Save", "Save qmcPACK input.", gtk.STOCK_SAVE)
        self.SaveAsAction = gtk.Action ("SaveAs", "Save _As", "Save qmcPACK input with a new name", gtk.STOCK_SAVE_AS)
        self.QuitAction   = gtk.Action("Quit", "_Quit", "Exit qmcGUI", gtk.STOCK_QUIT)
        FileGroup.add_action(FileAction)
        FileGroup.add_action(self.SaveAction)
        FileGroup.add_action(self.SaveAsAction)
        FileGroup.add_action(self.QuitAction)

        ViewGroup = gtk.ActionGroup("View")
        ViewAction = gtk.Action("View", "_View", "View operations", None)
        self.StructureAction = gtk.Action ("Structure", "_Structure", "View structure", None)
        ViewGroup.add_action(ViewAction)
        ViewGroup.add_action(self.StructureAction)

        HelpGroup  = gtk.ActionGroup("Help")
        HelpAction = gtk.Action("Help", "_Help", "User assistance", gtk.STOCK_HELP)
        self.AboutAction = gtk.Action ("About", "_About", "About qmcGUI", gtk.STOCK_ABOUT)
        HelpGroup.add_action (HelpAction);
        HelpGroup.add_action (self.AboutAction);
        
        self.insert_action_group (FileGroup, 0)
        self.insert_action_group (ViewGroup, -1)
        self.insert_action_group (HelpGroup, -1)
        descriptor = "<ui>                                    "\
                     "  <menubar name=\"MenuBar\">            "\
                     "    <menu action=\"File\">              "\
                     "      <menuitem action=\"Save\"/>       "\
                     "      <menuitem action=\"SaveAs\"/>     "\
                     "      <menuitem action=\"Quit\"/>       "\
                     "    </menu>                             "\
                     "    <menu action=\"View\">              "\
                     "      <menuitem action=\"Structure\"/>  "\
                     "    </menu>                             "\
                     "    <menu action=\"Help\">              "\
                     "      <menuitem action=\"About\"/>      "\
                     "    </menu>                             "\
                     "  </menubar>                            "\
                     "</ui>                                   "
        self.add_ui_from_string(descriptor)

class GUIAboutDialog(gtk.AboutDialog):
    def __init__(self):
        gtk.AboutDialog.__init__(self)
        self.set_name ("qmcPACK Input Builder")
        self.set_version("0.1")
        self.set_copyright("Copyright 2007, GNU Public License")
        self.set_website("http://cms.mcc.uiuc.edu/qmcpack/index.php/QMCPACK_Wiki_Home")
        self.set_authors(["Ken Esler (kesler@ciw.edu)"])
        

#############################################################
# function relativeto                                       #
# return a the name of file2 relative to the path of file1  #
# Example:                                                  #
#  file1 = "/home/kesler/BN/abc.xml"                        #
#  file2 = "/home/kesler/pseudo/B.xml"                      #
#  print relativeto (file1, file2) yields "../pseudo/B.xml" #
#############################################################

def relative2 (file1, file2):
    dir1 = os.path.dirname(file1)
    dir2 = os.path.dirname(file2)
    dirlist1 = dir1.split ('/')
    dirlist2 = dir2.split ('/')
    if (dirlist1[0] == ''):
        dirlist1.remove('')
    if (dirlist2[0] == ''):
        dirlist2.remove('')
    print "dirlist1 = " + repr(dirlist1)
    common = ""
    
    i = 0
    mindirs = min (len(dirlist1), len(dirlist2))
    while ((i < mindirs) and (dirlist1[i] == dirlist2[i])):
        common = common + '/' + dirlist1[i]
        i = i+1
    common = common 
    dirstrip  = dir1.replace(common,'',1)
    filestrip = file2.replace(common+'/','',1)
    striplist = dirstrip.split('/')
    rel = ""
    for d in striplist:
        if (d != ""):
            rel = rel + "../"
    rel = rel + filestrip
    return rel
    

def relativeto (file1, file2):
    dir = os.path.dirname (file2)
    common = os.path.commonprefix((file1, dir))
    dirstrip = dir.lstrip(common)
    filestrip = file2.lstrip(common)
    dirlist = dirstrip.split('/')
    rel = ""
    for d in dirlist:
        if (d != ""):
            rel = rel + "../"
    rel = rel + filestrip
    return rel


class GUI:
    def __init__(self):
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.set_title ('qmcPACK Input Builder')
        manager = MainManager()
        self.window.add_accel_group (manager.get_accel_group())
        manager.Setup()
        mainMenu = manager.get_widget("/MenuBar")
        mainBox = gtk.VBox()
        mainBox.pack_start (mainMenu, False, False)
        # Connect menu callbacks
        manager.SaveAsAction.connect("activate", self.save_as_callback)
        manager.SaveAction.connect("activate", self.save_callback)
        self.About = GUIAboutDialog()
        manager.AboutAction.connect("activate", self.about_callback)
        
        self.window.connect("delete_event", self.delete)
        notebook = gtk.Notebook()
        notebook.set_tab_pos (gtk.POS_LEFT)
        self.GeometryFrame = Geometry()
        self.WavefunctionFrame = Wavefunction(self.GeometryFrame)
        self.RunFrame = Run()
        notebook.append_page (self.GeometryFrame, gtk.Label("Geometry"))
        notebook.append_page (self.WavefunctionFrame, gtk.Label("Wave function"))
        notebook.append_page (self.RunFrame, gtk.Label("Run"))
        if (viewerOkay):
            self.Viewer = StructureViewer()
            notebook.append_page(self.Viewer, gtk.Label("Viewer"))
        mainBox.pack_start (notebook, False, False, 5)
        self.window.add (mainBox)
        # Setup dialogs
        buttons = (gtk.STOCK_CANCEL,gtk.RESPONSE_CANCEL,gtk.STOCK_SAVE,gtk.RESPONSE_OK)
        self.SaveAsDialog = gtk.FileChooserDialog("Save qmcPACK input as", \
                                                  self.window, gtk.FILE_CHOOSER_ACTION_SAVE, buttons)
        self.SaveAsDialog.hide()
        self.Filename = None
        self.window.show_all()

    def about_callback (self, action):
        self.About.run()
    
    def save_as_callback(self, action):
        if (self.SaveAsDialog.run() == gtk.RESPONSE_OK):
            self.Filename = self.SaveAsDialog.get_filename()
            self.write_qmcPACK(self.Filename)
            self.window.set_title (basename(self.Filename))
        self.SaveAsDialog.hide()
                

    def save_callback(self, action):
        if (self.Filename == None):
            self.save_as_callback(action)
        else:
            self.write_qmcPACK (self.Filename)

    def write_qmcPACK(self, filename):
        impl = getDOMImplementation()
        doc = impl.createDocument(None, "qmcPACK", None)
        topElem = doc.documentElement
        # Write out structural geometry elements
        systemElem  = doc.createElement("qmcsystem")
        systemElem.setAttribute("dim", "3")
        # Simulation cell information
        cellElem    = doc.createElement ("simulationcell")
        latticeElem = doc.createElement ("parameter")
        latticeElem.setAttribute("name", "lattice")
        latticeText = doc.createTextNode("\n"+self.GeometryFrame.get_lattice_text()+"      ")
        latticeElem.appendChild (latticeText)
        bcondsElem  = doc.createElement ("bconds")
        bcondsText  = doc.createTextNode(" " + self.GeometryFrame.get_periodic_text() + " ")
        bcondsElem.appendChild (bcondsText)
        cellElem.appendChild (latticeElem)
        cellElem.appendChild (bcondsElem)
        # Atom position information
        particleSetElem = doc.createElement ("particleset")
        particleSetElem.setAttribute ("name", "i")
        particleSetElem.setAttribute ("size", repr (self.GeometryFrame.get_num_atoms()))
        systemElem.appendChild(cellElem)

        # Ion types
        ionData = self.GeometryFrame.get_atom_data()
        bareIons = False
        pseudoIons = False
        for data in ionData:
            groupElem = doc.createElement ("group")
            groupElem.setAttribute("name", data[0])
            chargeElem = doc.createElement("parameter")
            chargeElem.setAttribute("name", "charge")
            chargeElem.appendChild(doc.createTextNode(repr(data[2])))
            valenceElem = doc.createElement("parameter")
            valenceElem.setAttribute("name", "valence")
            valenceElem.appendChild(doc.createTextNode(repr(data[2])))
            zElem = doc.createElement("parameter")
            zElem.setAttribute("name", "atomicnumber")
            zElem.appendChild(doc.createTextNode(repr(data[1])))
            groupElem.appendChild(chargeElem)
            groupElem.appendChild(valenceElem)
            groupElem.appendChild(zElem)
            particleSetElem.appendChild (groupElem)
            if ((data[3]=="Nonlocal PP") or (data[3] == "LocalPP")):
                pseudoIons = True
            else:
                bareIons = True
        # Add atom positions
        posElem = doc.createElement("attrib")
        posElem.setAttribute("name", "position")
        posElem.setAttribute("datatype", "posArray")
        posElem.setAttribute("condition", "1")
        positions = self.GeometryFrame.AtomPos.get_atom_positions()
        posText = "\n"
        for pos in positions:
            rowText = "    " + ("%12.6f %12.6f %12.6f" % (pos[0], pos[1], pos[2])) + "\n"
            posText = posText + rowText
        posText = posText + "      "
        posElem.appendChild(doc.createTextNode(posText))
        particleSetElem.appendChild(posElem)
        # Add atom types
        idElem = doc.createElement("attrib")
        idElem.setAttribute ("name", "ionid")
        idElem.setAttribute ("datatype", "stringArray")
        idText = "\n        "
        ids = self.GeometryFrame.AtomPos.get_atom_types()
        for id in ids:
            if (id != None):
                idText = idText + id + " "
        idText = idText + "\n      "
        idElem.appendChild(doc.createTextNode(idText))
        particleSetElem.appendChild(idElem)
        systemElem.appendChild(particleSetElem)

        # Add hamiltonian section
        hamElem = doc.createElement("hamiltonian")
        hamElem.setAttribute("name", "h0")
        hamElem.setAttribute("type", "generic")
        hamElem.setAttribute("targey", "e")
        # e-e interaction
        eeElem = doc.createElement("pairpot")
        eeElem.setAttribute("name", "ElecElec")
        eeElem.setAttribute("name", "coulomb")
        eeElem.setAttribute("source", "e")
        hamElem.appendChild(eeElem)
        # e-ion interaction
        # Pseudopotentials
        if (pseudoIons):
            pseudoElem = doc.createElement("pairpot")
            pseudoElem.setAttribute("type", "pseudo")
            pseudoElem.setAttribute("name", "PseudoPot")
            pseudoElem.setAttribute("source", "i")
            pseudoElem.setAttribute("wavefunction", "psi0")
            pseudoElem.setAttribute("format", "xml")
            for data in ionData:
                if (data[3] == "Nonlocal PP"):
                    pElem = doc.createElement ("pseudo")
                    pElem.setAttribute("elementType", data[0])
                    if (data[4] == None):
                        print "Warning: nonlocal pseudopotential file not set."
                    else:
                        pElem.setAttribute("href", \
                                           relative2(filename,data[4]))
                    pseudoElem.appendChild (pElem)
            hamElem.appendChild(pseudoElem)
        if (bareIons):
            bareElem = doc.createElement("pairpot")
            bareElem.setAttribute("name", "barIon")
            bareElem.setAttribute("type", "coulomb")
            hamElem.appendChild (bareElem)
        systemElem.appendChild(hamElem)
        topElem.appendChild (systemElem)

        #########################
        # Write run information #
        #########################
        for run in self.RunFrame.RunList:
            runtype = run.TypeCombo.get_active_text()
            runElem = doc.createElement ("qmc")
            runElem.setAttribute("method", runtype)
            runElem.setAttribute("target", "e")
            paramList = run.get_param_list()
            for param in paramList:
                paramElem = doc.createElement("parameter")
                paramElem.setAttribute("name", param[0])
                textElem = doc.createTextNode (param[1])
                paramElem.appendChild (textElem)
                runElem.appendChild (paramElem)
            topElem.appendChild (runElem)

        
        file = open (filename, "w")
        PrettyPrint (doc, file)
        file.close()
        

    def delete(self, widget, event=None):
        gtk.main_quit()
        return False

    def run(self):
        gtk.main()

if __name__ == "__main__":
    
    gui = GUI()
    gui.run()
    
