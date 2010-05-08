import pygtk
import gtk
import IO

class Orbitals(gtk.Frame):
    def h5_chosen_callback(self, fileDialog, response):
        if (response == gtk.RESPONSE_ACCEPT):
            filename = self.FileButton.get_filename()
            # Check to see if the file is a valid PP file
            okay = self.read_h5_file(filename)

    def read_eshdf (self, io):
        io.OpenSection('supercell')
        self.prim_vecs = io.ReadVar('primitive_vectors')
        io.CloseSection()
        io.OpenSection('atoms')

        io.CloseSection()

    def read_h5_file(self, filename):
        io = IO.IOSectionClass()
        if (not io.OpenFile(filename)):
            return False
        format = io.ReadVar('format')
        if (format == 'ES-HDF'):
            return self.read_eshdf (io)
        
        return False

    def __init__(self):
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
    def __init__(self):
        gtk.VBox.__init__(self)
        self.OrbitalsFrame = Orbitals()
        self.pack_start (self.OrbitalsFrame, False, False, 4)

        self.JastrowsFrame = Jastrows()
        self.pack_start (self.JastrowsFrame, False, False, 4)


        # self.Widgets.append (self.FileButton)
