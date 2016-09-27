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

class RunBlock(gtk.EventBox):
    def __init__(self):
        gtk.EventBox.__init__(self)
        self.Tips = gtk.Tooltips()
        #self.Tips.set_property("delay", 500)
        self.Frame = gtk.Frame()
        self.add (self.Frame)
        self.Frame.set_label ("Run #1")
        self.n = 1
        self.NumCols = 4

        self.MainTable = gtk.Table(1,self.NumCols)
        self.ButtonBox = gtk.VBox()
        self.ButtonBox.set_homogeneous(False)
        self.UpButton     = gtk.Button ("Move up")
        self.DownButton   = gtk.Button ("Move down")
        self.RemoveButton = gtk.Button ("Remove")
        self.ButtonBox.pack_start (self.UpButton,     False, False)
        self.ButtonBox.pack_start (self.RemoveButton, False, False)
        self.ButtonBox.pack_start (self.DownButton,   False, False)
        self.MainHBox  = gtk.HBox()
        self.MainHBox.pack_start (self.MainTable)
        self.MainHBox.pack_start (self.ButtonBox, False, False, padding=20)
        
        # Create all possible widgets for the table
        self.TypeBox = gtk.VBox()
        self.TypeBox.pack_start(gtk.Label("Run type"), False, False)
        self.TypeCombo = gtk.combo_box_new_text()
        self.TypeCombo.append_text("VMC")
        self.TypeCombo.append_text("DMC")
        self.TypeCombo.append_text("RMC")
        self.TypeCombo.append_text("Optimize")
        self.TypeCombo.connect ("changed", self.type_changed)
        self.TypeBox.pack_start (self.TypeCombo, False, False)
        self.col = 0
        self.row = 0
        self.WidgetList = []
        self.Frame.add (self.MainHBox)

        # Number of blocks
        self.BlocksBox = gtk.VBox()
        self.BlocksBox.pack_start (gtk.Label("Blocks"), False, False)
        self.BlocksButton = gtk.SpinButton(gtk.Adjustment(1.0, 1.0, 1.0e15, 1.0, 100.0))
        self.BlocksButton.set_digits(0)
        self.BlocksBox.pack_start (self.BlocksButton, False, False)

        # Steps per block
        self.StepsBox = gtk.VBox()
        self.StepsBox.pack_start (gtk.Label("Steps/block"), False, False)
        self.StepsButton = gtk.SpinButton(gtk.Adjustment(100.0, 1.0, 1.0e15, 1.0, 100.0))
        self.StepsButton.set_digits(0)
        self.StepsBox.pack_start (self.StepsButton, False, False)

        # Time step
        self.TimeStepBox = gtk.VBox()
        self.TimeStepBox.pack_start (gtk.Label("Time Step"), False, False)
        self.TimeStepButton = gtk.SpinButton(gtk.Adjustment(0.01, 0.0,100.0, 1.0e-4, 1.0e-2))
        self.TimeStepButton.set_digits(4)
        self.TimeStepBox.pack_start (self.TimeStepButton, False, False)

        # Walkers
        self.WalkersBox = gtk.VBox()
        self.WalkersLabel = gtk.Label("Walkers")
        self.WalkersBox.pack_start (self.WalkersLabel, False, False)
        self.WalkersButton = gtk.SpinButton\
                             (gtk.Adjustment(1.0, 1.0, 1e9, 1.0, 1.0e2))
        self.WalkersButton.set_digits(0)
        self.WalkersBox.pack_start (self.WalkersButton, False, False)

        # Walker record frequency
        self.RecordBox = gtk.EventBox()
        recordVBox = gtk.VBox()
        self.RecordBox.add (recordVBox)
        self.Tips.set_tip (self.RecordBox, "Specifies whether and how often the "+\
                           "walkers are recorded to the output file.")
        self.RecordLabel = gtk.Label ("Record Walkers")
        self.RecordWalkers = gtk.CheckButton ()
        self.RecordFreq = gtk.SpinButton\
                          (gtk.Adjustment(1.0, 1.0, 1.0e6, 1.0, 10.0))
        self.RecordFreq.set_digits(0)
        recordVBox.pack_start (self.RecordLabel, False, False)
        recordHBox = gtk.HBox()
        recordHBox.pack_start (self.RecordWalkers, False, False)
        recordHBox.pack_start (self.RecordFreq,    False, False)
        recordVBox.pack_start (recordHBox,     False, False)
        self.RecordWalkers.connect ("toggled", self.record_toggled)
        self.record_toggled(None)

        # Reference energy
        self.RefEnergyBox = gtk.VBox()
        self.RefEnergyLabel = gtk.Label("Reference E")
        self.RefEnergyBox.pack_start (self.RefEnergyLabel, False, False)
        self.RefEnergyButton = gtk.SpinButton(gtk.Adjustment(-0.5, -1e10, 1e10, 0.1, 1.0e0))
        self.RefEnergyButton.set_digits(6)
        self.RefEnergyBox.pack_start (self.RefEnergyButton, False, False)
        self.RefEnergyButton.set_width_chars(10)

        # Population control
        self.PopControlBox = gtk.VBox()
        self.PopControlLabel = gtk.Label("Pop Control")
        self.PopControlBox.pack_start (self.PopControlLabel, False, False)
        self.PopControlButton = gtk.SpinButton(gtk.Adjustment(50, 0.0, 1e10, 1.0, 1.0e1))
        self.PopControlButton.set_digits(0)
        self.PopControlBox.pack_start (self.PopControlButton, False, False)

        # Nonlocal Moves
        self.NonlocalBox = gtk.VBox()
        self.NonlocalBox.pack_start(gtk.Label("Nonlocal Moves"), False, False)
        self.NonlocalMoves = gtk.CheckButton()
        self.NonlocalMoves.set_active(False)
        NLHBox = gtk.HBox()
        NLHBox.pack_start (self.NonlocalMoves, True, False)
        self.NonlocalBox.pack_start (NLHBox, False, False)

        # Reconfig Moves
        self.ReconfigBox = gtk.VBox()
        self.ReconfigBox.pack_start(gtk.Label("Reconfiguration"), False, False)
        self.Reconfig = gtk.CheckButton()
        self.Reconfig.set_active(False)
        ReconfHBox = gtk.HBox()
        ReconfHBox.pack_start (self.Reconfig, True, False)
        self.ReconfigBox.pack_start (ReconfHBox, False, False)

        # Setup type-specific widget lists
        self.VMCList      = [self.TypeBox, self.BlocksBox, self.StepsBox, self.TimeStepBox,\
                             self.WalkersBox, self.RecordBox]
        self.DMCList      = [self.TypeBox, self.BlocksBox, self.StepsBox, self.TimeStepBox,\
                             self.WalkersBox, self.RefEnergyBox, self.PopControlBox,       \
                             self.NonlocalBox, self.ReconfigBox]
        self.RMCList      = [self.TypeBox, self.BlocksBox, self.StepsBox, self.TimeStepBox]
        self.OptimizeList = [self.TypeBox, self.BlocksBox, self.StepsBox, self.TimeStepBox]

        # Just add the type widget for now
        self.TypeCombo.set_active(0)
#        self.add_widget (self.TypeBox)
        self.show_all()


    def add_widget(self, widget):
        if (self.col == self.NumCols):
            self.col = 0
            self.row = self.row+1
            self.MainTable.resize(self.row+1, self.NumCols)
        self.WidgetList.append(widget)
        self.MainTable.attach (widget, self.col, self.col+1, self.row, self.row+1,\
                               xpadding=4, ypadding=4)
        self.col = self.col + 1
        self.show_all()


    def clear(self):
        for widget in self.WidgetList:
            self.MainTable.remove(widget)
        self.MainTable.resize(self.NumCols, 1)
        self.col = 0
        self.row = 0
        self.WidgetList = []
        
        
    def type_changed(self, widget):
        if (self.TypeCombo.get_active_text() == "VMC"):
            self.set_VMC()
        elif (self.TypeCombo.get_active_text() == "DMC"):
            self.set_DMC()
        elif (self.TypeCombo.get_active_text() == "RMC"):
            self.set_RMC()
        elif (self.TypeCombo.get_active_text() == "Optimize"):
            self.set_Optimize()

    def record_toggled(self, widget):
        self.RecordFreq.set_sensitive(self.RecordWalkers.get_active())

    def set_VMC(self):
        self.clear()
        for widget in self.VMCList:
            self.add_widget(widget)
            self.WalkersLabel.set_text ("Walkers")
        blue = gtk.gdk.Color(45000, 45000, 65535)
        self.modify_bg (gtk.STATE_NORMAL, blue)
        self.modify_bg (gtk.STATE_ACTIVE, blue)
        self.RecordBox.modify_bg (gtk.STATE_NORMAL, blue)
        self.RecordBox.modify_bg (gtk.STATE_ACTIVE, blue)

    def set_DMC(self):
        self.clear()
        for widget in self.DMCList:
            self.add_widget(widget)
        self.WalkersLabel.set_text ("Target Walkers")
        green = gtk.gdk.Color(45000, 65535, 45000)
        self.modify_bg (gtk.STATE_NORMAL, green)
        self.modify_bg (gtk.STATE_ACTIVE, green)

    def set_RMC(self):
        self.clear()
        for widget in self.RMCList:
            self.add_widget(widget)
        yellow = gtk.gdk.Color(65535, 65535, 45000)
        self.modify_bg (gtk.STATE_NORMAL, yellow)
        self.modify_bg (gtk.STATE_ACTIVE, yellow)

    def set_Optimize(self):
        self.clear()
        for widget in self.OptimizeList:
            self.add_widget(widget)
        red = gtk.gdk.Color(65535, 45000, 45000)
        self.modify_bg (gtk.STATE_NORMAL, red)
        self.modify_bg (gtk.STATE_ACTIVE, red)

    def set_run_number(self, n):
        self.Frame.set_label ("Run #" + repr(n))
        self.n = n

    def get_param_list(self):
        runtype = self.TypeCombo.get_active_text()
        paramList = []
        if (runtype == "VMC"):
            paramList.append(("blocks",         repr(self.BlocksButton.get_value())))
            paramList.append(("steps",          repr(self.StepsButton.get_value())))
            paramList.append(("timestep",       repr(self.TimeStepButton.get_value())))
            paramList.append(("walkers",        repr(self.WalkersButton.get_value())))
        elif (runtype == "DMC"):
            paramList.append(("blocks",         repr(self.BlocksButton.get_value())))
            paramList.append(("steps",          repr(self.StepsButton.get_value())))
            paramList.append(("timestep",       repr(self.TimeStepButton.get_value())))
            paramList.append(("target_walkers", repr(self.WalkersButton.get_value())))
            paramList.append(("ref_energy",     repr(self.RefEnergyButton.get_value())))
            paramList.append(("pop_control",    repr(self.PopControlButton.get_value())))
            reconf = "no"
            if (self.Reconfig.get_active()):
                reconf = "yes"
            nonlocal = "no"
            if (self.NonlocalMoves.get_active()):
                nonlocal = "yes"
            paramList.append(("reconfiguration", reconf))
            paramList.append(("nonlocalmoves", nonlocal))

        return paramList
            
        

class Run(gtk.Frame):
    def __init__(self):
        gtk.Frame.__init__(self, "Run")
        self.MainHBox = gtk.HBox()
        self.RunVBox = gtk.VBox()
        self.Scroll = gtk.ScrolledWindow()
        self.Scroll.add_with_viewport (self.RunVBox)
        self.Scroll.set_policy(gtk.POLICY_NEVER, gtk.POLICY_AUTOMATIC);
        self.MainHBox.pack_start (self.Scroll)
        AddBox = gtk.VBox()
        self.AddButton = gtk.Button ("Add Run")
        self.AddButton.connect("clicked", self.add_run)
        AddBox.pack_start (self.AddButton, False, False)
        self.MainHBox.pack_start (AddBox, False, False)
        self.RunList = []
        self.add (self.MainHBox)

    def add_run(self, widget):
        n = len(self.RunList) + 1
        runblock = RunBlock()
        runblock.set_run_number(n)
        runblock.UpButton.connect    ("clicked", self.move_up,    runblock)
        runblock.DownButton.connect  ("clicked", self.move_down,  runblock)
        runblock.RemoveButton.connect("clicked", self.remove_run, runblock)
        self.RunVBox.pack_start (runblock, False, False, 5)
        self.RunList.append(runblock)

    def clear_run_box(self):
        for r in self.RunList:
                self.RunVBox.remove(r)

    def fill_run_box(self):
        i = 1
        for r in self.RunList:
            r.set_run_number (i)
            self.RunVBox.pack_start(r, False, False, 5)
            i = i+1

    def move_up(self, widget, run):
        n = run.n - 1;
        if (n != 0):
            self.clear_run_box()
            tmp = self.RunList[n]
            self.RunList[n] = self.RunList[n-1]
            self.RunList[n-1] = tmp
            self.fill_run_box()

    def move_down(self, widget, run):
        n = run.n - 1;
        if (n != (len(self.RunList)-1)):
            self.clear_run_box()
            tmp = self.RunList[n]
            self.RunList[n] = self.RunList[n+1]
            self.RunList[n+1] = tmp
            self.fill_run_box()
            
    def remove_run(self, widget, run):
        self.clear_run_box()
        self.RunList.remove (run)
        self.fill_run_box()
        

    
        
