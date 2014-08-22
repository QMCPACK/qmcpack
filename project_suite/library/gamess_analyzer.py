
import os
from numpy import array,ndarray,abs
from generic import obj
from developer import DevBase
from debug import *
from simulation import SimulationAnalyzer
from gamess_input import GamessInput


class GamessAnalyzer(SimulationAnalyzer):
    def __init__(self,arg0=None):
        return

        self.path  = None
        self.input = None
        infile = None
        if isinstance(arg0,Simulation):
            sim = arg0
            infile = os.path.join(sim.locdir,sim.infile)
        else:
            infile = arg0
        #end if
        if infile!=None:
            self.path = os.path.dirname(infile)
            self.input = GamessInput(infile)
        #end if
    #end def __init__


    def analyze(self):
        None
    #end def analyze
#end class GamessAnalyzer
