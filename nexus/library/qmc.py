##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  qmc.py                                                            #
#    Intended for eventual unified descriptions of QMC input files.  #
#                                                                    #
#    Implementation incomplete.                                      #
#                                                                    # 
#====================================================================#


from project_base import Pobj
import qmcpack


class HamiltonianDescriptor(Pobj):
    None
#end class HamiltonianDescriptor

class WavefunctionDescriptor(Pobj):
    None
#end class WavefunctionDescriptor

class QMCDescriptor(Pobj):
    #general description of a qmc simulation
    #  information transferred to qmcpack or casino input
    None
#end class QMCDescriptor

