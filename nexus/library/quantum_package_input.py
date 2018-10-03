##################################################################
##  (c) Copyright 2018-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  quantum_package_input.py                                          #
#    Supports I/O for Quantum Package input data.                    #
#                                                                    #
#  Content summary:                                                  #
#    QuantumPackageInput                                             #
#      SimulationInput class for Quantum Package.                    #
#                                                                    #
#    generate_quantum_package_input                                  #
#      User-facing function to create arbitrary input.               #
#====================================================================#


from simulation import SimulationInput



class QuantumPackageInput(SimulationInput):
    def __init__(self,filepath=None):
        if filepath!=None:
            self.read(filepath)
        #end if
    #end def __init__


    def read_text(self,text,filepath=None):
        self.not_implemented()
    #end def read_text


    def write_text(self,filepath=None):
        self.not_implemented()
    #end def write_text


    def incorporate_system(self,system):
        self.not_implemented()
    #end def incorporate_system
#end class QuantumPackageInput



def generate_quantum_package_input(**kwargs):
    
    qpi = QuantumPackageInput()

    return qpi
#end def generate_quantum_package_input
