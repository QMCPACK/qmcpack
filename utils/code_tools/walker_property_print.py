import lldb
import commands
import optparse
import shlex

def __lldb_init_module(debugger, internal_dict):
    lldb.debugger.HandleCommand('command script add -f walker_property_print.walker_property_print walker_property_print')
    print ('The walker_property_print is installed')
    
def walker_property_print (valobj,internal_dict):
     """valobj: a Walker_t which you need the properties of
        internal_dict: an LLDB support object not to be used"""
     data = valobj.GetValueForExpressionPath('.Properties.X.X').GetPointeeData(0,9 * 8)
     error = lldb.SBError()
     values = [ data.GetDouble(error,i * 8) for i in range(9) ]
     names = ['LogPsi','SignPsi','UmbrellaWeight','R2Accepted','R2Proposed','DriftScale','AltEnergy','LocalEnergy','LocalPotential','NumProperties']
     str_out ="\nWalkerProperties:\n"
     for n,v in zip(names,values):
         str_out += "{} {:4.12f}\n".format(n,v)
     return "{}".format(str_out)
 
