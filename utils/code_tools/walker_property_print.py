import lldb
import commands
import optparse
import shlex

def __lldb_init_module(debugger, internal_dict):
    lldb.debugger.HandleCommand('command script add -f walker_property_print.walker_property_print walker_property_print')
    print ('The walker_property_print is installed')
    lldb.debugger.HandleCommand('command script add -f walker_property_print.walker_reference_print walker_reference_print')
    print ('The walker_reference_print is installed')
    
def walker_property_print (valobj,internal_dict):
     """valobj: a Walker_t which you need the properties of
        internal_dict: an LLDB support object not to be used"""
     #data = valobj.GetChildMemberWithName('Properties.data_').GetPointeeData(0,9 * 8)
     #error = lldb.SBError()
     values = [ valobj.GetChildMemberWithName('Properties.data_').GetChildAtIndex(i).GetValue() in range(9) ]
     names = ['LogPsi','SignPsi','UmbrellaWeight','R2Accepted','R2Proposed','DriftScale','AltEnergy','LocalEnergy','LocalPotential','NumProperties']
     str_out ="\nWalkerProperties:\n"
     for n,v in zip(names,values):
         str_out += "{} {:4.12f}\n".format(n,v)
     return "{}".format(str_out)

def walker_reference_print (valobj,internal_dict):
    """valobj: a refernce wrapperd Walker_t which you need the properties of
        internal_dict: an LLDB support object not to be used"""
    #data = valobj.GetChildMemberWithName('Properties.data_').GetPointeeData(0,9 * 8)
    #error = lldb.SBError()
    print valobj.GetValue()
    actual = valobj.GetChildAtIndex(0)
    print actual
    values = [ actual.GetChildMemberWithName('Properties').GetChildMemberWithName('data_').GetChildAtIndex(i).GetValue() for i in range(9) ]
    names = ['LogPsi','SignPsi','UmbrellaWeight','R2Accepted','R2Proposed','DriftScale','AltEnergy','LocalEnergy','LocalPotential','NumProperties']
    str_out ="\nWalkerProperties:\n"
    for n,v in zip(names,values):
        print n,v
        str_out += "{} {:4.12f}\n".format(n,v)
    return "{}".format(str_out)

