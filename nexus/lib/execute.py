##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  execute.py                                                        #
#    Support for control of local process execution.                 #
#                                                                    #
#  Content summary:                                                  #
#    execute                                                         #
#      Execute an arbitrary shell command within the user            #
#      environment of the local machine and wait for its completion. #
#                                                                    #
#====================================================================#


from subprocess import Popen,PIPE

from utilities import to_str

def execute(command,verbose=False,skip=False):
    out,err = '',''
    returncode = 0
    if skip:
        if verbose:
            print('Would have executed:\n  '+command)
        #end if
    else:
        if verbose:
            print('Executing:\n  '+command)
        #end if
        process = Popen(command,shell=True,stdout=PIPE,stderr=PIPE,close_fds=True)
        out,err = process.communicate()
        returncode = process.returncode
    #end if
    return to_str(out),to_str(err),returncode
#end def execute
