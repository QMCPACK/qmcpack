
from subprocess import Popen,PIPE

def execute(command,verbose=False,skip=False):
    out,err = '',''
    if skip:
        if verbose:
            print 'Would have executed:\n  '+command
        #end if
    else:
        if verbose:
            print 'Executing:\n  '+command
        #end if
        out,err = Popen(command,shell=True,stdout=PIPE,stderr=PIPE,close_fds=True).communicate()
    #end if
    return out,err
#end def execute
