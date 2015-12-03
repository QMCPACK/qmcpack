##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  debug.py                                                          #
#    User/developer interface to insert interactive breakpoints.     #
#    Useful for rapid development and debugging.                     #
#                                                                    #
#  Content summary:                                                  #
#    ci                                                              #
#      Function to enter interactive environment.                    #
#      Usage: ci(ls(),gs()).                                         #
#                                                                    #
#    ls                                                              #
#      Alias for standard locals function.                           #
#                                                                    #
#    gs                                                              #
#      Alias for standard globals function.                          #
#                                                                    #                                        
#====================================================================#


import code
import inspect

def ci(locs=None,globs=None):
    if locs is None or globs is None:
        cur_frame = inspect.currentframe()
        caller_frame = cur_frame.f_back
        locs  = caller_frame.f_locals
        globs = caller_frame.f_globals
    #end if
    code.interact(local=dict(globs,**locs))
#end def ci

ls = locals
gs = globals
interact = ci



