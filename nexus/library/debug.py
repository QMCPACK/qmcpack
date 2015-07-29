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

def ci(locs,globs):
    code.interact(local=dict(globs,**locs))
#end def ci

ls = locals
gs = globals
interact = ci



