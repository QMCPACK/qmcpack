##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


import code

def ci(locs,globs):
    code.interact(local=dict(globs,**locs))
#end def ci

ls = locals
gs = globals
interact = ci



