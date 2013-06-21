
import code

def ci(locs,globs):
    code.interact(local=dict(globs,**locs))
#end def ci

ls = locals
gs = globals
interact = ci



