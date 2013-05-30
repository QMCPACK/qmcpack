
import code

def ci(locs,globs):
    code.interact(local=dict(locs,**globs))
#end def ci

ls = locals
gs = globals



