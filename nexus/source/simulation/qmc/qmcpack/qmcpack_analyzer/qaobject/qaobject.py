
from generic import obj
from project_base import Pobj
from plotter import Plotter

class QAobject(Pobj):

    _global = obj()
    _global.dynamic_methods_objects=[]

    plotter = Plotter()

    opt_methods = set(['opt','linear','cslinear'])

    def __init__(self):
        return
    #end def __init__

    @staticmethod
    def condense_name(name):
        return name.strip().lower().replace(' ','_').replace('-','_').replace('__','_')
    #end def condense_name


    def _register_dynamic_methods(self):
        QAobject._global.dynamic_methods_objects.append(self)
        return
    #end def _register_dynamic_methods

    def _unlink_dynamic_methods(self):
        for o in QAobject._global.dynamic_methods_objects:
            o._unset_dynamic_methods()
        #end for
        return
    #end def _unlink_dynamic_methods

    def _relink_dynamic_methods(self):
        for o in QAobject._global.dynamic_methods_objects:
            o._reset_dynamic_methods()
        #end for
        return
    #end def _relink_dynamic_methods
#end class QAobject
