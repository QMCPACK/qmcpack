
from generic import obj
from developer import DevBase
from project_base import Pobj
from plotter import Plotter


class QAobj_base(Pobj):
    None
#end class QAobj_base


class QAobject(QAobj_base):

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


    _allowed_settings = set(['optimize'])
    _default_settings = obj(
        optimize = 'variance'
        )
    for var,val in _default_settings.iteritems():
        QAobj_base.__dict__[var] = val
    #end for

    @classmethod
    def settings(cls,**kwargs):
        vars = set(kwargs.keys())
        invalid = vars-cls._allowed_settings
        if len(invalid)>0:
            allowed = list(cls._allowed_settings)
            allowed.sort()
            invalid = list(invalid)
            invalid.sort()
            cls.class_error('attempted to set unknown variables\n  unknown variables: {0}\n  valid options are: {1}'.format(invalid,allowed))
        #end if
        for var,val in kwargs.iteritems():
            QAobj_base.__dict__[var] = val
        #end for
    #end settings
#end class QAobject



class Checks(Pobj):
    def __init__(self,label=''):
        self._label = label
        self._exclusions = set()
    #end def __init__

    def exclude(self,value):
        self._exclusions.add(value)
    #end def exclude

    def valid(self):
        valid = True
        for name,value in self.iteritems():
            if not (isinstance(name,str) and name.startswith('_')):
                if not value in self._exclusions:
                    valid = valid and value
                #end if
            #end if
        #end if
        self._valid = valid
        return valid
    #end def valid

    def write(self,pad=''):
        pad2 = pad+'  '
        if not '_valid' in self:
            self.valid()
        #end if
        valid = self._valid
        if valid:
            self.log(pad+self._label+' is valid')
        else:
            self.log(pad+self._label+' is invalid')
            for name,value in self.iteritems():
                if not (isinstance(name,str) and name.startswith('_')):
                    if value in self._exclusions:
                        self.log(pad2+name+' could not be checked')
                    elif value:
                        self.log(pad2+name+' is valid')
                    else:
                        self.log(pad2+name+' is invalid')
                    #end if
                #end if
            #end for
        #end if
    #end def write
#end class Checks
