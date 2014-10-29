

from generic import obj


class DevBase(obj):
    user_interface_data = obj()
    dev_instructions_data = obj()

    def not_implemented(self):
        #self.error('a base class function has not been implemented','Developer')
        self.error('a base class function has not been implemented',trace=True)
    #end def not_implemented

    @classmethod
    def set_user_interface(cls,class_variables=None,class_functions=None,member_variables=None,member_functions=None):
        if class_variables is None:
            class_variables = []
        if class_functions is None:
            class_functions = []
        if member_variables is None:
            member_variables = []
        if member_functions is None:
            member_functions = []
        ui = cls.user_interface_data
        ui.class_variables = class_variables
        ui.class_functions = class_functions
        ui.member_variables= member_variables
        ui.member_functions= member_functions
    #end def set_user_interface

    @classmethod
    def set_dev_instruction(cls,situation='writing a derived class',class_variables=None,class_functions=None,member_variables=None,member_functions=None):
        if class_variables is None:
            class_variables = []
        if class_functions is None:
            class_functions = []
        if member_variables is None:
            member_variables = []
        if member_functions is None:
            member_functions = []
        ins = obj()
        cls.dev_instructions_data[situation] = ins
        ins.class_variables = class_variables
        ins.class_functions = class_functions
        ins.member_variables= member_variables
        ins.member_functions= member_functions
    #end def set_dev_instruction
        
    @classmethod
    def write_dev_data(cls,s,dev_data,indent=''):
        p1 = '  '
        p2 = 2*p1
        p3 = 3*p1
        dpad = indent+p3
        for name,value in dev_data.iteritems():
            s+=indent+p1+name+'\n'
            if isinstance(value,list):
                for v in value:
                    s+=indent+p2+v+'\n'
                #end for
            else:
                for v,description in value.iteritems():
                    s+=indent+p2+v+'\n'
                    s+=dpad+description.replace('\n','\n'+dpad)
                #end for
            #end for
        #end for
        return s
    #end def write_dev_data

    @classmethod
    def user_interface(cls):
        s='User Interface for '+cls.__name__+' Class\n'
        s+= (len(s)-1)*'='+'\n'
        s+=cls.write_dev_data(cls.user_interface_data)
        return s
    #end def user_interface

    @classmethod
    def developer_instructions(cls):
        s='Developer Instructions for '+cls.__name__+' Class\n'
        s+= (len(s)-1)*'='+'\n'
        i1='  '
        i2=2*i1
        i3=3*i1
        for situation,instruction in cls.developer_instructions.iteritems():
            s+=i1+'situation: '+situation+'\n'
            s+=i2+'define the following variables and functions:'
            s+=cls.write_dev_data(instructions,i3)
        #end for
        return s
    #end def developer_instructions
#end class DevBase






class Void:
    void_items = dict()

    @classmethod
    def _unavailable(cls,self):
        sid = id(self)
        if sid in Void.void_items:
            module,item = Void.void_items[id(self)]
        else:
            module,item = None,None
        #end if
        if module is None and item is None:
            msg = 'encountered a void item from an unavailable module'
        elif module is None:
            msg = 'item '+str(item)+' is from an unavailable module'
        elif item is None:
            msg = 'encountered a void item from unavailable module '+str(module)+'  \nthis python module must be installed on your system to use this feature'
        else:
            msg = 'item '+str(item)+' is from unavailable module '+str(module)+'  \nthis python module must be installed on your system to use this feature'
        #end if
        DevBase.class_error(msg,'Void')
    #end def _unavailable

        
    @classmethod
    def _class_unavailable(cls):
        msg = 'encountered a void item from an unavailable module'
        DevBase.class_error(msg,'Void')
    #end def _class_unavailable
        


    def __init__(self,module=None,item=None):
        Void.void_items[id(self)] = module,item
    #end def __init__


    #list of magic functions taken from the following sources
    #  http://web.archive.org/web/20110131211638/http://diveintopython3.org/special-method-names.html
    #  http://www.ironpythoninaction.com/magic-methods.html


    #class methods
    @classmethod
    def __instancecheck__(cls,*args,**kwargs):
        Void._class_unavailable()
    def __subclasscheck__(cls,*args,**kwargs):
        Void._class_unavailable()
    def __subclasshook__(cls,*args,**kwargs):
        Void._class_unavailable()
    

    #member methods
    def __new__(self,*args,**kwargs):
        Void._unavailable(self)
    def __eq__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ne__(self,*args,**kwargs):
        Void._unavailable(self)
    def __lt__(self,*args,**kwargs):
        Void._unavailable(self)
    def __le__(self,*args,**kwargs):
        Void._unavailable(self)
    def __gt__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ge__(self,*args,**kwargs):
        Void._unavailable(self)
    def __nonzero__(self,*args,**kwargs):
        Void._unavailable(self)
    def __subclasses__(self,*args,**kwargs):
        Void._unavailable(self)
    def __call__(self,*args,**kwargs):
        Void._unavailable(self)
    def __hash__(self,*args,**kwargs):
        Void._unavailable(self)
    #def __del__(self,*args,**kwargs):
    #    Void._unavailable(self)
    def __dir__(self,*args,**kwargs):
        Void._unavailable(self)
    def __getitem__(self,*args,**kwargs):
        Void._unavailable(self)
    def __setitem__(self,*args,**kwargs):
        Void._unavailable(self)
    def __delitem__(self,*args,**kwargs):
        Void._unavailable(self)
    def __len__(self,*args,**kwargs):
        Void._unavailable(self)
    def __contains__(self,*args,**kwargs):
        Void._unavailable(self)
    def __iter__(self,*args,**kwargs):
        Void._unavailable(self)
    def __reversed__(self,*args,**kwargs):
        Void._unavailable(self)
    def __missing__(self,*args,**kwargs):
        Void._unavailable(self)
    def __length_hint__(self,*args,**kwargs):
        Void._unavailable(self)
    def __repr__(self,*args,**kwargs):
        Void._unavailable(self)
    def __str__(self,*args,**kwargs):
        Void._unavailable(self)
    def __unicode__(self,*args,**kwargs):
        Void._unavailable(self)
    def __getattr__(self,*args,**kwargs):
        Void._unavailable(self)
    def __setattr__(self,*args,**kwargs):
        Void._unavailable(self)
    def __delattr__(self,*args,**kwargs):
        Void._unavailable(self)
    def __getattribute__(self,*args,**kwargs):
        Void._unavailable(self)
    def __add__(self,*args,**kwargs):
        Void._unavailable(self)
    def __sub__(self,*args,**kwargs):
        Void._unavailable(self)
    def __mul__(self,*args,**kwargs):
        Void._unavailable(self)
    def __floordiv__(self,*args,**kwargs):
        Void._unavailable(self)
    def __div__(self,*args,**kwargs):
        Void._unavailable(self)
    def __truediv__(self,*args,**kwargs):
        Void._unavailable(self)
    def __mod__(self,*args,**kwargs):
        Void._unavailable(self)
    def __divmod__(self,*args,**kwargs):
        Void._unavailable(self)
    def __pow__(self,*args,**kwargs):
        Void._unavailable(self)
    def __lshift__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rshift__(self,*args,**kwargs):
        Void._unavailable(self)
    def __and__(self,*args,**kwargs):
        Void._unavailable(self)
    def __xor__(self,*args,**kwargs):
        Void._unavailable(self)
    def __or__(self,*args,**kwargs):
        Void._unavailable(self)
    def __neg__(self,*args,**kwargs):
        Void._unavailable(self)
    def __pos__(self,*args,**kwargs):
        Void._unavailable(self)
    def __abs__(self,*args,**kwargs):
        Void._unavailable(self)
    def __invert__(self,*args,**kwargs):
        Void._unavailable(self)
    def __complex__(self,*args,**kwargs):
        Void._unavailable(self)
    def __int__(self,*args,**kwargs):
        Void._unavailable(self)
    def __long__(self,*args,**kwargs):
        Void._unavailable(self)
    def __float__(self,*args,**kwargs):
        Void._unavailable(self)
    def __oct__(self,*args,**kwargs):
        Void._unavailable(self)
    def __hex__(self,*args,**kwargs):
        Void._unavailable(self)
    def __index__(self,*args,**kwargs):
        Void._unavailable(self)
    def __enter__(self,*args,**kwargs):
        Void._unavailable(self)
    def __exit__(self,*args,**kwargs):
        Void._unavailable(self)
    def __get__(self,*args,**kwargs):
        Void._unavailable(self)
    def __set__(self,*args,**kwargs):
        Void._unavailable(self)
    def __delete__(self,*args,**kwargs):
        Void._unavailable(self)
    def __doc__(self,*args,**kwargs):
        Void._unavailable(self)
    def __dict__(self,*args,**kwargs):
        Void._unavailable(self)
    def __slots__(self,*args,**kwargs):
        Void._unavailable(self)
    def __class__(self,*args,**kwargs):
        Void._unavailable(self)
    def __bases__(self,*args,**kwargs):
        Void._unavailable(self)
    def __name__(self,*args,**kwargs):
        Void._unavailable(self)
    def __all__(self,*args,**kwargs):
        Void._unavailable(self)
    def __file__(self,*args,**kwargs):
        Void._unavailable(self)
    def __module__(self,*args,**kwargs):
        Void._unavailable(self)
    #def __metaclass__(self,*args,**kwargs):
    #    Void._unavailable(self)
    def __import__(self,*args,**kwargs):
        Void._unavailable(self)
    def __radd__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rsub__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rmul__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rtruediv__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rfloordiv__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rmod__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rdivmod__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rpow__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rlshift__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rrshift__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rand__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rxor__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ror__(self,*args,**kwargs):
        Void._unavailable(self)
    def __iadd__(self,*args,**kwargs):
        Void._unavailable(self)
    def __isub__(self,*args,**kwargs):
        Void._unavailable(self)
    def __imul__(self,*args,**kwargs):
        Void._unavailable(self)
    def __itruediv__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ifloordiv__(self,*args,**kwargs):
        Void._unavailable(self)
    def __imod__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ipow__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ilshift__(self,*args,**kwargs):
        Void._unavailable(self)
    def __irshift__(self,*args,**kwargs):
        Void._unavailable(self)
    def __iand__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ixor__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ior__(self,*args,**kwargs):
        Void._unavailable(self)
    def __round__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ceil__(self,*args,**kwargs):
        Void._unavailable(self)
    def __floor__(self,*args,**kwargs):
        Void._unavailable(self)
    def __trunc__(self,*args,**kwargs):
        Void._unavailable(self)
    def __bool__(self,*args,**kwargs):
        Void._unavailable(self)
    def __copy__(self,*args,**kwargs):
        Void._unavailable(self)
    def __deepcopy__(self,*args,**kwargs):
        Void._unavailable(self)
    def __getstate__(self,*args,**kwargs):
        Void._unavailable(self)
    def __reduce__(self,*args,**kwargs):
        Void._unavailable(self)
    def __reduce_ex__(self,*args,**kwargs):
        Void._unavailable(self)
    def __getnewargs__(self,*args,**kwargs):
        Void._unavailable(self)
    def __setstate__(self,*args,**kwargs):
        Void._unavailable(self)
    def __hash__(self,*args,**kwargs):
        Void._unavailable(self)
    def __bytes__(self,*args,**kwargs):
        Void._unavailable(self)
    def __format__(self,*args,**kwargs):
        Void._unavailable(self)
    def __next__(self,*args,**kwargs):
        Void._unavailable(self)
#end class Void


def unavailable(module,*items):
    voids = []
    if len(items)==0:
        voids.append(Void(module))
    #end if
    for item in items:
        voids.append(Void(module,item))
    #end for
    if len(voids)==1:
        return voids[0]
    else:
        return voids
    #end if
#end def unavailable
