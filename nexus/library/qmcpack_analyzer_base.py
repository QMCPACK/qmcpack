##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  qmcpack_analyzer_base.py                                          #
#    Data object and analyzer base classes for QmcpackAnalyzer.      #
#    Maintains data global to these classes.                         #
#                                                                    #
#  Content summary:                                                  #
#    QAobject                                                        #
#      Base class for all QMCPACK analyzer components.               #
#      Exposes settings options to the user.                         #
#                                                                    #
#    Checks                                                          #
#      Class to assess overall validity based on stored results of   #
#      many checks (boolean values). Only use so far is to validate  #
#      the structure of Trace files. See qmcpack_method_analyzers.py.#
#                                                                    #
#    Plotter                                                         #
#      Wrapper class for mayavi visualization of isosurfaces and     #
#      surface slices. Previously used to visualize energy densities.#
#      See qmcpack_quantity_analyzers.py and spacegrid.py.           #
#                                                                    #
#    QAdata                                                          #
#      Represents stored data from QMCPACK's output files.           #
#      Classification marks it as a target for potential merging,    #
#      e.g. twist averaging.                                         #
#                                                                    #
#    QAHDFdata                                                       #
#      Specialization of QAdata for data from HDF files.             #
#                                                                    #
#    QAanalyzer                                                      #
#      Base class for analyzer classes. Analyzers load and analyze   #
#      data. Base class functionality includes recursive traversal   #
#      of nested analyzer object structures for loading and          #
#      analyzing data.                                               #
#                                                                    #
#====================================================================#


from numpy import minimum,resize
from generic import obj
from developer import DevBase
from hdfreader import HDFgroup
from debug import *




import numpy as np

class Plotter(DevBase):
    def __init__(self):
        self.initialized = False
        return
    #end def __init__

    def ensure_init(self):
        if not self.initialized:
            from enthought.mayavi import mlab
            from enthought.tvtk.api import tvtk
            self.mlab = mlab
            self.tvtk = tvtk

            self.show   = mlab.show
            self.plot3d = mlab.plot3d
            self.mesh   = mlab.mesh

            self.initialized = True
        #end if
    #end def ensure_init

    def isosurface(self,points,scalars,contours,dimensions,name='val'):
        self.ensure_init()
        mlab = self.mlab
        tvtk = self.tvtk
        sg=tvtk.StructuredGrid(dimensions=dimensions,points=points)
        sg.point_data.scalars = scalars
        sg.point_data.scalars.name = name
        d = mlab.pipeline.add_dataset(sg)
        iso = mlab.pipeline.iso_surface(d)
        if isinstance(contours,int):
            iso.contour.number_of_contours = contours
        elif isinstance(contours,list):
            iso.contour.auto_contours = False
            iso.contour.contours = contours
        else:
            self.error('isosurface contours must be an int or list\n  a '+str(type(contours))+' was provided instead')
        #end if
        return
    #end def isosurface

    def surface_slice(self,x,y,z,scalars,options=None):
        scale = 1.0
        opacity= 1.0
        if options!=None:
            if 'norm_height' in options:
                scale = options.norm_height/abs(scalars.max())
            if 'scale' in options:
                scale = options.scale
            if 'opacity' in options:
                opacity = options.opacity
        #end if
        self.ensure_init()
        from numerics import surface_normals
        self.mesh(x,y,z,opacity=.2)
        surfnorm = scale*surface_normals(x,y,z)
        xs=x.copy()
        ys=y.copy()
        zs=z.copy()
        xs[...] = x[...] + surfnorm[...,0]*scalars[...]
        ys[...] = y[...] + surfnorm[...,1]*scalars[...]
        zs[...] = z[...] + surfnorm[...,2]*scalars[...]
        self.mesh(xs,ys,zs,scalars=scalars,opacity=opacity)
        return
    #end def surface_slice
#end class Plotter



class QAobj_base(DevBase):
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
        #optimize = 'variance'
        optimize = 'lastcost'
        #optimize = 'energy_within_variance_tol'  # also ewvt
        )
    QAobj_base.class_set(**_default_settings)

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
        QAobj_base.class_set(**kwargs)
    #end settings
#end class QAobject



class Checks(DevBase):
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




class QAinformation(obj):
    None
#end class QAinformation



class QAdata(QAobject):
    def zero(self):
        for value in self:
            value[:] = 0
        #end for
        #self.sum()
    #end def zero

    def minsize(self,other):
        for name,value in self.iteritems():
            if name in other:
                self[name] = resize(value,minimum(value.shape,other[name].shape))
            else:
                self.error(name+' not found in minsize partner')
            #end if
        #end for
        #self.sum()
    #end def minsize

    def accumulate(self,other):
        for name,value in self.iteritems():
            if name in other:
                value += other[name][0:len(value)]
            else:
                self.error(name+' not found in accumulate partner')
            #end if
        #end for
        #self.sum()
    #end def accumulate

    def normalize(self,normalization):
        for value in self:
            value/=normalization
        #end for
        #self.sum()
    #end def normalize


    def sum(self):
        s = 0
        for value in self:
            s+=value.sum()
        #end for
        print '                sum = {0}'.format(s)
    #end def sum
#end class QAdata



class QAHDFdata(QAdata):
    def zero(self):
        for name,value in self.iteritems():
            if isinstance(value,HDFgroup):
                value.zero('value','value_squared')
            #end if
        #end for
    #end def zero

    def minsize(self,other):
        for name,value in self.iteritems():
            if isinstance(value,HDFgroup):
                if name in other and isinstance(other[name],HDFgroup):
                    value.minsize(other[name],'value','value_squared')
                else:
                    self.error(name+' not found in minsize partner')
                #end if
            #end if
        #end for
    #end def minsize

    def accumulate(self,other):
        for name,value in self.iteritems():
            if isinstance(value,HDFgroup):
                if name in other and isinstance(other[name],HDFgroup):
                    value.accumulate(other[name],'value','value_squared')
                else:
                    self.error(name+' not found in accumulate partner')
                #end if
            #end if
        #end for
    #end def accumulate

    def normalize(self,normalization):
        for value in self:
            if isinstance(value,HDFgroup):
                value.normalize(normalization,'value','value_squared')
            #end if
        #end for
    #end def normalize
#end class QAHDFdata




class QAanalyzer(QAobject):

    verbose_vlog = False

    capabilities = None
    request      = None
    run_info     = None
    method_info  = None

    opt_methods = set(['opt','linear','cslinear'])
    vmc_methods = set(['vmc'])
    dmc_methods = set(['dmc'])


    def __init__(self,nindent=0):
        self.info = QAinformation(
            initialized = False,
            data_loaded = False,
            analyzed    = False,
            failed      = False,
            nindent     = nindent
            )
        self.vlog('building '+self.__class__.__name__)
    #end def __init__

    def subindent(self):
        return self.info.nindent+1
    #end def indent

    def vlog(self,msg,n=0):
        if QAanalyzer.verbose_vlog:
            self.log(msg,n=self.info.nindent+n)
        #end if
    #end def vlog

    def reset_indicators(self,initialized=None,data_loaded=None,analyzed=None):
        if initialized!=None:
            self.info.initialized = initialized
        #end if
        if data_loaded!=None:
            self.info.data_loaded = data_loaded
        #end if
        if analyzed!=None:
            self.info.analyzed = analyzed
        #end if
    #end def reset_indicators

    def init_sub_analyzers(self):
        self.not_implemented()
    #end def init_sub_analyzers

    def load_data_local(self):
        None
    #end def load_data_local

    def remove_data_local(self):
        if 'data' in self:
            del self.data
        #end if
    #end def remove_data_local

    def analyze_local(self):
        None
    #end def analyze_local

    def set_global_info(self):
        None
    #end def set_global_info

    def unset_global_info(self):
        None
    #end def unset_global_info

    #def traverse(self,function,block_name=None,callpost=True,**kwargs):
    #    if not callpost:
    #        cls.__dict__[func_name](self,**kwargs)
    #    #end if
    #    if block_name is None or not self.info[block_name]: 
    #        for name,value in self.iteritems():
    #            if isinstance(value,QAanalyzer):
    #                value.traverse(value,func_name,block_name,callpost,**kwargs)
    #            elif isinstance(value,QAanalyzerCollection):
    #                for n,v in value.iteritems():
    #                    if isinstance(v,QAanalyzer):
    #                        v.traverse(v,func_name,block_name,callpost,**kwargs)
    #                    #end if
    #                #end for
    #            #end if
    #        #end for
    #    #end if
    #    if block_name!=None:
    #        self.info[block_name] = True
    #    #end if
    #    if callpost:
    #        cls.__dict__[func_name](self,**kwargs)
    #    #end if
    ##end def traverse

    def propagate_indicators(self,**kwargs):
        self.reset_indicators(**kwargs)
        for name,value in self.iteritems():
            if isinstance(value,QAanalyzer):
                value.propagate_indicators(**kwargs)
            elif isinstance(value,QAanalyzerCollection):
                for n,v in value.iteritems():
                    if isinstance(v,QAanalyzer):
                        v.propagate_indicators(**kwargs)
                    #end if
                #end for
            #end if
        #end for
    #end def propagate_indicators

    def load_data(self):
        if not self.info.data_loaded:
            self.vlog('loading '+self.__class__.__name__+' data',n=1)
            self.load_data_local()
            self.info.data_loaded = True
        #end if
        for name,value in self.iteritems():
            if isinstance(value,QAanalyzer):
                value.load_data()
            elif isinstance(value,QAanalyzerCollection):
                for n,v in value.iteritems():
                    if isinstance(v,QAanalyzer):
                        v.load_data()
                    #end if
                #end for
            #end if
        #end for
    #end def load_data

    def analyze(self,force=False):
        self.set_global_info()
        if not self.info.data_loaded:
            self.load_data_local()
            self.info.data_loaded = True
        #end if
        for name,value in self.iteritems():
            if isinstance(value,QAanalyzer):
                value.analyze(force)
            elif isinstance(value,QAanalyzerCollection):
                for n,v in value.iteritems():
                    if isinstance(v,QAanalyzer):
                        v.analyze(force)
                    #end if
                #end for
            #end if
        #end for
        if not self.info.analyzed or force:
            self.vlog('analyzing {0} data'.format(self.__class__.__name__),n=1)
            self.analyze_local()
            self.info.analyzed = True
        #end if
        self.unset_global_info()
    #end def analyze


    def remove_data(self):
        self.vlog('removing '+self.__class__.__name__+' data',n=1)
        names = list(self.keys())
        for name in names:
            if isinstance(self[name],QAdata):
                del self[name]
            #end if
        #end for                
        for name,value in self.iteritems():
            if isinstance(value,QAanalyzer):
                value.remove_data()
            elif isinstance(value,QAanalyzerCollection):
                for n,v in value.iteritems():
                    if isinstance(v,QAanalyzer):
                        v.remove_data()
                    #end if
                #end for
            #end if
        #end for
    #end def remove_data


    def zero_data(self):
        self.vlog('zeroing '+self.__class__.__name__+' data',n=1)
        for value in self:
            if isinstance(value,QAdata):
                value.zero()
            #end if
        #end if
        for name,value in self.iteritems():
            if isinstance(value,QAanalyzer):
                value.zero_data()
            elif isinstance(value,QAanalyzerCollection):
                for n,v in value.iteritems():
                    if isinstance(v,QAanalyzer):
                        v.zero_data()
                    #end if
                #end for
            #end if
        #end for
    #end def zero_data


    def minsize_data(self,other):
        self.vlog('minsizing '+self.__class__.__name__+' data',n=1)
        for name,value in self.iteritems():
            if isinstance(value,QAdata):
                if name in other and isinstance(other[name],value.__class__):
                    value.minsize(other[name])
                else:
                    self.error('data '+name+' not found in minsize_data partner')
                #end if
            #end if
        #end if
        for name,value in self.iteritems():
            if isinstance(value,QAanalyzer):
                if name in other and isinstance(other[name],value.__class__):
                    ovalue = other[name]
                else:
                    self.error('analyzer '+name+' not found in minsize_data partner')
                #end if
                value.minsize_data(ovalue)
            elif isinstance(value,QAanalyzerCollection):
                if name in other and isinstance(other[name],QAanalyzerCollection):
                    ovalue = other[name]
                else:
                    self.error('collection '+name+' not found in minsize_data partner')
                #end if
                for n,v in value.iteritems():
                    if isinstance(v,QAanalyzer):
                        if n in ovalue and isinstance(ovalue[n],v.__class__):
                            ov = ovalue[n]
                        else:
                            self.error('analyzer '+n+' not found in minsize_data partner collection '+name)
                        #end if
                        v.minsize_data(ov)
                    #end if
                #end for
            #end if
        #end for
    #end def minsize_data


    def accumulate_data(self,other):
        self.vlog('accumulating '+self.__class__.__name__+' data',n=1)
        for name,value in self.iteritems():
            if isinstance(value,QAdata):
                if name in other and isinstance(other[name],value.__class__):
                    value.accumulate(other[name])
                else:
                    self.error('data '+name+' not found in accumulate_data partner')
                #end if
            #end if
        #end if
        for name,value in self.iteritems():
            if isinstance(value,QAanalyzer):
                if name in other and isinstance(other[name],value.__class__):
                    ovalue = other[name]
                else:
                    self.error('analyzer '+name+' not found in accumulate_data partner')
                #end if
                value.accumulate_data(ovalue)
            elif isinstance(value,QAanalyzerCollection):
                if name in other and isinstance(other[name],QAanalyzerCollection):
                    ovalue = other[name]
                else:
                    self.error('collection '+name+' not found in accumulate_data partner')
                #end if
                for n,v in value.iteritems():
                    if isinstance(v,QAanalyzer):
                        if n in ovalue and isinstance(ovalue[n],v.__class__):
                            ov = ovalue[n]
                        else:
                            self.error('analyzer '+n+' not found in accumulate_data partner collection '+name)
                        #end if
                        v.accumulate_data(ov)
                    #end if
                #end for
            #end if
        #end for
    #end def accumulate_data


    def normalize_data(self,normalization):
        self.vlog('normalizing '+self.__class__.__name__+' data',n=1)
        for value in self:
            if isinstance(value,QAdata):
                value.normalize(normalization)
            #end if
        #end if
        for name,value in self.iteritems():
            if isinstance(value,QAanalyzer):
                value.normalize_data(normalization)
            elif isinstance(value,QAanalyzerCollection):
                for n,v in value.iteritems():
                    if isinstance(v,QAanalyzer):
                        v.normalize_data(normalization)
                    #end if
                #end for
            #end if
        #end for
    #end def normalize_data

#end class QAanalyzer



class QAanalyzerCollection(QAobject):
    None
#end class QAanalyzerCollection
