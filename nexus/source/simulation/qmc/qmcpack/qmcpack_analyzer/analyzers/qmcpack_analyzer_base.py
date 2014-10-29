
from numpy import minimum,resize
from generic import obj
from hdfreader import HDFgroup
from qaobject import QAobject


class QAinformation(obj):
    None
#end class QAinformation



class QAdata(QAobject):
    def zero(self):
        for value in self:
            value[:] = 0
        #end for
    #end def zero

    def minsize(self,other):
        for name,value in self.iteritems():
            if name in other:
                self[name] = resize(value,minimum(value.shape,other[name].shape))
            else:
                self.error(name+' not found in minsize partner')
            #end if
        #end for
    #end def minsize

    def accumulate(self,other):
        for name,value in self.iteritems():
            if name in other:
                value += other[name][0:len(value)]
            else:
                self.error(name+' not found in accumulate partner')
            #end if
        #end for
    #end def accumulate

    def normalize(self,normalization):
        for value in self:
            value/=normalization
        #end for
    #end def normalize
#end class QAdata



class QAHDFdata(QAdata):
    def zero(self):
        for value in self:
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

    capabilities = None
    request      = None
    run_info     = None
    method_info  = None

    opt_methods = set(['opt','linear','cslinear'])
    vmc_methods = set(['vmc'])
    dmc_methods = set(['dmc'])

    def __init__(self):
        self.info = QAinformation(
            initialized = False,
            data_loaded = False,
            analyzed    = False
            )
    #end def __init__

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

    def traverse(self,function,block_name=None,callpost=True,**kwargs):
        if not callpost:
            cls.__dict__[func_name](self,**kwargs)
        #end if
        if block_name is None or not self.info[block_name]: 
            for name,value in self.iteritems():
                if isinstance(value,QAanalyzer):
                    value.traverse(value,func_name,block_name,callpost,**kwargs)
                elif isinstance(value,QAanalyzerCollection):
                    for n,v in value.iteritems():
                        if isinstance(v,QAanalyzer):
                            v.traverse(v,func_name,block_name,callpost,**kwargs)
                        #end if
                    #end for
                #end if
            #end for
        #end if
        if block_name!=None:
            self.info[block_name] = True
        #end if
        if callpost:
            cls.__dict__[func_name](self,**kwargs)
        #end if
    #end def traverse

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
            self.analyze_local()
            self.info.analyzed = True
        #end if
        self.unset_global_info()
    #end def analyze


    def remove_data(self):
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
