##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  cascade.py                                                        #
#    Enables workflow/cascade manipulation without resorting to      #
#    user-side programming.  Users define a template workflow        #
#    and then perform common database/tree operations to create      #
#    complex and varied workflows.  For example, a node in a simple  #
#    workflow can be selected and any input variable can be          #
#    scanned over a range.  This will cause multiple independent     #
#    subcascades to be generated.  Supports graph visualization of   #
#    cascades/workflows.                                             #
#                                                                    #
#    Illustrative examples can be run by executing this file:        #
#      python cascade.py                                             #
#                                                                    #
#    Implementation is incomplete.                                   #
#                                                                    #                                        
#====================================================================#


import os
import string
from copy import deepcopy
from types import ClassType,TypeType
import tempfile
from numpy import ndarray

from generic import obj
from developer import DevBase,unavailable
from simulation import Simulation
from project_manager import ProjectManager
from debug import ci,ls,gs


try:
    from pydot import Dot,Node,Edge
    dot_unavailable = False
except ImportError:
    Dot,Edge,Node = unavailable('Dot','Edge','Node')
    dot_unavailable = True
#end try
try:
    import Image
    image_unavailable = False
except ImportError:
    Image = unavailable('Image')
    image_unavailable = True
#end try
try:
    import psutil
    psutil_unavailable = False
except ImportError:
    psutil = unavailable('psutil')
    psutil_unavailable = True
#end try


# these are general and should be used elsewhere
#   for present name resolution and access patterns, obj, list, and tuple will be used
#def dict_like(v):
#    return hasattr(v,'iteritems') and not isinstance(v,(ClassType,TypeType))
##end def dict_like
#
#
#def list_like(v):
#    return hasattr(v,'__iter__') and hasattr(v,'__len__') and not isinstance(v,(ClassType,TypeType))
##end def list_like

def obj_like(v):
    return isinstance(v,obj)
#end def obj_like


def list_like(v):
    return isinstance(v,(list,tuple))
#end def list_like


def func_like(v):
    return hasattr(v,'__call__') or isinstance(v,ClassType)
#end def func_like



class vdict(DevBase):
    def __init__(self,vtype):
        if not func_like(vtype):
            self.error('vdict can only be initialized with a class or function\n  you provided a {0}'.format(type(vtype)))
        #end if
        self._type = vtype
    #end def __init__


    def __call__(self,**data):
        self.set(**data)
        return self
    #end def __call__


    def _transform(self,data):
        vtype = self._type
        del self._type
        try:
            image = vtype(**self)
        except Exception,e:
            self.error('could not transform vdict into its underlying type\n  you provided type {0}\n  the following error was encountered while attempting to transform:\n  {1}\n  the function/type you provided may not accept one of the keyword arguments or there may be a bug in the underlying function code\n  please check and try again'.format(self.type.__class__.__name__,e))
        #end try
        self._type = vtype
        return image
    #end def _transform
#end class vdict




class link(DevBase):
    def __init__(self,variable_path):
        if not isinstance(variable_path,str):
            self.error('variable path must be a string\n  you provided {0}\n  which is of type {1}'.format(variable_path,variable_path.__class__.__name__))
        #end if
        self.variable_path = variable_path
    #end def __init__
#end class link




class PathGenerator(DevBase):
    def __init__(self,variables,name=None,path='',paths=None,vpaths=None,links=None,parent=None):
        if paths is None:
            paths = []
            self.paths = paths
        #end if
        if vpaths is None:
            vpaths = []
            self.vpaths = paths
        #end if
        if links is None:
            links = []
            self.links = links
        #end if
        if obj_like(variables):
            for n,value in variables.iteritems():
                if len(path)==0:
                    full_name = n
                else:
                    full_name = path+'.'+n
                #end if
                full_path = full_name
                PathGenerator(value,full_name,full_path,paths,vpaths,links,variables)
            #end for
        elif list_like(variables):
            i=0
            for value in variables:
                full_name = path+'[{0}]'.format(i)
                full_path = full_name
                PathGenerator(value,full_name,full_path,paths,vpaths,links,variables)
                i+=1
            #end for
        #end if
        if name!=None:
            paths.append(name)
            if isinstance(parent,vdict) and not isinstance(variables,vdict):
                vpaths.append(name)
            #end if
            if isinstance(variables,link):
                links.append((name,variables.variable_path))
            #end if
        #end if
    #end def __init__

    
    def get_paths(self):
        return set(self.paths),set(self.vpaths),self.links
    #end def get_paths
#end class PathGenerator




class CascadeElement(DevBase):

    graph_node_settings = obj(
        default = obj(
            fillcolor = 'cyan',
            fontcolor = 'black'
            ),
        blocked = obj(
            fillcolor = 'black',
            fontcolor = 'white'
            ),
        selected = obj(
            fillcolor = 'red',
            fontcolor = 'black'
            ),
        clustered = obj(
            fillcolor = 'green',
            fontcolor = 'black'
            ),
        branched = obj(
            fillcolor = 'orange',
            fontcolor = 'black'
            )
        )

    state_variables = set('name eid cid'.split())

    nelements = 0
    @classmethod
    def checkout_id(cls):
        eid = cls.nelements
        cls.nelements+=1
        return eid
    #end def checkout_id


    def __init__(self,name,sim_generator):
        self.state = obj(
            name = name,
            eid  = self.checkout_id(),
            cid  = -1
            )
        self.variables      = None
        self.variable_paths = None
        self.allowed_paths  = None
        self.sim_generator  = sim_generator
        self.generate       = func_like(sim_generator)
        self.generated      = False
        self.dependencies   = obj()
        self.dependents     = obj()
        self.sim            = None
        self.clone          = None
        self.variations     = obj()
        self.graph_node     = Node(self.state.eid,style='filled',shape='Mrecord')
        self.graph_edges    = obj()
        self.status         = None
        self.reset_status()
        if not self.generate and not isinstance(sim_generator,Simulation):
            self.error('cascade element {0} must be a Simulation object or a function that generates one\n  you provided a {1}'.format(name,sim_generator.__class__.__name__))
        #end if
    #end def __init__


    def reset_status(self):
        self.blocked    = False
        self.selected   = False
        self.clustered  = False
        self.cluster_id = -1
        self.branched   = False
        self.update_status()
    #end def reset_status


    def update_status(self,graph=False,s=None,v=None,e=None):
        if self.branched:
            self.status = 'branched'
        elif self.clustered:
            self.status = 'clustered'
        elif self.selected:
            self.status = 'selected'
        elif self.blocked:
            self.status = 'blocked'
        else:
            self.status = 'default'
        #end if
        if graph:
            self.update_graph_elements(s=s,v=v,e=e)
        #end if
    #end def update_status

        
    def update_graph_elements(self,s=None,v=None,e=None):
        if isinstance(s,str):
            s = [s]
        #end if
        if isinstance(v,str):
            v = [v]
        #end if
        if isinstance(e,str):
            e = [e]
        #end if
        nset = self.graph_node_settings[self.status]
        gn = self.graph_node
        gn.set_fillcolor(nset.fillcolor)
        gn.set_fontcolor(nset.fontcolor)
        name = self.state.name
        ladd = ''
        if s!=None:
            ladd += ' {0}={1}'.format(s[0],self.state[s[0]])
            for var in s[1:]:
                ladd += ', {0}={1}'.format(var,self.state[var])
            #end for
        #end if
        if e!=None:
            ladd += ' {0}={1}'.format(e[0],self[e[0]])
            for var in e[1:]:
                ladd += ', {0}={1}'.format(var,self[var])
            #end for
        #end if
        if self.status=='clustered':
            ladd = ' {0}'.format(self.cluster_id)+ladd
        #end if
        label = name + ladd
        if len(self.variations)>0:
            vnames = list(self.variations.keys())
            vnames.sort()
            ls = '<<table border="0" cellborder="0"><tr><td align="center">{0}</td></tr>'.format(label)
            for vn in vnames:
                vs = self.variations[vn]
                ls += '<tr><td align="center" port="{0}"><font point-size="10">{1}</font></td></tr>'.format(vn,vs)
            #end for
            ls+='</table>>'
            label = ls
        #end if
        gn.set_label(label)
    #end def update_graph_elements


    def acquire_variables(self,variables):#,variable_paths,allowed_paths):
        self.variables      = deepcopy(variables)
        #self.variable_paths = deepcopy(variable_paths)
        #self.allowed_paths  = set(allowed_paths)
    #end def acquire_variables


    def depends(self,other,*quantities):
        oid = other.state.eid
        quantities = set(quantities)
        edges = obj()
        for quantity in quantities:
            edges[quantity] = Edge(other.graph_node, self.graph_node,
                        label=quantity, fontsize='10.0')
        #end for
        if oid in self.dependencies:
            dep = self.dependencies[oid]
            for quantity in quantities:
                dep.quantities.add(quantity)
                dep.edges[quantity] = edges[quantity]
            #end for
        else:
            dep = obj(element=other,quantities=quantities,edges=edges)
            self.dependencies[oid] = dep
        #end if
        other.dependents[self.state.eid] = self
    #end def depends


    def add_variation(self,vname,variation):
        if not isinstance(variation,str):
            self.error('expected variation string, received {0}'.format(variation.__class__.__name__))
        #end if
        self.variations[vname] = variation
    #end def add_variation


    def select(self):
        self.selected = True
    #end def select

            
    def unselect(self):
        self.selected = False
    #end def unselect


    def cluster(self):
        self.clustered = True
    #end def cluster

        
    def uncluster(self):
        self.clustered  = False
        self.cluster_id = -1
    #end def uncluster


    def branch(self):
        self.branched = True
    #end def branch


    def unbranch(self):
        self.branched = False
    #end def unbranch


    def sever_links(self):
        self.sever_dependency_links()
        self.sever_dependent_links()
    #end def sever_links

        
    def sever_dependency_links(self):
        eid = self.state.eid
        for dep in self.dependencies:
            other = dep.element
            del other.dependents[eid]
        #end for
        self.dependencies.clear()
    #end def sever_dependency_links


    def sever_dependent_links(self):
        eid = self.state.eid
        for other in self.dependents:
            del other.dependencies[eid]
        #end for
        self.dependents.clear()
    #end def sever_dependent_links


    def sever_dependency_link(self,other):
        oid = other.state.eid
        eid = self.state.eid
        if not oid in self.dependencies:
            self.error('cannot sever dependency link of {0}\n  cascade element {0} {1} does not depend on cascade element {2} {3}'.format(self.state.name,eid,other.state.name,oid))
        #end if
        del other.dependents[eid]
        del self.dependencies[oid]
    #end def sever_dependency_link
        
        
    def sever_dependent_link(self,other):
        oid = other.state.eid
        eid = self.state.eid
        if not oid in self.dependents:
            self.error('cannot sever dependent link of {3}\n  cascade element {0} {1} does not depend on cascade element {2} {3}'.format(other.state.name,oid,self.state.name,eid))
        #end if
        del other.dependencies[eid]
        del self.dependents[oid]
    #end def sever_dependent_link


    def sever_up_links(self):
        self.sever_dependency_links()
    #end def sever_up_links


    def sever_down_links(self):
        self.sever_dependent_links()
    #end def sever_down_links


    def attempt_selection(self,sel):
        s = self.state
        v = self.variables
        try:
            exec('selected = {0}'.format(sel))
        except Exception,e:
            self.error('  selection check failed\n  selection string provided: {0}\n  this string may not be formatted correctly\n  please revise the selection string and try again\n  the following exception was caught and may be helpful to you:\n{1}'.format(sel,e))
        #end try
        if not isinstance(selected,bool):
            self.error('  selection check failed\n  selection string provided: {0}\n  evaluation of the selection string did nor result in a boolean value\n  this string may not be formatted correctly\n  please revise the selection string and try again'.format(sel))
        #end if
        return selected
    #end def attempt_selection


    def block(self):
        if not self.blocked:
            self.blocked = True
            self.update_status()
            if self.generated:
                self.sim.block = True
            #end if
            for dep in self.dependents:
                dep.block()
            #end for
        #end if
    #end def block


    def unblock(self):
        if self.blocked:
            self.blocked = False
            if self.generated:
                self.sim.block = False
            #end if
            for dep in self.dependents:
                dep.unblock()
            #end for
        #end if
    #end def unblock
      

    def spread_cid(self,cid,cascade,visited=None):
        if visited is None:
            visited = set()
        #end if
        eid = self.state.eid
        self.state.cid = cid
        cascade[eid] = self
        visited.add(eid)
        for dep in self.dependents:
            if not dep.state.eid in visited:
                dep.spread_cid(cid,cascade,visited)
            #end if
        #end for
        for odep in self.dependencies:
            dep = odep.element
            if not dep.state.eid in visited:
                dep.spread_cid(cid,cascade,visited)
            #end if
        #end for
    #end def spread_cid
      

    def spread_cluster_id(self,cluster_id,cluster,visited=None):
        if visited is None:
            visited = set()
        #end if
        eid = self.state.eid
        if self.clustered and self.cluster_id==-1:
            self.cluster_id = cluster_id
            cluster[eid] = self
        #end if
        visited.add(eid)
        for dep in self.dependents:
            if not dep.state.eid in visited and dep.clustered and dep.cluster_id==-1:
                dep.spread_cluster_id(cluster_id,cluster,visited)
            #end if
        #end for
        for odep in self.dependencies:
            dep = odep.element
            if not dep.state.eid in visited and dep.clustered and dep.cluster_id==-1:
                dep.spread_cluster_id(cluster_id,cluster,visited)
            #end if
        #end for
    #end def spread_cluster_id


    def collect_dependencies(self,coll,exclude_self=False):
        eid = self.state.eid
        if not eid in coll or exclude_self:
            if not exclude_self:
                coll[eid] = self
            #end if
            for dep in self.dependencies:
                elem = dep.element
                elem.collect_dependencies(coll)
            #end for
        #end if
    #end def collect_dependents
        

    def collect_dependents(self,coll,exclude_self=False):
        eid = self.state.eid
        if not eid in coll or exclude_self:
            if not exclude_self:
                coll[eid] = self
            #end if
            for elem in self.dependents:
                elem.collect_dependents(coll)
            #end for
        #end if
    #end def collect_dependents


    def collect_up(self,coll,exclude_self=False):
        self.collect_dependencies(coll,exclude_self)
    #end def collect_up


    def collect_down(self,coll,exclude_self=False):
        self.collect_dependents(coll,exclude_self)
    #end def collect_down


    def mark_up(self,var,value,exclude_self=False,state=False):
        if not exclude_self:
            if state:
                vset = self.state
            else:
                vset = self
            #end if
            vset[var] = value
        #end if
        for dep in self.dependencies:
            dep.element.mark_up(var,value)
        #end for
    #end def mark_up


    def mark_down(self,var,value,exclude_self=False,state=False):
        if not exclude_self:
            if state:
                vset = self.state
            else:
                vset = self
            #end if
            vset[var] = value
        #end if
        for element in self.dependents:
            element.mark_down(var,value)
        #end for
    #end def mark_down


    def make_clone(self):
        clone = CascadeElement(self.state.name,self.sim_generator)
        clone.acquire_variables(self.variables)
        clone.variations = deepcopy(self.variations)
        clone.state.cid = self.state.cid
        clone.blocked  = self.blocked
        clone.selected = self.selected
        self.clone = clone
    #end def make_clone


    def link_clone(self):
        clone = self.clone
        if clone is None:
            self.error('attempted to link non-existent clone')
        #end if
        for dep in self.dependencies:
            other = dep.element
            quantities = dep.quantities
            if other.clone!=None:
                clone.depends(other.clone,*quantities)
            else:
                clone.depends(other,*quantities)
            #end if
        #end for
        for other in self.dependents:
            if other.clone is None:
                self.error('attempted to link to non-existent dependent clone')
            #end if
            dep = other.dependencies[self.state.eid]
            quantities = dep.quantities
            other.clone.depends(clone,*quantities)
        #end for
    #end def link_clone


    def yield_clone(self,lst):
        clone = self.clone
        if clone is None:
            self.error('attempted to yield non-existent clone')
        #end if
        lst.append(clone)
        self.clone = None
    #end def yield_clone


    def make_sim(self):
        if not self.generated:
            if self.generate:
                inputs = self.variables[self.state.name]
                try:
                    self.sim = self.sim_generator(**inputs)
                except Exception,e:
                    self.error('generating simulations failed\n  cascade element {0} named {1} encountered an error when attempting to generate its corresponding simulation\n  some of the variables provided may not match the function you provided (in the elements argument of Cascade) to generate the simulation\n  please check the compatibility of the input variables with the simulation generating function and try again\n  the variables listed below and error message encountered may be of some help\n  variables provided:\n{2}\n  error message encountered:\n{3}'.format(self.state.eid,self.state.name,self.variables,e))
                #end if
                if not isinstance(self.sim,Simulation):
                    self.error('generating simulations failed\n  simulation generating function provided in the elements argument of Cascade for simulations of type {0} did not produce a Simulation object\n  object type produced: {1}\n  please check the function provided for {0} and try again'.format(self.state.name,self.sim.__class__.__name__))
            else:
                self.sim = self.sim_generator
            #end if
            self.generated = True
        #end if
    #end def make_sim


    def link_sim(self):
        deps = []
        for dep in self.dependencies:
            dlist = [dep.element]
            dlist.extend(dep.quantities)
            deps.append(dlist)
        #end for
        if len(deps)>0:
            self.sim.depends(deps)
        #end if
    #end def link_sim


    def yield_sim(self,lst):
        lst.append(self.sim)
    #end def yield_sim
#end class CascadeElement




class CascadeElementCollection(DevBase):
    def __init__(self,*elements):
        self.add_elements(*elements)
    #end def __init__

    def add_elements(self,*elements):
        for c in elements:
            self[c.state.eid] = c
        #end for
    #end def add_elements

    def remove_elements(self,*elements):
        for c in elements:
            del self[c.state.eid]
        #end for
    #end def remove_elements

    def collect_selected(self,*elements):
        for c in elements:
            if c.selected:
                self[c.state.eid] = c
            #end if
        #end for
    #end def collect_selected

    def disjoint(self,fragment=False):
        ids = set()
        for c in self:
            for dep in c.dependencies:
                ids.add(dep.element.state.eid)
            #end for
            for element in c.dependents:
                ids.add(element.state.eid)
            #end for
        #end for
        if fragment:
            disjoint = len(set(self.keys())-ids)>0
        else:
            disjoint = not len(ids-set(self.keys()))>0
        #end if
        return disjoint
    #end def disjoint

    def operate(self,op,*args,**kwargs):
        for c in self:
            op(c,*args,**kwargs)
        #end for
    #end def operate

    def select(self):           self.operate(CascadeElement.select)
    def unselect(self):         self.operate(CascadeElement.unselect)
    def cluster(self):          self.operate(CascadeElement.cluster)
    def uncluster(self):        self.operate(CascadeElement.uncluster)
    def branch(self):           self.operate(CascadeElement.branch)
    def unbranch(self):         self.operate(CascadeElement.unbranch)
    def make_clone(self):       self.operate(CascadeElement.make_clone)
    def link_clone(self):       self.operate(CascadeElement.link_clone)
    def yield_clone(self,coll): self.operate(CascadeElement.yield_clone,coll)
    def block(self):            self.operate(CascadeElement.block)
    def unblock(self):          self.operate(CascadeElement.unblock)
    def sever_links(self):      self.operate(CascadeElement.sever_links)
    def make_sim(self):         self.operate(CascadeElement.make_sim)
    def link_sim(self):         self.operate(CascadeElement.link_sim)
    def yield_sim(self,lst):    self.operate(CascadeElement.yield_sim,lst)
#end class CascadeElementCollection




class SubCascades(DevBase):
    def __init__(self,*elements):
        all = CascadeElementCollection(*elements)

        cascades = obj()
        ids = list(all.keys())
        cid = 0
        while len(ids)>0:
            cascade = CascadeElementCollection()
            nucleus = all[ids[0]]
            nucleus.spread_cid(cid,cascade)
            for c in cascade:
                ids.remove(c.state.eid)
            #end for
            cascades[cid] = cascade
            cid+=1
        #end while

        self.cascades = cascades
        self.all = all
    #end def __init__


    def add_cluster(self,cluster):
        if not isinstance(cluster,CascadeElementCollection):
            self.error('could not add sub-cascade\n  sub-cascade must be of type CascadeElementCollection\n  you provided a {0}'.format(cascade.__class__.__name__))
        #end if
        if len(cluster)==0:
            self.error('attempted to cluster with no elements to SubCascades\n  this is probably a mistake')
        #end if
        eid = cluster.keys()[0]
        ce = cluster[eid]
        if cluster.disjoint():
            cid = max(self.cascades.keys())+1
            ce.spread_cid(cid,cluster)
            self.cascades[cid] = cluster
        else:
            cid = ce.state.cid
            if not cid in self.cascades:
                self.error('attempted to add unknown cluster to SubCascades\n  the cluster is not disjoint, but the cascade id (cid) of the cluster does not match any existing id')
            #end if
            self.cascades[cid].add_elements(*cluster)
        #end if
        self.all.add_elements(*cluster)
    #end def add_cluster


    def delete(self,*elements):
        #elem_to_delete = obj()
        #for element in elements:
        #    element.collect_dependents(elem_to_delete)
        ##end for
        #for eid,element in elem_to_delete.iteritems():
        for element in elements:
            eid = element.state.eid
            cid = element.state.cid
            element.sever_links()
            element.clear()
            del self.cascades[cid][eid]
            del self.all[eid]
        #end for
    #end def delete
#end class SubCascades




class Cascade(DevBase):


    project_manager = ProjectManager()


    def __init__(self,description=None,elements=None,dependencies=None,variables=None,postpone=False):
        if description!=None and not isinstance(description,str):
            self.error('Cascade description must be a string\n  you provided a {0} with value {1}\n  please provide a description string and try again'.format(description.__class__.__name__,description))
        #end if
        self.set(
            description       = description,
            initialized      = False,
            elements         = obj(),
            dependencies     = obj(),
            variables        = obj(),
            variable_paths   = obj(),
            selector_by_name = obj(),
            cascade          = None,  # set to SubCascades in init_cascade
            simulations      = None,  # set in generate_simulations
            selection        = CascadeElementCollection(),
            generated_simulations = False
            )
        if elements!=None:
            self.add_elements(elements)
        #end if
        if dependencies!=None:
            self.add_dependencies(dependencies)
        #end if
        if variables!=None:
            self.add_variables(variables)
        #end if
        if not postpone:
            self.initialize()
        #end if
    #end def __init__


    def add_elements(self,elements,dependencies=None,variables=None,overwrite=False):
        if self.initialized:
            self.error('attempted to add elements after initialization\n  please note that initialization automatically occurs upon creation of a Cascade unless input argument postpone=True')
        #end if
        if not obj_like(elements):
            self.error('elements must be of type obj or similar')
        #end if
        if not overwrite:
            for name,elem in elements.iteritems():
                if name in self.elements:
                    self.error('attempted to overwrite element {0}\n  if this is the behavior you want, set overwrite to True'.format(name))
                else:
                    self.elements[name] = CascadeElement(name,elem)
                #end if
            #end for
        else:
            for name,elem in elements.iteritems():
                self.elements[name] = CascadeElement(name,elem)
            #end for
        #end if
        for name in elements.keys():
            self.selector_by_name[name] = (False,None)
        #end for
        if dependencies!=None:
            self.add_dependencies(dependencies)
        #end if
        if variables!=None:
            self.add_variables(variables,overwrite)
        #end if
    #end def add_elements


    def add_dependencies(self,dependencies):
        if self.initialized:
            self.error('attempted to add dependencies after initialization\n  please note that initialization automatically occurs upon creation of a Cascade unless input argument postpone=True')
        #end if
        if not obj_like(dependencies):
            self.error('dependencies must be of type obj or similar')
        #end if
        for name,deps_in in dependencies.iteritems():
            if not name in self.elements:
                self.error('cannot add dependencies to element {0}, it does not exist'.format(name))
            elif not obj_like(deps_in):
                self.error('dependencies for element {0} must be of type obj or similar\n  you provided a {2}'.format(name,deps_in.__class__.__name__))
            else:
                if not name in self.dependencies:
                    deps = obj()
                    self.dependencies[name] = deps
                else:
                    deps = self.dependencies[name]
                #end if
                for quantity,simtag in deps_in.iteritems():
                    if not quantity in deps:
                        deps[quantity] = simtag
                    else:
                        self.error('attempted to add dependency {0} to simulation {1}, but it is already present\n  {2} is currently obtained from {3}\n  attempted to request {4} from {5}'.format(quantity,name,quantity,deps[quantity],quantity,simtag))
                    #end if
                #end for
            #end if
        #end for
    #end def add_dependencies


    def add_variables(self,variables,overwrite=False):
        if self.initialized:
            self.error('attempted to add variables after initialization\n  please note that initialization automatically occurs upon creation of a Cascade unless input argument postpone=True')
        #end if
        if not obj_like(variables):
            self.error('variables must be of type obj or similar')
        #end if
        if not overwrite:
            for name,varset in variables.iteritems():
                if name in self.variables:
                    self.error('attempted to overwrite variable {0}\n  if this is the behavior you want, set overwrite to True'.format(name))
                else:
                    self.variables[name] = deepcopy(varset)
                #end if
            #end for
        else:
            for name,varset in variables.iteritems():
                self.variables[name] = deepcopy(varset)
            #end for
        #end if
    #end def add_variables


    def initialize(self):
        if self.initialized:
            self.error('attempted to re-initialize after initialization')
        #end if
        self.elem_names = set(self.elements.keys())
        self.init_variables()
        self.init_cascade()
        self.initialized=True
    #end def initialize


    def init_variables(self,depth=2):
        if self.initialized:
            self.error('attempted to call init_variables after initialization')
        #end if
        elem_names = self.elem_names
        velems = set(self.variables.keys())
        not_present = elem_names-velems
        for name in not_present:
            if isinstance(self.elements[name].sim_generator,Simulation):
                not_present.remove(name)
            #end if
        #end for
        if len(not_present)>0:
            self.error('cascade cannot be constructed\n  the following simulations  are elements of the cascade, but variables have not been provided for them:\n  {0}\n  please provide these variables and try again'.format(list(not_present)))
        #end if

        # make variable paths
        pg = PathGenerator(self.variables)
        npaths,vpaths,links = pg.get_paths()
        # make a second pass if there are links
        if len(links)>0:
            for path,lpath in links:
                if not lpath in npaths:
                    self.error('variable link cannot be resolved\n  link found at: {0}\n  requested destination: {1}\n  this destination does not exist\n  please revise the link path provided to variables and try again'.format(path,lpath))
                #end if
                v = self.variables
                try:
                    exec('val=v.{1}; v.{0} = deepcopy(val)'.format(path,lpath))
                except Exception,e:
                    self.error('variable link resolution failed\n  this is a bug, please contact the developer')
                #end try
                if not isinstance(val,(obj,list)):
                    self.error('links can only point to list or obj types\n  link found at: {0}\n  requested destination: {1}\n  this destination is of type {2}\n  please ensure that the target of the link is a list or obj and try again'.format(path,lpath,val.__class__.__name__))
                #end if
            #end for
            pg = PathGenerator(self.variables)
            npaths,vpaths,links2 = pg.get_paths()
            if len(links2)>0:
                unresolved = ''
                for path,lpath in links2:
                    unresolved+='\n    {0} -> {1}'.format(path,lpath)
                #end for
                self.error('some variable links have not been resolved\n  this can be caused by multi-level or circular links\n  please note that only simple links are allowed\n  (a simple link points at a variable that contains no links)\n  please revise your links and try again\n  the list of unresolved links provided below may be of some help:{0}'.format(unresolved))
            #end if
        #end if

        # separate paths for each element
        vpaths = obj()
        for ename in self.elem_names:
            epaths = obj()
            vpaths[ename] = epaths
            enamed = ename+'.'
            enameb = ename+'['
            for path in npaths:
                if path.startswith(enamed) or path==ename or path.startswith(enameb):
                    epaths[path] = path
                #end if
            #end for
        #end for

        # make path shortcuts
        drange = range(depth)
        for ename,epaths in vpaths.iteritems():
            names   = []
            name_count = dict()
            for path in epaths:
                tokens = path.split('.')
                l = len(tokens)
                for d in drange:
                    if len(tokens)>d:
                        name = tokens[l-1-d]
                        for n in tokens[l-d:]:
                            name+='.'+n
                        #end for
                        names.append((l,name,path))
                        idx = l,name
                        if not idx in name_count:
                            name_count[idx] = 1
                        else:
                            name_count[idx] += 1
                        #end if
                    #end if
                #end for
            #end for
            names.sort()
            for l,name,path in names:
                if name_count[l,name]==1 and not name in epaths:
                    epaths[name] = path
                #end if
            #end for
        #end for

        # resolve links
        for ename,epaths in vpaths.iteritems():
            enamed = ename+'.'
            enameb = ename+'['
            for path,name in links:
                if path.startswith(enamed) or path==ename or path.startswith(enameb):
                    if not name in epaths:
                        for n,p in all_paths.iteritems():
                            if p.startswith(name):
                                epaths[n] = p.replace(name,path)
                            #end if
                        #end for
                        epaths[name] = path
                    #end if
                #end if
            #end for
        #end for

        # make non-namespace prefixed versions of paths
        for ename,epaths in vpaths.iteritems():
            enamed = ename+'.'
            enamed_len = len(enamed)
            names = set(epaths.keys())
            for name in names:
                path = epaths[name]
                if path==ename or name==ename:
                    del epaths[name]
                #end if
            #end for
            names = set(epaths.keys())
            for name in names:
                path = epaths[name]
                new = False
                if name.startswith(enamed):
                    new_name = name[enamed_len:]
                    new = True
                else:
                    new_name = name
                #end if
                if path.startswith(enamed):
                    new_path = path[enamed_len:]
                    new = True
                #end if
                if new:
                    del epaths[name]
                    if not new_name in epaths:
                        epaths[new_name] = new_path
                    #end if
                #end if
            #end for
        #end for

        # collect valid variable names
        vvars = []
        apaths = obj()
        for ename,epaths in vpaths.iteritems():
            vset = set(epaths.keys())
            vvars.extend(vset)
            apaths[ename] = vset
        #end for
        self.variable_paths  = vpaths
        self.valid_variables = set(vvars)
        self.allowed_paths   = apaths

        # distribute variable sets to cascade elements
        for ename,element in self.elements.iteritems():
            vvals = self.variables[ename]
            #vpaths = self.variable_paths[ename]
            #apaths = self.allowed_paths[ename]
            #element.acquire_variables(vvals,vpaths,apaths)
            element.acquire_variables(vvals)
        #end for
    #end def init_variables


    def init_cascade(self):
        if self.initialized:
            self.error('attempted to call init_cascade after initialization')
        #end if
        elem_names = self.elem_names
        dependent_keys = set(self.dependencies.keys())
        dependency_keys = []
        for deps in self.dependencies:
            dependency_keys.extend(deps.values())
        #end for
        dependency_keys = set(dependency_keys)
        not_present = dependent_keys-elem_names
        if len(not_present)>0:
            self.error('cascade cannot be constructed\n  the following simulations  are not elements of the cascade, but they depend on simulations in the cascade:\n  {0}\n  please provide these cascade elements and try again'.format(list(not_present)))
        #end if
        not_present = dependency_keys-elem_names
        if len(not_present)>0:
            self.error('cascade cannot be constructed\n  the following simulations  are not elements of the cascade, but simulations in the cascade depend on them:\n  {0}\n  please provide these cascade elements and try again'.format(list(not_present)))
        #end if
        for name,deps in self.dependencies.iteritems():
            dependent = self.elements[name]
            for quantity,depname in deps.iteritems():
                dependency = self.elements[depname]
                dependent.depends(dependency,quantity)
            #end for
        #end for
        self.cascade = SubCascades(*self.elements)
    #end def init_cascade


    def generate_simulations(self):
        if self.generated_simulations:
            self.error('attempted to generate simulations, but simulations have already been generated')
        #end if
        self.simulations = []
        all = self.cascade.all
        all.make_sim()
        all.link_sim()
        all.yield_sim(self.simulations)
        self.generated_simulations = True
    #end def generate_simulations


    def run(self):
        if not self.generated_simulations:
            self.generate_simulations()
        #end if
        self.project_manager.add_simulations(self.simulations)
        self.project_manager.run_project()
    #end def run


    def graph(self,label=None,label_append=None,show=True,pause=True,savefile=None,format='png',s=None,v=None):
        if dot_unavailable:
            self.log('  graph requested but the Dot/pydot module is not present\n  you will need to install Dot/pydot if you want to use this feature')
            return
        #end if
        if label!=None and not isinstance(label,str):
            self.error('label provided for graph must be a string\n  you provided a {0} with value {1}\n  please provide a label string and try again'.format(label.__class__.__name__,label))
        #end if
        if label is None:
            label = self.description
        #end if
        if isinstance(label_append,str):
            label += label_append
        #end if

        graph = Dot(graph_type='digraph')
        if label!=None:
            graph.set_label(label)
            graph.set_labelloc('t')
        #end if
        for c in self.cascade.all:
            c.update_status(graph=True,s=s,v=v)
        #end for
        for c in self.cascade.all:
            graph.add_node(c.graph_node)
        #end for
        for c in self.cascade.all:
            for dep in c.dependencies:
                for edge in dep.edges:
                    graph.add_edge(edge)
                #end for
            #end for
        #end for
        ext = '.'+format
        if savefile is None:
            fout = tempfile.NamedTemporaryFile(suffix=ext)
            savefile = fout.name
        elif not isinstance(savefile,str):
            self.error('savefile must be a string\n  you provided a {0} with value {1}\n  please provide a valid filename string and try again'.format(type(savefile),savefile))
        elif not savefile.endswith(ext):
            savefile += ext
            fout = open(savefile,'rb')
        #end if
        graph.write(savefile,format=format)
        if show:
            if image_unavailable:
                self.log('  cannot view cascade graph: module Image is not present\n  you will have to install Image if you want to use this feature')
                return
            #end if
            if pause:
                if psutil_unavailable:
                    self.log('  cannot pause viewing of cascade graph image: module psutil is not present\n  you will have to install psutil if you want to use this feature')
                else:
                    pids = set()
                    for proc in psutil.process_iter():
                        if proc.name=='display':
                            pids.add(proc.pid)
                        #end if
                    #end for
                #end if
            #end if
            image = Image.open(savefile)
            image.show()
            if pause and not psutil_unavailable:
                raw_input('Press enter to close cascade diagram: ')
                for proc in psutil.process_iter():
                    if proc.name=='display' and not proc.pid in pids:
                        proc.kill()
                    #end if
                #end for
            #end if
        #end if
    #end def graph


    name_characters = set(string.letters+string.digits+'_')
    path_characters = set(string.letters+string.digits+'_.[]')
    def select(self,selector,direction=None,cont=False,continued=False,graph=False,**graph_args):
        continued = continued or cont
        sel = selector
        if not isinstance(sel,str):
            self.error('argument to select must be a string\n  you provided a {0}\n  please change the input to a selection string and try again'.format(type(sel)))
        #end if
        if sel=='all':
            self.unselect()
            self.selection.add_elements(*self.cascade.all)
            self.selection.select()
            if graph:
                self.graph(**graph_args)
            #end if
            return
        elif sel=='none':
            self.unselect()
            if graph:
                self.graph(**graph_args)
            #end if
            return
        #end if
        #parse the selection string
        nc = self.name_characters
        pc = self.path_characters
        ns = len(sel)
        vset = obj(s=[],v=[])
        for c,paths in vset.iteritems():
            i=0
            while i<ns:
                if sel[i]==c:
                    if i==0 and ns>1 and sel[1]=='.' or i<ns and not sel[i-1] in nc and sel[i+1]=='.':
                        i+=2
                        s = i
                        while i<ns and sel[i] in pc:
                            i+=1
                        #end while
                        paths.append(sel[s:i])
                    else:
                        i+=1
                    #end if
                else:
                    i+=1
                #end if
            #end while
        #end for
        vset.s = set(vset.s)
        vset.v = set(vset.v)
        if len(vset.s)>0:
            svars = CascadeElement.state_variables
            dne = vset.s-svars
            if len(dne)>0:
                dne = list(dne)
                dne.sort()
                svars = list(svars)
                svars.sort()
                options=''
                for svar in svars:
                    options+='\n    '+svar
                #end for
                self.error('encountered invalid selection string\n  selection string provided: {0}\n  the following requested state variables (s.*) do not exist:\n    {1}\n  valid options are:{2}'.format(sel,list(dne),options))
            #end if
        #end if
        if len(vset.v)>0:
            vvars = self.valid_variables
            dne = vset.v-vvars
            if len(dne)>0:
                dne = list(dne)
                dne.sort()
                vvars = list(vvars)
                vvars.sort()
                options=''
                for vvar in vvars:
                    options+='\n    '+vvar
                #end for
                self.error('encountered invalid selection string\n  selection string provided: {0}\n  the following requested variables (v.*) do not exist:\n    {1}\n  valid options are:{2}'.format(sel,list(dne),options))
            #end if
        #end if
        selector_by_name = self.selector_by_name
        for name,vpaths in self.variable_paths.iteritems():
            apaths = self.allowed_paths[name]
            attempt_select = vset.v<=apaths
            if attempt_select:
                seln = str(sel)
                for valias in vset.v:
                    vpath = vpaths[valias]
                    if valias!=vpath:
                        seln = seln.replace(valias,vpath)
                    #end if
                #end for
            else:
                seln = None
            #end if
            selector_by_name[name] = (attempt_select,seln)
        #end for
                    
        selection = []
        if continued:
            subset = self.selection
        else:
            subset = self.cascade.all
        #end if
        self.cascade.all.unselect()
        for c in subset:
            attempt_select,seln = selector_by_name[c.state.name]
            if attempt_select:
                print 
                selected = c.attempt_selection(seln)
                if selected:
                    selection.append(c)
                #end if
            #end if
        #end for
        if len(selection)==0:
            self.error('no cascade elements selected\n  selection string did not match any cascade elements\n  selection string provided: {0}\n  please revise your selection request and try again'.format(sel))
        #end if
        self.selection.clear()
        if direction!=None:
            self.expand_selection(selection,self.selection,direction,exclude=True,location='select')
        #end if
        self.selection.add_elements(*selection)
        self.selection.select()
        if graph:
            self.graph(label_append=' select',**graph_args)
        #end if
    #end def select


    def unselect(self):
        self.cascade.all.unselect()
        self.selection.clear()
    #end def unselect


    def expand_selection(self,selection,expanded_selection,direction,exclude=False,location='expand_selection'):
        if direction=='up':
            for c in selection:
                c.collect_up(expanded_selection,exclude_self=exclude)
            #end for
        elif direction=='down':
            for c in selection:
                c.collect_down(expanded_selection,exclude_self=exclude)
            #end for
        elif direction=='updown':
            for c in selection:
                c.collect_up(expanded_selection,exclude_self=exclude)
                c.collect_down(expanded_selection,exclude_self=True)
            #end for
        else:
            self.error('invalid direction supplied to {0}\n  you provided: {1}\n  valid options are: up, down, updown'.format(location,direction))
        #end if
    #end def expand_selection


    def check_selection(self,opname):
        if len(self.selection)==0:
            self.error('operation {0} cannot be performed\n  no cascade elements are currently selected\n  please use the select operation to choose cascade elements prior to using operation {0}'.format(opname))
        #end if
    #end def check_selection
        

    def check_generated_simulations(self,opname):
        if self.generated_simulations:
            self.error('operation {0} cannot be performed\n  operation {0} modifies cascade structure and no operations of this type are permitted after simulation objects have been generated\n  note that simulations are generated by the following operations: generate_simulations'.format(opname))
        #end if
    #end def check_generated_simulations


    def block(self):
        self.check_selection('block')
        self.selection.block()
    #end def block
        

    def unblock(self):
        self.check_selection('unblock')
        self.selection.unblock()
    #end def block


    def delete(self):
        self.check_generated_simulations('delete')
        self.check_selection('delete')
        self.cascade.delete(*self.selection)
        self.unselect()
    #end def delete


    branch_types = (list,tuple,ndarray)
    def branch(self,variables,values,graph_clusters=False,graph_branches=False,**graph_args):
        # check for input errors
        self.check_generated_simulations('branch')
        self.check_selection('branch')
        if not isinstance(values,self.branch_types):
            self.error('branch failed\n  values provided to branch must be in a list, tuple, or array\n  you provided a {0} for values\n  please assemble branch values in a list or array and try again'.format(values.__class__.__name__))
        #end if
        if len(values)==0:
            self.error('branch failed\n  list of values provided to branch must contain one or more elements\n  please revise value list and try again')
        #end if
        if isinstance(variables,str):
            variables = [variables]
        #end if
        if not isinstance(variables,self.branch_types):
            self.error('branch failed\n  variables provided to branch must be in a list, tuple, or array\n  you provided a {0} for variables\n  please assemble branch variables in a list try again'.format(variables.__class__.__name__))
        #end if 
        for var in variables:
            if not isinstance(var,str):
                self.error('branch failed\n  variables provided to branch must be names (strings)\n  you provided a {0} with value {1}\n  please provide a valid string instead and try again'.format(var.__class__.__name__,var))
            elif not var in self.valid_variables:
                names = ''
                for name in self.valid_variables.keys():
                    names += '\n    '+name
                #end for
                self.error('branch failed\n  variable provided to branch is not a valid name\n  variable name provided: {0}\n  please select a valid variable name and try again\n  the list of valid names provided below may be of some help:'.format(var,names))
            #end if
        #end for
        nvars = len(variables)
        for iv in range(len(values)):
            vals = values[iv]
            if nvars>1 and len(vals)!=nvars:
                self.error('branch failed\n  number of values provided at index {0} of values list does not match the number of variables\n  number of variables: {1}\n  number of values: {2}\n  index in value list\n  variable names: {3}\n  values provided: {4}\n  please provide {1} values instead of {2} at value list index {0} and try again'.format(iv,len(variables),len(vals),variables,vals))
            #end if
        #end for

        # now that the input is (hopefully) clean, proceed with the branch
        # identify which elements in the selection depend on variables requested
        varying = obj()
        for var in variables:
            varying[var] = obj()
        #end for
        for c in self.selection:
            name = c.state.name
            if name in self.allowed_paths:
                allowed = self.allowed_paths[name]
                for var in variables:
                    if var in allowed:
                        varying[var][c.state.eid] = c
                    #end if
                #end for
            #end if
        #end for
        all_varying = []
        all_ids = set()
        for var,elems in varying.iteritems():
            for eid,c in elems.iteritems():
                if not eid in all_ids:
                    all_ids.add(eid)
                    all_varying.append(c)
                #end if
            #end for
        #end for
        del all_ids
        for c in all_varying:
            c.mark_down('clustered',True)
        #end for
        ids = set()
        for c in self.cascade.all:
            if c.clustered:
                ids.add(c.state.eid)
            #end if
        #end for

        # find all elements downstream of any varying element and group them into clusters
        clusters = obj()
        cluster_id = 0
        for nucleus in all_varying:
            if len(ids)>0:
                cluster = CascadeElementCollection()
                nucleus.spread_cluster_id(cluster_id,cluster)
                if len(cluster)>0:
                    for c in cluster:
                        ids.remove(c.state.eid)
                    #end for
                    clusters[cluster_id] = cluster
                    cluster_id+=1
                #end if
            #end if
        #end for
        if len(ids)>0:
            self.error('branch failed\n  not all clusters were identified\n  this is a bug, please contact the developer')
        #end if

        # graph the clusters, if requested
        if graph_clusters:
            self.graph(label_append=' clusters',**graph_args)
        #end if

        # branch the clusters and assign a new cascade id if they are disjoint
        nbranches = len(values)
        vset = obj()
        for iv in range(nbranches):
            # pack the values
            vals = values[iv]
            if nvars==1:
                vset[variables[0]] = vals
            else:
                for ivar in range(nvars):
                    vset[variables[ivar]] = vals[ivar]
                #end for
            #end if
            # assign values into the varying elements
            for var,val in vset.iteritems():
                for c in varying[var]:
                    v = c.variables
                    try:
                        exec('v.{0} = val'.format(var))
                    except Exception,e:
                        self.error('branch failed\n  the following exception was encountered during variable assignment:\n{0}\n  this is a bug, please contact the developer\n'.format(e))
                    #end try
                    if isinstance(val,(int,float,bool,str)):
                        vval = val
                    else:
                        vval = '*'
                    #end if
                    c.add_variation(var,'{0}={1}'.format(var,vval))
                #end for
            #end for
                    
            # the last set of variables is assigned only, no branch is created
            if iv==nbranches-1:
                break
            #end if

            # create a branch of each cluster
            for cluster in clusters:
                branch = CascadeElementCollection()
                cluster.make_clone()
                cluster.link_clone()
                cluster.yield_clone(branch)
                self.cascade.add_cluster(branch)
                branch.branch() # set the branched flag of the branch
            #end for
        #end for
        for cluster in clusters:
            cluster.branch() # set the branched flag of the branch
        #end for

        # graph branches, if requested
        if graph_branches:
            self.graph(label_append=' branches',**graph_args)
        #end if

        # update the selection to account for the branches
        self.selection.clear()
        self.selection.collect_selected(*self.cascade.all)

        # dissolve the clusters
        self.cascade.all.uncluster()
        self.cascade.all.unbranch()
    #end def branch
#end class Cascade




if __name__=='__main__':
    from project import generate_physical_system
    from project import generate_pwscf,generate_pw2qmcpack
    from project import generate_qmcpack,linear,loop,vmc
    from project import Job



    tests = obj(
        sim_like     = 0,
        multi_branch = 0,
        relax        = 1
        )



    if tests.sim_like:
        system = vdict(generate_physical_system)
        job    = vdict(Job)

        c = Cascade(
            elements = obj(
                scf = generate_pwscf,
                p2q = generate_pw2qmcpack,
                vmc = generate_qmcpack
                ),

            dependencies = obj(
                p2q = obj( orbitals='scf' ),
                vmc = obj( orbitals='p2q' )
                ),

            variables = obj(
                system = system(
                    structure = 'dimer',
                    cell      = 'big',
                    constants = 5.6
                    ),

                scf = obj(
                    pretend = True,
                    a = 1,
                    b = 2,
                    hubbard_u = obj( Cu = 3.5 ),
                    pseudos = ['Cu.ncpp','O.xml'],
                    system  = link('system'),
                    job     = job(nodes=4,hours=5,app_name='pw.x')
                    ),

                p2q = obj(
                    param = 5,
                    job = job(cores=1,hours=1)
                    ),

                vmc = obj(
                    pseudos   = ['Cu.xml','O.xml'],
                    system    = link('system'),
                    job       = job(nodes=20, threads=8,
                                    hours=2,  minutes=30,
                                    app_name='qmcapp'),
                    qmc_calcs = [
                        loop(max = 4,
                             qmc = linear(
                                energy   = 1.0,
                                timestep = 0.5,
                                samples  = 16000
                                )
                             ),
                        loop(max = 4,
                             qmc = linear(
                                unreweightedvariance = 1.0,
                                timestep = 0.5,
                                samples  = 16000
                                )
                             ),
                        vmc(
                            blocks   = 200,
                            steps    = 10,
                            timestep = 0.5
                            )
                        ]
                    )

                )

            )


        #print c.variable_paths
        #exit()

        demo = 0

        if demo:
            print
            print 'untouched graph'
            c.select('none')
            c.graph()

            print
            print 'p2q selected'
            c.select('v.param==5')
            c.graph()

            print
            print 'p2q blocked'
            c.block()
            c.select('none')
            c.graph()

            print
            print 'p2q unblocked'
            c.select('v.param==5')
            c.unblock()
            c.select('none')
            c.graph()

            print
            print 'p2q deleted'
            c.select('v.param==5')
            c.delete()
            c.graph()
        #end if

        c.select('v.param==5')
        c.branch('param',[0,1,2],graph_clusters=0)

        exit()



        c.elements.p2q.select()
        #c.elements.p2q.clustered=True
        c.elements.p2q.block()

        c.graph()
    #end if


    if tests.multi_branch:
        def f(x):
            print x
        #end def

        c = Cascade(
            description = 'Multi-Branch Example',
            elements = obj(
                a = f,
                b = f,
                c = f,
                d = f, 
                e = f,
                f = f,
                g = f,
                h = f,
                i = f,
                j = f,
                k = f,
                l = f,
                m = f,
                n = f,
                o = f
                ),
            dependencies = obj(
                c = obj(a='a'),
                d = obj(a='a',b='b'),
                e = obj(b='b'),
                f = obj(c='c'),
                g = obj(d='d'),
                h = obj(d='d'),
                i = obj(e='e'),
                j = obj(e='e'),
                k = obj(f='f'),
                l = obj(f='f'),
                m = obj(h='h'),
                n = obj(i='i'),
                o = obj(i='i')
                ),
            variables = obj(
                a = obj(a=1,s=obj(s11=obj(v1=1,v2=2),s12=obj(v3=3,a=6))),
                b = obj(b=1),
                c = obj(c=1),
                d = obj(d=1,x=1),
                e = obj(e=1),
                f = obj(f=1),
                g = obj(g=1),
                h = obj(h=1,y=0),
                i = obj(i=1,x=1),
                j = obj(j=1),
                k = obj(k=1),
                l = obj(l=1,x=1),
                m = obj(m=1,x=1),
                n = obj(n=1),
                o = obj(o=1),
                )
            )

        #print
        #print 'bare graph'
        #c.graph()
        #
        #print
        #print 'selecting x'
        #c.select('v.x==1')
        #c.graph()


        print
        print 'branching on x'
        #c.graph()
        #c.select('all')
        c.select('s.name=="d"',direction='down',graph=1)
        #c.branch('x',[1,2,3],graph_clusters=1,graph_branches=1,pause=1)
        #c.branch('y',[1,2],graph_clusters=1,graph_branches=1,pause=1)
        c.branch(('x','y'),[(1,3),(2,2),(3,1)],graph_clusters=1,graph_branches=1,pause=1)
        c.graph()


    #end if


    if tests.relax:
        c = Cascade(
            description = 'Relax Example',
            elements = obj(
                relax = generate_pwscf,
                scf   = generate_pwscf,
                p2q   = generate_pw2qmcpack,
                opt   = generate_qmcpack,
                dmc   = generate_qmcpack
                ),
            dependencies = obj(
                scf = obj(structure='relax'),
                p2q = obj(orbitals='scf'),
                opt = obj(structure='relax',orbitals='p2q'),
                dmc = obj(structure='relax',orbitals='p2q',jastrow='opt')
                ),
            variables = obj(
                relax = obj(s='s1'),
                scf   = obj(func='lda'),
                p2q   = obj(cores=1),
                opt   = obj(samples=16000),
                dmc   = obj(timestep=.01)
                )
            )

        #print c.variable_paths
        #exit()

        #c.graph()


        c.select('all')
        c.branch('s',['s1','s2'],graph_clusters=1,graph_branches=1)

        c.branch('func',['lda','pbe'],graph_clusters=1,graph_branches=1)



        c.select('v.func=="lda"',direction='down',graph=0)
        c.select('s.name=="opt"',cont=1,graph=1)
        c.delete()
        c.graph(s=['eid'])
        
    #end if

# next steps
#   check that branch w/ selection works
#     (select only middle branch of Multi-Branch example)
#   check that multivariate branching works
#   test Cascade.run()
#     do this with a minimal and unbranched real workflow w/ generate only
#   test select/branch in realistic inputs scenario
#     real workflow
#     branch and select a few times on nested/linked data
#     check that generate only works
#     check that runs actually complete
#     check that further cascade modifying calls fail
#     perhaps should expose generate_simulations to user so run() can be commented
#       in or out w/o affecting code below 
#       (will still need to load existing simulations... how to do this?)
#   evaluate vdicts when make_sim is called
#     init_variables needs to store a data structure spanning all element types
#       that contains paths to vdicts
#     the data structure should be ordered by depth so that high-depth vdicts
#       are appropriately evaluated first
#     the appropriate subset of paths needs to be provided to make_sim so that
#       it can transform all vdicts to their underlying types in the right order
#  implement show_variables
#  test delete
#  implement merge
#  implement chain
#  implement link
#

#end if
