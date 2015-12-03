##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  generator.py                                                      #
#    Generates database of objects with simple select and vary       #
#    operations.                                                     #
#                                                                    #
#  Content summary:                                                  #
#    Generator                                                       #
#      Main class that manages and expands the object database.      #
#                                                                    #                                        
#====================================================================#


#! /usr/bin/env python

from types import FunctionType,StringType,DictType
from developer import DevBase


class Gobj(DevBase):
    None
#end class Gobj


class Entry(Gobj):
    None
#end class Entry


class Generator(Gobj):
    def __init__(g,database=None):
        g.vars = Gobj()
        g.database = [] #a collection of entry objects, which contain fields
        g.selection = None #this will be a Generator when created by select
        g.parent = None
        g.alias_paths = None #mapping of entry variables to target object paths (a=obj.obj1.v)
        g.modes = Gobj()
        g.modes.exit_on_fail = True
        if database is None:
            database=[Entry()]
        #end if
        for e in database:
            g.database.append(e)
        #end for
    #end def __init__
        

    def exit(g):
        if g.modes.exit_on_fail:
            exit()
        #end if
    #end def exit


    def error(g,message):
        print 'Generator Error: '+message+'\n'
    #end def error

        
    def setvars(g,*args,**assignments):
        if len(args)>0 and type(args[0])==DictType:
            assign=args[0]
        else:
            assign = assignments
        #end if
        for k,v in assign.iteritems():
            g.vars[k]=v
        #end for
        return g
    #end def set


    def delvars(g,*varlist):
        for varname in varlist:
            del g.vars[varname]
        #end for
        return g
    #end def delvars


    def showvars(g):
        print 'Variables:'
        print g.vars
        return g
    #end def showvars        


    def select(g,inpt='skipped',**equalities):
        inp = inpt
        if g.selection!=None:
            del g.selection
        #end if
        if inp=='skipped':
            inp=''
            for k,v in equalities.iteritems():
                if type(v)==StringType:
                    val = '"'+v+'"'
                else:
                    val = str(v)
                #end if
                inp += 'e.'+k+'=='+val+' and '
            #end for
            inp = inp[:-5]
        #end if
        inptype = type(inp)
        if inp is None or inp=='none':
            e = Entry()
            selection = [e]
            g.expand(selection)
        #end if
        elif inp=='all':
            selection = g.database
        else:
            selection = []
            if inptype==FunctionType:
                f=inp
                for e in g.database:
                    if f(e):
                        selection.append(e)
                    #end if
                #end for
            elif inptype==StringType:
                rule = inp
                v = g.vars
                for e in g.database:
                    exec 'ein = '+rule
                    if ein:
                        selection.append(e)
                    #end if
                #end for
            else:
                g.error('Invalid input type for function select: '+inptype)
                g.error('  valid types are: '+str(Types.function)+' '+str(Types.string))
                g.exit()
            #end if
        #end if
        g.selection = Generator(selection)
        g.selection.parent = g
        return g
    #end def select


    def check_select(g):
        if g.selection is None:
            g.select('all')
        #end if
    #end def check_select


    def set(g,*actions,**assignments):
        g.check_select()
        if len(actions)>0:
            v = g.vars
            for e in g.selection.database:
                for a in actions:
                    exec a
                #end for
            #end for                    
        else:                
            for e in g.selection.database:
                for k,val in assignments.iteritems():
                    e[k]=val
                #end for
            #end for
        #end if
        return g
    #end def set


    def delete(g,*varlist):
        g.check_select()
        for e in g.selection.database:
            for varname in varlist:
                del e[varname]
            #end for
        #end for
        return g
    #end def delete


    def generate(g,**covariables):
        g.check_select()
        first = True
        same_length = True
        for k,v in covariables.iteritems():
            if first:
                clen = len(v)
                first = False
            #end if
            same_length = same_length and len(v)==clen
        #end for
        if not same_length:
            g.error('Arguments to generate (covariables) must have the same length.')
            g.exit()
        #end if
        if clen>0:
            for e in g.selection.database:
                for k,v in covariables.iteritems():
                    e[k]=v[0]
                #end for
            #end for
            expansion = []
            for i in range(1,clen):
                for eo in g.selection.database:
                    e = eo.copy()
                    for k,v in covariables.iteritems():
                        e[k]=v[i]
                    #end for
                    expansion.append(e)
                #end for
            #end for
            g.selection.expand(expansion)
        #end if
        return g
    #end def generate


    def expand(g,expansion):
        g.database.extend(expansion)
        if g.parent!=None:
            g.parent.expand(expansion)
        #end if
    #end def expand


    def show(g,*args):
        i=-1
        if len(args)==0:
            for e in g.selection.database:
                i+=1
                print 'Entry '+str(i)+':'
                print e
            #end for
        else:
            for e in g.selection.database:
                i+=1
                print 'Entry '+str(i)+':'
                for a in args:
                    print '  '+a+' = '+str(e[a])
                #end for
                print
            #end for            
        #end if
        return g
    #end def show           


    def generate_objects(g,obj,alias_paths=None):
        if alias_paths!=None:
            g.alias_paths={}
            for k,v in alias_paths.iteritems():
                g.alias_paths[k] = v.split('.')
            #end for
        #end if
        if g.alias_paths!=None:
            objects = []
            for e in g.selection.database:
                o = obj.copy()
                for k,vpath in g.alias_paths.iteritems():
                    v = o
                    for i in range(len(vpath)-1):
                        v = v[vpath[i]]
                    #end for
                    vloc = vpath[-1]
                    v[vloc] = e[k]
                #end for
            #end for
        #end if
        return g
    #end def generate_objects            
#end class Generator
