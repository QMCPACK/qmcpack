##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  qmcpack_variations.py                                             #
#    Originally intended to support QMCPACK benchmark runs that      #
#    require input variations.  Used for a short time, now obsolete. #
#                                                                    #
#====================================================================#


import os
from datetime import datetime
# jtk library imports
from generic import obj                                   # generic.py
from project import Qmcpack,QmcpackInput,QmcpackAnalyzer  # project.py
from project import Job,nexus_core
from developer import DevBase                             # developer.py
from debug import *                                       # debug.py



class QmcpackVariations(DevBase):
    simulations = []
    machine = None
    build_path = None
    app_loc    = None

    build_vars  = set('revision compiler wftype'.split())
    job_vars    = set('nodes cores processes threads hours minutes queue options'.split())
    inflag_vars = set('async_swap'.split())


    @classmethod
    def settings(cls,
                 machine       = None,
                 build_path    = None,
                 app_loc       = 'bin/qmcapp',
                 generate_only = False):
        cls.machine       = machine
        cls.build_path    = build_path
        cls.app_loc       = app_loc
        cls.generate_only = generate_only
    #end def settings

    
    def __init__(self,name,source,variations,
                 build_order        = None,
                 variable_order     = None,
                 build_delimiter    = '_',
                 variable_delimiter = '__',
                 machine            = None,
                 build_path         = None,
                 app_loc            = None):
        if not isinstance(name,str):
            self.error('name must be a string\n  you provided: '+str([name]))
        #end if
        if build_order is None:
            self.error('build order must be provided')
        elif isinstance(build_order,str):
            build_order = build_order.split()
        #end if
        if variable_order is None:
            self.error('variable order must be provided')
        elif isinstance(variable_order,str):
            variable_order = variable_order.split()
        #end if
        if machine!=None:
            self.machine = machine
        #end if
        if build_path!=None:
            self.build_path = build_path
        #end if
        if app_loc!=None:
            self.app_loc = app_loc
        #end if
        self.set(
            name               = name,
            source             = source,
            variations         = variations,
            build_order        = build_order,
            variable_order     = variable_order,
            build_delimiter    = str(build_delimiter),
            variable_delimiter = str(variable_delimiter),
            source_info        = obj(),
            sims               = obj()
            )
        self.get_source()
        self.make_variations()
    #end def __init__

    
    def get_source(self):
        source = self.source
        if isinstance(source,str) and os.path.exists(source):
            self.get_source_from_directory()
        else:
            self.error('getting source from {0} has not yet been implemented'.format(source.__class__.__name__))
        #end if
    #end def get_source


    def get_source_from_directory(self):
        source = self.source
        if os.path.isfile(source):
            filepath = source
            path,file = os.path.split(filepath)
        else:
            path = source
            file = 'qmc.xml'
            filepath = os.path.join(path,file)
        #end if

        # read the input file template
        qi = QmcpackInput(filepath)

        # find the wavefunction hdf5 file, if present
        wf_h5file = None
        ds = qi.get('determinantset')
        if ds!=None and 'href' in ds:
            wf_h5file = ds.href
        #end if

        # find the pseudopotential files, if present
        pp_files = qi.get_pp_files()

        # create a representation of the physical system 
        # from the input file template
        system = qi.return_system()

        # place an xyz file of the structure in the source directory
        if system!=None:
            system.structure.write_xyz(os.path.join(path,'structure.xyz'))
        #end if

        self.source_info.set(
            type      = 'template',
            input     = qi,
            wf_h5file = wf_h5file,
            pp_files  = pp_files,
            system    = system,
            path      = path,
            infilepath= filepath,
            infile    = file
            )
    #end def get_source_from_directory

    
    def make_variations(self):
        source      = self.source
        source_info = self.source_info
        build_order = self.build_order
        variable_order = self.variable_order
        build_delimiter = self.build_delimiter
        variable_delimiter = self.variable_delimiter

        boset = set(build_order)
        voset = set(variable_order)

        invalid = boset - self.build_vars
        if len(invalid)>0:
            self.error('variables provided in build_order are invalid\n  invalid entries: '+str(list(invalid))+'\b  valid options are: '+str(list(self.build_vars)))
        #end if

        #assume template input file source for now
        qpinput,wf_h5file,pp_files,system,source_path = self.source_info.tuple('input','wf_h5file','pp_files','system','path')

        #get the current date and time
        now = datetime.now()
        date_time = datetime.strftime(now,'%y%m%d_%H%M%S')


        errors = False

        for variation in self.variations:
            error = False
            vars = set(variation.keys())
            build_vars  = vars & self.build_vars
            job_vars    = vars & self.job_vars
            inflag_vars = vars & self.inflag_vars
            input_vars  = vars-build_vars-job_vars-inflag_vars

            invars = inflag_vars | input_vars

            if build_vars<boset:
                self.error('not enough build variables in variation\n  build variables required: '+str(build_order),exit=False,trace=False)
                error = True
            #end if
            if invars<voset:
                self.error('not enough input variables in variation\n  input variables required: '+str(variable_order),exit=False,trace=False)
                error = True
            #end if

            # determine the build directory
            build_list = []
            build_dir  = ''
            for var in build_order:
                if var in variation:
                    val = variation[var]
                    build_list.append(val)
                    if var=='revision' and isinstance(val,int):
                        val = 'r'+str(val)
                    #end if
                    build_dir += str(val)+build_delimiter
                else:
                    self.error('variable '+str(var)+' in build_order is not in variation',exit=False,trace=False)
                    error = True
                #end if
            #end for
            if len(build_list)>0:
                build_dir = build_dir[:-len(build_delimiter)]
            #end if

            # determine the input var directory
            var_list = []
            var_dir  = ''
            for var in variable_order:
                if var in variation:
                    val = variation[var]
                    var_list.append(val)
                    var_dir += '{0}-{1}{2}'.format(var,val,variable_delimiter)
                else:
                    self.error('variable '+str(var)+' in variable_order is not in variation',exit=False,trace=False)
                    error = True
                #end if
            #end for
            if len(var_list)>0:
                var_dir = var_dir[:-len(variable_delimiter)]
            #end if

            if error:
                self.error('  variation provided:\n'+str(variation),exit=False,trace=False)
            #end if

            # make the job object 
            jin = variation.obj(*job_vars)
            if len(inflag_vars)>0:
                inflags = ''
                for var in inflag_vars:
                    val = variation[var]
                    if val:
                        inflags+=' --'+var
                    #end if
                #end for
                jin.app_flags = inflags
            #end if
            job = Job(**jin)

            # set the simulation path
            path = os.path.join(self.name,self.machine,var_dir,build_dir+build_delimiter+date_time)

            # fill out the template input file
            #  assume simple and unique find and replace variations for now
            input = qpinput.copy()
            assignments = variation.obj(*input_vars)
            input.assign(**assignments)

            #  add the relative path location of the wavefunction file
            runpath = os.path.join(nexus_core.local_directory,nexus_core.runs,path)

            wftmp  = os.path.join(source_path,wf_h5file)
            wfpath = os.path.relpath(wftmp,runpath)

            ds = input.get('determinantset')
            if ds!=None:
                ds.href = wfpath
            #end if

            # check that the build exists
            app_loc = os.path.join(self.build_path,build_dir,self.app_loc)
            if not os.path.exists(app_loc) and not nexus_core.generate_only:
                print '    Error: no qmcapp at '+app_loc
                error = True
            #end if
            job.app_name = app_loc

            # locate pp files
            files = []
            if pp_files!=None:
                files = list(pp_files)
            #end if

            # make the simulation object
            qmc = Qmcpack(
                identifier = 'qmc',
                path       = path,
                job        = job,
                system     = system,
                input      = input,
                files      = files 
                )
            self.sims[tuple(var_list+build_list)] = qmc

            # add it to the global simulation list
            self.simulations.append(qmc)

            errors = errors or error
        #end for

        if errors:
            self.error('errors encountered, see above')
        #end if
    #end def make_variations

#end class QmcpackVariations
