##################################################################
##  (c) Copyright 2019-  by Jaron T. Krogel                     ##
##################################################################


"""
Track versions of dependencies supported by Nexus.
"""

import sys
from platform import python_version,python_version_tuple
import datetime
import importlib


nexus_version = 2,1,0
"""
Current Nexus version.
"""

python_supported = 'python3'
"""
Current Python family supported.
"""

# Require Python 3
if python_version_tuple()<('3','0','0'):
    print('\nNexus is compatible only with Python 3.\n  You attempted to run with Python version {}.\n  Please rerun with Python 3.\n'.format(python_version()))
    sys.exit(1)
#end if

years_supported = 2
"""
Policy for how many years back Nexus will extend support to dependencies.
"""

support_cutoff_date = (datetime.datetime.now()-datetime.timedelta(days=365*years_supported)).date()
"""
Cutoff date for support.
"""

current_date = datetime.datetime.now().date()
"""
Current date.
"""


required_dependencies = set([python_supported,'numpy'])
"""
Required dependencies for Nexus.
"""


currently_supported = [
    ('python3'    , (3,  6,  0) ),
    ('numpy'      , (1, 13,  1) ),
    ('scipy'      , (0, 19,  1) ), # optional
    ('h5py'       , (2,  7,  1) ), # optional
    ('matplotlib' , (2,  0,  2) ), # optional
    ('pydot'      , (1,  2,  3) ), # optional
    ('spglib'     , (1,  9,  9) ), # optional
    ('seekpath'   , (1,  4,  0) ), # optional
    ('pycifrw'    , (4,  3,  0) ), # optional
    ('cif2cell'   , (1,  2, 10) ), # optional
    ]
"""
Currently supported versions of Nexus dependencies.
"""


import_aliases  = dict(
    pycifrw = 'CifFile',
    )

raw_version_data = dict(
    # python 3 releases
    #   3.3  https://www.python.org/dev/peps/pep-0398/
    #   3.4  https://www.python.org/dev/peps/pep-0429/
    #   3.5  https://www.python.org/dev/peps/pep-0478/
    #   3.6  https://www.python.org/dev/peps/pep-0494/
    #   3.7  https://www.python.org/dev/peps/pep-0537/
    #   3.8  https://www.python.org/dev/peps/pep-0569/
    #   3.9  https://www.python.org/dev/peps/pep-0596/
    #   3.10 https://www.python.org/dev/peps/pep-0619/
    python3 = '''
        3.3.0  2012-09-29
        3.3.1  2013-04-06
        3.3.2  2013-05-13
        3.3.3  2013-11-16
        3.3.4  2014-02-09
        3.4.0  2014-03-16
        3.4.1  2014-05-18
        3.4.2  2014-10-06
        3.4.3  2015-02-25
        3.5.0  2015-09-13
        3.5.1  2015-12-06
        3.5.2  2016-06-26
        3.6.0  2016-12-23
        3.6.1  2017-03-21
        3.6.2  2017-07-17
        3.6.3  2017-10-03
        3.6.4  2017-12-19
        3.6.5  2018-03-28
        3.7.0  2018-06-27
        3.7.1  2018-10-20
        3.7.2  2018-12-24
        3.7.3  2019-03-25
        3.7.4  2019-07-08
        3.8.0  2019-10-21
        3.9.0  2020-06-08
        3.10.0 2021-10-04
        ''',
    # numpy releases
    #   https://github.com/numpy/numpy/releases
    numpy = '''
        1.10.0  2015-10-05
        1.10.1  2015-10-12
        1.10.2  2015-12-14
        1.10.3  2016-01-06
        1.10.4  2016-01-06
        1.11.0  2016-03-27
        1.11.1  2016-06-25
        1.11.2  2016-10-03
        1.11.3  2016-12-18
        1.12.0  2017-01-15
        1.12.1  2017-03-18
        1.13.0  2017-06-07
        1.13.1  2017-07-06
        1.13.2  2017-09-27
        1.13.3  2017-09-29
        1.14.0  2018-01-06
        1.14.1  2018-02-20
        1.14.2  2018-03-12
        1.14.3  2018-04-28
        1.14.4  2018-06-06
        1.14.5  2018-06-12
        1.15.0  2018-07-23
        1.15.1  2018-08-21
        1.15.2  2018-09-23
        1.15.3  2018-10-22
        1.15.4  2018-11-04
        1.16.0  2019-01-13
        1.16.1  2019-01-31
        1.16.2  2019-02-26
        1.16.3  2019-04-21
        1.16.4  2019-05-28
        1.17.0  2019-07-26
        1.17.1  2019-08-27
        ''',
    # scipy releases
    #   https://github.com/scipy/scipy/releases
    scipy = '''
        0.16.0  2015-07-24
        0.16.1  2015-10-24
        0.17.0  2016-01-23
        0.17.1  2016-05-12
        0.18.0  2016-07-25
        0.18.1  2016-09-19
        0.19.0  2017-03-09
        0.19.1  2017-06-23
        1.0.0   2017-10-25
        1.0.1   2018-03-24
        1.1.0   2018-05-05
        1.2.0   2018-12-17
        1.2.1   2019-02-09
        1.3.0   2019-05-17
        1.3.1   2019-08-08
        ''',
    # h5py releases
    # https://pypi.org/project/h5py/#history
    h5py = '''
        2.4.0  2015-01-05
        2.5.0  2015-04-08
        2.6.0  2016-04-08
        2.7.0  2017-03-18
        2.7.1  2017-09-01
        2.8.0  2018-05-13
        2.9.0  2018-12-19
        ''',
    # matplotlib releases
    #  https://pypi.org/project/matplotlib/#history
    matplotlib = '''
        1.5.0  2015-10-29
        1.5.1  2016-01-10
        1.5.2  2016-08-18
        1.5.3  2016-09-09
        2.0.0  2017-01-17
        2.0.1  2017-05-02
        2.0.2  2017-05-10
        2.1.0  2017-10-07
        2.1.1  2017-12-11
        2.2.0  2018-03-06
        2.2.2  2018-03-17
        2.2.3  2018-08-11
        3.0.0  2018-09-18
        3.0.1  2018-10-25
        3.0.2  2018-11-10
        3.0.3  2019-02-28
        3.1.0  2019-05-18
        3.1.1  2019-07-02
        ''',
    # pydot
    #   https://pypi.org/project/pydot/#history
    pydot = '''
        1.0.28  2012-01-02
        1.0.29  2016-05-17
        1.1.0   2016-05-23
        1.2.0   2016-07-01
        1.2.1   2016-07-01
        1.2.2   2016-07-01
        1.2.3   2016-10-06
        1.2.4   2017-12-25
        1.3.0   2018-11-19
        1.4.0   2018-12-01
        1.4.1   2018-12-12
        ''',
    # spglib releases
    #   https://pypi.org/project/spglib/#history
    spglib = '''
        1.6.0   2014-05-20
        1.8.3   2015-12-15
        1.9.3   2016-05-11
        1.9.5   2016-09-15
        1.9.6   2016-10-17
        1.9.7   2016-10-19
        1.9.8   2016-11-01
        1.9.9   2016-12-14
        1.9.10  2017-10-02
        1.10.0  2017-10-21
        1.10.1  2017-10-27
        1.10.2  2017-12-12
        1.10.3  2018-01-13
        1.10.4  2018-08-01
        1.11.0  2018-11-08
        1.11.1  2018-11-12
        1.11.2  2018-12-07
        1.12.0  2019-01-29
        1.12.1  2019-02-01
        1.12.2  2019-02-05
        1.13.0  2019-07-02
        1.14.0  2019-07-30
        1.14.1  2019-07-30
        ''',
    # pycifrw releases
    #   https://pypi.org/project/PyCifRW/#history
    pycifrw = '''
        4.1    2015-01-03
        4.2    2016-03-01
        4.3    2017-02-27
        4.4    2018-02-12
        4.4.1  2019-02-27
        ''',
    # cif2cell releases
    #   https://pypi.org/project/cif2cell/#history
    #   https://sourceforge.net/projects/cif2cell/files/
    cif2cell = '''
        1.2.2   2014-08-26
        1.2.7   2015-05-18
        1.2.10  2016-01-19
        2.0.0   2019-10-07
        ''',
    # seekpath releases
    #   https://pypi.org/project/seekpath/#history
    seekpath = '''
        1.0.0  2016-09-29
        1.0.1  2016-09-29
        1.0.2  2016-10-10
        1.1.0  2016-10-26
        1.1.1  2016-11-21
        1.2.0  2016-12-19
        1.3.0  2017-01-20
        1.4.0  2017-04-04
        1.5.0  2017-09-21
        1.6.0  2017-10-13
        1.7.0  2017-11-30
        1.7.1  2017-11-30
        1.8.0  2017-12-22
        1.8.1  2018-02-12
        1.8.2  2018-07-17
        1.8.3  2018-10-04
        1.8.4  2018-10-04
        1.9.0  2019-07-25
        1.9.1  2019-07-29
        1.9.2  2019-07-30
        1.9.3  2019-09-03
        ''',
    )
"""
Version and date information for Nexus dependencies.
"""


class VersionError(Exception):
    """
    (`Internal API`) Exception unique to versions module.
    """
    None
#end class VersionError


def version_error(msg):
    """
    (`Internal API`) Print an error message and raise VersionError.

    Parameters
    ----------
    msg : `str`
        The error message to be printed.
    """
    print('\nVersion error:\n  '+msg.replace('\n','\n  ')+'\nexiting.\n')
    raise(VersionError)
#end def version_error


def time_ago(date):
    """
    (`Internal API`) Compute time in years from current date.

    Parameters
    ----------
    date : `datetime.date`
        A date prior to today's date.

    Returns
    -------
    years_ago : `float or int`
        Number of years separating today's date from date of interest.
        Floating point value includes months as fraction.
    """
    delta = current_date-date
    return float(delta.days)/365
#end def time_ago



def process_version(*version):
    """
    (`Internal API`) Process version number from a range of inputted formats.

    Parameters
    ----------
    *version : `int or str or tuple(str) or tuple(int)
        Version number as string (e.g. `"1.0.0"`) tuple ( `('1','0','0')` or 
        `(1,0,0)` or `(1,0)`, etc.) or int (e.g. `1`). 

    Returns
    -------
    version : `tuple(int)`
        Version tuple with at least 3 entries.
    """
    if len(version)==1:
        version = version[0]
    #end if
    if isinstance(version,str):
        version = version.replace('+','').replace('-','')
        version = version.split('.')
    elif isinstance(version,int):
        version = (version,)
    #end if
    version = tuple(map(int,version))
    while len(version)<3:
        version += (0,)
    #end while
    return version
#end def process_version


def version_to_string(version):
    """
    (`Internal API`) Convert version tuple to string.
    """
    vs = '{}.{}.{}'.format(*version)
    return vs
#end def version_to_string



class Versions(object):
    """
    Handles version information for Nexus dependencies.

    Used to produce printed summaries of currently supported versions 
    for Nexus dependencies as well as versions supported by the current 
    age policy.  Also used to print summarized information about versions 
    of dependencies detected on the current machine and recommendations 
    for Nexus users regarding the state of the their installation.

    The main data items for consumption by this class are 
    `support_cutoff_date`, `required_dependencies`, `currently_supported`, 
    and `raw_version_data`.

    This class is not intended for direct use by Nexus users.  See the 
    three main API functions `check_versions`, `current_versions`, and 
    `policy_versions` below.

    This class follows a singleton pattern.
    """

    instances = []

    def __init__(self):
        self.initialize()
        if len(Versions.instances)>0:
            self.error('Only one Versions object can be instantiated!')
        #end if
        Versions.instances.append(self)
    #end def __init__


    def initialize(self):
        """
        (`Internal API`)  Parse version data and determine status of 
        module dependencies on the current machine via attempted import.
        """

        self.__dict__.clear()

        # parse raw version data
        version_data  = dict()
        dependencies  = set()
        full_names    = dict()
        for name,vd_string in raw_version_data.items():
            vd = dict()
            n = 0
            for line in vd_string.strip().splitlines():
                version_string,date_string = line.split()
                version = process_version(version_string)
                date = datetime.date(*map(int,date_string.split('-')))
                vd[version] = dict(
                    version_string = version_string,
                    version        = version,
                    date_string    = date_string,
                    date           = date,
                    )
                n+=1
            #end for
            full_names[name] = name.lower()
            name = name.lower()
            version_data[name] = vd
            if not name.startswith('python'):
                dependencies.add(name)
            #end if
        #end for
        dependencies.add(python_supported)

        cur_supp = dict()
        ordered_dependencies = []
        for name,version in currently_supported:
            ordered_dependencies.append(name)
            cur_supp[name] = version
        #end for

        optional_dependencies = dependencies-required_dependencies

        # store general dependency data
        self.full_names            = full_names
        self.version_data          = version_data
        self.ordered_dependencies  = ordered_dependencies
        self.dependencies          = dependencies
        self.required_dependencies = required_dependencies
        self.optional_dependencies = optional_dependencies
        self.currently_supported   = cur_supp

        # determine which dependencies are available
        dependency_available = dict()
        dependency_version   = dict()
        dependency_supported = dict()
        for name in self.ordered_dependencies:
            version   = None
            supported = None
            if not name.startswith('python'):
                module = None
                try:
                    module_name = name
                    if name in import_aliases:
                        module_name = import_aliases[name]
                    #end if
                    module = importlib.import_module(module_name)
                    available = True
                except:
                    available = False
                #end try
                if available:
                    if '__version__' in module.__dict__:
                        supported_version = self.currently_supported[name]
                        version = process_version(module.__version__)
                        supported = version >= supported_version
                    #end if
                #end if
            else:
                python_major_supported = int(python_supported[-1])
                version = process_version(python_version())
                available = version[0]==python_major_supported
                if available:
                    supported_version = self.currently_supported[name]
                    supported = version >= supported_version
                #end if
            #end if
            dependency_available[name] = available
            dependency_version[name]   = version
            dependency_supported[name] = supported
        #end for

        # store data regarding dependency availability and support
        self.dependency_available = dependency_available
        self.dependency_version   = dependency_version
        self.dependency_supported = dependency_supported
    #end def initalize


    def error(self,msg):
        """
        (`Internal API`) Raise an error.
        """
        version_error(msg)
    #end def error


    def is_dependency(self,name):
        """
        (`Internal API`) Return whether a name corresponds to a Nexus dependency.
        """
        return name.lower() in self.dependencies
    #end def is_dependency


    def check_dependency(self,name):
        """
        (`Internal API`) Check whether a name corresponds to a Nexus dependency, and raise an error if it doesn't. 
        """
        name = name.lower()
        if name not in self.dependencies:
            self.error('"{}" is not a dependency of Nexus'.format(name))
        #end if
        return name
    #end def check_dependency


    def available(self,name):
        """
        (`Internal API`) Return whether a Nexus dependency is available 
        on the current machine.
        """
        name = self.check_dependency(name)
        return self.dependency_available[name]
    #end def available


    def supported(self,name,version=None):
        """
        (`Internal API`) Return whether a Nexus dependency on the current 
        machine has a version that is currently supported.
        """
        name = self.check_dependency(name)
        if version is not None:
            version = process_version(version)
            supported_version = self.currently_supported[name]
            return version>=supported_version
        else:
            if not self.dependency_available[name]:
                return False
            else:
                return self.dependency_supported[name]
            #end if
        #end if
    #end def supported


    def policy_supported_version(self,name=None):
        """
        (`Internal API`) Determine versions of Nexus dependencies that 
        comply with the age policy relative to today's date.
        """
        if name is not None:
            name = name.lower()
            if name not in self.dependencies:
                self.error('"{}" is not a dependency of Nexus'.format(name))
            #end if
            vdata = self.version_data[name]
            versions = list(reversed(sorted(vdata.keys())))
            sv = None
            for version in versions:
                vd = vdata[version]
                sv = version
                if vd['date'] <= support_cutoff_date:
                    break
                #end if
            #end for
            return sv
        else:
            supported_versions = dict()
            for name in self.dependencies:
                sv = self.policy_supported_version(name)
                supported_versions[name] = sv
            #end for
            return supported_versions
        #end if
    #end def policy_supported_version


    def write_supplied_versions(self,supplied_versions,age=True,opt_req=False):
        """
        (`Internal API`) Write information about a set of dependencies.
        """
        write_age = age
        s = ''
        for name in self.ordered_dependencies:
            if name in supplied_versions:
                version = supplied_versions[name]
                vd = self.version_data[name][version]
                date_str = vd['date_string']
                date     = vd['date']
                age = time_ago(date)
                version = version_to_string(version)
                if write_age:
                    s += '  {:<12} = {:<8}  (dated: {} , age: {:2.1f} years )\n'.format(name,version,date_str,age)
                elif opt_req:
                    if name in self.required_dependencies:
                        s += '  {:<12} = {:<8}   (required)\n'.format(name,version)
                    else:
                        s += '  {:<12} = {:<8}   (optional)\n'.format(name,version)
                    #end if
                else:
                    s += '  {:<12} = {:<8}\n'.format(name,version)
                #end if
            #end if
        #end for
        return s
    #end def write_supplied_versions


    def write_current_versions(self,age=True,opt_req=False,header=None):
        """
        (`Internal API`) Write information about currently supported versions of Nexus dependencies. 
        """
        supported_versions = self.currently_supported
        if header is None:
            s = '\nNexus dependencies currently supported:\n'.format(years_supported)
        else:
            s = header
        #end if
        s += self.write_supplied_versions(supported_versions,age=age,opt_req=opt_req)
        return s
    #end def write_current_versions


    def write_policy_versions(self):
        """
        (`Internal API`) Write information about versions of Nexus dependencies current with the age policy as of today.
        """
        supported_versions = self.policy_supported_version()
        s = '\nNexus dependencies for {} year policy:\n'.format(years_supported)
        s += self.write_supplied_versions(supported_versions)
        return s
    #end def write_policy_versions


    def write_available_versions(self,status=True,opt_req=False):
        """
        (`Internal API`) Write information about versions of Nexus dependencies that are available on the current machine.
        """
        available_versions = self.dependency_version
        s = '\nNexus dependencies available on current machine:\n'.format(years_supported)
        for name in self.ordered_dependencies:
            if self.dependency_available[name]:
                version = available_versions[name]
                supported_version = self.currently_supported[name]
                if version is None:
                    version = '(unknown)'
                    support = '(unknown)'
                else:
                    if version >= supported_version:
                        support = 'supported'
                    else:
                        support = 'unsupported'
                    #end if
                    version = version_to_string(version)
                #end if
                supported_version = version_to_string(supported_version)
                if status:
                    s += '  {:<12} = {:<9}   status: {:<11}   oldest supported: {:<8}\n'.format(name,version,support,supported_version)
                elif opt_req:
                    if name in self.required_dependencies:
                        s += '  {:<12} = {:<9}   (required)\n'.format(name,version)
                    else:
                        s += '  {:<12} = {:<9}   (optional)\n'.format(name,version)
                    #end if
                else:
                    s += '  {:<12} = {:<9}\n'.format(name,version)
                #end if
            #end if
        #end for
        return s
    #end def write_available_versions


    def print_current_versions(self):
        """
        (`Internal API`) Print information about currently supported versions of Nexus dependencies. 
        """
        print(self.write_current_versions())
    #end def print_current_versions


    def print_policy_versions(self):
        """
        (`Internal API`) Print information about versions of Nexus dependencies current with the age policy as of today.
        """
        print(self.write_policy_versions())
    #end def print_policy_versions


    def print_available_versions(self,status=True):
        """
        (`Internal API`) Print information about versions of Nexus dependencies that are available on the current machine.
        """
        print(self.write_available_versions(status=status))
    #end def print_available_versions


    def check(self,write=True,exit=False,full=False,n=0,pad='  '):
        """
        (`Internal API`) Check whether all required Nexus dependencies are present on the current machine.  Optionally write detailed information for the Nexus user about the status of their installation.
        """
        serr = ''
        header ='\nChecking for Nexus dependencies on the current machine...\n'
        s = versions.write_available_versions(status=False,opt_req=True)
        wcv = '\nNexus dependencies recommended for full functionality:\n'
        s += versions.write_current_versions(age=False,opt_req=True,header=wcv)
        available_dependencies = set()
        unavailable_dependencies = set()
        req_supported    = []
        opt_supported    = []
        req_unsupported  = []
        opt_unsupported  = []
        req_unknown      = []
        opt_unknown      = []
        req_missing      = []
        opt_missing      = []
        for name in self.ordered_dependencies:
            required = name in self.required_dependencies
            if self.available(name):
                supported = self.supported(name)
                unknown = supported == None
                if unknown:
                    if required:
                        req_unknown.append(name)
                    else:
                        opt_unknown.append(name)
                    #end if
                elif not supported:
                    if required:
                        req_unsupported.append(name)
                    else:
                        opt_unsupported.append(name)
                    #end if
                else:
                    if required:
                        req_supported.append(name)
                    else:
                        opt_supported.append(name)
                    #end if
                #end if
            else:
                if required:
                    req_missing.append(name)
                else:
                    opt_missing.append(name)
                #end if
            #end if
        #end for

        all_req_supported = set(req_supported)==self.required_dependencies
        all_opt_supported = set(opt_supported)==self.optional_dependencies
        all_supported = all_req_supported and all_opt_supported

        #req_missing.append('numpy')
        if len(req_missing)>0:
            serr = 'Some required Nexus dependencies are missing.\nMissing dependencies: {}\nPlease check your Python installation.'.format(req_missing)
            s += '\nRequired dependencies are missing:\n'
            for name in req_missing:
                s += '  {} is missing.  Install {} or greater.\n'.format(name,version_to_string(self.currently_supported[name]))
            #end for
            s += '\nNexus will not work.\n  Please install the missing dependencies above.\n'
        else:
            if all_supported:
                s += '\nAll Nexus dependencies are met.\n  Both core workflow and optional features should work.\n'
            elif all_req_supported:
                s += '\nAll required Nexus dependencies are met.\n  Core workflow features should work.\n  Some optional features may not.\n  See below for more information.\n'
            else:
                s += '\nRequired dependencies are present, but some are unsupported.\n  Core workflow features may still work.\n  Please install updated versions if problems are encountered.\n'
            #end if
            if not all_req_supported:
                s += '\nRequired dependencies in need of user check or update:\n'
                for name in req_unknown:
                    s += '  {} version is unknown.  Check for {} or greater.\n'.format(name,version_to_string(self.currently_supported[name]))
                #end for
                for name in req_unsupported:
                    s += '  {} version {} is outdated.  Update to {} or greater.\n'.format(name,version_to_string(self.dependency_version[name]),version_to_string(self.currently_supported[name]))
                #end for
            #end if
            if not all_opt_supported:
                s += '\nSome optional dependencies are missing or merit an update.\n  These modules are not needed for core workflow operation.\n  Optional features related to outdated modules may still work.\n  Please install updated versions if problems are encountered.\n'
                if len(opt_missing)>0:
                    s += '\nOptional dependencies that are missing:\n'
                    for name in opt_missing:
                        s += '  {:<10} is missing.  Install {} or greater.\n'.format(name,version_to_string(self.currently_supported[name]))
                    #end for
                #end if
                if len(opt_unknown)>0 or len(opt_unsupported)>0:
                    s += '\nOptional dependencies benefitting from user check or update:\n'
                    for name in opt_unknown:
                        s += '  {:<10} version is unknown.  Check for {} or greater.\n'.format(name,version_to_string(self.currently_supported[name]))
                    #end for
                    for name in opt_unsupported:
                        s += '  {:<10} version {} is outdated.  Update to {} or greater.\n'.format(name,version_to_string(self.dependency_version[name]),version_to_string(self.currently_supported[name]))
                    #end for
                #end if
            #end if
        #end if
        s = n*pad+header+s.replace('\n','\n'+(n+1)*pad)+'\n'
        if write:
            print(s)
        #end if
        error = len(serr)>0
        if error and exit:
            self.error(serr)
        #end if
        if not full:
            return error
        else:
            return error,s,serr
        #end if
    #end def check
#end class Versions



# store availability of various required/optional dependencies
numpy_available      = False
scipy_available      = False
h5py_available       = False
matplotlib_available = False
pydot_available      = False
spglib_available     = False
pycifrw_available    = False
seekpath_available   = False
cif2cell_available   = False

numpy_supported      = False
scipy_supported      = False
h5py_supported       = False
matplotlib_supported = False
pydot_supported      = False
spglib_supported     = False
pycifrw_supported    = False
seekpath_supported   = False
cif2cell_supported   = False

try: # versioning info is never worth failure
    versions = Versions()

    numpy_available      = versions.available('numpy')
    scipy_available      = versions.available('scipy')
    h5py_available       = versions.available('h5py')
    matplotlib_available = versions.available('matplotlib')
    pydot_available      = versions.available('pydot')
    spglib_available     = versions.available('spglib')
    pycifrw_available    = versions.available('pycifrw')
    seekpath_available   = versions.available('seekpath')
    cif2cell_available   = versions.available('cif2cell')

    numpy_supported      = versions.supported('numpy')
    scipy_supported      = versions.supported('scipy')
    h5py_supported       = versions.supported('h5py')
    matplotlib_supported = versions.supported('matplotlib')
    pydot_supported      = versions.supported('pydot')
    spglib_supported     = versions.supported('spglib')
    pycifrw_supported    = versions.supported('pycifrw')
    seekpath_supported   = versions.supported('seekpath')
    cif2cell_supported   = versions.supported('cif2cell')
except:
    versions = None
#end try


def check_versions(write=True,exit=False):
    """
    (`User API`)  Print detailed information about the status of Nexus 
    dependencies on the current machine.
    """
    if versions is None:
        print('\nProblem encountered in Nexus version checking code.\nThis is a developer error.\nPlease run the Nexus tests via the nxs-test script for more details and report this issue to the developers.')
        error = True
    else:
        error = versions.check(write=write,exit=exit)
    #end if
    return error
#end def check_versions


def current_versions():
    """
    (`User API`)  Print information about currently supported versions of 
    Nexus dependencies.
    """
    if versions is None:
        print('\nProblem encountered in Nexus version checking code.\nThis is a developer error.\nPlease run the Nexus tests via the nxs-test script for more details and report this issue to the developers.')
    else:
        versions.print_current_versions()
    #end if
#end def current_versions


def policy_versions():
    """
    (`User API`)  Print information about versions of Nexus dependencies 
    that meet the age policy for today's date.
    """
    if versions is None:
        print('\nProblem encountered in Nexus version checking code.\nThis is a developer error.\nPlease run the Nexus tests via the nxs-test script for more details and report this issue to the developers.')
    else:
        versions.print_policy_versions()
    #end if
#end def policy_versions
