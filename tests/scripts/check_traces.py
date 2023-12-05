#! /usr/bin/env python3

# Performs consistency checks between traces.h5 and scalar.dat/stat.h5/dmc.dat files
# Jaron Krogel/ORNL

# Type the following to view documentation for command line inputs:
#   >check_traces.py -h

# For usage examples, type:
#   >check_traces.py -x

# check_traces.py packages obj and HDFreader classes from Nexus.
#   Note that h5py is required (which depends on numpy).
#   This script is compatible with both Python 2 and 3.

import os
import sys
from copy import deepcopy
import numpy as np
import h5py
from optparse import OptionParser

# Returns failure error code to OS.
# Explicitly prints 'fail' after an optional message.
def test_fail():
    print('\n\nTest status: fail')
    sys.exit(1)
#end def test_fail


# Returns success error code to OS.
# Explicitly prints 'pass' after an optional message.
def test_pass():
    print('\n\nTest status: pass')
    sys.exit(0)
#end def test_pass



######################################################################
# from generic.py
######################################################################

logfile = sys.stdout

def log(msg,n=0,indent='  '):
    if not isinstance(msg,str):
        msg = str(msg)
    #end if
    if n>0:
        indent = n*indent
        msg=indent+msg.replace('\n','\n'+indent)
    #end if
    logfile.write(msg+'\n')
#end def log


def error(msg,header=None,n=0):
    post_header=' error:'
    if header is None:
        header = post_header.lstrip()
    else:
        header += post_header
    #end if
    log('\n  '+header,n=n)
    log(msg.rstrip(),n=n)
    log('  exiting.\n')
    test_fail()
#end def error



class object_interface(object):
    _logfile = sys.stdout

    def __len__(self):
        return len(self.__dict__)
    #end def __len__

    def __contains__(self,name):
        return name in self.__dict__
    #end def

    def __getitem__(self,name):
        return self.__dict__[name]
    #end def __getitem__

    def __setitem__(self,name,value):
        self.__dict__[name]=value
    #end def __setitem__

    def __delitem__(self,name):
        del self.__dict__[name]
    #end def __delitem__

    def __iter__(self):
        for item in self.__dict__:
            yield self.__dict__[item]
        #end for
    #end def __iter__

    def __repr__(self):
        s=''
        for k in sorted(self.keys()):
            if not isinstance(k,str) or k[0]!='_':
                v=self.__dict__[k]
                if hasattr(v,'__class__'):
                    s+='  {0:<20}  {1:<20}\n'.format(str(k),v.__class__.__name__)
                else:
                    s+='  {0:<20}  {1:<20}\n'.format(str(k),type(v))
                #end if
            #end if
        #end for
        return s
    #end def __repr__

    def __str__(self,nindent=1):
        pad = '  '
        npad = nindent*pad
        s=''
        normal = []
        qable  = []
        for k,v in self.items():
            if not isinstance(k,str) or k[0]!='_':
                if isinstance(v,object_interface):
                    qable.append(k)
                else:
                    normal.append(k)
                #end if
            #end if
        #end for
        normal = sorted(normal)
        qable  = sorted(qable)
        indent = npad+18*' '
        for k in normal:
            v = self[k]
            vstr = str(v).replace('\n','\n'+indent)
            s+=npad+'{0:<15} = '.format(str(k))+vstr+'\n'
        #end for
        for k in qable:
            v = self[k]
            s+=npad+str(k)+'\n'
            s+=v.__str__(nindent+1)
            if isinstance(k,str):
                s+=npad+'end '+k+'\n'
            #end if
        #end for
        return s
    #end def __str__

    def keys(self):
        return self.__dict__.keys()
    #end def keys

    def values(self):
        return self.__dict__.values()
    #end def values

    def items(self):
        return self.__dict__.items()
    #end def items

    def copy(self):
        return deepcopy(self)
    #end def copy

    def clear(self):
        self.__dict__.clear()
    #end def clear

    def log(self,*items,**kwargs):
        log(*items,**kwargs)
    #end def log

    def error(self,message,header=None,n=0):
        if header is None:
            header = self.__class__.__name__
        #end if
        error(message,header,n=n)
    #end def error
#end class object_interface



class obj(object_interface):
    def __init__(self,*vars,**kwargs):
        for var in vars:
            if isinstance(var,(dict,object_interface)):
                for k,v in var.items():
                    self[k] = v
                #end for
            else:
                self[var] = None
            #end if
        #end for
        for k,v in kwargs.items():
            self[k] = v
        #end for
    #end def __init__

    def append(self,value):
        self[len(self)] = value
    #end def append

    def first(self):
        return self[min(self.keys())]
    #end def first
#end class obj

######################################################################
# end from generic.py
######################################################################



######################################################################
# from developer.py
######################################################################

class DevBase(obj):
    None
#end class DevBase

######################################################################
# end from developer.py
######################################################################




######################################################################
# from hdfreader.py
######################################################################
import keyword
from inspect import getmembers

class HDFglobals(DevBase):
    view = False
#end class HDFglobals


class HDFgroup(DevBase):
    def _escape_name(self,name):
        if name in self._escape_names:
            name=name+'_'
        #end if
        return name
    #end def escape_name

    def _set_parent(self,parent):
        self._parent=parent
        return
    #end def set_parent

    def _add_dataset(self,name,dataset):
        self._datasets[name]=dataset
        return 
    #end def add_dataset

    def _add_group(self,name,group):
        group._name=name
        self._groups[name]=group
        return 
    #end def add_group

    def _contains_group(self,name):
        return name in self._groups.keys()
    #end def _contains_group

    def _contains_dataset(self,name):
        return name in self._datasets.keys()
    #end def _contains_dataset


    def __init__(self):
        self._name=''
        self._parent=None
        self._groups={};
        self._datasets={};
        self._group_counts={}

        self._escape_names=None
        self._escape_names=set(dict(getmembers(self)).keys()) | set(keyword.kwlist)
        return
    #end def __init__


    def _remove_hidden(self,deep=True):
        if '_parent' in self:
            del self._parent
        #end if
        if deep:
            for name,value in self.items():
                if isinstance(value,HDFgroup):
                    value._remove_hidden()
                #end if
            #end for
        #end if
        for name in list(self.keys()):
            if name[0]=='_':
                del self[name]
            #end if
        #end for
    #end def _remove_hidden


    # read in all data views (h5py datasets) into arrays
    #   useful for converting a single group read in view form to full arrays
    def read_arrays(self):
        self._remove_hidden()
        for k,v in self.items():
            if isinstance(v,HDFgroup):
                v.read_arrays()
            else:
                self[k] = np.array(v)
            #end if
        #end for
    #end def read_arrays


    def get_keys(self):
        if '_groups' in self:
            keys = list(self._groups.keys())
        else:
            keys = list(self.keys())
        #end if
        return keys
    #end def get_keys
#end class HDFgroup




class HDFreader(DevBase):
    
    def __init__(self,fpath,verbose=False,view=False):
        
        HDFglobals.view = view

        if verbose:
            print('  Initializing HDFreader')
        #end if

        self.fpath=fpath
        if verbose:
            print('    loading h5 file')
        #end if

        try:
            self.hdf = h5py.File(fpath,'r')
        except IOError:
            self._success = False
            self.hdf = obj(obj=obj())
        else:
            self._success = True
        #end if

        if verbose:
            print('    converting h5 file to dynamic object')
        #end if

        #convert the hdf 'dict' into a dynamic object
        self.nlevels=1
        self.ilevel=0
        #  Set the current hdf group
        self.obj = HDFgroup()
        self.cur=[self.obj]
        self.hcur=[self.hdf]

        if self._success:
            cur   = self.cur[self.ilevel]
            hcur  = self.hcur[self.ilevel]
            for kr,v in hcur.items():
                k=cur._escape_name(kr)
                if isinstance(v, h5py.Dataset):
                    self.add_dataset(cur,k,v)
                elif isinstance(v, h5py.Group):
                    self.add_group(hcur,cur,k,v)
                else:
                    self.error('encountered invalid type: '+str(type(v)))
                #end if
            #end for
        #end if

        if verbose:
            print('  end HDFreader Initialization')
        #end if

        return
    #end def __init__


    def increment_level(self):
        self.ilevel+=1
        self.nlevels = max(self.ilevel+1,self.nlevels)
        if self.ilevel+1==self.nlevels:
            self.cur.append(None)
            self.hcur.append(None)
        #end if
        self.pad = self.ilevel*'  '
        return
    #end def increment_level

    def decrement_level(self):
        self.ilevel-=1
        self.pad = self.ilevel*'  '
        return
    #end def decrement_level

    def add_dataset(self,cur,k,v):
        if not HDFglobals.view:
            cur[k] = np.array(v)
        else:
            cur[k] = v
        #end if
        cur._add_dataset(k,cur[k])
        return
    #end def add_dataset

    def add_group(self,hcur,cur,k,v):
        cur[k] = HDFgroup()
        cur._add_group(k,cur[k])
        cur._groups[k]._parent = cur
        self.increment_level()
        self.cur[self.ilevel]  = cur._groups[k]
        self.hcur[self.ilevel] = hcur[k]

        cur   = self.cur[self.ilevel]
        hcur  = self.hcur[self.ilevel]
        for kr,v in hcur.items():
            k=cur._escape_name(kr)
            if isinstance(v, h5py.Dataset):
                self.add_dataset(cur,k,v)
            elif isinstance(v, h5py.Group):
                self.add_group(hcur,cur,k,v)
            #end if
        #end for
    #end def add_group
#end class HDFreader

######################################################################
# end from hdfreader.py
######################################################################



# Represents QMCPACK data file.
#   Used to read scalar.dat, stat.h5, dmc.dat, traces.h5
class DataFile(DevBase):

    aliases = None

    def __init__(self,filepath=None,quantities=None):
        self.data       = None
        self.filepath   = filepath
        self.quantities = None
        if quantities is not None:
            self.quantities = list(quantities)
        #end if
        if filepath is not None:
            self.read(filepath)
            if self.aliases is not None:
                for name,alias in self.aliases.items():
                    if name in self.data:
                        self.data[alias] = self.data[name]
                        del self.data[name]
                    #end if
                #end for
            #end if
            if quantities is not None:
                missing = set(quantities)-set(self.data.keys())
                if len(missing)>0:
                    self.error('some quantities are missing from file "{}"\nquantities present: {}\nquantities missing: {}'.format(self.filepath,sorted(self.data.keys()),sorted(missing)))
                #end if
            #end if
        #end if
    #end def __init__

    def read(self,filepath):
        None
    #end def read
#end class DataFile


# Used to parse scalar.dat and dmc.dat files
class DatFile(DataFile):
    def read(self,filepath):
        lt = np.loadtxt(filepath)
        if len(lt.shape)==1:
            lt.shape = (1,len(lt))
        #end if

        data = lt[:,1:].transpose()

        fobj = open(filepath,'r')
        variables = fobj.readline().split()[2:]
        fobj.close()

        self.data = obj()
        for i,vname in enumerate(variables):
            self.data[vname]=data[i,:]
        #end for
    #end def read
#end class DatFile


# Parses scalar.dat
class ScalarDatFile(DatFile):
    aliases = obj(BlockWeight='Weight')
#end class ScalarDat


# Parses dmc.dat
class DmcDatFile(DatFile):
    None
#end class DmcDatFile



# Parses scalar data from stat.h5
class ScalarHDFFile(DataFile):
    def read(self,filepath):
        hr = HDFreader(filepath)
        if not hr._success:
            self.error('hdf file read failed\nfile path: '+filepath)
        #end if
        data = hr.obj
        data._remove_hidden()

        self.data = obj()
        for name,value in data.items():
            if 'value' in value:
                self.data[name] = value.value.flatten()
            #end if
        #end for
    #end def read
#end class ScalarHDFFile



# Parses and organizes data from traces.h5.
#   QMCPACK writes one traces.h5 for each MPI task.
#   At every MC step, data from each walker is written to this file.
class TracesFileHDF(DataFile):
    def __init__(self,filepath=None,blocks=None):
        self.info = obj(
            blocks              = blocks,
            particle_sums_valid = None,
            )
        DataFile.__init__(self,filepath)
    #end def __init__


    def read(self,filepath=None):
        # Open the traces.h5 file
        hr = HDFreader(filepath)
        if not hr._success:
            self.error('hdf file read failed\nfile path: '+filepath)
        #end if
        hdf = hr.obj
        hdf._remove_hidden()

        # Translate from flat table structure to nested data structure.
        #   Do this for top level "int_data" and "real_data" HDF groups
        for name,buffer in hdf.items():
            self.init_trace(name,buffer)
        #end for

        # Sum trace data over walkers into per-step and per-block totals
        self.accumulate_scalars()
    #end def read


    # Translate from serialized traces table format to fully dimensioned 
    #   data resolved by physical quantity.
    def init_trace(self,name,fbuffer):
        trace = obj()
        if 'traces' in fbuffer:
            # Wide data table
            #   Each row corresponds to a single step of a single walker.
            #   Each row contains serialized scalar (e.g. LocalEnergy) 
            #   and array (e.g. electron coordinates) data. 
            ftrace = fbuffer.traces

            # Number of rows (walkers*steps for single mpi task)
            nrows = len(ftrace)

            # Serialization layout of each row is stored in "layout".
            #   Layout is separated into a few potential domains:
            #     scalar domain  : scalar quantities such as LocalEnergy
            #                      domain name is "scalars"
            #     electron domain: array quantities dimensioned like number of electrons (e.g. electron positions)
            #                      domain name follows particleset (often "e")
            #     ion domain     : array quantities dimensioned like number of ions
            #                      domain name follows particleset (often "ion0")
            for dname,fdomain in fbuffer.layout.items():
                domain = obj()
                # Get start and end row indices for each quantity
                for qname,fquantity in fdomain.items():
                    q = obj()
                    for vname,value in fquantity.items():
                        q[vname] = value
                    #end for

                    # extract per quantity data across all walkers and steps
                    quantity = ftrace[:,q.row_start:q.row_end]

                    # reshape from serialized to multidimensional data for the quantity
                    if q.unit_size==1:
                        shape = [nrows]+list(fquantity.shape[0:q.dimension])
                    else:
                        shape = [nrows]+list(fquantity.shape[0:q.dimension])+[q.unit_size]
                    #end if
                    quantity.shape = tuple(shape)
                    #if len(fquantity.shape)==q.dimension:
                    #    quantity.shape = tuple([nrows]+list(fquantity.shape))
                    ##end if
                    domain[qname] = quantity
                #end for
                trace[dname] = domain
            #end for
        else:
            self.error('traces are missing in file "{}"'.format(self.filepath))
        #end if
        # rename "int_data" and "real_data" as "int_traces" and "real_traces"
        self[name.replace('data','traces')] = trace
    #end def init_trace


    # Perform internal consistency check between per-walker single 
    #   particle energies and per-walker total energies.
    def check_particle_sums(self,tol):
        t = self.real_traces

        # Determine quantities present as "scalars" (total values) and also per-particle
        scalar_names = set(t.scalars.keys())
        other_names = []
        for dname,domain in t.items():
            if dname!='scalars':
                other_names.extend(domain.keys())
            #end if
        #end for
        other_names = set(other_names)
        sum_names = scalar_names & other_names

        # For each quantity, determine if the sum over particles matches the total
        same = True
        for qname in sum_names:
            # Get total value for each quantity
            q = t.scalars[qname]

            # Perform the sum over particles
            qs = 0*q
            for dname,domain in t.items():
                if dname!='scalars' and qname in domain:
                    tqs = domain[qname].sum(1)
                    if len(tqs.shape)==1:
                        qs[:,0] += tqs
                    else:
                        qs[:,0] += tqs[:,0]
                    #end if
                #end if
            #end for

            # Compare total and summed quantities
            qsame = (abs(q-qs)<tol).all()
            if qsame:
                log('{:<16} matches'.format(qname),n=3)
            else:
                log('{:<16} does not match'.format(qname),n=3)
            #end if
            same = same and qsame
        #end for
        self.info.particle_sums_valid = same
        return self.info.particle_sums_valid
    #end def check_particle_sums


    # Sum trace data over walkers into per-step and per-block totals
    def accumulate_scalars(self):
        # Get block and step information for the qmc method
        blocks = self.info.blocks
        if blocks is None:
            self.scalars_by_step  = None
            self.scalars_by_block = None
            return
        #end if

        # Get real and int valued trace data
        tr = self.real_traces
        ti = self.int_traces

        # Names shared by traces and scalar files
        scalar_names = set(tr.scalars.keys())

        # Walker step and weight traces
        st = ti.scalars.step
        wt = tr.scalars.weight
        if len(st)!=len(wt):
            self.error('weight and steps traces have different lengths')
        #end if

        # Compute number of steps and steps per block
        steps = st.max()+1
        steps_per_block = steps//blocks

        # Accumulate weights into steps and blocks
        ws   = np.zeros((steps,))
        wb   = np.zeros((blocks,))
        #ws2  = np.zeros((steps,))
        for t in range(len(wt)): # accumulate over walker population per step
            ws[st[t]] += wt[t]
            #ws2[st[t]] += 1.0 # scalar.dat/stat.h5 set weights to 1.0 in dmc
        #end for
        s = 0
        for b in range(blocks): # accumulate over steps in a block
            wb[b] = ws[s:s+steps_per_block].sum()
            #wb[b] = ws2[s:s+steps_per_block].sum()
            s+=steps_per_block
        #end for

        # Accumulate walker population into steps
        ps  = np.zeros((steps,))
        for t in range(len(wt)):
            ps[st[t]] += 1
        #end for

        # Accumulate quantities into steps and blocks
        #   These are the values directly comparable with data in 
        #   scalar.dat, stat.h5, and dmc.dat.
        scalars_by_step  = obj(Weight=ws,NumOfWalkers=ps)
        scalars_by_block = obj(Weight=wb)
        qs   = np.zeros((steps,))
        qb   = np.zeros((blocks,))
        qs2  = np.zeros((steps,))
        quantities = set(tr.scalars.keys())
        quantities.remove('weight')
        for qname in quantities:
            qt = tr.scalars[qname]
            if len(qt)!=len(wt):
                self.error('quantity {0} trace is not commensurate with weight and steps traces'.format(qname))
            #end if
            qs[:] = 0
            #qs2[:] = 0
            for t in range(len(wt)):
                qs[st[t]] += wt[t]*qt[t]
                #qs2[st[t]] += 1.0*qt[t]
            #end for
            qb[:] = 0
            s=0
            for b in range(blocks):
                qb[b] = qs[s:s+steps_per_block].sum()
                #qb[b] = qs2[s:s+steps_per_block].sum()
                s+=steps_per_block
            #end for
            qb = qb/wb
            qs = qs/ws
            scalars_by_step[qname]  = qs.copy()
            scalars_by_block[qname] = qb.copy()
        #end for
        self.scalars_by_step  = scalars_by_step
        self.scalars_by_block = scalars_by_block
    #end def accumulate_scalars
#end class TracesFileHDF



# Aggregates data from the full collection of traces.h5 files for a 
#   single series (e.g. VMC == series 0) and compares aggregated trace
#   values to data in scalar.dat, stat.h5, and dmc.dat.
class TracesAnalyzer(DevBase):

    # Read data from scalar.dat, stat.h5, dmc.dat and all traces.h5 for the series
    def __init__(self,options):

        self.particle_sums_valid = None
        self.scalar_dat_valid    = None
        self.stat_h5_valid       = None
        self.dmc_dat_valid       = None

        self.failed = False

        self.options = options.copy()

        prefix = options.prefix
        series = options.series
        method = options.method
        mpi    = options.mpi
        pseudo = options.pseudo
        path   = options.path

        # Determine the quantities to check
        dmc_dat_quants = ['Weight','LocalEnergy','NumOfWalkers']
        scalar_quants = ['LocalEnergy','Kinetic','LocalPotential',
                         'ElecElec','IonIon']
        if not pseudo:
            scalar_quants.append('ElecIon')
        else:
            scalar_quants.extend(['LocalECP','NonLocalECP'])
        #end if
        scalar_dat_quants = scalar_quants+['Weight']
        stat_h5_quants  = scalar_quants
        if self.options.quantities is not None:
            for qlist in (scalar_dat_quants,stat_h5_quants,dmc_dat_quants):
                old = list(qlist)
                del qlist[:]
                for q in old:
                    if q in self.options.quantities or q=='Weight':
                        qlist.append(q)
                    #end if
                #end for
            #end for
        #end if

        # Make paths to scalar, stat, dmc and traces files
        prefix = prefix+'.s'+str(series).zfill(3)

        scalar_file = os.path.join(path,prefix+'.scalar.dat')
        stat_file   = os.path.join(path,prefix+'.stat.h5')
        dmc_file    = os.path.join(path,prefix+'.dmc.dat')

        trace_files = []
        if mpi==1:
            tf = os.path.join(path,prefix+'.traces.h5')
            trace_files.append(tf)
        else:
            for n in range(mpi):
                tf = os.path.join(path,prefix+'.p'+str(n).zfill(3)+'.traces.h5')
                trace_files.append(tf)
            #end for
        #end if

        # Check that all output files exist
        files = [scalar_file,stat_file]
        if method=='dmc':
            files.append(dmc_file)
        #end if
        files.extend(trace_files)
        for filepath in files:
            if not os.path.exists(filepath):
                self.error('filepath {} does not exist'.format(filepath))
            #end if
        #end for

        # Load scalar, stat, dmc, and traces files

        # Load scalar.dat
        self.scalar_dat = ScalarDatFile(scalar_file,scalar_dat_quants)

        # Load stat.h5
        self.stat_h5 = ScalarHDFFile(stat_file,stat_h5_quants)

        # Load dmc.dat
        self.dmc_dat = None
        if method=='dmc':
            self.dmc_dat = DmcDatFile(dmc_file,dmc_dat_quants)
        #end if

        # Load all traces.h5
        self.data = obj()
        blocks = len(self.scalar_dat.data.first())
        for filepath in sorted(trace_files):
            trace_file = TracesFileHDF(filepath,blocks)
            self.data.append(trace_file)
        #end for
        assert(len(self.data)==mpi)

    #end def __init__


    # Check that per-particle quantities sum to total/scalar quantities 
    #   in each traces.h5 file separately.
    def check_particle_sums(self,tol):
        same = True
        for trace_file in self.data:
            log('Checking traces file: {}'.format(os.path.basename(trace_file.filepath)),n=2)
            same &= trace_file.check_particle_sums(tol=tol)
        #end for
        self.particle_sums_valid = same
        self.failed |= not same
        self.pass_fail(same,n=2)
        return same
    #end def check_particle_sums


    # Check aggregated traces data against scalar.dat
    def check_scalar_dat(self,tol):
        valid = self.check_scalar_file('scalar.dat',self.scalar_dat,tol)
        self.scalar_dat_valid = valid
        self.pass_fail(valid,n=2)
        return valid
    #end def check_scalar_dat


    # Check aggregated traces data against stat.h5
    def check_stat_h5(self,tol):
        valid = self.check_scalar_file('stat.h5',self.stat_h5,tol)
        self.stat_h5_valid = valid
        self.pass_fail(valid,n=2)
        return valid
    #end def check_stat_h5


    # Shared checking implementation for scalar.dat and stat.h5
    def check_scalar_file(self,file_type,scalar_file,tol):

        # Check that expected quantities are present
        qnames = scalar_file.quantities
        trace_names = set(self.data[0].scalars_by_block.keys())
        missing = set(qnames)-trace_names
        if len(missing)>0:
            self.error('{} file check failed for series {}\ntraces file is missing some quantities\nquantities present: {}\nquantities missing: {}'.format(file_type,self.options.series,sorted(trace_names),sorted(missing)))
        #end if

        scalars_valid  = True
        scalars        = scalar_file.data

        # Sum traces data for each quantity across all traces.h5 files
        summed_scalars = obj()
        for qname in qnames:
            summed_scalars[qname] = np.zeros(scalars[qname].shape)
        #end for
        wtot = np.zeros(summed_scalars.first().shape)
        for trace_file in self.data:
            # scalar.dat/stat.h5 are resolved per block
            w = trace_file.scalars_by_block.Weight
            wtot += w
            for qname in qnames:
                q = trace_file.scalars_by_block[qname]
                summed_scalars[qname] += w*q
            #end for
        #end for

        # Compare summed trace data against scalar.dat/stat.h5 values
        for qname in qnames:
            qscalar = scalars[qname]
            qb = summed_scalars[qname]/wtot
            match = abs(qb-qscalar)<tol
            all_match = match.all()
            self.log('{:<16} {}/{} blocks match'.format(qname,match.sum(),len(match)),n=2)
            if not all_match:
                for b,(m,qfile,qtrace) in enumerate(zip(match,qscalar,qb)):
                    if not m:
                        log('{:>3}  {: 16.12f}  {: 16.12f}  {: 16.12f}'.format(b,qfile,qtrace,qfile-qtrace),n=3)
                    #end if
                #end for
            #end if
            scalars_valid &= all_match
        #end for

        self.failed |= not scalars_valid

        return scalars_valid
    #end def check_scalar_file


    # Check aggregated traces data against dmc.dat
    def check_dmc_dat(self,tol):

        # Some DMC steps are excluded due to a known bug in QMCPACK weights
        dmc_steps_exclude = self.options.dmc_steps_exclude

        # Check that expected quantities are present
        dmc_file = self.dmc_dat
        qnames = dmc_file.quantities
        trace_names = set(self.data[0].scalars_by_step.keys())
        missing = set(qnames)-trace_names
        if len(missing)>0:
            self.error('dmc.dat check failed for series {}\ntraces file is missing some quantities\nquantities present: {}\nquantities missing: {}'.format(self.options.series,sorted(trace_names),sorted(missing)))
        #end if
        weighted = set(['LocalEnergy'])

        dmc_valid = True
        dmc       = dmc_file.data

        # Sum traces data for each quantity across all traces.h5 files
        summed_scalars = obj()
        for qname in qnames:
            summed_scalars[qname] = np.zeros(dmc[qname].shape)
        #end for
        wtot = np.zeros(summed_scalars.first().shape)
        for trace_file in self.data:
            # dmc.dat is resolved per step
            w = trace_file.scalars_by_step.Weight
            wtot += w
            for qname in qnames:
                q = trace_file.scalars_by_step[qname]
                if qname in weighted:
                    summed_scalars[qname] += w*q
                else:
                    summed_scalars[qname] += q
                #end if
            #end for
        #end for

        # Compare summed trace data against dmc.dat values
        for qname in qnames:
            qdmc = dmc[qname]
            if qname in weighted:
                qb = summed_scalars[qname]/wtot
            else:
                qb = summed_scalars[qname]
            #end if
            match = abs(qb-qdmc)<tol
            all_match = match.all()
            self.log('{:<16} {}/{} steps match'.format(qname,match.sum(),len(match)),n=2)
            if not all_match:
                for s,(m,qfile,qtrace) in enumerate(zip(match,qdmc,qb)):
                    if not m:
                        log('{:>3}  {: 16.12f}  {: 16.12f}  {: 16.12f}'.format(s,qfile,qtrace,qfile-qtrace),n=3)
                    #end if
                #end for
            #end if
            if dmc_steps_exclude>0:
                all_match = match[dmc_steps_exclude:].all()
            #end if
            dmc_valid &= all_match
        #end for

        if dmc_steps_exclude>0:
            log('\nExcluding first {} DMC steps from match check.'.format(dmc_steps_exclude),n=2)
        #end if

        self.dmc_dat_valid = dmc_valid
        self.pass_fail(dmc_valid,n=2)

        self.failed |= not dmc_valid

        return dmc_valid
    #end def check_dmc_dat


    # Print a brief message about pass/fail status of a subtest
    def pass_fail(self,passed,n):
        if passed:
            self.log('\nCheck passed!',n=n)
        else:
            self.log('\nCheck failed!',n=n)
        #end if
    #end def pass_fail
#end class TracesAnalyzer



examples = '''

===================================================================
Example 1: QMCPACK VMC/DMC with scalar-only traces, single mpi task
===================================================================

Contents of QMCPACK input file (qmc.in.xml):
--------------------------------------------
<simulation>
   <project id="qmc" series="0">
      ...
   </project>

   <qmcsystem>
      ...
   </qmcsystem>

   # write traces files w/o array info (scalars only)
   <traces array="no"/>  

   # vmc run, scalars written to stat.h5 
   <qmc method="vmc" move="pbyp">
      <estimator name="LocalEnergy" hdf5="yes"/>
      ...
   </qmc>

   # dmc run, scalars written to stat.h5 
   <qmc method="dmc" move="pbyp" checkpoint="-1">
     <estimator name="LocalEnergy" hdf5="yes"/>
     ...
   </qmc>

</simulation>

QMCPACK execution:
------------------
export OMP_NUM_THREADS=1
mpirun -np qmcpack qmc.in.xml

QMCPACK output files:
---------------------
qmc.s000.qmc.xml
qmc.s000.scalar.dat
qmc.s000.stat.h5
qmc.s000.traces.h5
qmc.s001.cont.xml
qmc.s001.dmc.dat
qmc.s001.qmc.xml
qmc.s001.scalar.dat
qmc.s001.stat.h5
qmc.s001.traces.h5

check_traces.py usage:
----------------------
check_traces.py -p qmc -s '0,1' -m 'vmc,dmc' --dmc_steps_exclude=1

Either execute check_traces.py in /your/path/to/qmcpack/output as above, or do:

check_traces.py -p qmc -s '0,1' -m 'vmc,dmc' --dmc_steps_exclude=1 /your/path/to/qmcpack/output


====================================================================
Example 2: QMCPACK VMC/DMC, selective scalar traces, single mpi task
====================================================================

Contents of QMCPACK input file (qmc.in.xml):
--------------------------------------------
<simulation>
   <project id="qmc" series="0">
      ...
   </project>

   <qmcsystem>
      ...
   </qmcsystem>

   # write traces files w/o array info (scalars only)
   <traces array="no">
      <scalar_traces>
         Kinetic ElecElec  # only write traces of Kinetic and ElecElec
      </scalar_traces>
   </traces>

   # vmc run, scalars written to stat.h5 
   <qmc method="vmc" move="pbyp">
      <estimator name="LocalEnergy" hdf5="yes"/>
      ...
   </qmc>

   # dmc run, scalars written to stat.h5 
   <qmc method="dmc" move="pbyp" checkpoint="-1">
     <estimator name="LocalEnergy" hdf5="yes"/>
     ...
   </qmc>

</simulation>

QMCPACK execution:
------------------
export OMP_NUM_THREADS=1
mpirun -np qmcpack qmc.in.xml

QMCPACK output files:
---------------------
qmc.s000.qmc.xml
qmc.s000.scalar.dat
qmc.s000.stat.h5
qmc.s000.traces.h5
qmc.s001.cont.xml
qmc.s001.dmc.dat
qmc.s001.qmc.xml
qmc.s001.scalar.dat
qmc.s001.stat.h5
qmc.s001.traces.h5

check_traces.py usage:
----------------------
check_traces.py -p qmc -s '0,1' -m 'vmc,dmc' -q 'Kinetic,ElecElec' --dmc_steps_exclude=1


===================================================================
Example 3: QMCPACK VMC/DMC, selective array traces, single mpi task
===================================================================

Contents of QMCPACK input file (qmc.in.xml):
--------------------------------------------
<simulation>
   <project id="qmc" series="0">
      ...
   </project>

   <qmcsystem>
      ...
   </qmcsystem>

   # write traces files w/ all scalar info and select array info
   <traces>
      <array_traces>
         position Kinetic  # write per-electron positions and kinetic energies
      </array_traces>
   </traces>

   # vmc run, scalars written to stat.h5 
   <qmc method="vmc" move="pbyp">
      <estimator name="LocalEnergy" hdf5="yes"/>
      ...
   </qmc>

   # dmc run, scalars written to stat.h5 
   <qmc method="dmc" move="pbyp" checkpoint="-1">
     <estimator name="LocalEnergy" hdf5="yes"/>
     ...
   </qmc>

</simulation>

QMCPACK execution:
------------------
export OMP_NUM_THREADS=1
mpirun -np qmcpack qmc.in.xml

QMCPACK output files:
---------------------
qmc.s000.qmc.xml
qmc.s000.scalar.dat
qmc.s000.stat.h5
qmc.s000.traces.h5
qmc.s001.cont.xml
qmc.s001.dmc.dat
qmc.s001.qmc.xml
qmc.s001.scalar.dat
qmc.s001.stat.h5
qmc.s001.traces.h5

check_traces.py usage:
----------------------
check_traces.py -p qmc -s '0,1' -m 'vmc,dmc' --psum --dmc_steps_exclude=1

'''



if __name__=='__main__':

    # Define command line options
    usage = '''usage: %prog [options] [path]'''
    parser = OptionParser(usage=usage,add_help_option=False,version='%prog 0.1')

    parser.add_option('-h','--help',dest='help',
                      action='store_true',default=False,
                      help='Print help information and exit (default=%default).'
                      )
    parser.add_option('-x','--examples',dest='examples',
                      action='store_true',default=False,
                      help='Print usage examples and exit (default=%default).'
                      )
    parser.add_option('-p','--prefix',dest='prefix',
                      default='qmc',
                      help='Series number(s) to check (default=%default).'
                      )
    parser.add_option('-s','--series',dest='series',
                      default='None',
                      help='Series number(s) to check (default=%default).'
                      )
    parser.add_option('-m','--methods',dest='methods',
                      default='None',
                      help='QMC method for each series.  Can be "vmc" or "dmc" for each series (default=%default).'
                      )
    parser.add_option('-q','--quantities',dest='quantities',
                      default='default',
                      help='QMC method for each series.  Can be "vmc" or "dmc" for each series (default=%default).'
                      )
    parser.add_option('-n','--mpi',dest='mpi',
                      default='1',
                      help='Number of MPI tasks in the original QMCPACK run.  This is also the number of traces.h5 files produced by a single VMC or DMC section (default=%default).'
                      )
    parser.add_option('--psum',dest='particle_sum',
                      action='store_true',default=False,
                      help='Check sums of single particle energies (default=%default).'
                      )
    parser.add_option('--pseudo',dest='pseudo',
                      action='store_true',default=True,
                      help='QMC calculation used pseudopotentials (default=%default).'
                      )
    parser.add_option('--dmc_steps_exclude',dest='dmc_steps_exclude',
                      default='0',
                      help='Exclude a number of DMC steps from being checked.  This option is temporary and will be removed once a bug in the DMC weights for the first step is fixed (default=%default).'
                      )
    parser.add_option('--tol',dest='tolerance',
                      default='1e-8',
                      help='Tolerance to check (default=%default).'
                      )

    opt,paths = parser.parse_args()

    options = obj(**opt.__dict__)

    # Process command line options
    if options.help:
        log('\n'+parser.format_help().strip()+'\n')
        sys.exit(0)
    #end if

    if options.examples:
        log(examples)
        sys.exit(0)
    #end if

    tol = float(options.tolerance)

    if len(paths)==0:
        options.path = './'
    elif len(paths)==1:
        options.path = paths[0]
    else:
        error('Only a single path is accepted as input.\nPaths provided:\n{}'.format(paths))
    #end if
    if not os.path.exists(options.path):
        error('Path to QMCPACK run does not exist.\nPath provided: {}'.format(options.path))
    elif os.path.isfile(options.path):
        error('Path to QMCPACK run is actually a file.\nOnly directory paths are accepted.\nPath provided: {}'.format(options.path))
    #end if

    def process_list(slist,type=lambda x: x):
        tokens = slist.strip('"').strip("'").replace(',',' ').split()
        lst = [type(t) for t in tokens]
        return lst
    #end def process_list

    options.series  = process_list(options.series,int)
    options.methods = process_list(options.methods)
    options.mpi     = int(options.mpi)
    options.dmc_steps_exclude = int(options.dmc_steps_exclude)
    if options.quantities=='default':
        options.quantities = None
    else:
        options.quantities = process_list(options.quantities)
    #end if

    if len(options.series)!=len(options.methods):
        error('"series" and "methods" must match in length.')
    #end if
    valid_methods = ['vmc','dmc']
    invalid = set(options.methods)-set(valid_methods)
    if len(invalid)>0:
        error('Invalid entries given for "methods".\nValid options are: {}\nYou provided: {}'.format(valid_methods,sorted(invalid)))
    #end if

    valid_quantities = '''
        LocalEnergy Kinetic LocalPotential ElecElec IonIon ElecIon LocalECP NonLocalECP
        '''.split()
    if options.quantities is not None:
        invalid = set(options.quantities)-set(valid_quantities)
        if len(invalid)>0:
            error('Invalid entries given for "quantities".\nValid options are: {}\nYou provided: {}'.format(valid_quantities,sorted(invalid)))
        #end if
    #end if


    # Parse all files across all requested series and compare traces vs scalar/stat/dmc

    log('\nChecking match between traces and scalar/stat/dmc files\n')

    log('\nOptions provided:\n'+str(options).rstrip())

    failed = False

    series_in  = options.series
    methods_in = options.methods

    del options.series
    del options.methods

    # Loop over vmc/dmc series
    for series,method in zip(series_in,methods_in):

        options.series = series
        options.method = method

        log('\n\nChecking series {} method={}'.format(series,method))

        # Read scalar.dat, stat.h5, dmc.dat, and *traces.h5 for the series
        ta = TracesAnalyzer(options)

        # Check traces data against scalar/stat/dmc files
        if method=='vmc':
            if options.particle_sum:
                log('\nChecking sums of single particle energies',n=1)
                ta.check_particle_sums(tol)
            #end if

            log('\nChecking scalar.dat',n=1)
            ta.check_scalar_dat(tol)

            log('\nChecking stat.h5',n=1)
            ta.check_stat_h5(tol)

        elif method=='dmc':
            if options.particle_sum:
                log('\nChecking sums of single particle energies',n=1)
                ta.check_particle_sums(tol)
            #end if

            log('\nSkipping checks of scalar.dat and stat.h5',n=1)
            log('Statistics for these files are currently computed\n'
                'after branching. Since traces are written before \n'
                'branching, these files cannot be reconstructed \n'
                'from the traces.',n=2)

            log('\nChecking dmc.dat',n=1)
            ta.check_dmc_dat(tol)
        #end if

        failed |= ta.failed
    #end for

    # Print final pass/fail message
    if failed:
        test_fail()
    else:
        test_pass()
    #end if
#end if
