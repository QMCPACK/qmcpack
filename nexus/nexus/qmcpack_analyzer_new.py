
import os
import subprocess
import copy
import numpy as np

from nexus.developer import DevBase,obj
from nexus.qmcpack_input import QmcpackInput


class QmcpackScalarInfo(DevBase):
    def __init__(self):
        self.aliases = obj(
            LocalEnergy    = 'E'  ,
            Kinetic        = 'T'  ,
            LocalPotential = 'V'  ,
            ElecElec       = 'Vee',
            ElecIon        = 'Vei',
            IonIon         = 'Vii',
            LocalECP       = 'Vl' ,
            NonLocalECP    = 'Vnl',
            Variance       = 'Ev' ,
            LocalEnergy_sq = 'E2' ,
            BlockWeight    = 'bw' ,
            BlockCPU       = 'bc' ,
            AcceptRatio    = 'ar' ,
            Efficiency     = 'eff',
            TotalTime      = 'tt' ,
            TotalSamples   = 'ts' ,
            DiffEff        = 'de' ,
            Weight         = 'w'  ,
            NumOfWalkers   = 'nw' ,
            LivingFraction = 'lf' ,
            AvgSentWalkers = 'asw',
            )

        self.inv_aliases = obj(
            E   = 'LocalEnergy'   ,
            T   = 'Kinetic'       ,
            V   = 'LocalPotential',
            Vee = 'ElecElec'      ,
            Vei = 'ElecIon'       ,
            Vii = 'IonIon'        ,
            Vl  = 'LocalECP'      ,
            Vnl = 'NonLocalECP'   ,
            Ev  = 'Variance'      ,
            E2  = 'LocalEnergy_sq',
            bw  = 'BlockWeight'   ,
            bc  = 'BlockCPU'      ,
            ar  = 'AcceptRatio'   ,
            eff = 'Efficiency'    ,
            tt  = 'TotalTime'     ,
            ts  = 'TotalSamples'  ,
            de  = 'DiffEff'       ,
            w   = 'Weight'        ,
            nw  = 'NumOfWalkers'  ,
            lf  = 'LivingFraction',
            asw = 'AvgSentWalkers',
            )

        self.nonenergy = {
            'BlockWeight','BlockCPU','AcceptRatio','Efficiency',
            'TotalTime','TotalSamples','DiffEff','Weight',
            'NumOfWalkers','LivingFraction','AvgSentWalkers'}

        self.integer = {'TotalSamples', 'NumOfWalkers'}

        self.constant = {'IonIon','KEcorr','MPC'}

        self.analyze = set(self.aliases.keys())
    #end def __init__
#end class QmcpackScalarInfo
scalar_info = QmcpackScalarInfo()

    


class ReadScalarIssues(DevBase):
    issues = (
        'no_file',          # file does not exist
        'empty_file',       # file is completely empty
        'bad_header',       # header is malformed
        'no_data',          # file contains no data
        'nan_vals',         # some data values are NaN
        'unparsable_vals',  # some data values couldn't be read
        'uneven_cols',      # not all rows had the same length
        'incomplete',       # file does not contain all expected data
        'no_usable_vals',   # file has no usable data
        )

    def __init__(self):
        for issue in self.issues:
            self[issue] = False
    #end def __init__

    def add(self,issue):
        assert issue in self.issues
        self[issue] = True
    #end def add

    def issue_set(self):
        return {issue for issue in self.keys() if self[issue]}
    #end def issue

    def failed(self):
        issues = self.issue_set()
        issues -= {'nan_vals','incomplete'}
        return len(issues)>0
    #end def failed
#end class ReadScalarIssues



def read_scalar_file(filepath,
                     issues       = False,
                     add_variance = False,
                     remove_index = False,
                     trim_nan     = True,
                     nrows        = None,
                     dict_type    = dict,
                     ):
    '''
    Robustly read scalar.dat or dmc.dat files

    If an unhandled exception is raised, this code needs fixing.
    '''
    import subprocess
    import warnings
    assert isinstance(filepath,str)
    assert isinstance(issues,bool)
    ret_issues = issues
    data   = dict_type()
    issues = ReadScalarIssues()
    if not ret_issues:
        ret = data
    else:
        ret = data,issues
    # check if file exists
    if not os.path.exists(filepath):
        issues.add('no_file')
        return ret
    # parse header
    f = open(filepath,'r')
    var_names        = None
    malformed_header = False
    for line in f:
        line = line.strip()
        if len(line)==0:
            continue
        if line.startswith('#'):
            tokens = line.split()[1:]
            if len(tokens)==0:
                malformed_header = True
                break
            else:
                var_names = tokens
        else:
            malformed_header = True
        break
    f.close()
    if malformed_header:
        issues.add('bad_header')
    elif var_names is None:
        issues.add('empty_file')
    if len(issues.issue_set())>0:
        return ret
    # attempt data parse via loadtxt
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        try:
            lt = np.loadtxt(filepath)
            parsed = True
        except:
            parsed = False
    if parsed:
        if len(lt.shape)==1:
            lt.shape = (1,len(lt))
        data_arr = lt.T
        if len(data_arr)==0:
            issues.add('no_data')
            return ret
        for n,var in enumerate(var_names):
            data[var] = data_arr[n,:].ravel()
    # robust fallback data parsing
    if not parsed:
        nvars = len(var_names)
        no_data         = True
        uneven_cols     = False
        unparsable_vals = False
        data_tokens = []
        f = open(filepath,'r')
        for line in f:
            line = line.strip()
            if len(line)==0 or line[0]=='#':
                continue
            no_data = False
            tokens = line.split()
            if len(tokens)!=nvars:
                uneven_cols = True
            vals = []
            for t in tokens:
                try:
                    vals.append(float(t))
                except:
                    unparsable_vals = True
            if uneven_cols or unparsable_vals:
                break
            data_tokens.extend(vals)
        f.close()
        if no_data:
            issues.add('no_data')
        if unparsable_vals:
            issues.add('unparsable_vals')
        if uneven_cols:
            issues.add('uneven_cols')
        if len(data_tokens)>0:
            assert len(data_tokens)%nvars==0
            data_arr = np.array(data_tokens,dtype=float)
            data_arr = np.reshape(data_arr,(len(data_arr)//nvars,nvars))
            data_arr = data_arr.T
            for n,var in enumerate(var_names):
                data[var] = data_arr[n,:].ravel()
            parsed = True
        else:
            return ret
    # adjust entries
    if add_variance and 'LocalEnergy_sq' in data and 'LocalEnergy' in data:
        data['Variance'] = data['LocalEnergy_sq'] - data['LocalEnergy']**2
        del data['LocalEnergy_sq']
    if remove_index and 'index' in data:
        del data['index']
    # trim nan's
    if parsed and trim_nan:
        trim_lens = []
        for k,d in data.items():
            nan_entries = ~np.isfinite(d)
            inan_first = nan_entries.argmax()
            if not np.isfinite(d[inan_first]):
                trim_lens.append(inan_first)
        if len(trim_lens)>0:
            issues.add('nan_vals')
            ntrim = max(trim_lens)
            for k,d in data.items():
                data[k] = d[:ntrim]
    # check for usable data
    nusable = max([len(v) for v in data.values()])
    if nusable==0:
        issues.add('no_usable_vals')
    # check if all expected rows are present
    if nrows is not None:
        wc_out = subprocess.check_output(['wc', '-l', filepath])
        line_count = int(wc_out.split()[0])
        if line_count<nrows+1:
            issues.add('incomplete')
    return ret
#end def read_scalar_file



def qmcpack_analyzer_outfiles(qmc,prefix,series,group_index=None):
    assert qmc in {'opt','vmc','dmc'}
    ss     = 's'+str(series).zfill(3)
    if group_index is None:
        prefix = f'{prefix}.{ss}.'
    else:
        gs = 'g'+str(group_index).zfill(3)
        prefix = f'{prefix}.{gs}.{ss}.'        
    if qmc=='vmc':
        postfixes = ['scalar.dat']
    elif qmc=='dmc':
        postfixes = ['scalar.dat','dmc.dat']
    elif qmc=='opt':
        postfixes = ['scalar.dat','opt.xml','vp.h5']
    outfiles = tuple(prefix+pf for pf in postfixes)
    return outfiles
#end def qmcpack_analyzer_outfiles



class QmcpackInputInfo(DevBase):
    '''
    Collects prefix and qmc run info from a QMCPACK input file

    Lists expected scalar/dmc data and opt param files for the run
    '''

    def __init__(self,filepath):
        self.filepath     = filepath
        self.qmc_type     = None
        self.prefix       = None
        self.group_index  = None
        self.series_start = 0
        self.has_twist    = False
        self.qmc_info     = None
        self.read(filepath)
    #end def __init__

    def read(self,filepath):
        if not os.path.exists(filepath):
            self.error(f'provided qmcpack input file does not exist.\nFilepath: {filepath}')
        # parse filepath
        filename = os.path.split(filepath)[1]
        ftokens = filename.split('.')
        group_index = None
        for t in ftokens:
            if t.startswith('g'):
                try:
                    gi = int(t[1:])
                except:
                    gi = None
                if gi is not None:
                    group_index = gi
                    break
        self.group_index=group_index
        # parse input
        qi = QmcpackInput(filepath)
        qi.pluralize()
        # prefix, series
        project = qi.get('project')
        self.prefix = project.id
        if 'series' in project:
            self.series_start = project.series
        # twist
        twistnum = qi.get('twistnum')
        twist    = qi.get('twist')
        if twistnum is not None or twist is not None:
            self.has_twist = True
        # qmc method info
        shr = dict(blocks=1,timestep=0.)
        defaults = dict(
            vmc = dict(qmc='vmc',warmupsteps=0  ,steps=1,**shr),
            dmc = dict(qmc='dmc',warmupsteps=200,steps=1,**shr),
            opt = dict(qmc='opt',warmupsteps=0  ,**shr))
        opt_methods = set(['linear'])
        series = self.series_start
        def qinfo_from_qmc(qmc,prefix,series):
            if qmc.method in opt_methods:
                qmc.method = 'opt'
            qinfo = obj(**defaults[qmc.method])
            for k,v in qmc.items():
                if k in qinfo:
                    qinfo[k] = v
            qinfo.series = series
            qinfo.outfiles = qmcpack_analyzer_outfiles(qinfo.qmc,prefix,series,self.group_index)
            return qinfo
        qmc_info = obj()
        for qmc in qi.simulation.calculations:
            if 'max' in qmc: #loop
                loop_max = qmc.max
                opt_calcs = qmc.unroll()
                assert len(opt_calcs)==loop_max
                for n in range(loop_max):
                    qmc_info[series] = qinfo_from_qmc(opt_calcs[n],self.prefix,series)
                    series += 1
            else:
                qmc_info[series] = qinfo_from_qmc(qmc,self.prefix,series)
                series += 1
        qmc_types = set([q.qmc for q in qmc_info.values()])
        if 'dmc' in qmc_types:
            self.qmc_type = 'dmc'
        elif 'opt' in qmc_types:
            self.qmc_type = 'opt'
        elif 'vmc' in qmc_types:
            self.qmc_type = 'vmc'
        else:
            raise RuntimeError('impossible branch')
        self.qmc_info = qmc_info
    #end def read
#end class QmcpackInputInfo
