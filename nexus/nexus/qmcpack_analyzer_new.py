
import os
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
    """Record distinct conditions encountered while reading scalar data.

    Each issue is represented by a boolean attribute initialized to ``False``.
    Issues are added with :meth:`add`, and the names of active issues can be
    obtained with :meth:`issue_set`.

    Attributes
    ----------
    no_file : bool
        The requested path does not exist.  No attempt is made to read a
        header or data.
    empty_file : bool
        The file contains no nonblank line from which a header can be read.
        A file containing a valid header but no data rows instead receives
        ``no_data``.
    bad_header : bool
        The first nonblank line is not a comment containing at least one
        column name.
    bad_col_count : bool
        After excluding a tolerated short final numeric row and corrupt
        trailing rows, all rows have the same number of columns, but that
        number differs from the number of header columns.  This differs from
        ``uneven_cols``, which indicates inconsistent row widths.
    no_data : bool
        A valid header is present, but there are no nonblank, noncomment data
        rows.
    nan_vals : bool
        At least one recognized scalar column contains a non-finite value.
        The index and ignored trailing columns do not participate in this
        check.  If ``trim_nan=True``, every affected row is removed from the
        returned data; the issue is still recorded.  Such rows count toward
        ``nrows`` if they are otherwise complete and numeric.
    unparsable_vals : bool
        A row containing a token that cannot be converted to a float occurs
        before a later complete numeric row.  This specifically identifies
        corruption within the data.  Unparsable rows after the final complete
        row produce ``corrupt_end`` instead.
    uneven_cols : bool
        Rows outside the tolerated corrupt end do not all have the same
        token count.  Both numeric and internally unparsable rows participate
        in this check.  A single final numeric row with too few columns is
        treated as a partial write and is excluded, while a row with too many
        columns is not excluded.
    incomplete : bool
        The number of complete numeric rows differs from the requested
        ``nrows``.  A complete row has exactly the header width and every
        token converts to a float; NaNs and infinities count as converted
        values.  Counting occurs before recognized-column NaN filtering.
    no_usable_vals : bool
        Data rows exist, but no usable contiguous data prefix can be returned,
        or filtering removes every row in that prefix.
    corrupt_end : bool
        One or more rows after the final complete numeric row contain tokens
        that cannot be converted to floats.  This is a recoverable condition:
        the valid prefix is returned and the corrupt rows are excluded from
        column consistency checks.  A short partial numeric row may occur
        before or after such garbage without changing the classification; it
        is not itself corrupt end data.

    Notes
    -----
    ``failed()`` treats ``nan_vals``, ``incomplete``, and ``corrupt_end`` as
    recoverable conditions.  Every other active issue is considered a parsing
    failure.  When corruption occurs in the middle of a file, only the
    contiguous valid prefix preceding it is returned, even though later rows
    are inspected for issue classification and ``nrows`` validation.
    """

    issues = (
        'no_file',          # file does not exist
        'empty_file',       # file is completely empty
        'bad_header',       # header is malformed
        'bad_col_count',    # header and data column counts do not match
        'no_data',          # file contains no data
        'nan_vals',         # some data values are NaN
        'unparsable_vals',  # some data values couldn't be read
        'uneven_cols',      # not all rows had the same length
        'incomplete',       # file does not contain all expected data
        'no_usable_vals',   # file has no usable data
        'corrupt_end',      # trailing rows contain nonnumeric garbage
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
        issues -= {'nan_vals','incomplete','corrupt_end'}
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
    # Parse every data row before classifying corruption and column counts.
    nvars = len(var_names)
    rows = []
    f = open(filepath,'r')
    for line in f:
        line = line.strip()
        if len(line)==0 or line[0]=='#':
            continue
        tokens = line.split()
        vals = []
        parsable = True
        for t in tokens:
            try:
                vals.append(float(t))
            except:
                parsable = False
        rows.append(obj(
            ncols    = len(tokens),
            parsable = parsable,
            complete = parsable and len(tokens)==nvars,
            values   = vals,
            ))
    f.close()
    if len(rows)==0:
        issues.add('no_data')
        if nrows is not None and nrows!=0:
            issues.add('incomplete')
        return ret

    # Unparsable rows after the final complete row are recoverable end
    # corruption, even if a short partial numeric row follows them.
    complete_indices = [n for n,row in enumerate(rows) if row.complete]
    if len(complete_indices)>0:
        last_complete = complete_indices[-1]
    else:
        last_complete = -1
    corrupt_end_rows = {
        n for n,row in enumerate(rows)
        if not row.parsable and n>last_complete
        }
    if len(corrupt_end_rows)>0:
        issues.add('corrupt_end')

    # A short final numeric row represents a write still in progress.  It is
    # not used for column consistency checks or as complete scalar data.
    parsable_indices = [n for n,row in enumerate(rows) if row.parsable]
    partial_final_row = None
    if len(parsable_indices)>0:
        n = parsable_indices[-1]
        if rows[n].ncols<nvars:
            partial_final_row = n
    column_rows = [
        row for n,row in enumerate(rows)
        if n not in corrupt_end_rows and n!=partial_final_row
        ]

    col_counts = {row.ncols for row in column_rows}
    if len(col_counts)>1:
        issues.add('uneven_cols')
    elif (len(col_counts)==1 and next(iter(col_counts))!=nvars and
          any(row.parsable for row in column_rows)):
        issues.add('bad_col_count')

    # Unparsable data is fatal only if complete numeric data resumes later.
    if any(not row.parsable and n<last_complete for n,row in enumerate(rows)):
        issues.add('unparsable_vals')

    # Only return the contiguous complete prefix.  Later complete rows are
    # still counted for nrows and used to detect middle corruption above.
    data_rows = []
    for row in rows:
        if not row.complete:
            break
        data_rows.append(row.values)

    complete_rows = sum(row.complete for row in rows)
    if nrows is not None and complete_rows!=nrows:
        issues.add('incomplete')

    if len(data_rows)>0:
        data_arr = np.asarray(data_rows,dtype=float).T
    else:
        data_arr = None

    # Scalar quantities must be contiguous and directly follow index.
    # Ignore all columns beginning with the first unrecognized quantity.
    scalar_names = scalar_info.analyze | scalar_info.constant
    selected_cols = []
    first_scalar = 0
    if var_names[0]=='index':
        selected_cols.append((0,'index'))
        first_scalar = 1
    for n in range(first_scalar,nvars):
        var = var_names[n]
        if var not in scalar_names:
            break
        selected_cols.append((n,var))
    if data_arr is not None:
        for n,var in selected_cols:
            data[var] = data_arr[n,:].ravel()
    # adjust entries
    if add_variance and 'LocalEnergy_sq' in data and 'LocalEnergy' in data:
        data['Variance'] = data['LocalEnergy_sq'] - data['LocalEnergy']**2
        del data['LocalEnergy_sq']
    if remove_index and 'index' in data:
        del data['index']
    # detect non-finite scalar values and optionally remove affected rows
    scalar_data = [d for k,d in data.items() if k in scalar_names]
    if len(scalar_data)>0:
        usable_rows = np.ones(len(scalar_data[0]),dtype=bool)
        for d in scalar_data:
            usable_rows &= np.isfinite(d)
        if not usable_rows.all():
            issues.add('nan_vals')
            if trim_nan:
                for k,d in data.items():
                    data[k] = d[usable_rows]
    # check for usable data
    if len(data)==0:
        nusable = 0
    else:
        nusable = max([len(v) for v in data.values()])
    if nusable==0:
        issues.add('no_usable_vals')
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
        if project is None:
            self.prefix = 'default_project'
        else:
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
        method_types = dict(
            opt          = 'opt',
            linear       = 'opt',
            cslinear     = 'opt',
            linear_batch = 'opt',
            vmc          = 'vmc',
            vmc_batch    = 'vmc',
            dmc          = 'dmc',
            dmc_batch    = 'dmc',
            )
        series = self.series_start
        def qinfo_from_qmc(qmc,prefix,series):
            method = qmc.method
            if method not in method_types:
                raise ValueError(f'unrecognized qmc method: {method}')
            qmc_type = method_types[method]
            qinfo = obj(**defaults[qmc_type])
            for k,v in qmc.items():
                if k in qinfo:
                    qinfo[k] = v
            qinfo.series = series
            qinfo.outfiles = qmcpack_analyzer_outfiles(qinfo.qmc,prefix,series,self.group_index)
            return qinfo
        qmc_info = obj()
        calculations = qi.get('calculations')
        if calculations is None:
            calculations = ()
        for qmc in calculations:
            if 'max' in qmc: #loop
                for calc in qmc.unroll():
                    qmc_info[series] = qinfo_from_qmc(calc,self.prefix,series)
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
        self.qmc_info = qmc_info
    #end def read
#end class QmcpackInputInfo
