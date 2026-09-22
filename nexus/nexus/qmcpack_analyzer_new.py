##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#python standard library imports
import os
from numbers import Integral

#custom library imports
from .developer import DevBase,obj
from .qmcpack_input import QIxml,QmcpackInput,collection,loop,project,simulation


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

        self.inv_aliases = obj()
        for k,v in self.aliases.items():
            self.inv_aliases[v] = k

        self.nonenergy = {
            'BlockWeight','BlockCPU','AcceptRatio','Efficiency',
            'TotalTime','TotalSamples','DiffEff','Weight',
            'NumOfWalkers','LivingFraction','AvgSentWalkers',
            }

        self.integer = {'TotalSamples','NumOfWalkers'}

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
        The first trusted data row has a uniform but incorrect width: its
        number of tokens differs from the number of header columns and no
        earlier data row established the header width.  Parsing stops at this
        row.  This differs from ``uneven_cols``, where an earlier row already
        had the expected width.  A tolerated short final numeric row and
        trailing corrupt text do not produce this issue.
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
        An exact-width row containing a token that cannot be converted to a
        float occurs before a later complete numeric row.  Such a row is
        skipped, does not count toward ``nrows``, and does not stop extraction.
        A wrong-width unparsable row also produces this issue if complete data
        occurs later, but that row stops extraction because its width is
        wrong.  Unparsable rows after the final trusted complete row produce
        ``corrupt_end`` instead.
    uneven_cols : bool
        A row has a different token count from the header after an earlier
        row established the expected width.  Parsing stops at the first such
        row, and that row and all following rows are excluded from extraction
        and ``nrows`` validation.  Both numeric and nonnumeric rows can expose
        this condition.  A short final numeric row is instead treated as a
        partial write and does not produce this issue; an overlong final
        numeric row does.
    incomplete : bool
        The number of complete numeric rows differs from the requested
        ``nrows``.  A complete row has exactly the header width and every
        token converts to a float; NaNs and infinities count as converted
        values.  Counting occurs before recognized-column NaN filtering.
    nrows_unchecked : bool
        No expected ``nrows`` value was supplied, so row-count completeness
        was not validated.  This is recoverable for ordinary parsing, but
        causes :meth:`complete` to return ``False`` because completeness
        cannot be established.
    no_usable_vals : bool
        No fully parsable, exact-width row before the first width error has
        usable recognized scalar data.  With ``trim_nan=True``, a row is usable
        only when all recognized scalar values are finite; with
        ``trim_nan=False``, parsed NaNs and infinities are allowed.  The index
        is not scalar data, except that it counts as usable when it is the only
        column in the file.
    corrupt_end : bool
        One or more unparsable rows occur after the final trusted complete
        numeric row.  These rows are trailing garbage rather than middle
        corruption.  A short final numeric row may occur before or after such
        garbage without becoming corrupt itself.  Trailing corruption is
        recoverable and does not by itself prevent complete trusted data from
        being returned.

    Notes
    -----
    ``failed()`` treats ``nan_vals``, ``incomplete``, ``nrows_unchecked``, and
    ``corrupt_end`` as recoverable conditions.  Every other active issue is
    considered a parsing failure.  When corruption occurs in the middle of a
    file, later exact-width complete rows are returned after the corrupt row.
    In contrast, the first wrong-width row is a trust boundary: that row and
    all later rows are excluded from extraction and ``nrows`` validation.
    Later rows may be inspected only to distinguish middle corruption from a
    corrupt end and a final partial numeric row.
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
        'nrows_unchecked',  # no expected row count was provided
        'no_usable_vals',   # file has no usable data
        'corrupt_end',      # trailing rows contain nonnumeric garbage
        )

    def __init__(self,*,nrows_checked=False):
        for issue in self.issues:
            self[issue] = False
        if not nrows_checked:
            self.add('nrows_unchecked')
    #end def __init__

    def add(self,issue):
        if issue not in self.issues:
            raise ValueError(f'unrecognized scalar-read issue: {issue}')
        self[issue] = True
    #end def add

    def issue_set(self):
        return {issue for issue in self.keys() if self[issue]}
    #end def issue_set

    def failed(self):
        issues = self.issue_set()
        issues -= {'nan_vals','incomplete','nrows_unchecked','corrupt_end'}
        return len(issues)>0
    #end def failed

    def complete(self,*,allow_nan=False):
        """Return whether the requested number of usable rows is present.

        Parameters
        ----------
        allow_nan : bool, optional
            If ``False``, recognized non-finite scalar values prevent the
            result from being complete.  If ``True``, rows containing such
            values count toward completeness.  The default is ``False``.

        Returns
        -------
        complete : bool
            ``True`` only when ``nrows`` was supplied, exactly that many
            complete numeric rows were present, and no fatal parsing issue was
            encountered.  Recoverable trailing corruption does not prevent
            completeness.  ``nan_vals`` prevents completeness unless
            ``allow_nan=True``.
        """
        if not isinstance(allow_nan,bool):
            raise TypeError('allow_nan must be a bool')
        incomplete_issues = {
            'no_file',
            'empty_file',
            'bad_header',
            'bad_col_count',
            'no_data',
            'unparsable_vals',
            'uneven_cols',
            'incomplete',
            'nrows_unchecked',
            'no_usable_vals',
            }
        if any(self[issue] for issue in incomplete_issues):
            return False
        if self.nan_vals and not allow_nan:
            return False
        return True
    #end def complete
#end class ReadScalarIssues



def read_scalar_file(
    filepath,
    *,
    issues       = False,
    add_variance = False,
    remove_index = False,
    trim_nan     = True,
    nrows        = None,
    dict_type    = dict,
    ):
    """Robustly read a QMCPACK scalar.dat or dmc.dat file.

    Parameters
    ----------
    filepath : str
        Path to the scalar data file.
    issues : bool, optional
        If ``True``, return a :class:`ReadScalarIssues` object with the data.
    add_variance : bool, optional
        If ``True``, construct ``Variance`` from ``LocalEnergy_sq`` and
        ``LocalEnergy`` when both quantities are present.
    remove_index : bool, optional
        If ``True``, omit the index column from the returned data.
    trim_nan : bool, optional
        If ``True``, remove rows containing non-finite recognized scalar data.
    nrows : int, optional
        Expected number of complete, parsable data rows.
    dict_type : type, optional
        Mapping type used to hold returned scalar arrays.

    Returns
    -------
    data : mapping
        Parsed recognized scalar columns.
    issues : ReadScalarIssues
        Issue information, returned only when ``issues=True``.
    """
    import numpy as np

    if not isinstance(filepath,str):
        raise TypeError('filepath must be a str')
    if not isinstance(issues,bool):
        raise TypeError('issues must be a bool')
    ret_issues = issues
    data   = dict_type()
    issues = ReadScalarIssues(nrows_checked=nrows is not None)
    if not ret_issues:
        ret = data
    else:
        ret = data,issues
    # check if file exists
    if not os.path.exists(filepath):
        issues.add('no_file')
        return ret
    # parse header
    var_names        = None
    malformed_header = False
    with open(filepath,'r') as fobj:
        for line in fobj:
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
    if malformed_header:
        issues.add('bad_header')
    elif var_names is None:
        issues.add('empty_file')
    if issues.no_file or issues.empty_file or issues.bad_header:
        return ret
    # Read every row so that a width boundary can be classified without using
    # any data beyond it for extraction or nrows validation.
    nvars = len(var_names)
    rows = []
    with open(filepath,'r') as fobj:
        for line in fobj:
            line = line.strip()
            if len(line)==0 or line[0]=='#':
                continue
            tokens = line.split()
            vals = []
            parsable = True
            for token in tokens:
                try:
                    vals.append(float(token))
                except ValueError:
                    parsable = False
            rows.append(obj(
                ncols    = len(tokens),
                parsable = parsable,
                complete = parsable and len(tokens)==nvars,
                values   = vals,
                ))
    if len(rows)==0:
        issues.add('no_data')
        if nrows is not None and nrows!=0:
            issues.add('incomplete')
        return ret

    wrong_width = next(
        (n for n,row in enumerate(rows) if row.ncols!=nvars), len(rows)
        )
    trusted_rows = rows[:wrong_width]

    # Determine whether the width boundary is a tolerated partial final write,
    # recoverable trailing garbage, or an actual column-count error.  Later
    # rows inform this classification but remain untrusted.
    if wrong_width<len(rows):
        boundary_row = rows[wrong_width]
        later_complete = any(row.complete for row in rows[wrong_width+1:])
        later_parsable = any(row.parsable for row in rows[wrong_width+1:])
        partial_final = (
            boundary_row.parsable and boundary_row.ncols<nvars and
            not later_parsable
            )
        corrupt_boundary = not boundary_row.parsable and not later_complete
        if corrupt_boundary:
            issues.add('corrupt_end')
        elif partial_final:
            if any(not row.parsable for row in rows[wrong_width+1:]):
                issues.add('corrupt_end')
        else:
            if wrong_width==0:
                issues.add('bad_col_count')
            else:
                issues.add('uneven_cols')
            if not boundary_row.parsable and later_complete:
                issues.add('unparsable_vals')

    # Exact-width unparsable rows do not stop parsing.  They are middle
    # corruption when complete data resumes, and corrupt-end data otherwise.
    trusted_complete = [n for n,row in enumerate(trusted_rows) if row.complete]
    last_trusted_complete = trusted_complete[-1] if trusted_complete else -1
    if any(not row.parsable and n<last_trusted_complete
           for n,row in enumerate(trusted_rows)):
        issues.add('unparsable_vals')
    if any(not row.parsable and n>last_trusted_complete
           for n,row in enumerate(trusted_rows)):
        issues.add('corrupt_end')

    # Skip exact-width nonnumeric rows, but retain all complete rows up to the
    # first width error.
    data_rows = [row.values for row in trusted_rows if row.complete]
    complete_rows = len(data_rows)
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
    # Detect non-finite recognized scalar values and optionally remove their
    # rows.  Index and ignored columns do not participate in this check.
    scalar_data = [d for k,d in data.items() if k in scalar_names]
    ndata_rows = len(data_rows)
    if len(scalar_data)>0:
        finite_rows = np.ones(ndata_rows,dtype=bool)
        for d in scalar_data:
            finite_rows &= np.isfinite(d)
        if not finite_rows.all():
            issues.add('nan_vals')
            if trim_nan:
                for k,d in data.items():
                    data[k] = d[finite_rows]
        if trim_nan:
            nusable = int(finite_rows.sum())
        else:
            nusable = ndata_rows
    elif var_names==['index']:
        # Index is normally metadata, but an index-only file has no other
        # possible payload and its complete rows therefore count as usable.
        nusable = ndata_rows
    else:
        nusable = 0
    if nusable==0:
        issues.add('no_usable_vals')
    return ret
#end def read_scalar_file



def qmcpack_analyzer_outfiles(qmc,prefix,series,group_index=None):
    """Return expected output filenames for one QMCPACK calculation."""
    if qmc not in {'opt','vmc','dmc'}:
        raise ValueError(f'unrecognized qmc type: {qmc}')
    series_label = 's'+str(series).zfill(3)
    if group_index is None:
        prefix = f'{prefix}.{series_label}.'
    else:
        group_label = 'g'+str(group_index).zfill(3)
        prefix = f'{prefix}.{group_label}.{series_label}.'
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
    """Collect analysis information from a QMCPACK input file.

    Information extraction is best-effort.  Apart from errors constructing
    :class:`QmcpackInput` itself, unavailable or malformed input data do not
    raise exceptions.  The attributes supported by this class retain their
    initial value of ``None`` when the corresponding information cannot be
    determined.  In this event, ``incomplete`` is set to ``True``.

    Parameters
    ----------
    filepath : str or os.PathLike
        Path to the QMCPACK input file.

    Attributes
    ----------
    filepath : str or os.PathLike
        Input path supplied at construction.
    qmc_type : {``'opt'``, ``'vmc'``, ``'dmc'``} or None
        Major calculation type.  Mixed sequences are classified in the
        precedence order DMC, optimization, then VMC.
    prefix : str or None
        Project identifier used as the output-file prefix.
    group_index : int or None
        Group index obtained from a ``gNNN`` component of the input filename.
    series_start : int or None
        Initial QMCPACK series number.  A calculation sequence without an
        explicit project series starts at zero.
    has_twist : bool or None
        Whether ``twistnum`` or ``twist`` is present in the simulation input.
    qmc_info : obj or None
        Mapping from integer series numbers to per-calculation ``obj``
        records.  Each record contains ``series``, normalized ``qmc`` type,
        ``warmupsteps``, ``blocks``, ``timestep``, and ``outfiles``.  VMC and
        DMC records also contain ``steps``.  Input values override method
        defaults.  ``outfiles`` is a tuple of expected output filenames, or
        ``None`` when no project prefix is available.  Loop calculations are
        expanded in execution order.  If any QMC section is malformed, the
        complete mapping is unavailable and remains ``None``.
    incomplete : bool
        Whether any other attribute initialized by the constructor is
        ``None`` after reading.  Thus an optional value that is absent, such
        as ``group_index`` for an ungrouped input, also marks the information
        as incomplete.
    """

    def __init__(self,filepath):
        self.filepath     = filepath
        self.qmc_type     = None
        self.prefix       = None
        self.group_index  = None
        self.series_start = None
        self.has_twist    = None
        self.qmc_info     = None
        self.incomplete   = False
        self.read(filepath)
        information = (
            self.filepath,
            self.qmc_type,
            self.prefix,
            self.group_index,
            self.series_start,
            self.has_twist,
            self.qmc_info,
            )
        self.incomplete = any(value is None for value in information)
    #end def __init__

    def read(self,filepath):
        # Construction can still raise if the file itself cannot be read as a
        # QMCPACK input.  All information extraction after this point is
        # best-effort and leaves affected attributes at their None defaults.
        qi = QmcpackInput(filepath)

        # parse filepath
        filename = os.path.split(str(filepath))[1]
        ftokens  = filename.split('.')
        for token in ftokens:
            if token.startswith('g') and token[1:].isdigit():
                self.group_index = int(token[1:])
                break

        # QmcpackInput can also represent individual input elements.  Only a
        # full simulation contains the information collected here.
        if not hasattr(qi,'__contains__') or 'simulation' not in qi:
            return
        sim = qi['simulation']
        if not isinstance(sim,simulation):
            return

        # prefix, series
        series_query_failed = False
        qproject = sim['project'] if 'project' in sim else None
        if qproject is not None and not isinstance(qproject,project):
            series_query_failed = True
        elif qproject is not None:
            if 'id' in qproject and isinstance(qproject.id,str):
                self.prefix = qproject.id
            if 'series' in qproject:
                series_start = qproject.series
                if isinstance(series_start,Integral):
                    self.series_start = int(series_start)
                else:
                    series_query_failed = True

        # twist
        twistnum,twist = sim.get(('twistnum','twist'))
        self.has_twist = twistnum is not None or twist is not None

        # qmc method info
        qmc_info = obj()
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

        # extract calculation list
        if 'calculations' in sim:
            calculations = sim.calculations
        elif 'qmc' in sim:
            calculations = sim.qmc
        else:
            calculations = None
        if calculations is None:
            calculations = ()
        elif isinstance(calculations,collection):
            calculations = tuple(calculations)
        elif isinstance(calculations,QIxml):
            calculations = (calculations,)
        else:
            return
        if len(calculations)>0 and series_query_failed:
            return
        elif len(calculations)>0 and self.series_start is None:
            self.series_start = 0

        def qmc_from_loop(qmc_loop):
            """Expand a QMCPACK loop element."""
            if 'max' not in qmc_loop:
                return None
            loop_count = qmc_loop.max
            if (not isinstance(loop_count,Integral) or
                isinstance(loop_count,bool) or loop_count<0):
                return None
            if 'calculations' in qmc_loop:
                loop_qmc = qmc_loop.calculations
            elif 'qmc' in qmc_loop:
                loop_qmc = (qmc_loop.qmc,)
            else:
                loop_qmc = ()
            if isinstance(loop_qmc,collection):
                loop_qmc = tuple(loop_qmc)
            elif isinstance(loop_qmc,QIxml):
                loop_qmc = (loop_qmc,)
            elif not isinstance(loop_qmc,(tuple,list)):
                return None
            return tuple(loop_qmc)*int(loop_count)
        #end def qmc_from_loop

        def qinfo_from_qmc(qmc,prefix,series):
            """Extract basic information from a QMC input section."""
            if not isinstance(qmc,QIxml) or 'method' not in qmc:
                return None
            method = qmc.method
            if not isinstance(method,str) or method not in method_types:
                return None
            qmc_type = method_types[method]
            qinfo = obj(**defaults[qmc_type])
            for k,v in qmc.items():
                if k in qinfo:
                    qinfo[k] = v
            qinfo.series = series
            if prefix is None:
                qinfo.outfiles = None
            else:
                qinfo.outfiles = qmcpack_analyzer_outfiles(
                    qinfo.qmc,prefix,series,self.group_index
                    )
            return qinfo
        #end def qinfo_from_qmc

        # extract info from calculations
        series = self.series_start
        qmc_info_valid = True
        for calculation in calculations:
            if isinstance(calculation,loop):
                qmc_calculations = qmc_from_loop(calculation)
                if qmc_calculations is None:
                    qmc_info_valid = False
                    break
            else:
                qmc_calculations = (calculation,)
            for qmc in qmc_calculations:
                qinfo = qinfo_from_qmc(qmc,self.prefix,series)
                if qinfo is None:
                    qmc_info_valid = False
                    break
                qmc_info[series] = qinfo
                series += 1
            if not qmc_info_valid:
                break
        if not qmc_info_valid:
            self.qmc_info = None
            self.qmc_type = None
            return
        # determine the major run type
        qmc_types = {q.qmc for q in qmc_info.values()}
        if 'dmc' in qmc_types:
            self.qmc_type = 'dmc'
        elif 'opt' in qmc_types:
            self.qmc_type = 'opt'
        elif 'vmc' in qmc_types:
            self.qmc_type = 'vmc'
        self.qmc_info = qmc_info
    #end def read
#end class QmcpackInputInfo
