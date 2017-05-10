
import argparse
from collections import namedtuple, defaultdict, OrderedDict
import difflib
import glob
import json
import os.path
import StringIO
import sys

# See https://gcc.gnu.org/onlinedocs/gcc/Invoking-Gcov.html or the gcov man page
# for information on the gcov file format.

FunctionSummary = namedtuple('FunctionSummary', ['name','called_count','returned','blocks_exe'])
# name - function name (mangled if C++)
# called_count - number of calls.  Zero if uncovered.
# returned - percentage of calls returned - usually 0% or 100%
# blocks_exe - percentage of basic blocks executed

CallInfo = namedtuple('CallInfo', ['type','returned'])
# Basic block exits due to calls

BranchInfo = namedtuple('BranchInfo', ['type','taken','is_executed','is_fallthrough','is_throw'])
# Basic block exits due to branches

LineInfo = namedtuple('LineInfo', ['count','lineno','src'])
# count - number of times the line was executed.
#         Also '-' for non-executable code or '#####' or '=====' for unexecuted code
# lineno - line number in source file
# src - source line contents

# --- Return from read_gcov
GcovFile = namedtuple('GcovFile', ['fname', 'tags', 'line_info', 'function_info', 'call_branch_info', 'branch_info','func_ranges'])
# fname - gcov file name
# tags - dict mapping tag name to value (tags['Source'] is the original source file)
# line_info - ordered dict mapping line num to LineInfo
# function_info - dict mapping line number to list of FunctionSummary
# call_branch_info -
# branch_info -
# func_ranges -  FunctionRange object that maps line numbers to function name

# --- Return from summarize_all_files ---
CoverageSummary = namedtuple('CoverageSummary', ['covered_files','uncovered_files','total_files','covered_functions','uncovered_functions', 'total_functions'])



# map line #'s ranges to function
class FunctionRange(object):
  def __init__(self):
    self.ranges = []

  def set_range(self, func, start, end):
    self.ranges.append( (start,end,func))

  def find_func(self, line):
    funcs = []
    for start,end,func in self.ranges:
      if line >= start and line <= end:
        funcs.append(func)
    return funcs


def get_keyword_value(kw, pos, parts):
    """
        Match expected keyword to a position in an array.
        The array 'parts' is most likely generated from splitting a line.
        Returns the item in parts after the keyword.
        Some keywords can have spaces, which means they occupy multiple slots in the
        parts array.
    """

    if len(parts) < pos:
        print 'Not enough parts in line for kw %s'%kw
        print '   line parts: %s'%parts
        return ''

    kws = kw.split()
    offset = len(kws)
    for i,kw in enumerate(kws):
        if parts[pos+i] != kw:
            print 'error, keyword expected:%s, found: %s'%(kw, parts[pos+i])
            print '   line parts: %s'%parts
            return ''

    if len(parts) < pos+offset:
        print 'Not enough parts to line for kw value %s'%kw
        print '   line parts: %s'%parts
        return ''

    return parts[pos+offset]



def get_percent(s):
    '''Convert strings of the form '<x>%' to integer'''

    if len(s) == 0:
        print 'get_percent empty string'
        return 0

    if not s.endswith('%'):
        print 'not percent: %s'%s
        return int(s)

    if len(s) == 1:
        return 0

    return int(s[:-1])


def read_gcov(fname):
    with open(fname, 'r') as f:
        line_info = OrderedDict()
        func_info = dict()
        func_location = defaultdict(list)  # which line does the function info follow?
        #call_branch_location = defaultdict(list)
        branch_location = defaultdict(list)
        tags = OrderedDict()
        lines_with_info = set()
        current_src_line = 0
        current_function = None
        func_ranges = FunctionRange()
        prev_function_start_line = 0
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('function'):
                parts = line.split()
                func_name = get_keyword_value('function', 0, parts)
                called_count = int(get_keyword_value('called', 2, parts))
                ret_percent = get_percent(get_keyword_value('returned', 4, parts))
                block_exe_percent = get_percent(get_keyword_value('blocks executed', 6, parts))

                prev_function_end_line = current_src_line
                if current_function:
                  # Empty functions, like static initializers
                  tmp_end = prev_function_end_line
                  if prev_function_end_line < prev_function_start_line:
                    tmp_end = prev_function_start_line
                  func_ranges.set_range(current_function, prev_function_start_line, tmp_end)

                current_function = func_name
                prev_function_start_line = current_src_line + 1

                fs = FunctionSummary(func_name, called_count, ret_percent, block_exe_percent)
                func_info[func_name] = fs
                func_location[current_src_line].append(fs)
                lines_with_info.add(current_src_line)
                continue

            if line.startswith('call'):
                parts = line.split()
                exed = 0
                if 'returned' in line:
                    exed = get_percent(get_keyword_value('returned', 2, parts))
                elif 'never executed' not in line:
                    print "Unexpected value in call branch info: %s"%parts

                lines_with_info.add(current_src_line)
                branch_location[current_src_line].append(CallInfo('call', exed))

                continue

            if line.startswith('branch'):
                parts = line.split()
                taken = 0
                is_executed = True
                is_fallthrough = 'fallthrough' in line
                is_throw = 'throw' in line
                if 'taken' in line:
                    taken = get_percent(get_keyword_value('taken', 2, parts))
                elif 'never executed' in line:
                    is_executed = False
                else:
                    print "Unexpected value in branch info: %s"%parts

                lines_with_info.add(current_src_line)
                branch_location[current_src_line].append(BranchInfo('branch', taken, is_executed, is_fallthrough, is_throw))
                continue

            parts = line.split(':')
            exe_count = parts[0]

#  Lines in .gcov files:
#  Normal tag line (marked by line number is zero)
#        -:    0:Source:sample.cpp
#  Error message in tags
#        -:    0:Source is newer than graph
#  Normal line
#        1:   26:  int j = 0;

            if len(parts) > 1:
                line_no = int(parts[1])
                if line_no == 0:
                    if parts[2] == 'Source is newer than graph':
                      tag = 'error'
                      value = parts[2]
                    else:
                      tag = parts[2]
                      value = parts[3]
                    tags[tag] = value
                else:
                    parts = line.split(':',2)
                    current_src_line = line_no
                    line_info[line_no] = LineInfo(exe_count, line_no, parts[2])
            else:
                print 'Unhandled line: ',parts

    # Last function in the file
    if current_function:
      func_ranges.set_range(current_function, prev_function_start_line, current_src_line+1)

    return GcovFile(fname, tags, line_info, func_location, None, branch_location, func_ranges)

def write_gcov(gcov, f=sys.stdout):
    for t,v in gcov.tags.iteritems():
        if t == 'error':
          print >>f, '%8s-:%5d:%s'%(' ', 0, v)
        else:
          print >>f, '%8s-:%5d:%s:%s'%(' ', 0, t, v)

    for line_no, line_info in gcov.line_info.iteritems():
        assert(line_no == line_info.lineno)
        print >>f, '%9s:%5d:%s'%(line_info.count, line_no, line_info.src)
        #if line_no in gcov.call_branch_info:
        #    for i,cb in enumerate(gcov.call_branch_info[line_no]):
        #        if cb.returned == 0:
        #            print 'call    %d never executed'%i
        #        else:
        #            print 'call    %d returned %d%%'%(i, cb.returned)

        if gcov.branch_info and line_no in gcov.branch_info:
            for i,br in enumerate(gcov.branch_info[line_no]):
                if br.type == 'call':
                    if br.returned == 0:
                        print >>f, 'call    %d never executed'%i
                    else:
                        print >>f, 'call    %d returned %d%%'%(i, br.returned)
                elif br.type == 'branch':
                    if not br.is_executed:
                        print >>f, 'branch  %d never executed'%i
                    else:
                        ft = ''
                        if br.is_fallthrough:
                            ft = ' (fallthrough)'
                        if br.is_throw:
                            ft = ' (throw)'
                        print >>f, 'branch  %d taken %d%%%s'%(i, br.taken, ft)
                else:
                    print 'unknown branch info type: ',br.type

        if gcov.function_info and line_no in gcov.function_info:
            for f1 in gcov.function_info[line_no]:
                print >>f, 'function %s called %d returned %d%% blocks executed %d%%'%(f1.name, f1.called_count, f1.returned, f1.blocks_exe)

def read_gcov_files_glob(target_dir):
    flist = glob.glob(os.path.join(target_dir,'*.gcov'))
    return read_gcov_files(flist)


def read_gcov_files(flist):
    all_files = dict()
    for fname in flist:
        #func_info, tags, line_info = read_gcov(fname)
        gcov = read_gcov(fname)
        source = gcov.tags['Source']
        all_files[source] = gcov

    return all_files


def summarize_all_files(all_files):
    covered_file = set()
    uncovered_file = set()
    covered_func = set()
    uncovered_func = set()
    total_files = 0
    total_func = 0
    for fname, val in all_files.iteritems():
        total_files += 1
        file_touched = False
        for f_line_no,func_infos in val.function_info.iteritems():
            for func_info in func_infos:
                total_func += 1
                if func_info.called_count == 0:
                    if func_info.name not in covered_func:
                        uncovered_func.add(func_info.name)
                else:
                    if func_info.name in uncovered_func:
                        uncovered_func.discard(func_info.name)

                    covered_func.add(func_info.name)
                    file_touched = True

        if file_touched:
            covered_file.add(fname)
        else:
            uncovered_file.add(fname)

    return CoverageSummary(covered_file, uncovered_file, total_files, covered_func, uncovered_func, total_func)

def print_summary(cov_sum):
    print 'Total files: %d Touched files: %d  Uncovered files:%d'% \
            (cov_sum.total_files, len(cov_sum.covered_files), len(cov_sum.uncovered_files))
    print 'Total func: %d  Touched func: %d   Uncovered func:%d'% \
            (cov_sum.total_functions, len(cov_sum.covered_functions), len(cov_sum.uncovered_functions))


def export_as_json(all_files, export_name='export.json'):
    """
        Write JSON in same format as llvm-cov export.
        Only works for function summaries.
        No source ranges.
    """
    root = dict()
    root['version'] = '2.0.0'
    root['type'] = 'convert.gcov.json.export'

    file_list = []
    func_list = []
    for fname, gcov in all_files.iteritems():
        file_entry = dict()
        file_entry['filename'] = fname
        file_entry['expansions'] = []
        file_entry['segments'] = []

        count = 0
        covered = 0
        func_entry = dict()
        for func, called in gcov.function_info.iteritems():
            func_entry['name'] = func
            func_entry['count'] = called
            func_entry['filenames'] = [fname]
            func_entry['regions'] = []
            count += 1
            if called > 0:
                covered += 1

        percent_covered = 0
        if count > 0:
            percent_covered = 100*covered/count
        func_summary = {'count':count,'covered':covered, 'percent':percent_covered}

        file_entry['summary'] = {'functions':func_summary}
        file_list.append(file_entry)

        func_list.append(func_entry)

    root['data'] = [{'files':file_list, 'functions':func_list}]
    fout = open(export_name, 'w')
    json.dump(root, fout, indent=2)
    fout.close()

if __name__ == '__main__':

    desc = """
    Read gcov files.

    The --action (-a) option specifies what action to take:
      summarize (s) - print summary of covered/uncovered lines and functions (default)
      write (w)     - write gcov file to <input_base>.out.gcov
      test (t)      - compare initial gcov file with file after read/write to ensure idempotency
      export (e)    - export to json file (unfinished)
      dump (d)      - Print sections of file to stdout
    """
    parser = argparse.ArgumentParser(description=desc,formatter_class=argparse.RawDescriptionHelpFormatter)

    summarize_action = ['summarize','s']
    write_action = ['write','w']
    export_action = ['export','e']
    dump_action = ['dump','d']
    test_action = ['test','t']
    actions = summarize_action + write_action + export_action + dump_action + test_action
    parser.add_argument('input_files',nargs='*')
    parser.add_argument('-a','--action',default='summarize',choices=actions)
    parser.add_argument('-o','--output')
    args = parser.parse_args()

    print 'input files = ', args.input_files
    if args.action in summarize_action:
        all_files = read_gcov_files(args.input_files)
        cov_sum = summarize_all_files(all_files)
        print_summary(cov_sum)

    if args.action in export_action:
        export_file = 'export.json'
        if args.output:
          export_file = args.output
        all_files = read_gcov_files(args.input_files)
        export_as_json(all_files, export_file)

    if args.action in write_action:
        output = '.out'
        if args.output:
          output = args.output

        for input_file in args.input_files:
          gcov = read_gcov(input_file)
          output_name = input_file.replace('.gcov',output + '.gcov')
          if input_file == output_name:
            print 'Would overwrite input file, skipping:',input_file
          else:
            write_gcov(gcov, open(output_name, 'w'))

    if args.action in test_action:
        for input_file in args.input_files:
          gcov = read_gcov(input_file)
          output_buf = StringIO.StringIO()
          write_gcov(gcov, output_buf)
          #output_buf.write('extra')  # how to test the diff
          output_buf.seek(0)
          output_name = input_file.replace('.gcov','out.gcov')
          input_buf = open(input_file,'r')
          diff = difflib.context_diff(input_buf.readlines(), output_buf.readlines(),input_file,output_name)
          for d in diff:
            print d,

    if args.action in dump_action:
        fname = args.input_files[0]
        gcov = read_gcov(fname)
        print '** File tags **'
        for k,v in gcov.tags.iteritems():
            print k,v
        print '** Function info **'
        for f,c in gcov.function_info.iteritems():
            print f,c
        print '** line info **'
        for lineno, c in gcov.line_info.iteritems():
            print lineno, c

