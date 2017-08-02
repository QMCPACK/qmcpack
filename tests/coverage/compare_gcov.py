
import argparse
from collections import OrderedDict, namedtuple, defaultdict
import glob
import os
import sys
import demangle
import read_gcov
import merge_gcov

def get_gcov_files(dir_name):
  """
    Get list of *.gcov files in a directory.  Returns list of bare filenames.
  """
  fs = glob.glob(os.path.join(dir_name,'*.gcov'))
  fs = [os.path.basename(f) for f in fs]
  return fs


def keep_file(gcov):
  """
      Some files should not be considered in coverage.
      Files to skip include system files (in /usr), the unit tests themselves
      (in tests/ directories) and the unit test infrastructure (in external_codes/)
  """
  source_file = gcov.tags['Source']

  path_elements = source_file.split('/')

  # Only looking at specific depths of directories. This might need
  # to be expanded, for example if unit test directories contain
  # subdirectories.
  try:
    if path_elements[-2] == 'tests':
      #print 'Unit test, skipping'
      return False
    if path_elements[-3] == 'external_codes':
      #print 'External code, skipping'
      return False
  except IndexError:
    pass

  if source_file.startswith('/usr/'):
    return False

  return True

def remove_unwanted_file(gcov, fname, dir_base):
  if keep_file(gcov):
    if get_total_covered_lines(gcov).covered > 0:
      return True


  if fname.endswith('.gcov'):
    os.unlink(os.path.join(dir_base,fname))
  else:
    print 'Filter error, attempting to remove a non-gcov file: ',fname
  return False

# There are two locations for the source filename associated with each gcov
# file.  First, the actual filename, without the '.gcov' suffix. If the -p
# option to gcov is used, the full path is contained in the filename, with
# the directory separator replaced with '#'.
# The second location is the 'Source:' tag inside the gcov file.


def read_and_filter_gcov_files(fnames, directory):
  new_fnames = set()
  gcov_map = dict()
  for fname in fnames:
    gcov = read_gcov.read_gcov(os.path.join(directory, fname))
    keep = remove_unwanted_file(gcov, fname, directory)
    if keep:
      source_file = gcov.tags['Source']
      new_fnames.add(fname)
      gcov_map[fname] = gcov

  return new_fnames, gcov_map

def merge_gcov_files_in_dir(fnames, directory, output_dir=None, src_prefix=None):
  to_merge = defaultdict(list)
  for fname in fnames:
    # Files from 'gcov -l' have '##' to separate the two parts of the path
    if '##' in fname:
      original_name, src_name = fname.split('##')
      to_merge[src_name].append(fname)
    else:
      to_merge[fname].append(fname)

  if output_dir is None:
    output_dir = ''

  for output_fname, input_fnames in to_merge.iteritems():
      inputs = [os.path.join(directory, fname) for fname in input_fnames]
      merge_gcov.merge_gcov_files(inputs, os.path.join(output_dir,output_fname),
                                  src_prefix_to_add=src_prefix)


def compare_gcov_dirs(dir_base, dir_unit):
  base = set(get_gcov_files(dir_base))
  unit = set(get_gcov_files(dir_unit))

  base_names, base_gcov_map = read_and_filter_gcov_files(base, dir_base)
  unit_names, unit_gcov_map = read_and_filter_gcov_files(unit, dir_unit)

  both = unit.intersection(base_names)
  only_base = base_names.difference(unit_names)
  only_unit = unit_names.difference(base_names)

  print 'Files in both: ',len(both)
  print 'Files only in base: ',len(only_base)
  print 'Files only in unit: ',len(only_unit)

  return both, only_base, only_unit, base_gcov_map, unit_gcov_map


def compare_gcov_files(both, only_base, only_unit, base_gcov_map, unit_gcov_map, dir_unit, dir_diff):
  for fname in both:
    gcov_base = base_gcov_map[fname]
    gcov_unit = unit_gcov_map[fname]
    gcov_diff = compute_gcov_diff(gcov_base, gcov_unit)
    out_fname = os.path.join(dir_diff, fname)
    if gcov_diff:
      read_gcov.write_gcov(gcov_diff, open(out_fname, 'w'))

  # completely uncovered in unit tests
  # Assign all to unit tests
  # Assign to diff only if some coverage in base
  for fname in only_base:
    gcov_base = base_gcov_map[fname]
    handle_uncovered(gcov_base, fname, dir_diff)
    handle_uncovered(gcov_base, fname, dir_unit, always_copy=True)

def handle_uncovered(gcov_base, fname, dir_diff, always_copy=False):
  # Need to see if there are any covered lines in the file
  file_coverage = get_total_covered_lines(gcov_base)
  if always_copy or file_coverage.covered > 0:
    gcov_diff = mark_as_uncovered(gcov_base)
    out_fname = os.path.join(dir_diff, fname)
    if gcov_diff:
      read_gcov.write_gcov(gcov_diff, open(out_fname, 'w'))

def mark_as_uncovered(gcov_base):
  diff_line_info = OrderedDict()
  for line in gcov_base.line_info.iterkeys():
    base_line = gcov_base.line_info[line]
    Uncov_norm= '#####'
    Nocode = '-'

    diff_count = base_line.count
    try:
      base_count = int(base_line.count)
      diff_count = Uncov_norm
    except ValueError:
      if base_line.count == Uncov_norm:
        diff_count = Nocode


    diff_line_info[line] = read_gcov.LineInfo(diff_count, line, base_line.src)

  return read_gcov.GcovFile(gcov_base.fname, gcov_base.tags, diff_line_info, None, None, None, gcov_base.func_ranges)


class FileCoverage:
  def __init__(self, covered=0, uncovered=0, total=0):
    self.covered = covered
    self.uncovered = uncovered
    self.total = total

  def __add__(self, o):
    return FileCoverage(self.covered + o.covered, self.uncovered + o.uncovered, self.total + o.total)


class CompareCoverage:
  def __init__(self, base=FileCoverage(), unit=FileCoverage(), rel=FileCoverage()):
    self.base = base
    self.unit = unit
    self.rel = rel

  def __add__(self, o):
    return CompareCoverage(self.base + o.base, self.unit + o.unit, self.rel + o.rel)


def get_total_covered_lines(gcov):
  if not keep_file(gcov):
    return FileCoverage(0, 0, 0)

  total_covered_lines = 0
  total_uncovered_lines = 0
  total_lines = 0

  Nocode = '-'

  for line in gcov.line_info.iterkeys():
    base_line = gcov.line_info[line]


    base_count = 0
    try:
      base_count = int(base_line.count)
    except ValueError:
      pass

    use_line = True
    if gcov.func_ranges:
      funcs = gcov.func_ranges.find_func(line)
      for func_name in funcs:
        if func_name:
          if func_name.startswith("_GLOBAL__"):
            use_line = False
            base_count = 0
          if 'static_initialization_and_destruction' in func_name:
            use_line = False
            base_count = 0

    if use_line:
      if base_line != Nocode:
        total_lines += 1

      if base_count == 0:
        total_uncovered_lines += 1
      else:
        total_covered_lines += 1

  return FileCoverage(total_covered_lines, total_uncovered_lines, total_lines)


def compute_gcov_diff(gcov_base, gcov_unit, print_diff=False, print_diff_summary=False, coverage_stats=None):
  diff_line_info = OrderedDict()
  if print_diff or print_diff_summary:
    print 'file ',gcov_base.tags['Source']

  total_base_count = 0

  base_totals_total = 0
  base_totals_covered = 0
  base_totals_uncovered = 0

  unit_totals_total = 0
  unit_totals_covered = 0
  unit_totals_uncovered = 0
  unit_totals_rel_covered = 0
  unit_totals_rel_uncovered = 0

  for line in gcov_base.line_info.iterkeys():
    if line not in gcov_unit.line_info:
      print 'error, line not present: %d ,line=%s'%(line, gcov_base.line_info[line])
    base_line = gcov_base.line_info[line]
    unit_line = gcov_unit.line_info[line]

    Uncov_norm = '#####'
    Uncov_exp  = '====='
    Uncov = [Uncov_norm, Uncov_exp]
    Nocode = '-'
    diff_count = None
    base_count = 0
    try:
      base_count = int(base_line.count)
    except ValueError:
      pass

    unit_count = 0
    try:
      unit_count = int(unit_line.count)
    except ValueError:
      pass

    # Skip some compiler-added functions
    #  Assuming base and unit are the same
    if gcov_base.func_ranges:
      funcs = gcov_base.func_ranges.find_func(line)
      for func_name in funcs:
        if func_name:
          if func_name.startswith("_GLOBAL__"):
            base_count = 0
            unit_count = 0
          if 'static_initialization_and_destruction' in func_name:
            base_count = 0
            unit_count = 0

    total_base_count += base_count


    if base_line.count != Nocode:
      base_totals_total += 1
    if base_line.count in Uncov:
      base_totals_uncovered += 1
    if base_count > 0:
      base_totals_covered += 1

    if unit_line.count != Nocode:
      unit_totals_total += 1
    if unit_line.count in Uncov:
      unit_totals_uncovered += 1
    if unit_count > 0:
      unit_totals_covered += 1

    line_has_diff = False

    #  base_line.count is the string value for the count
    #  base_count is the converted integer value
    #     will be 0 for No code or uncovered

    if base_line.count == Nocode and unit_line.count == Nocode:
      diff_count = Nocode

    # doesn't work well if one side is nocode and the other has count 0
    #elif base_line.count == Nocode and unit_line.count != Nocode:
    #  diff_count = Nocode
    #  line_has_diff = True
    #elif base_line.count != Nocode and unit_line.count == Nocode:
    #  diff_count = Nocode
    #  line_has_diff = True

    elif base_line.count == Nocode and unit_count == 0:
      diff_count = Nocode
    elif base_count == 0 and unit_line.count == Nocode:
      diff_count = Nocode

    elif base_line.count in Uncov and unit_line.count in Uncov:
      #diff_count = 9999   # special count to indicate uncovered in base?
      diff_count = Nocode  # Not sure of the right solution

    elif base_count != 0 and unit_count == 0:
      diff_count = Uncov_norm
      unit_totals_rel_uncovered += 1
      line_has_diff = True
    elif base_count == 0 and unit_count != 0:
      diff_count = Nocode
      line_has_diff = True
    elif base_count != 0 and unit_count != 0:
      diff_count = unit_count
      unit_totals_rel_covered += 1
    elif base_count == 0 and unit_count == 0:
      diff_count = 0

    if print_diff and line_has_diff:
      if gcov_base.line_info[line].src.strip() != gcov_unit.line_info[line].src.strip():
        print 'line diff, base: ',gcov_base.line_info[line]
        print '     unit: ',gcov_unit.line_info[line]
      print '%9s  %9s %6d : %s'%(base_line.count, unit_line.count, line, gcov_unit.line_info[line].src)
    if diff_count is None:
      print 'Unhandled case: ',line,diff_count,base_line.count,unit_line.count

    diff_line_info[line] = read_gcov.LineInfo(diff_count, line, base_line.src)

  if print_diff_summary:
    print 'Base        Unit'
    print '%4d/%-4d   %4d/%-4d  covered lines'%(base_totals_covered, base_totals_total, unit_totals_covered, unit_totals_total)
    print '%4d/%-4d   %4d/%-4d  uncovered lines'%(base_totals_uncovered, base_totals_total, unit_totals_uncovered, unit_totals_total)
    print '            %4d/%-4d  uncovered lines relative to base'%(unit_totals_rel_uncovered, (unit_totals_rel_uncovered+unit_totals_rel_covered))


  if coverage_stats is not None:
    base_coverage = FileCoverage(base_totals_covered, base_totals_uncovered, base_totals_total)
    unit_coverage = FileCoverage(unit_totals_covered, unit_totals_uncovered, unit_totals_total)
    rel_coverage = FileCoverage(unit_totals_rel_covered, unit_totals_rel_uncovered, unit_totals_rel_covered + unit_totals_rel_uncovered)

    cc = CompareCoverage(base_coverage, unit_coverage, rel_coverage)
    coverage_stats[gcov_base.tags['Source']] = cc

  if total_base_count == 0:
    return None

  return read_gcov.GcovFile(gcov_base.fname, gcov_base.tags, diff_line_info, None, None, None, gcov_base.func_ranges)


FunctionCoverageInfo = namedtuple('FunctionCoverageInfo',['uncovered','total','names'])

def compute_gcov_func_diff(gcov_base, gcov_unit, print_diff=False, print_diff_summary=False, func_coverage_stats=None):
  if print_diff or print_diff_summary:
    print 'file ',gcov_base.tags['Source']

  uncovered_funcs = []
  nfunc = 0
  for func_line in gcov_base.function_info.keys():
    base_func = gcov_base.function_info[func_line]
    unit_func = gcov_unit.function_info[func_line]

    nfunc += len(base_func)

    if len(unit_func) == 0 :
      for idx in range(len(base_func)):
        uncovered_funcs.append(base_func[idx].name)
      # function summary not even printed.  Seems to be an issue with templates
      continue
    for idx in range(min(len(base_func), len(unit_func))):
      if base_func[idx].called_count > 0 and unit_func[idx].called_count == 0:
        uncovered_funcs.append(gcov_base.function_info[func_line][idx].name)

  nfunc_uncovered = len(uncovered_funcs)
  if func_coverage_stats is not None:
    func_coverage_stats[gcov_base.tags['Source']] = FunctionCoverageInfo(nfunc_uncovered,nfunc,uncovered_funcs)


def print_function_coverage(func_coverage):
  print 'Function coverage'
  for name,fc in func_coverage.iteritems():
    print name,fc.total,fc.uncovered
    for func in demangle.demangle(fc.names):
      print '  ',func

def by_uncovered(x,y):
  if x[1].uncovered == y[1].uncovered:
    return 0;
  if x[1].uncovered > y[1].uncovered:
    return -1
  return 1


def print_coverage_summary(coverage_stats):
  sorted_keys = sorted(coverage_stats.iteritems(), cmp=by_uncovered)

  for source,stats in sorted_keys:
    percent = 0

    if stats.total > 0:
      percent = 100.0*stats.covered/stats.total
    print source,'%4d/%-4d'%(stats.covered,stats.total),stats.uncovered,'%.2f'%percent

def summarize_coverage_summary(coverage_stats):
    by_dirs = defaultdict(FileCoverage)
    for source, stats in coverage_stats.iteritems():
      src = source
      while True:
        dirname = os.path.dirname(src)
        if stats.total > 0:
          #by_dirs[dirname] += stats
          by_dirs[dirname] += stats
        if dirname == '':
          break
        if os.path.basename(src) in ['src','/','build']:
          break
        if src == '/':
          break
        src = dirname

    print '\n by Dir \n'
    ordered = OrderedDict()
    for k in sorted(by_dirs.keys()):
      ordered[k] = by_dirs[k]
    print_coverage_summary(ordered)

def compare_gcov(fname, dir_base, dir_unit, dir_diff=None):
  gcov_base = read_gcov.read_gcov(os.path.join(dir_base, fname))
  gcov_unit = read_gcov.read_gcov(os.path.join(dir_unit, fname))
  gcov_diff = compute_gcov_diff(gcov_base, gcov_unit)
  if dir_diff and gcov_diff:
    out_fname = os.path.join(dir_diff, fname)
    read_gcov.write_gcov(gcov_diff, open(out_fname, 'w'))

def print_gcov_diff(both, only_base, only_unit, base_gcov_map, unit_gcov_map):
  coverage_summary = dict()
  func_coverage_summary = dict()
  for fname in both:
    gcov_base = base_gcov_map[fname]
    gcov_unit = unit_gcov_map[fname]
    compute_gcov_diff(gcov_base, gcov_unit, print_diff_summary=False, print_diff=False, coverage_stats=coverage_summary)

    #print 'Function Info'
    compute_gcov_func_diff(gcov_base, gcov_unit, print_diff=False, func_coverage_stats=func_coverage_summary)

  print '\nCompletely uncovered files (relative to base)\n'
  for fname in only_base:
    gcov_base = base_gcov_map[fname]
    base_coverage = get_total_covered_lines(gcov_base)
    if base_coverage.covered > 0:
      unit_coverage = FileCoverage(0, 0, 0)
      rel_coverage = FileCoverage(0, base_coverage.covered, base_coverage.covered)
      src_name = gcov_base.tags['Source']
      if True:
        print src_name,'%4d/%-4d'%(base_coverage.covered, base_coverage.total)

      coverage_summary[src_name] = CompareCoverage(base_coverage, unit_coverage, rel_coverage)

  rel_coverage_summary = {src_name:x.rel for src_name,x in coverage_summary.iteritems()}
  print_coverage_summary(rel_coverage_summary)
  summarize_coverage_summary(rel_coverage_summary)


  func_coverage_summary2 = {src_name:FileCoverage(x.total-x.uncovered, x.uncovered, x.total) for src_name,x in func_coverage_summary.iteritems()}
  print_function_coverage(func_coverage_summary)
  summarize_coverage_summary(func_coverage_summary2)



if __name__ ==  '__main__':
  parser = argparse.ArgumentParser(description="Compare GCOV files")

  # Compare gcov files in base and unit directories and write results to output directory
  # This is used to generate the gcov files for the CDash report
  compare_action = ['compare','c']

  # Just delete unwanted files from directory
  process_action = ['process','p']

  # Compare gcov files in base and unit directories and print results
  diff_action = ['diff','d']

  # Merge files with same source (from gcov -l)
  merge_action = ['merge','m']

  actions = compare_action + process_action + diff_action + merge_action

  parser.add_argument('-a','--action',default='compare',choices=actions)
  parser.add_argument('--base-dir',
                      required=True,
                      help="Directory containing the input base gcov files")
  parser.add_argument('--unit-dir',
                      help="Directory containing the input unit test (target) gcov files")
  parser.add_argument('--output-dir',
                      help="Directory to write the difference gcov files")
  parser.add_argument('-f','--file',
                      help="Limit analysis to single file")
  parser.add_argument('-p','--prefix',
                      help="Source prefix to add when merging files")
  args = parser.parse_args()

  if args.action in compare_action:
    if not args.unit_dir:
      print '--unit-dir required for compare'
      sys.exit(1)

    if not args.output_dir:
      print '--output-dir required for compare'
      sys.exit(1)

    both, only_base, only_unit, base_gcov_map, unit_gcov_map = compare_gcov_dirs(args.base_dir, args.unit_dir)
    compare_gcov_files(both, only_base, only_unit, base_gcov_map, unit_gcov_map, args.unit_dir, args.output_dir)

  if args.action in process_action:
    base = set(get_gcov_files(args.base_dir))
    read_and_filter_gcov_files(base, args.base_dir)

  if args.action in merge_action:
    base = set(get_gcov_files(args.base_dir))
    merge_gcov_files_in_dir(base, args.base_dir, args.output_dir, args.prefix)

  if args.action in diff_action:
    if not args.unit_dir:
      print '--unit-dir required for diff'
      sys.exit(1)

    if args.file:
      compare_gcov(args.file, args.base_dir, args.unit_dir)
    else:
      both, only_base, only_unit, base_gcov_map, unit_gcov_map = compare_gcov_dirs(args.base_dir, args.unit_dir)

      print_gcov_diff(both, only_base, only_unit, base_gcov_map, unit_gcov_map)
