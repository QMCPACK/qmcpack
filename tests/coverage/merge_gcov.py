
import read_gcov
import argparse
from collections import OrderedDict
import os.path
import sys

def merge_gcov_files(fnames, output_fname, src_prefix_to_add=None):
  gcovs = [read_gcov.read_gcov(f) for f in fnames]
  source_names = [g.tags['Source'] for g in gcovs]
  for source in source_names[1:]:
    if source_names[0] != source:
      print 'Source file names do not match',source_names[0],source
      return

  gcov_out = merge_gcovs(gcovs)

  # The gcov '--source-prefix' options strips a directory from the file name and the
  #  'Source' tag.   Restore it to the tag so gcovr can find the source file.
  if src_prefix_to_add is not None:
    src_filename = gcov_out.tags['Source']
    gcov_out.tags['Source'] = os.path.join(src_prefix_to_add, src_filename)

  out_f = sys.stdout
  if output_fname:
    out_f = open(output_fname, 'w')
  read_gcov.write_gcov(gcov_out, out_f)


def merge_gcovs(gcovs):
  output_line_info = OrderedDict()

  for line in gcovs[0].line_info.iterkeys():
      lines = [g.line_info[line] for g in gcovs]

      Uncov_norm= '#####'
      Uncov_throw= '====='
      Uncov = [Uncov_norm, Uncov_throw]
      Nocode = '-'

      output_count = None

      lines_count = []
      for lc in lines:
        line_count = 0
        try:
          line_count = int(lc.count)
        except:
          pass
        lines_count.append(line_count)

      any_lines_covered = reduce(lambda x,y:x or y, [a > 0 for a in lines_count])
      any_lines_uncovered = reduce(lambda x,y:x or y, [a.count in Uncov for a in lines])
      all_lines_nocode = reduce(lambda x,y:x and y, [a.count == Nocode for a in lines])

      if any_lines_covered:
        output_count = sum(lines_count)
      elif any_lines_uncovered:
        output_count = Uncov_norm
      elif all_lines_nocode:
        output_count = Nocode

      #if line1_count > 0 or line2_count > 0:

      #  output_count = line1_count + line2_count
      #elif line1.count == Uncov_norm or line2.count == Uncov_norm:
      #  output_count = Uncov_norm
      #elif line1.count == Nocode and line2.count == Nocode:
      #  output_count = Nocode

      if output_count is None:
        print 'Unhandled situation:'
        for idx,line in enumerate(lines):
          print '  line %d: '%idx,line

      output_line_info[line] = read_gcov.LineInfo(str(output_count), line, lines[0].src)


  return read_gcov.GcovFile(gcovs[0].fname, gcovs[0].tags, output_line_info, None, None, None, gcovs[0].func_ranges)


# this function was replaced by the more general merge_gcov_files that can
#  merge an arbitrary number of files.  However, it's easier to understand
# the logic for two files, so this function is left here.
def merge_two_gcov_files(fname1, fname2, output_fname):
  gcov1 = read_gcov.read_gcov(fname1)
  gcov2 = read_gcov.read_gcov(fname2)

  source_file1 = gcov1.tags['Source']
  source_file2 = gcov2.tags['Source']
  if source_file1 != source_file2:
    print 'Source files do not match',source_file1,source_file2
    return

  gcov_out = merge_two_gcovs(gcov1, gcov2)

  out_f = sys.stdout
  if output_fname:
    out_f = open(output_fname, 'w')
  read_gcov.write_gcov(gcov_out, out_f)


def merge_two_gcovs(gcov1, gcov2):
  output_line_info = OrderedDict()

  for line in gcov1.line_info.iterkeys():
      line1 = gcov1.line_info[line]
      line2 = gcov2.line_info[line]

      Uncov_norm= '#####'
      Nocode = '-'

      output_count = None

      line1_count = 0
      try:
        line1_count = int(line1.count)
      except:
        pass

      line2_count = 0
      try:
        line2_count = int(line2.count)
      except:
        pass

      if line1_count > 0 or line2_count > 0:
        output_count = line1_count + line2_count
      elif line1.count == Uncov_norm or line2.count == Uncov_norm:
        output_count = Uncov_norm
      elif line1.count == Nocode and line2.count == Nocode:
        output_count = Nocode

      if output_count is None:
        print 'Unhandled situation:'
        print '  line1: ',line1
        print '  line2: ',line2

      output_line_info[line] = read_gcov.LineInfo(str(output_count), line, line1.src)

  return read_gcov.GcovFile(gcov1.fname, gcov1.tags, output_line_info, None, None, None, gcov1.func_ranges)


if __name__ == '__main__':
  desc = """
  Merge different gcov files corresponding to the same source
  """

  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument('input_files',nargs='*')
  parser.add_argument('-o','--output')

  args = parser.parse_args()

  #merge_two_gcov_files(args.input_files[0], args.input_files[1], args.output)
  merge_gcov_files(args.input_files, args.output)
