
from subprocess import *
import sys


def demangle(names):
  """Run c++filt on a list of mangled C++ names.  Returns a list of unmangled names."""
  cmd = ['c++filt']

  p = Popen(cmd, stdin=PIPE,stdout=PIPE,stderr=PIPE)
  inp = '\n'.join(names)
  stdout, stderr = p.communicate(inp)
  if stderr:
    print 'stderr:',stderr

  lines = stdout.split('\n')
  return lines

if __name__ == '__main__':
  plain = demangle(sys.argv[1:])
  print plain
