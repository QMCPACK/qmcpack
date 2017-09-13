#!/usr/bin/env python

if __name__ == '__main__':
  import os
  cwd = os.getcwd()
  fname = os.path.join(cwd,'../../../Testing/Temporary/LastTest.log')

  from mmap import mmap
  with open(fname,'r+') as f:
    mm = mmap(f.fileno(),0)
  # end with

  fsuccess = mm.find('Finite difference test: PASS') != -1
  rsuccess = mm.find('Ratio test: PASS') != -1

  if fsuccess and rsuccess:
    exit(0)
  else:
    exit(1)

# end __main__
