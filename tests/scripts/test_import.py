#! /usr/bin/env python3

import sys

# Test for existence of module in python install

if len(sys.argv) < 2:
  print('Usage: test_import.py <module>')
  sys.exit(1)

module_to_test = sys.argv[1]
#print('testing ',module_to_test)
try:
  mod = __import__(module_to_test)
  print('True')
except ImportError:
  print('False')
