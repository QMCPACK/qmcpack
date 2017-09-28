#from __future__ import print_function

import argparse
import os
import sys

def run_one_nexus_test(nexus_file):
  okay = True
  g = dict()

  g['override_generate_only_setting'] = True
  try:
    execfile(nexus_file, g)
  except Exception as e:
    print('Error ',type(e),e)
    okay = False
  except SystemExit:
    print('system exit called')
    okay = False

  if okay:
    print("  pass")
  else:
    print("  FAIL")

  return okay


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Test Nexus examples')
  parser.add_argument('test_dir',
                      help='Directory of test to run')
  parser.add_argument('--nexus-file',
                     help='Name of Nexus file for this example')
  parser.add_argument('--nexus-path',
                      help='Path to Nexus library (for setting PYTHONPATH)')
  args = parser.parse_args()

  test_dir = args.test_dir
  if not os.path.exists(test_dir):
    print("Test not found: ", test_dir)
    sys.exit(1)

  curr_dir = os.getcwd()
  os.chdir(test_dir)


  nexus_path = args.nexus_path
  if nexus_path:
    if not os.path.exists(nexus_path):
      print("Nexus path not found: ", nexus_path)
      sys.exit(1)
    sys.path.append(nexus_path)

  nexus_file = args.nexus_file
  if not os.path.exists(nexus_file):
      print("Nexus file not found: ", nexus_file)
      sys.exit(1)

  ret = run_one_nexus_test(args.nexus_file)

  os.chdir(curr_dir)

  if ret:
    sys.exit(0)
  sys.exit(1)

