from __future__ import print_function

import xml.etree.ElementTree as ET
import os.path
import sys
import os

# Read timing information from the .info.xml file and output highlights
# in a form that can be read by CDash.


def get_timer_from_list(timers, timer_names):
  for timer_name in timer_names:
    for timer in timers:
      if timer.find('name').text == timer_name:
        return timer
  return None


def get_incl_time(timers, timer_names):
  timer = get_timer_from_list(timers, timer_names)
  time_incl = 0.0
  if timer is not None:
    time_incl = float(timer.find('time_incl').text)
  return time_incl


def get_performance_info(info_fname):
  tree = ET.parse(info_fname)
  timing = tree.find('timing')
  timers = timing.findall(".//timer")
  # Alternative with XPath syntax to find a particular timer
  # vmc_timers = timing.findall(".//timer[name='VMCSingleOMP']")

  vmc_time = get_incl_time(timers, ['VMCSingleOMP', 'VMCcuda'])
  dmc_time = get_incl_time(timers, ['DMCOMP', 'DMCcuda'])

  return {'VMC Time': vmc_time, 'DMC Time': dmc_time}


# Output on stdout in the right format will be included in the CDash information.
# Format found in this email:
#   https://cmake.org/pipermail/cmake/2010-April/036574.html
def print_for_cdash(vals):
  for k, v in vals.iteritems():
    print('<DartMeasurement name="%s" type="numeric/double">%g</DartMeasurement>' % (k, v))


def get_info_file(fname):
  """
     Construct the info file name from the project id.
     Read the project id from the qmcpack input file.
  """

  if fname.endswith('.info.xml'):
    return fname

  info_fname = ''
  try:
    tree = ET.parse(fname)
  except IOError as e:
    print('Assuming xml input file, unable to open:',fname)
    print('  Error ',e)
    return None
  node = tree.find('project')
  if node is not None:
    base = node.attrib.get('id')
    if base:
      info_fname = base + '.info.xml'
  else:
    print("project node in XML file node found")
  path = os.path.dirname(fname)
  return os.path.join(path, info_fname)


if __name__ == '__main__':
  if len(sys.argv) < 1:
    print('Usage: process_perf.py <*.info.xml file>|<xml input file>')
    sys.exit(1)

  fname = sys.argv[1]
  if not os.path.exists(fname):
    print('Input file not found: ', fname)
    sys.exit(1)

  info_fname = get_info_file(fname)
  if not info_fname:
    print('Info file not extracted from: ', fname)
    sys.exit(1)

  if not os.path.exists(info_fname):
    print('Info file does not exist: %s'%info_fname)
    print('  Current directory: %s'%os.getcwd())
    sys.exit(1)
  vals = get_performance_info(info_fname)
  print_for_cdash(vals)
