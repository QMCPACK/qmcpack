#! /usr/bin/env python3
from __future__ import print_function

# Rewrite parts of a QMCPACK input file
#   walkers
#   project id

import argparse
import xml.etree.ElementTree as ET

# utilities
def add_or_change_parameter(parent_node, param_name, new_value):
  nodes = parent_node.findall("./parameter[@name='%s']"%param_name)
  if nodes:
    nodes[0].text = new_value
  else:
    n = ET.SubElement(parent_node, "parameter", {"name":param_name})
    n.text = new_value

def change_parameter(parent_node, param_name, new_value):
  nodes = parent_node.findall("./parameter[@name='%s']"%param_name)
  if nodes:
    nodes[0].text = new_value

def add_or_change_attribute(node, attrib_name, new_value):
  node.attrib[attrib_name] = new_value


# utilities
def change_number_of_walkers(tree, new_number_of_walkers):
  qmc_nodes = tree.findall(".//qmc")
  for qmc in qmc_nodes:
    change_parameter(qmc, 'walkers', str(new_number_of_walkers))

def use_NLPP_algorithm(tree):
    nodes = tree.findall(".//hamiltonian/pairpot[@type='pseudo']")
    if len(nodes) == 0:
        print('NLPP not found')
        return
    for node in nodes:
      add_or_change_attribute(node, 'algorithm', 'batched')

def use_delayed_update(tree, delay):
    nodes = tree.findall(".//wavefunction/determinantset/slaterdeterminant")
    if len(nodes) == 0:
        print('slaterdeterminant not found')
        return
    for node in nodes:
      add_or_change_attribute(node, 'delay_rank', str(delay))

def change_jastrow(tree, j3_tree):
    wf_nodes = tree.findall('.//wavefunction')
    if len(wf_nodes) != 1:
        print('No wavefunction nodes found, or more than one')
        return
    wf_node = wf_nodes[0]
    jastrow_nodes = wf_node.findall("./jastrow")

    for j_node in jastrow_nodes:
        wf_node.remove(j_node)

    j3_nodes = j3_tree.findall('.//jastrow')
    for j3_node in j3_nodes:
        wf_node.append(j3_node)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Adjust QMCPACK input files")
  parser.add_argument('input_file',
                      help="Input XML file")
  parser.add_argument('-b', '--batch', action='store_true',
                      help="Use the batched algorithm for non-local pseudopotential evaluation")
  parser.add_argument('-d', '--delay',
                      help="Use the delayed update for determinants")
  parser.add_argument('-i', '--inplace', action='store_true',
                      help="Edit files in place")
  parser.add_argument('-j', '--J123',
                      help="Use one, two and three body Jastrow factors")
  parser.add_argument('-o', '--output',
                      help="Ouput XML file")
  parser.add_argument('-w', '--walker',
                      help="Change the number of walkers")

  args = parser.parse_args()
  fname_in = args.input_file

  tree = ET.parse(fname_in)

  if args.delay:
    use_delayed_update(tree, args.delay)

  if args.batch:
    use_NLPP_algorithm(tree)

  if args.J123:
    j3_tree = ET.parse(args.J123)
    change_jastrow(tree, j3_tree)

  if args.walker:
    change_number_of_walkers(tree,args.walker)

  if args.output:
    tree.write(args.output)
  else:
    if args.inplace:
      tree.write(fname_in)
    else:
      print("Nothing has been changed. Use -i or -o options for in-place change or output to a file.")

