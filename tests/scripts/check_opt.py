#! /usr/bin/env python3

from __future__ import print_function

from optparse import OptionParser
import xml.etree.ElementTree as ET
import sys
import os


def exit_fail(msg=None):
    if msg!=None:
        print(msg)
    #end if
    print('Test status: fail')
    exit(1)
#end def exit_fail


def exit_pass(msg=None):
    if msg!=None:
        print(msg)
    #end if
    print('Test status: pass')
    exit(0)
#end def exit_pass


# Open the XML file and return coefficient values
def get_opt_coeff(file):
    tree = ET.parse(file)
    root = tree.getroot()
    j_coeff = []
    bf_coeff = []
    for elem in root.findall('wavefunction/jastrow/correlation/'):
        for i in elem.text.split():
         j_coeff.append(float(i))

    for elem in root.findall('wavefunction/determinantset/backflow/transformation/correlation/'):
        for i in elem.text.split():
         bf_coeff.append(float(i))

    return j_coeff,bf_coeff
#end def get_opt_coeff


passfail = {True:'pass',False:'fail'}
def run_opt_test(options):


    prefix_file = options.prefix+'.s'+str(options.series).zfill(3)+'.opt.xml'    

    # Set default for the reference
    ref_file = './qmc-ref/{0}'.format(prefix_file) if options.ref is None else options.ref

    if not os.path.exists(prefix_file):
       exit_fail("Test not found:" + prefix_file)       

    if not os.path.exists(ref_file):
       exit_fail("Reference not found:" + ref_file)

    j_output, bf_output = get_opt_coeff(prefix_file)
    j_reference, bf_reference = get_opt_coeff(ref_file)

    if len(j_output) != len(j_reference):
       exit_fail('Number of coefficient in test({0}) does not match with the reference({1})'.format(len(j_output),len(j_reference)))

    success = True
    j_tolerance = 1e-06
    bf_tolerance = 1e-05
    j_deviation = []
    bf_deviation = []
   
    for i in range(len(j_output)):
       j_deviation.append(abs(float(j_output[i])-float(j_reference[i])))
       quant_success = j_deviation[i] <= j_tolerance 
       if quant_success is False:
           success &= quant_success
       #end if
    #end for
  
    msg='\n  Testing Series: {0}\n'.format(options.series)
    msg+='   reference Jastrow coefficients   : {0}\n'.format(j_reference)
    msg+='   computed  Jastrow coefficients   : {0}\n'.format(j_output)
    msg+='   pass tolerance           : {0: 12.6f}\n'.format(j_tolerance) 
    msg+='   deviation from reference : {0}\n'.format(j_deviation)
    msg+='   status of this test      :   {0}\n'.format(passfail[success])

    if bf_output or bf_reference:
      if len(bf_output) != len(bf_reference):
        exit_fail('Number of coefficient in test({0}) does not match with the reference({1})'.format(len(bf_output),len(bf_reference)))
      
      for i in range(len(bf_output)):
       bf_deviation.append(abs(float(bf_output[i])-float(bf_reference[i])))
       quant_success = bf_deviation[i] <= bf_tolerance
       if quant_success is False:
           success &= quant_success
       #end if
      #end for 

      msg+='\n  Testing Series: {0}\n'.format(options.series)
      msg+='   reference Backflow coefficients   : {0}\n'.format(bf_reference)
      msg+='   computed  Backflow coefficients   : {0}\n'.format(bf_output)
      msg+='   pass tolerance           : {0: 12.6f}\n'.format(bf_tolerance)
      msg+='   deviation from reference : {0}\n'.format(bf_deviation)
      msg+='   status of this test      :   {0}\n'.format(passfail[success])


    return success, msg

#end def run_opt_test


if __name__ == '__main__':

    parser = OptionParser(
       usage='usage: %prog [options]',
       add_help_option=False,
       )
    parser.add_option('-h','--help',dest='help',
                     action='store_true',default=False,
                     help='Print help information and exit (default=%default).'
                     )
    parser.add_option('-p','--prefix',dest='prefix',
                     default='qmc',
                     help='Prefix for output files (default=%default).'
                     )
    parser.add_option('-s','--series',dest='series',
                     default='0',
                     help='Output series to analyze (default=%default).'
                     )
    parser.add_option('-r','--ref',dest='ref',
                     help='Reference to check output files (default=./qmc-ref/$PREFIX).'
                     )
  
    options,files_in = parser.parse_args()

    if options.help:
     print('\n'+parser.format_help().strip())
     exit()
   #end if

    success,msg = run_opt_test(options)

    if success:
       exit_pass(msg)
    else:
       exit_fail(msg)

#end if
