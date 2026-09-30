################################################################################
## This file is distributed under the University of Illinois#NCSA Open Source
## License.  See LICENSE file in top directory for details.
##
## Copyright (c) 2019 QMCPACK developers.
##
## File developed by:
## Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
##
## File created by:
## Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
################################################################################

import os
import sys
import hashlib
import subprocess
import re

from os.path import join, getsize

help_flag = re.compile(r'(-h|--help|-\?)')
for arg in sys.argv:
    if help_flag.match(arg):
        print ("""Formats all C/C++/CUDA files in a directory recursively with clang-format 
until no change is observed.
usage: python full_clang-format.py <directory>
directory is taken to be . if not supplied""")
        exit()
if(sys.argv[1]):
    format_root = sys.argv[1]
else:
    format_root = os.getcwd()

#first time through
formatted_hashes = {}
code_suffixes = re.compile(r'.*\.(cpp|c|cu|h|hpp)$')
for root, dirs, files in os.walk(format_root):
    for name in files:
        if code_suffixes.match(name) != None:
            path_name = join(root,name)
            print(path_name)
            with open(path_name, 'r') as f_for_hash:
                sha256 = hashlib.sha256()
                contents = f_for_hash.read()
                sha256.update(str.encode(contents))
                formatted_hashes[path_name] = sha256.digest()
            subprocess.run(["clang-format", "-i", path_name])
    #this prevents mucking about in the .git directories
    if '.git' in dirs:
        dirs.remove('.git')

#now we have the list of formatted files
#we use this to keep formatting until there are no changes
try:
    #popitem will throw a KeyError when we finally drain the file list
    while True:
        (path_name,hash) = formatted_hashes.popitem()
        print (path_name)
        subprocess.run(["clang-format", "-i", path_name])
        with open(path_name, 'r') as f_for_hash:
            sha256 = hashlib.sha256()
            contents = f_for_hash.read()
            sha256.update(str.encode(contents))
            digest = sha256.digest()
            if hash != digest:
                formatted_hashes[path_name] = digest
except KeyError as ke:
    print(ke)
    print("Done formatting")

                
