#!/usr/bin/env spack python

import spack.config
# spack uses this as well so it must be there
import argparse

parser = argparse.ArgumentParser(description='For chosen spec of gcc add gfortran as fortran compiler for all clang compiler configs')
parser.add_argument('gcc_spec', type=string, help='spack spec to get gfortran from')
args = parser.parse_args()

conf_scope = spack.config.ConfigScope('user/linux', os.path.join(os.environ['HOME'],'.spack','linux'))
compilers_config = conf_scope.get_section_filename('compilers')

gcc_configs = [ c for c in conf_scope.get_section('compilers')['compilers'] if args.gcc_spec) in c['compiler']['spec'] ]
if len(gcc_configs) != 1:
    raise NameError('gcc version must evaluate to a single gcc spec')
gcc = gcc_configs[0]

llvm_configs = [ c for c in conf_scope.get_section('compilers')['compilers'] if "clang" in c['compiler']['spec'] ]
for lconf in llvm_configs:
    lconf['compiler']['paths']['f77'] = gcc['compiler']['paths']['f77']
    lconf['compiler']['paths']['fc'] = gcc['compiler']['paths']['fc']
    lconf['compiler']['modules'].append(args.gcc_spec)

conf_scope.write_section('compilers')
