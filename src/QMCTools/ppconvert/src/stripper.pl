#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory 
#//
#// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory 
#//////////////////////////////////////////////////////////////////////////////////////



#!/usr/bin/perl -w
use strict;

my $line;
while (my $ln = <STDIN>) {
    $ln =~ s/^\s//;
    $ln =~ s/common\/blitz.*?\.h \\?\n?//g;
    $ln =~ s/common\/blitz.*?\.cc \\?\n?//g;
#    $ln =~ s/\/home.*?einspline.*?\.h \\?\n?//g;
    $ln =~ s/\/home\/lshulen\/include.*?\.h\s?\\?\n?//g;
    $line .= $ln   
}

print "$line\n";
