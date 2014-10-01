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
