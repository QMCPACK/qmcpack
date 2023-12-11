#!/usr/bin/env perl
use strict;
use File::Copy;

my $dir = `pwd`;
chomp($dir);
print "Changing all scripts to use configure script at: $dir/setup-qmc-conf.pl\n"; 
my $newline = 'my %config = do \"' . "$dir/setup-qmc-conf.pl" .'\";';

opendir DIR, ".";
my @files = grep { $_ ne '.' && $_ ne '..'} readdir DIR;
closedir DIR;

foreach my $file (@files) {
    if ($file =~ /\.pl$/ && !($file eq "ConfigureScripts.pl") ) {
	my $sedline = 'sed \'/my %config = do/c ' . $newline . "\' $file";
#	print $sedline . "\n";
	`$sedline > $file.tmp`;
	move ("$file.tmp",  "$file");
        `chmod 755 $file`;
    }

}
			    
