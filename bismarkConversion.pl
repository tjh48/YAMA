#!/usr/bin/env/perl

#bismarkConversion 0.0.1 - converter for bismark_methylation_extractor outputs to a format liked by the R package segmentSeq

# Copyright (C) 2009 Thomas J. Hardcastle <tjh48@cam.ac.uk>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  You should have received a copy of 
#    the GNU General Public License along with this program. If not, see 
#    <http://www.gnu.org/licenses/>.


use List::Util qw(min);
use File::Path 'rmtree';

use strict ;
use warnings;
use Getopt::Long;
$| = 1;

my $logfile;
my $logFH;
my $splittingLength = 10000000;

sub printlog {
    my $message = pop(@_);
    print $message;
#    open(my $logFH, ">>", "$logfile") or die "cannot open > $logfile: $!";
    print {$logFH} $message;
#    close($logFH);
}

sub dielog {
    my $message = pop(@_);
#    open(my $logFH, ">>", "$logfile") or die "cannot open > $logfile: $!";
    print {$logFH} $message;
#    close($logFH);
    die $message;
}



my $result_dir ;
my $result_name;
my $localSortFlag;
my $sort_parallel = 1;
my $sort_memory="50%";
my $bismarkName;
my $data_dir;

# Get ARGS
GetOptions(
    'bismarkName:s'    => \$bismarkName,
    'd:s'    => \$data_dir,
    'o:s'    => \$result_dir,
    'n:s'    => \$result_name,
    'p:i'    => \$sort_parallel,
    'S:s'    => \$sort_memory
) ;

$result_dir = "." unless $result_dir;
$data_dir = "." unless $data_dir;
die "\n No input file specified.\n" unless $bismarkName;

$result_name = "$bismarkName.segSeq" unless $result_name;
$result_name =~ s/.*\///;

my @bismarkOutData = (["CpG", "z"], ["CHG", "x"], ["CHH", "h"]);

foreach my $bisrow (0..@bismarkOutData-1) {
    my $peacepipe = "sed \'1,1d\' $data_dir/$bismarkOutData[$bisrow][0]_$bismarkName | cut -f 2- | sort -k2,2n -k3,3n -k4,4 | uniq -c | 
awk 'BEGIN {strand = \"\"; chr = \"\"; base = 0; umeth = 0; meth = 0}
{if(\$3 == chr && \$4 == base) {if(\$5 == \"$bismarkOutData[$bisrow][1]\") {umeth += \$1} else {meth += \$1}} else {if(chr != \"\") {print chr \"\t\" base \"\t*\t\" meth \"\t\" umeth}; strand = \$2; chr = \$3; base = \$4; if(\$5 == \"$bismarkOutData[$bisrow][1]\") {umeth = \$1; meth = 0} else {meth = \$1; umeth = 0}}} END {print chr \"\t\" base \"\t*\t\" meth \"\t\" umeth}' > $result_dir/$bismarkOutData[$bisrow][0]_$result_name\n";
system("$peacepipe");
}




