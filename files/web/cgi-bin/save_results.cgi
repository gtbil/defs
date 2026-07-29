#!/usr/bin/env perl

#######################################################################################################
# Copyright Notice and Disclaimer for BatchPrimer3
#
# Copyright (c) 2007, Depatrtment of Plant Sciences, University of California,
# and Genomics and Gene Discorvery Unit, USDA-ARS-WRRC.
#
# Authors: Dr. Frank M. You
#
#
# Redistribution and use in source and binary forms of this script, with or without
# modification,are permitted provided that the following conditions are met:
#
#   1. Redistributions must reproduce the above copyright notice, this list of conditions and the
#      following disclaimer in the documentation and/or other materials provided with the distribution.
#      Redistributions of source code must also reproduce this information in the source code itself.
#   2. If the program is modified, redistributions must include a notice (in the same places as above)
#      indicating that the redistributed program is not identical to the version distributed by Whitehead
#      Institute.
#   3. All advertising materials mentioning features or use of this software must display the following
#      acknowledgment: This product includes software developed by Department of Plant Sciences, UC Davis
#      and Genomics and Gene Discorvery Unit, USDA-ARS-WRRC.
#   4. The name of the UC Davis and USDA-ARS-WRRC may not be used to endorse or promote products derived
#      from this software without specific prior written permission. 
#######################################################################################################

use CGI;
use strict;
use warnings;

my $cgi = new CGI;

my $file_name = $cgi->param("FILE");
my $type = $cgi->param("TYPE");
if (!open (TABLE_FILE, "<$file_name")) {
    print "can\'t create the file $file_name. Please check if your directory is correct. " .
    "Probably you need to change the mode of the directory: chmod 777 ";
}

my $results = "";
foreach (<TABLE_FILE>) {
  $results = $results.$_;
}

my $export_file = 'results_table.xls';
if ($type eq "txt") {
    $export_file = 'results_table.txt';
} elsif ($type eq "excel") {
    $export_file = 'results_table.xls';
}

if ($type eq "txt") {
    # Send the text to a file output over http
    print "Content-type: application/octet-stream\n";
} elsif ($type eq "excel") {
    # Send the text to a file output over http
    print "Content-type: application/excel\n";
}

print "Content-Disposition: attachment; filename=\"".$export_file."\"\n\n";
print $results;
exit;

