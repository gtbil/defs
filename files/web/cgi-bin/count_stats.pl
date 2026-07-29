#!/usr/local/bin/perl -w

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

use strict;
use warnings;
use CGI;
use Fcntl qw ( :flock);


my $q = CGI->new();
#print $q->start_html();
print $q->header();
# always get back to user after perl code
my ($users, $total_access) = &get_statistics(); 


#print "<script language=\"JavaScript\" type=\"text/javascripti\">";
print "document.write('<b>&nbsp;&nbsp;Users: $users  <br>');";
print "document.write('&nbsp;&nbsp;Times of access: $total_access</b>')";
#print "</script>";
#print $q->end_html();

# end of the main program

sub get_statistics {
    my $log_file = "batchprimer3_log.txt";

	my $count = 0;
    my %hash = ();
 	
    if (-e $log_file) {
        open (FILE, "<$log_file") or die ("Can't open $!");
        flock(FILE, LOCK_EX) or die ("Can't get exclusive lock: $!");
	
        while (my $line = <FILE>) {
           chomp($line);
           next if $line =~ /^s*$/;
           $count++;
           my @arr = split(/\t/, $line);
           $hash{$arr[1]} = "a";
        }
        flock (FILE, LOCK_UN) or die ("Can't  unlock file $log_file");
        close(FILE);
    }
    
    my $users = scalar(keys(%hash));   
    return ($users, $count); 
}
