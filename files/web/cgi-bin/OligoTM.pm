#!/usr/bin/perl -w

#######################################################################################################
# Copyright Notice and Disclaimer for BatchPrimer3
#
# Copyright (c) 2007, Depatrtment of Plant Sciences, University of California,
# and Genomics and Gene Discorvery Unit, USDA-ARS-WRRC.
#
# Authors: Dr. Frank M. You
#
# This program is implemented based on OligoTM.c program from Primer3 package. 
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

package OligoTM;

use strict;

our $MAX_NN_TM_LENGTH = 36;
our $OLIGOTM_ERROR = -999999;
our %S;
our %H;
our %G;

$S{"AA"} = 240;
$S{"AC"} = 173;
$S{"AG"} = 208;
$S{"AT"} = 239;
$S{"AN"} = 215;
$S{"CA"} = 129;
$S{"CC"} = 266;
$S{"CG"} = 278;
$S{"CT"} = 208;
$S{"CN"} = 220;
$S{"GA"} = 135;
$S{"GC"} = 267;
$S{"GG"} = 266;
$S{"GT"} = 173;
$S{"GN"} = 210;
$S{"TA"} = 169;
$S{"TC"} = 135;
$S{"TG"} = 129;
$S{"TT"} = 240;
$S{"TN"} = 168;
$S{"NA"} = 168;
$S{"NC"} = 210;
$S{"NG"} = 220;
$S{"NT"} = 215;
$S{"NN"} = 203;

$H{"AA"} =  91;
$H{"AC"} =  65;
$H{"AG"} =  78;
$H{"AT"} =  86;
$H{"AN"} =  80;
$H{"CA"} =  58;
$H{"CC"} =  110;
$H{"CG"} =  119;
$H{"CT"} =  78;
$H{"CN"} =  91;
$H{"GA"} =  56;
$H{"GC"} =  111;
$H{"GG"} =  110;
$H{"GT"} =  65;
$H{"GN"} =  85;
$H{"TA"} =  60;
$H{"TC"} =  56;
$H{"TG"} =  58;
$H{"TT"} =  91;
$H{"TN"} =  66;
$H{"NA"} =  66;
$H{"NC"} =  85;
$H{"NG"} =  91;
$H{"NT"} =  80;
$H{"NN"} =  80;

$G{"AA"} =  1900;
$G{"AC"} =  1300;
$G{"AG"} =  1600;
$G{"AT"} =  1500;
$G{"AN"} =  1575;
$G{"CA"} =  1900;
$G{"CC"} =  3100;
$G{"CG"} =  3600;
$G{"CT"} =  1600;
$G{"CN"} =  2550;
$G{"GA"} =  1600;
$G{"GC"} =  3100;
$G{"GG"} =  3100;
$G{"GT"} =  1300;
$G{"GN"} =  2275;
$G{"TA"} =  900;
$G{"TC"} =  1600;
$G{"TG"} =  1900;
$G{"TT"} =  1900;
$G{"TN"} =  1575;
$G{"NA"} =  1575;
$G{"NC"} =  2275;
$G{"NG"} =  2550;
$G{"NT"} =  1575;
$G{"NN"} =  1994;


sub new {
    my $type = shift;
    my $class = ref($type) || $type;
    
    my $self = {
        tm => 0,
        oligo_dg => 0,
        end_oligo_dg => 0
    };
    
    bless ($self, $class);
    return $self;
}

sub oligo_tm {
    my $self = shift;
    my $seq = shift;
    my $dna_conc = shift;  
    my $salt_conc = shift;  
    my $tm = tm ($seq, $dna_conc, $salt_conc);
    $self->{tm} = $tm;
    return $self->{tm};
}

# The melting temperature is calculated using the formula based
# on the nearest neighbour thermodynamic theory
sub tm {
    my $s = shift;         # sequence
    my $dna_conc = shift;  # "c" is the nucleic acid molar concentration
    		    	   # (determined empirically, W.Rychlik et.al.1990),
                           # (default value is 0.2 µM for Allawi and SantaLucia thermodynamic parameters),
    my $salt_conc = shift;  #salt molar concentration (default value is 50 mM).
    
    my $i = 0;
    my $j = 108;
    my $dh = 0.0;
    my $ds = 0.0;

    if (length($s) < 2) {
    	return $OLIGOTM_ERROR;
    }
    
    # convert the degeneration code
    $s = degenerateIUBcode($s);
 
    for( my $k = 0; $k < length($s) - 1; $k++) {
 	my $s1 = uc(substr($s, $k, 2));
    	$i += $H{$s1};
        $j += $S{$s1};
    }
    $dh = $i * (-100);
    $ds = $j * (-0.1);
    return ($dh / ($ds + 1.987 * log($dna_conc /4000000000.0)) - 273.15) +
    	16.6 * (log($salt_conc / 1000)  / log(10));
}

#Calculate the primer stability (Delta G)
sub oligo_dg {
    my $self = shift;
    my $seq = shift;
   
    my $i = 0;
    if(length($seq)< 2) {
        return $OLIGOTM_ERROR;
    }
    my $len = length($seq);
    for(my $j = 0; $j < $len - 1; $j++) {
        my $s1 = uc(substr($seq, $j, 2));
        $i += $G{$s1};
    }

    my $oligo_dg = $i / 1000.0;

    $self->{oligo_dg} = $oligo_dg;
    return $self->{oligo_dg};
}

#Calculate the 3' end stability (Delta G)
sub end_oligo_dg {
    my ($self, $seq, $end_len) = @_;

    my $len = length($seq);
    my $end_oligo_dg = $len >= $end_len ? oligo_dg(substr($seq, $len - $end_len)) : oligo_dg($seq);
    
    $self->{end_oligo_dg} = $end_oligo_dg;
    return $self->{end_oligo_dg};
    
}

sub degenerateIUBcode {
    my $s = shift;
    $s =~ tr/MRWSYK/AAACCG/;  # here use only one nucleotide to represent the alleles. It's not good but no better solution
    return $s;
}


return 1;

