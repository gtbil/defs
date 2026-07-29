#!/usr/bin/perl -w

#######################################################################################################
# Copyright Notice and Disclaimer for BatchPrimer3
#
# Copyright (c) 2007, Depatrtment of Plant Sciences, University of California,
# and Genomics and Gene Discorvery Unit, USDA-ARS-WRRC.
#
# Authors: Dr. Frank M. You
#
#
# Redistribution and use in source and binary forms of the BStats script, with or without
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

package BStats;

use strict;
my $nc = 0;   
my $average = 0;
my $variance = 0;
my $xlarge = 0;
my $xsmall = 0;
my $sd = 0;
my $cv = 0;
my @x;



sub new {
    my $type = shift;
    my $class = ref($type) || $type;
    
    my $data = shift;
    @x = @{$data};
    $nc = scalar(@x);
    
    &do_calculation();

    my $self = {
        nc => $nc,
        average => $average,
        variance => $variance,
        xlarge => $xlarge,
        xsmall => $xsmall,
        sd => $sd,
        cv => $cv       
    };
    
    
    bless ($self, $class);
    return $self;
}


sub size {
    my $self = shift;
    return $self->{nc};
}

sub do_calculation {
    return if $nc == 0;
    
    my $sum = 0.0;
    my $sum2 = 0.0;
    $xlarge = $x[0];
    $xsmall = $x[0];
    
    for (my $i = 1; $i < $nc; $i++) {
        if ($xlarge < $x[$i]) {
            $xlarge = $x[$i];
        }
        if ($xsmall > $x[$i]) {
            $xsmall = $x[$i];
        }
        $sum += $x[$i];
        $sum2 += $x[$i] * $x[$i];
    }
       
    $average = $sum / $nc;
    $variance = ($sum2 - $sum * $sum /$nc)/ ($nc - 1);
    $sd = sqrt($variance);
    if ($average != 0) {
        $cv = $sd / $average * 100;
    } else {
        $cv = 0;
    }

    $average = int($average * 100 + 0.5)/100.0;
    $variance = int($variance * 100 + 0.5)/100.0;
    $sd = int($sd * 100 + 0.5)/100.0;
    $cv = int($cv * 100 + 0.5)/100.0;
    

}
   
sub average {
    my $self = shift;
    return $self->{average};
}
    
sub variance {
    my $self = shift;
    return $self->{variance};
}

sub cv {
    my $self = shift;
    return $self->{cv};
}
    
sub max {
    my $self = shift;
    return $self->{xlarge};
}
    
sub min {
    my $self = shift;
    return $self->{xsmall};
}

sub sd {
    my $self = shift;
    return $self->{sd};
}


return 1;

