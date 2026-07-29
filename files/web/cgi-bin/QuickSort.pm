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

package QuickSort;

use strict;

my @dataArray;  #array reference
my @indexArray = ();

sub new {
    my $type = shift;
    my $class = ref($type) || $type;
    
    my $self = {};
    
    bless ($self, $class);
    return $self;
}


    
#create a QuickSort object and sort the user's iinteger array and return sorted integer object
#array and index array with original order

sub sort {
    my $self = shift;
    my $array = shift;
    
    #kepp original order in the index Array
    for (my $i = 0; $i < @$array; $i++) {
        $indexArray[$i] = $i;
        $dataArray[$i] = @$array[$i]; 
    }
    &quicksort(0, scalar(@dataArray)-1);
    return (\@dataArray, \@indexArray);
} 

# recursive quicksort that breaks array up into sub
# arrays and sorts each of them.
sub quicksort {
    my ($p, $r) = @_;
    if ($p < $r) {
        my $q = &partition($p,$r);
        $q-- if ( $q == $r);
        &quicksort($p, $q);
        &quicksort($q+1, $r);
    } 
} 

# Partition by splitting this chunk to sort in two and
# get all big elements on one side of the pivot and all
# the small elements on the other.
sub partition {
    my ($lo, $hi) = @_;
    my $pivot = $dataArray[$lo];
    while (1) {
        while ( $dataArray[$hi] >= $pivot  && $lo < $hi) {
            $hi--;
        }
        while ( $dataArray[$lo] < $pivot  && $lo < $hi ) {
            $lo++;
        }
        
        if ( $lo < $hi ) {
            # exchange objects on either side of the pivot
            my $element = $dataArray[$lo];
            my $index = $indexArray[$lo];
            $dataArray[$lo] = $dataArray[$hi];
            $indexArray[$lo] = $indexArray[$hi];
            $dataArray[$hi] = $element;
            $indexArray[$hi] = $index;
        }
        else {
            return $hi;
        } 
    }
}
    

# check if user's array is already sorted
sub isAlreadySorted() {
    for (my $i = 1; $i < @dataArray; $i++ ) {
        if ($dataArray[$i] < $dataArray[$i-1]) {
            return 0;
        }
    }
    return 1;
} 

return 1;

