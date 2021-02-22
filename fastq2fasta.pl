#!/usr/bin/env perl
use strict;
while(<>){
    if (m/^@(.*)$/){
	my $l1 = ">$1\n";
	my $l2 = <>;
	my $dummy = <>;
	if ($dummy =~/^\+/){print $l1,$l2;}
	$dummy = <>;
    }
}
