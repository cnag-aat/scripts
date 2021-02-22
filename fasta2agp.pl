#!/usr/bin/env perl

use strict;
my $c=0;
open IN, "fastalength $ARGV[0]  | ";
while(<IN>){
chomp;
my @F=split;
	print(join("\t",($F[1],1,$F[0],++$c,"W",$F[1],1,$F[0],"+")),"\n");
}	
close IN;
