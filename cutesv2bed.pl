#!/usr/bin/env perl

while(<>){
	my @F=split;
	if(m/END=(\d+)/){
		my $end = $1; 
		print join("\t",($F[0],$F[1],$end,$F[2])),"\n";
	}
}
