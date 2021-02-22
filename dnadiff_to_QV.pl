#!/usr/bin/env perl 
use strict;
my $label = shift(@ARGV);
my $totalbases = 0;
my $indels = 0;
my $snps = 0;
while(my $line = <>){ #read in the tab-delimited text file of staffing hours on standard input or giving file as an argument
	if($line =~ /^AlignedBases\s+(\d+)/){
		$totalbases = $1;
	}
	if($line =~ /^TotalIndels\s+(\d+)/){
		$indels = $1;
	}
	if($line =~ /^TotalSNPs\s+(\d+)/){
		$snps = $1;
	}
}
my $qv = -10 * (log10(($snps + $indels)/$totalbases));
print "$qv\t$label\n";

sub log10 {
    my $n = shift;
    return log($n)/log(10);
  }