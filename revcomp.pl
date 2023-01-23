#!/usr/bin/env perl
open (IN,"FastaToTbl $ARGV[0] |");
open OUT, "|TblToFasta";
while(<IN>){
	chomp;
    my ($id,$seq) = split;
	print OUT "$id\t";
    	$seq=~tr/acgtACGT/tgcaTGCA/;
    	my $r=reverse($seq);
	print OUT "$r\n";
    
}
close IN;
close OUT;
