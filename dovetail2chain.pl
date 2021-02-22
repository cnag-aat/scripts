#!/usr/bin/env perl

#DOVETAIL FORMAT
#tName_dt	qName	qStart	qEnd	orientation	tStart	tEnd

#CHHAIN FORMAT
#chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
#aln_ungapped_block_length

use Getopt::Long;
my $length_file=0;
my $df=0;
my $lookup=0;
GetOptions('lookup|lu=s'=>\$lookup,
	   'dt=s'=>\$dt,
	   'lengths|len=s'=>\$length_file
	  );
#print STDERR "$lookup $length_file $dt\n" and exit();
open(my $lenFh,"<$length_file");

my %len=();
while(<$lenFh>){
  chomp;
  my ($l,$n)=split;
  $len{$n}=$l;
}
close $lenFh;

open(my $luFh,"<$lookup");
my %scaffold=();
while(<$luFh>){
  chomp;
  my ($new,$old)=split;
  $scaffold{$old}=$new;
}
close $luFh;
my $id = 1;
open(my $dtFh,"<$dt");
while(<$dtFh>){
  chomp;
  my ($qName_dt,$tName,$tStart,$tEnd,$orientation,$qStart,$qEnd)=split;
  my $tSize=$tEnd;
  my $qName=$scaffold{$qName_dt};
  my $qSize=$len{$qName};
  print(join(" ",("chain",$tSize,$tName,$tSize,"+",$tStart,$tEnd,$qName,$qSize,$orientation,$qStart,$qEnd,$id++)),"\n");
  print "$tSize\n\n";
}
close $dtFh;
