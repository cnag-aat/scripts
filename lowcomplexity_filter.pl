#!/usr/bin/env perl
 
use Pod::Usage;
use Getopt::Long;
use File::Basename;
 
$this_program = basename($0);
my $in = 0;
my $thresh = 0.5;
my $result = GetOptions(
			'f|fasta=s' => \$in,
			't|threshold=f'   => \$thresh
		       );
die "Must supply fasta for filtering using -f <in.fa>\n" if !$in;
$base = fileparse($in, qw(\.fasta \.fa));
my %lowcomp_len;
print STDERR "Running sdust $in...\n";
open SD, "sdust $in |";
while(<SD>){
  chomp;
  my ($name,$start,$end)=split;
  $lowcomp_len{$name}+=($end-$start);
}
close SD;
print STDERR "Filtering $in according to sdust results...\nOK reads are written to $base.hicomp.fasta\nLow-complexity reads in $base.lowcomp.fasta...\n";
open SEQ,"sed 's/\s+/====/g' $in | FastaToTbl | ";
open OK, "| TblToFasta |sed 's/====/ /g' > $base.hicomp.fasta";
open LOW, "| TblToFasta |sed 's/====/ /g' > $base.lowcomp.fasta";
while (<SEQ>) {
  chomp;
  my($seqname,$seq)=split;
  if (exists($lowcomp_len{$seqname}) && ($lowcomp_len{$seqname}/length($seq) > $thresh)){
    print LOW "$_\n";
  }else{
    print OK "$_\n";
  }
}
close OK;
close LOW;
close SEQ;
