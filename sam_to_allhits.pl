#!/usr/bin/env perl
 
use Pod::Usage;
use Getopt::Long;
use File::Basename;
 
$this_program = basename($0);
my $minmatch = 150;
$minqual = 0;
my $result = GetOptions(
	'i|input=s' => \$input,
	'o|output=s' => \$output,
	'b|buffer=i' => \$printlines,
	'm|minq=i'   => \$minqual,
	'u|uniq=s'   => \$uniq,
	'reads|r=s' => \$readfile,
	'm|minmatch=i' => \$minmatch
);

my %qseqlen;
open (Q,"fastalength $readfile |");
while(<Q>){
  chomp;
  my ($l,$r) =split;
  $qseqlen{$r}=$l;
}
close Q;
my %refseqlen;
while (<>) {
  if (m/SQ.*SN:(\S+).*LN:(\d+)/){
    $refseqlen{$1}=$2;
  }
	next if ($_=~/^\@/);
	$t++;
	print STDERR " $t lines parsed [$chk blocks].\r" unless ($t % 50000);
	my @samfields = split /\t/;
	my $strand = strand($samfields[1]);
	my $contig = $samfields[2];
	next if $contig eq "*";
	my $read   = $samfields[0];
	my $from   = $samfields[3];
	my $mapq = $samfields[4];
	next if $mapq < $minqual;
	my $refalnlen    = 0;
	my $matches = 0;
	my $insert = 0;
	my $del = 0;
	my $soft = 0;
	my $hard = 0;
	my $intron = 0;
	if ($samfields[4] < $minq) {
		$lowqual++;
	}
	while ($samfields[5]=~/(\d+)M/g) { $matches += $1; }
	while ($samfields[5]=~/(\d+)I/g) { $insert += $1; }
	while ($samfields[5]=~/(\d+)D/g) { $del += $1; }
	while ($samfields[5]=~/(\d+)S/g) { $soft += $1; }
	while ($samfields[5]=~/(\d+)H/g) { $hard += $1; }
	while ($samfields[5]=~/(\d+)N/g) { $intron += $1; }
	$refalnlen = $matches + $del + $intron;
	$to = $from + $refalnlen -1;
	my $qalnlen = $matches + $insert;
	my $readlen = $matches + $insert + $soft +$hard;
	# M+I+S=ALNQLEN
	# M
	my $lefthard = 0;
	my $righthard = 0;
	my $leftsoft = 0;
	my $rightsoft = 0;
	my $qstart = 0;
	my $qend = $readlen;
	if ($samfields[5]=~s/^(\d+)H//){$lefthard=$1;}
	if ($samfields[5]=~/^(\d+)S/){$leftsoft=$1;}
	if ($samfields[5]=~s/(\d+)H$//){$righthard=$1;}
	if ($samfields[5]=~/(\d+)S$/){$rightsoft=$1;}
	if ($strand eq "+"){
	  $qstart = 1 + $lefthard +$leftsoft;
	  $qend = $readlen - $righthard -$rightsoft;
	}else{
	  $qstart = 1 + $righthard + $rightsoft;
	  $qend = $readlen - $lefthard - $leftsoft;
	}
  next if $refalnlen < $minmatch;
	print("$read\t$qseqlen{$read}\t$qstart\t$qend\t+\t$contig\t$refseqlen{$contig}\t$from\t$to\t$strand\t$mapq\t");printf("%0.2f\n",(1-(($insert+$del)/(($qalnlen+$refalnlen)/2))));
}

sub strand{
  my $flag = shift;
  my $strand = '+';
  if ($flag & 0x10) {
	  $strand = '-';
  }
  return $strand;		
}
