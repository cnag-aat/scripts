#!/usr/bin/env perl
use strict;
use lib "/home/devel/talioto/myperlmods";
use lib "/home/devel/talioto/myperlmods/Bio";
use SeqOp;
use Getopt::Long;
use File::Basename qw( fileparse );
use Bio::DB::Fasta;
use Bio::Seq;
use Bio::SeqIO; 
# Arguments: FASTA, AGP

my $agp_file = 0;
my $fasta = 0;
my $nums=0;
my %chr;
my $sp = 0;
my $ssp = 0;
GetOptions(
	   'fasta|f:s'        => \$fasta,
	   'agp|a:s'         =>\$agp_file,
	  );
die "Usage: $0 -fasta <FASTA> -agp <AGP FILE>" if !($agp_file && $fasta);
if (!-e "$fasta.index") {
  `indexGenome.pl -f $fasta`;
}
my $db = Bio::DB::Fasta->new($fasta);

my $nums=0;
my @ids = $db->ids;
$nums=scalar @ids;
$sp = length($nums);
$ssp = $sp;
my $seq_out = Bio::SeqIO->new('-fh' => \*STDOUT,'-format' => 'fasta');
my $lastid=0;
my $last_seq='';


open AGP, "<$agp_file";
my $base = "SALSA";
my $ssid = 0;
my $sscount = 1;
my $scount = 1;
while (<AGP>) {
  chomp;
  next if m/^#/;
  next if m/^\s/;
  my @F = split /\s+/;

  #$lastid=$F[0] unless $lastid;
  if ($F[0] ne $lastid) {
    $ssid = sprintf("%s"."_ss%0$ssp"."i",$base,$sscount++);
    $lastid=$F[0];
    $last_seq='';
    $scount = 1;
  }

  if ($F[4] !~ m/(N|U)/i) {
    my ($start,$stop) = $F[8] ne '-'?($F[6], $F[7]):($F[7], $F[6]);
    my $seqid = sprintf("%ss%0$sp"."d",$ssid,$scount++);
    #print STDERR  "$lastid\n";	
    print_seq($seqid,$db->seq($F[5],$start,$stop));	
    print STDERR (join "\t",($F[0],$F[5],$seqid,$F[6],$F[7],$F[8])),"\n";
  }
} 
close AGP;

sub print_seq{
  my($id,$seq)=@_;	
  $seq=~s/\s+//g;
  my $seqobj = Bio::Seq->new( -display_id => "$id", -seq => $seq);
  $seq_out->write_seq($seqobj);
}

  sub lensort
  {
    $b->{len} <=> $a->{len}
  }

