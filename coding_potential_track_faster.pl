#!/usr/bin/env perl
use strict;
use Getopt::Long;
use lib "/project/devel/aateam/perlmods";
use lib "/apps/BIOPERL/1.6.1/lib/perl5";
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;
use File::Basename qw( fileparse );
use SeqOp;

# Need to add GFF2 FITS (geneid-like) output
my $fa= 0;
my $verbose = 0;
my $geneid_param = "/home/devel/talioto/param/human.090903.param";
my $frame = 1;
my $strand = "+";
my $window = 60;
GetOptions(
           'seq:s'          => \$fa,
	   'param:s'        => \$geneid_param,
	   'frame|f:s'       => \$frame,
	   'strand:s'       => \$strand,
	   'w:s'            => \$window
	   );
if ($strand =~/(-|r|R|minus)/){$strand = "-";}else{$strand = "+";}
# open the first file: table of codon usage frequencies

open(PARAM,"<$geneid_param") or die "$0: the file $geneid_param can not be opened: $!\n";
if($frame!~/[123]/){die "-frame must be 1, 2, or 3";}
# load the frequencies of codon usage into the hash table %pcodons
# this hash is indexed by a triplet of nucleotides or codon

my %llhex;
while (my $line = <PARAM>) {
  next if $line!~/Markov_Transition_probability_matrix/;
  while(my $entry = <PARAM>){
    last if $entry!~/\w/;
    last if $entry=~/^#/;
    chomp $entry;
    my($hex,$i,$frame,$lls)=split /\s+/,$entry;
    $llhex{$hex}{$frame}=$lls;
  }     
}
close(PARAM);

my $db = Bio::DB::Fasta->new($fa);

my @seqids = sort $db->ids;
my $start=int($window/2)+($frame-1)+1;

foreach my $id (@seqids){
  my $chrlen =$db->length($id);
  my $rem = $chrlen % 3;
  my $sobj = $db->get_Seq_by_id($id);
  my $seq = $sobj->seq;
  print "fixedStep  chrom=$id  start=$start  step=3\n";
  if($strand eq '+'){
    for (my $base = $frame;$base<($chrlen-$window+1);$base+=3){
      #my $dna = SeqOp::get_seq_BioDBFasta($db,$id,$base,$base+$window-1,'+');
      #my $dna = $sobj->subseq($base,$base+$window-1);
      my $dna = substr($seq,$base-1,$window);
      print compute_cp($dna,0),"\n";
    }
  }else{
    for (my $base = $frame;$base<($chrlen-$window+1);$base+=3){
      #my $dna = SeqOp::get_seq_BioDBFasta($db,$id,$base,$base+$window-1,'+');
      #my $dna = $sobj->subseq($base,$base+$window-1);
      #my $rdna = SeqOp::get_seq_BioDBFasta($db,$id,$base,$base+$window-1,'-');
      my $dna = substr($seq,$base-1,$window);
      my $rdna = reverse($dna);
      $rdna =~ tr/ACGTacgt/TGCATGCA/;
      print compute_cp($rdna,0),"\n";
    }
  }
  
}
sub compute_cp {
  my $seq   = shift;
  my $ucseq = uc($seq);
  #print STDERR $ucseq,"\n";
  my $frame = shift;
  my $f = $frame;
  my $score = 0;
  for(my $p=0;$p<(length($ucseq)-5);$p++){
    my $h=substr($ucseq,$p,6);
    #print STDERR "$h\t$f\t$llhex{$h}{$f}\n";
    $score+=$llhex{$h}{$f};
    $f++;
    if($f==3){$f=0;}
  }
  return $score;
}

