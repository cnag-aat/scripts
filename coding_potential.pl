#!/usr/bin/env perl
use strict;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;

# Need to add GFF2 FITS (geneid-like) output
my $fa= 0;
my $verbose = 0;
my $geneid_param = "/home/groups/assembly/talioto/repos/geneid/param/human.101007.scoring.param";
my $frame = 1;
my $strand = "+";
my $start = 1;
my $end = $start + 3;
my $seqname = '';

GetOptions(
'fasta:s'    => \$fa,
'seq:s'      => \$seqname,
'param:s'    => \$geneid_param,
'frame|f:s'  => \$frame,
'strand:s'   => \$strand,
'start:s'    => \$start,
'end:s'      => \$end
);
if ($strand =~/(-|r|R|minus)/){$strand = "-";}else{$strand = "+";}
# open the first file: table of codon usage frequencies

open(PARAM,"<$geneid_param") or die "$0: the file $geneid_param can not be opened: $!\n";
if($frame!~/[123]/){die "-frame must be 1, 2, or 3";}
# load the frequencies of codon usage into the hash table %pcodons
# this hash is indexed by a triplet of nucleotides or codon

#print STDERR "Frame: $frame\n";
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
my $chrlen =$db->length($seqname);
my $rem = $chrlen % 3;
my $sobj = $db->get_Seq_by_id($seqname);
my $seq = $sobj->seq;

if($strand eq '+'){
    my $dna = substr($seq,$start-2+$frame,$end-$start+1);
    print(">$seqname:",$start-1+$frame,"-",$end," cp:",compute_cp($dna,0),"\n$dna\n");
}else{
    my $dna = substr($seq,$start-1,$end-$start+2-$frame);
    my $rdna = reverse($dna);
    $rdna =~ tr/ACGTacgt/TGCATGCA/;
    print(">$seqname:",$end-$frame,"-",$start," cp:",compute_cp($rdna,0),"\n$rdna\n");
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
        if($f == 3){$f=0;}
    }
    return $score;
}

