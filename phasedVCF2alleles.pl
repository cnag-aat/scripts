#!/usr/bin/env perl

use Getopt::Long;
use File::Basename;
use warnings;
use strict;
use Getopt::Long 'HelpMessage';

=head1 NAME

Convert Phased VCF (Whatshap) + GFF + REF FASTA => GENE/TRANSCRIPT/CDS FASTA per haplotype.

=head1 SYNOPSIS

  --vcf,-v        Input VCF (required)
  --ref,-r,-f     Ref genomic FASTA (required)
  --gff,-g        GFF file of a single gene (required)


=head1 VERSION

1.0

=cut

my $name = "sample";
my $outname = "";
GetOptions(
	   'v|vcf:s'  => \ my $vcf,
	   'ref|r:s'    => \ my $ref,
	   'gff|g:s'      => \ my $gff,
	   'name|n:s'      => \ $name,
  	   'help'       => sub { HelpMessage(0) },
	  )or HelpMessage(1);
# die unless we got the mandatory argument
HelpMessage(1) unless ($vcf && $ref && $gff);
die "$gff not found" if !-e $gff;
die "$vcf not found" if !-e $vcf;
die "$ref not found" if !-e $ref;

open GFF,"<$gff" or die "Couldn't open $gff";
my $seqout1 = "";
my $seqout2 = "";
my $strand = "+";
while(my $l = <GFF>){
  if ($l=~m/gene_name=([^;]+)/){$outname = "$name.$1";}
  my @f = split "\t",$l;
  if($f[6] eq "-"){$strand = "-";}
  if($f[2] eq "CDS"){
    my $seq1 = `samtools faidx $ref $f[0]:$f[3]-$f[4] | bcftools consensus -H 1pIu $vcf | FastaToTbl | cut -f 2 -d" "`;
    chomp $seq1; $seqout1.=$seq1;
    my $seq2 = `samtools faidx $ref $f[0]:$f[3]-$f[4] | bcftools consensus -H 2pIu $vcf | FastaToTbl | cut -f 2 -d" "`;
    chomp $seq2; $seqout2.=$seq2;
  
  }
  if ($strand eq "+"){
    open (SEQ1,"|TblToFasta >$outname.hap1.fasta");
    print SEQ1 "$outname.hap1\t$seqout1\n";
    open (SEQ1,"|TblToFasta >$outname.hap2.fasta");
    print SEQ1 "$outname.hap2\t$seqout2\n";
  }else{
    open (SEQ1,"|TblToFasta | revcomp.pl >$outname.hap1.fasta");
    print SEQ1 "$outname.hap1\t$seqout1\n";
    open (SEQ1,"|TblToFasta | revcomp.pl >$outname.hap2.fasta");
    print SEQ1 "$outname.hap2\t$seqout2\n";
  }
}
