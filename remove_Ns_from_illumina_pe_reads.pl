#!/usr/bin/env perl

#Tyler Alioto 2022-05-19
#Fernando Cruz 2022-05-20

# Usage: /home/devel/fcruz/bin/scripts/remove_Ns_from_old_pe_reads.pl read1.fastq.gz read2.fastq.gz


# NOTEs:
# I do not use a regex just a split buy space and keep the first two elements SRA and read id
# Illumina PE format assumed incorporates two space @SRR6071634.1 HWI-D00572:122:C6PG2ANXX:4:1101:1237:1948 length=125 
# It should be modified by keeping just one element in the id to work for cases like this: @A00500:270:H7YGVDSX2:1:1101:8757:1000 1:N:0:GTTCCTCA

use strict;
use File::Basename;


my $read1 = shift;
my $read2 = shift;

open (R1,"zcat -f $read1 |");
open (R2,"zcat -f $read2 |");
my ($base1,$path1,$ext1) = fileparse($read1,qw(\.1\.fastq.gz \.1\.fastq));
my ($base2,$path2,$ext2) = fileparse($read2,qw(\.2\.fastq.gz \.2\.fastq));
open OUT1,"|gzip -c > $base1.Nfree.1.fastq.gz"; #added Nfree to the suffix extension
open OUT2,"|gzip -c > $base2.Nfree.2.fastq.gz"; #added Nfree to the suffix extension
my  $total_pairs_removed=0;
while (my $r1_header = <R1>){
  my $r1_seq = <R1>;
  my $r1_spacer = <R1>;
  my $r1_quality = <R1>;
  my $r2_header = <R2>;
  my $r2_seq = <R2>;
  my $r2_spacer = <R2>;
  my $r2_quality = <R2>;
  
  my @data1=split/\s/, $r1_header;#Read 1
  my $r1_id = $data1[0]." ".$data1[1];
  
  my @data2=split/\s/, $r2_header;#Read 2
  my $r2_id = $data2[0]." ".$data2[1];#Read 2

#  print "screening pair R1 $r1_id R2 $r2_id\n";
  die "$r2_id != $r1_id! Check if files are paired properly!\n" if $r2_id ne $r1_id;
  if ($r1_seq=~/N/ or $r2_seq=~/N/){
    $total_pairs_removed++;
  }else{
    print OUT1 $r1_header;
    print OUT1 $r1_seq;
    print OUT1 $r1_spacer;
    print OUT1 $r1_quality;
    print OUT2 $r2_header;
    print OUT2 $r2_seq;
    print OUT2 $r2_spacer;
    print OUT2 $r2_quality;
  }
}

close R1;
close R2;
close OUT1;
close OUT2;
print STDERR "Removed $total_pairs_removed pairs.\n";
