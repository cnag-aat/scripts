#!/usr/bin/env perl
use strict;
#use warnings;
use Getopt::Long;

my ($gfffile, $ids_list,$add, $discard);

GetOptions(
           'gff:s'           => \$gfffile, 
           'l:s'                => \$ids_list,  
           'a'              => \$add,   #if this option is given the lines with the id are going to be kept
           'b'               => \$discard, #If this option is given, the lines with the id are goin to be discarded
           );
           
open IDS,"< $ids_list" || die "cannot open $ids_list";
my %ids;
while(<IDS>){
  chomp;
  $ids{$_}++;
}
close IDS;

open GFF, "<", "$gfffile" || die "Cannot open $gfffile";

my @F8;
my $gene;
if ($add) {

   while (<GFF>) {
      chomp $_;
      my @line = split /\s+/, $_;  
      if ($line[2] eq 'gene') {
        if ($line[8] =~ m/ID=([^;]+)/) {$gene = $1;}
        elsif ($line[8] =~ m/ID=/) {$gene = $';}
     #    @F8 = split /=/, $line[8];
     #    $gene = $F8[1];
    #  print "$gene\n";
        print "$_\n" if defined($ids{"$gene"});
      }
      elsif ($line[2] =~ m/transcript|mRNA/) {
        if ($line[8] =~ m/Parent=([^;]+)/) {$gene = $1;}
        print "$_\n" if defined ($ids{"$gene"});
      }
      else {
        print "$_\n" if defined($ids{"$gene"}); 
      }
   }
}

if ($discard) {

   while (<GFF>) {
      chomp $_;
      my @line = split /\s+/, $_;  
      if ($line[2] eq 'gene') {
        if ($line[8] =~ m/ID=([^;]+)/) {$gene = $1;}
        elsif ($line[8] =~ m/ID=/) {$gene = $';}
        print "$_\n" unless exists($ids{$gene});
      } 
      elsif ($line[2] =~ m/transcript|mRNA/) {
        if ($line[8] =~ m/Parent=([^;]+)/) {$gene = $1;}
        print "$_\n" unless exists ($ids{$gene});
      }
      else{
        print "$_\n" unless exists  ($ids{$gene});
      }
  }
}

close GFF;
