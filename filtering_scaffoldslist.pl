#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

my ($gfffile, $ids_list,$add, $discard);

GetOptions(
           'gff:s'           => \$gfffile, 
           'l:s'             => \$ids_list,  
           'a'               => \$add,   #if this option is given the lines with the id are going to be kept
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
if ($add) {
  while (<GFF>) {
    chomp $_;
    next if ($_ =~ m/^$/);
    my @line = split /\t/, $_;  
    print "$_\n" if (exists $ids{$line[0]});
  }
}

elsif ($discard){
  while (<GFF>) {
    chomp $_;
    next if ($_ =~ m/^$/);
    my @line = split /\t/, $_;  
    print "$_\n" if (!exists $ids{$line[0]});
  }
}

close GFF;
