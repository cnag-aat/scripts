#!/usr/bin/env perl
use strict;
use Data::Dumper;
my $chr=0;
#chr_Pp04	10712648	10780842	213	W	pdulcis7_s0411	1	111761	-
#chr_Pp04	10854592	10864591	218	N	10000	scaffold	yes	peach:Pp04
#chr_Pp04	10864592	10918757	219	W	pdulcis7_s0183	209193	271603	+
my @components=();
my $chr = shift @ARGV;
while(<>){
  next if m/^#/;
  chomp;
  my @f=split;
  next if $f[4]!~m/(N|W)/;
  push @components, \@f;
}

# OK, now print reversed chain with strands switched, too
my @revcomponents = reverse(@components);
my $start=1;
my $end=0;
my $counter=1;
foreach my $comp (@revcomponents) {
  $comp->[0]=$chr if $chr;
  $comp->[3]=$counter++;
  if ($comp->[4] eq 'W'){
    if($comp->[8] eq '+'){
      $comp->[8] = '-';
    }else{
      $comp->[8] = '+';
    }
    $comp->[1]=$end+1;
    my $length = $comp->[7]-$comp->[6]+1;
    $comp->[2]=$end+$length;
  }elsif($comp->[4] eq 'N'){
    $comp->[1]=$end+1;
    my $length = $comp->[5];
    $comp->[2]=$end+$length;
  }
  print join("\t",@{$comp}),"\n";
  $end=$comp->[2];
}
