#!/usr/bin/env perl
use strict;
use Data::Dumper;
my $chr=0;
#chr_Pp04	10712648	10780842	213	W	pdulcis7_s0411	1	111761	-
#chr_Pp04	10854592	10864591	218	N	10000	scaffold	yes	peach:Pp04
#chr_Pp04	10864592	10918757	219	W	pdulcis7_s0183	209193	271603	+
my @components=();
#my $chr = shift @ARGV;
while(<>){
  next if m/^#/;
  chomp;
  my @f=split;
  next if $f[4]!~m/(N|W)/;
  push @components, \@f;
}

# OK, now print reversed chain with strands switched, too
#my @revcomponents = reverse(@components);
my $gap;
my $start=1;
my $end=0;
my $counter=1;
my $prevtype = 0;
my $prevchr = 0;
foreach my $comp (@components) {
  $comp->[3]=$counter++;
  if ($comp->[4] eq 'W'){
    #check if same chr and no previous gap, insert 10000bp gap
    if (($comp->[0] eq $prevchr)&&($prevtype eq 'W')){
	$gap->[0]=$comp->[0];
 	$gap->[1]=$end+1;
	$gap->[2]=$end+10000;$end=$gap->[2];
	$gap->[3]=$comp->[3];$comp->[3]=$counter++;
	$gap->[5]=10000;
	print join("\t",@{$gap}),"\n";	
    }
    $comp->[1]=$end+1;
    my $length = $comp->[7]-$comp->[6]+1;
    $comp->[2]=$end+$length;
    $prevtype='W';
  }elsif($comp->[4] eq 'N'){
    $gap = $comp;
    $comp->[1]=$end+1;
    my $length = $comp->[5];
    $comp->[2]=$end+$length;
    $prevtype='N';
  }
  print join("\t",@{$comp}),"\n";
  $prevchr=$comp->[0];
  $end=$comp->[2];
}
