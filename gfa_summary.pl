#!/usr/bin/env perl

#from input like this:
#S 1 <seq> LN:i:5098144 dp:f:1.0
#L	1	+	1	+	0M
  
#output something like this:
#1 length=5101183 depth=1.00x circular=true
use strict;
my $file = shift @ARGV;
open (IN,"<$file") or die "Can't open $file!\n";
my %seqs=();
my %lengths=();
while (<IN>){
  chomp;
  my @f = split;
  if($f[0] eq 'S'){ #Segment
    #parse length and depth from tags
    #print "$f[1] $f[3] $f[4]\n";
    my $len = 0;
    my $dep = 0;
    if($f[3]=~m/LN:i:(\d+)/){$len=$1;}
    if($f[4]=~m/dp:f:(\S+)/){$dep=$1;}
    $seqs{$f[1]}->{length}=$len;
    $seqs{$f[1]}->{depth}=sprintf("%.2f", $dep);
    $lengths{$f[1]}=$len;
  }elsif($f[0] eq 'L'){ #Link
    #circular if linked to itself
    if($f[1] eq $f[3]){
      $seqs{$f[1]}->{circular}=1;
    }
  }
}
foreach my $s (sort { $lengths{$b} <=> $lengths{$a} } keys %lengths) {
  print join(" ",
	     ($s,
	      "length=".$seqs{$s}->{length},
	      "depth=".$seqs{$s}->{depth}."x",
	      $seqs{$s}->{circular}?"circular=true":"circular=false")
	    ),"\n";
      }
