#!/usr/bin/env perl
use strict;

my $basename = shift @ARGV;
open FF,">$basename.FF.dist";
open FR,">$basename.FR.dist";
open RF,">$basename.RF.dist";
my $incompat=0;
my $diff=0;
my $total=0;
my $fr=0;
my $rf=0;
my $ff=0;
while(my $l1=<>){
  $total++;
    my $l2 = <>;
    chomp $l1;
    chomp $l2;
    my @f1 = split "\t",$l1;
    my @f2 = split "\t",$l2;
    my $strata1= $f1[3];
    my $strata2= $f2[3];
    next if $strata1 !~ /^(0:)*1(:.*)?$/;
    next if $strata2 !~ /^(0:)*1(:.*)?$/;
    next if $f1[4] eq '-';
    next if $f2[4] eq '-';
    my ($seq1,$strand1,$pos1,$matchstring1)=split ":",$f1[4];
    my ($seq2,$strand2,$pos2,$matchstring2)=split ":",$f2[4];

    if($seq1 ne $seq2){
      $diff++;
    }else{
    my $s1start = $pos1;
    my $s2start = $pos2;
    my $s1end = 0;
    my $s2end = 0;
    my $tab1start = 0;
    my $tab1end = 0;
    my $tab2start = 0;
    my $tab2end = 0;
    # now parse matchstring to get end coord in reference
    my $length1 = 0;
    my $length2 = 0;
    while($matchstring1=~/(([ACGT])|(>\d+(\+|-))|(\d+))/g){
	my $m = $1;
	if ($m=~/[ACGT]/){
	    $length1++;
	}elsif($m=~/^\d/){
	    $length1+=$m;
	}else{
	    $m=/>(\d+)(.*)/;
	    #if ($2 eq '+'){ #insertion in read, so do not count this length
	    if ($2 eq '+'){ #deletion in read, so add this to length of match
		$length1+=$1;
	    }
	}
    }
    while($matchstring2=~/(([ACGT])|(>\d+(\+|-))|(\d+))/g){
	my $m = $1;
	if ($m=~/[ACGT]/){
	    $length2++;
	}elsif($m=~/^\d/){
	    $length2+=$m;
	}else{
	    $m=/>(\d+)(.*)/;
	    #if ($2 eq '+'){ #insertion in read, so do not count this length
	    if ($2 eq '+'){ #deletion in read, so add this to length of match
		$length2+=$1;
	    }
	}
    }
    $s1end = $s1start + $length1 - 1;
    $s2end = $s2start + $length2 - 1;
    if ($strand1 eq '+' && $strand2 eq '-'){
      if($s1start < $s2end){
	print FR $s2end-$s1start+1,"\n";
	$fr++;
      }else{
	print RF $s1end-$s2start+1,"\n";
	#print STDERR join("\t",($strand1, $s1start,$s1end,$strand2,$s2start,$s2end, $s1end-$s2start+1,$l1,$l2)),"\n";
	$rf++;
      }
    }elsif ($strand2 eq '+' && $strand1 eq '-'){
      if($s2start < $s1end){
	print FR $s1end-$s2start+1,"\n";
	$fr++;
      }else{
	print RF $s2end-$s1start+1,"\n";
	#print STDERR join("\t",($strand1, $s1start,$s1end,$strand2,$s2start,$s2end,$s2end-$s1start+1,$l1,$l2)),"\n";
	$rf++;
      }
    }else{
      if ($strand1 eq '+'){
	if ($s2end>$s1start){
	  print FF $s2end-$s1start+1,"\n";
	  $ff++;
	}else{
	  $incompat++;
	}
      }else{
	if ($s1end>$s2start){
	  print FF $s1end-$s2start+1,"\n";
	  $ff++;
	}else{
	  $incompat++;
	}
      }
    }
  }

}
close FF;
close FR;
close RF;

`plotDist.R $basename.RF.dist $basename.RF`;
`plotDist.R $basename.FF.dist $basename.FF`;
`plotDist.R $basename.FR.dist $basename.FR`;

print STDERR "FR:\t$fr\nRF:\t$rf\nFF:\t$ff\ndiff:\t$diff\nincompatible pairs:\t$incompat\n";
