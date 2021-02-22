#!/usr/bin/env perl
use strict;
use Data::Dumper;
my $chr=0;
#chr_Pp04	10712648	10780842	213	W	pdulcis7_s0411	1	111761	-
#chr_Pp04	10854592	10864591	218	N	10000	scaffold	yes	peach:Pp04
#chr_Pp04	10864592	10918757	219	W	pdulcis7_s0183	209193	271603	+
my $prev_contig=0;
my $prev_gap=0;
my $last=0;
my $prev_component=undef;
my @ends=();
my @firsts=();
while(<>){
  next if m/^#/;
  chomp;
  my @f=split;
  next if $f[4]!~m/(N|W)/;
  if ($f[0] ne $chr){ #first contig
    if($last){push @ends,$last;}
    $chr=$f[0];
    if($f[4] eq "W"){ #contig/scaff
      my $comp = {"name"=>$f[5],"contig"=>\@f};
      push(@firsts,$comp);
      $last = $comp;
      $prev_contig=$comp;
      $prev_component="contig";
    }else{ #gap
      #this should NEVER happen
      print join("\t",@f),"\n";;
      die "WTF: chromosomes should not start with a gap\n";
      
    }
  }else{ #not first component
    if($f[4] eq "W"){ #contig/scaff
      if($prev_component ne "gap"){
	if($last){push @ends,$last;}
	my $comp = {"name"=>$f[5],"contig"=>\@f};
	push(@firsts,$comp);
	$last = $comp;
	$prev_contig=$comp;
      }else{
	my $comp = {"name"=>$f[5],"contig"=>\@f,"prev"=>$prev_contig,"gap"=>$prev_gap};
	$prev_contig->{next}=$comp;
	$last = $comp;
	$prev_contig = $comp;
      }
      $prev_component="contig";
    }elsif ($f[4] eq "N"){ #gap
      $prev_gap=\@f;
      $prev_component = "gap";
    }
  }
}
if($last){push @ends,$last;}

foreach my $first (@firsts) {

  my $s=$first;
  print STDERR "Starting at ",$s->{name},"\n";
  while (defined $s->{next}->{name}) {
    #next if !defined($s->{next}->{name});
    my $excised = 0;
    if (defined $s->{next}->{next}) {
      if (($s->{next}->{next}->{name} eq $s->{name})
	&&($s->{contig}->[8] eq $s->{next}->{next}->{contig}->[8])) {
	if ($s->{contig}->[8] eq "+") {
	  if ($s->{contig}->[7] == $s->{next}->{next}->{contig}->[6]  - 1) {
	    #remove contig in middle
	    $s->{contig}->[7]=$s->{next}->{next}->{contig}->[7];
	    if (defined $s->{next}->{next}->{next}){
	      $s->{next} = $s->{next}->{next}->{next};
	    }else{
	      $s->{next} = undef;
	    }
	    $excised=1;
	  }
	} else { #minus strand
	  if ($s->{contig}->[6] == $s->{next}->{next}->{contig}->[7]  + 1) {
	    #remove contig in middle
	    $s->{contig}->[6]=$s->{next}->{next}->{contig}->[6];
	    if (defined $s->{next}->{next}->{next}){
	      $s->{next} = $s->{next}->{next}->{next};
	    }else{
	      $s->{next} = undef;
	    }
	    $excised=1;
	  }
	}

      } elsif ((defined $s->{next}->{next}->{next} && $s->{next}->{next}->{next}->{name} eq $s->{name})
	      &&($s->{contig}->[8] eq $s->{next}->{next}->{next}->{contig}->[8])) {
	if ($s->{contig}->[8] eq "+") {
	  if ($s->{contig}->[7] == $s->{next}->{next}->{next}->{contig}->[6]  - 1) {
	    #remove contig in middle
	    $s->{contig}->[7]=$s->{next}->{next}->{next}->{contig}->[7];
	    if (defined $s->{next}->{next}->{next}->{next}){
	      $s->{next} = $s->{next}->{next}->{next}->{next};
	    }else{
	      $s->{next} = undef;
	    }
	    $excised=1;
	  }
	} else { #minus strand
	  if ($s->{contig}->[6] == $s->{next}->{next}->{next}->{contig}->[7]  + 1) {
	    #remove contig in middle
	    $s->{contig}->[6]=$s->{next}->{next}->{next}->{contig}->[6];
	    if (defined $s->{next}->{next}->{next}->{next}){
	      $s->{next} = $s->{next}->{next}->{next}->{next};
	    }else{
	      $s->{next} = undef;
	    }
	    $excised=1;
	  }
	}
      }elsif ((defined $s->{next}->{next}->{next}->{next} && $s->{next}->{next}->{next}->{next}->{name} eq $s->{name})
	      &&($s->{contig}->[8] eq $s->{next}->{next}->{next}->{next}->{contig}->[8])){
	if ($s->{contig}->[8] eq "+") {
	  if ($s->{contig}->[7] == $s->{next}->{next}->{next}->{next}->{contig}->[6]  - 1) {
	    #remove contig in middle
	    $s->{contig}->[7]=$s->{next}->{next}->{next}->{next}->{contig}->[7];
	    $excised=1;
	    if(defined $s->{next}->{next}->{next}->{next}->{next}){
	      $s->{next} = $s->{next}->{next}->{next}->{next}->{next};
	      print STDERR $s->{name},"\n";
	    }else{
	      $s->{next}=undef;
	    }
	  }
	} else { #minus strand
	  if ($s->{contig}->[6] == $s->{next}->{next}->{next}->{next}->{contig}->[7]  + 1) {
	    #remove contig in middle
	    $s->{contig}->[6]=$s->{next}->{next}->{next}->{next}->{contig}->[6];
	    $excised=1;
	    if(defined $s->{next}->{next}->{next}->{next}->{next}){
	      $s->{next} = $s->{next}->{next}->{next}->{next}->{next};
	      print STDERR $s->{name},"\n";
	    }else{
	      $s->{next}=undef;
	    }
	  }
	}
      }
    }
    if(!$excised){
      $s=$s->{next};
    }
  }   
}

# OK, now print the corrected chain out
my $c = 1; #let's make each chain unique
foreach my $first (@firsts) {
  my $s=$first;
  my $chr = $s->{contig}->[0] . "_$c";
  my $start=1;
  my $end=0;
  my $counter=1;
  $s->{contig}->[0]=$chr;
  $s->{contig}->[3]=$counter;
  $s->{contig}->[1]=$end+1;
  my $length = $s->{contig}->[7]-$s->{contig}->[6]+1;
  $s->{contig}->[2]=$end+$length;
  $end=$end+$length;
  print join("\t",@{$s->{contig}}),"\n";
  while (defined $s->{next}->{name}) {
    #print STDERR $s->{name},"\t",$s->{next}->{name},"\n";print STDERR Data::Dumper->Dump([$s->{next}]) if $s->{name} eq "pdulcis7_s0305";
    $s->{next}->{gap}->[0]=$chr;
    $s->{next}->{gap}->[3]=++$counter;
    $s->{next}->{gap}->[1]=$end+1;
    my $length = $s->{next}->{gap}->[5];
    $s->{next}->{gap}->[2]=$end+$length;
    $end=$end+$length;
    print join("\t",@{$s->{next}->{gap}}),"\n";
    $s->{next}->{contig}->[0]=$chr;
    $s->{next}->{contig}->[3]=++$counter;
    $s->{next}->{contig}->[1]=$end+1;
    my $length = $s->{next}->{contig}->[7]-$s->{next}->{contig}->[6]+1;
    $s->{next}->{contig}->[2]=$end+$length;
    $end=$end+$length;
    print join("\t",@{$s->{next}->{contig}}),"\n";
    $s=$s->{next};
  }
  $c++;
}
