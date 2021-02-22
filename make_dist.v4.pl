#!/usr/bin/env perl
use strict;
use Getopt::Long;
my $I = 0;
my $T = 4;
my $e = 0.04;
my $m = 0.02;
my $name = $$;
GetOptions(
           'T:s' =>\$T,
	   'I:s' => \$I,
           'e:s' => \$e,
	   'm:s' => \$m,
	   'base|root|name:s' =>\$name
          );
die "$0 -I index.gem [-T <threads> -e <edit distance> -m <mismatch rate>]" if !$I;
my $cmd = "/apps/DEVEL/gemtools/1.7/bin/gem-mapper -I $I -q offset-33 -T $T --fast-mapping=0 -m $m -e $e --granularity 100000 --unique-mapping";
print STDERR "Mapping with command: $cmd\n";
open (OUT,"|$cmd | gem2dist.pl $name");
while(<>){print OUT;}




