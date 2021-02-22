#! /usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $in;
my $out;
my $result = GetOptions ("inFile=s" => \$in,
			 "outFile=s"=> \$out
    );

my %coverage;
my $count = 0;
my @pids;
my @head;
my $maxCov = 1;

# Run took 13 hours and 125 mb
# Open the in and out files

if($in =~ /\.gz$/){open(IN, "zcat $in |") or die "Could not open the inFile: \"$in\"\n";} else{open(IN, "<$in") or die " Could not open the inFile: \"$in\"\n";}
print "The inFile: \"$in\" was successfully opened!\n";
open(OUT, ">$out") or die "Could not open the outFile: \"$out\"\n";
print "The outFile: \"$out\" was successfully opened!\n";

# Look for the header

my $head = <IN>;
chomp $head;
if($head =~ /^chrom/){@head = split("\t", $head); print "Header of the inFile: \"$in\" found!\n";} else{die "The inFile: \"$in\" did not contain a valid header starting with \"chrom\"\n";}

# Create subhashes for each center

my $i = 3;
while($i < @head)
{
    $coverage{$i}{head} = $head[$i];
    $i++;
}

my $check = 1000000;
my $test = 0;
while(<IN>)
{
    chomp;
    my @line = split("\t", $_);
    $count += $line[2] - $line[1];
    $i = 3;
    while($i < @head)
    {
	
	if(defined $line[$i]){$coverage{$i}{$line[$i]} += $line[2] - $line[1];}else{print STDERR join("\t", @line), "\n"; die "$i element did not exist\n";}
	if(defined $line[$i] && $line[$i] > $maxCov){$maxCov = $line[$i];}
	$i++;
    }
    if($test == $check){print "$check lines read...\n";$check += 1000000;}
#last if($test > 100);
    $test++;
}
close IN;

my $headout = "#coverage";
$i = 3;
while($i < @head)
{
    $headout .= "\t".$coverage{$i}{head};
    $i++;
}
$headout .= "\n";
print OUT $headout;

my $j = 0;
while($j <= $maxCov)
{
    my $line = $j;
    $i = 3;
    while($i < @head)
    {
	if(exists $coverage{$i}{$j}){$line .= "\t".$coverage{$i}{$j};}
	else{$line .= "\t0";}
	$i++;
    }
    $line .= "\n";
    print OUT $line;
    $j++;
}
print "Program successfully finished!\nTotal bases = $count\nMaximum coverage = $maxCov\n";

# $i = 0;
# foreach(@head)
# {
# if(exists $coverage{$i}{head}){print $coverage{$i}{head}, "\t$i\n";}
# $i++;
# }


close OUT;
