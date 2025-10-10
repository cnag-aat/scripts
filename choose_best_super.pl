#!/usr/bin/env perl
use strict;
use Getopt::Long;
use File::Copy "cp";

my $gaps1 = '';
my $gaps2 = '';
my $in1 = '';
my $in2 = '';
my $out1 = '';
my $out2 = '';
my $pref = "hap2";
my $numsupers=23;

GetOptions(
'g1:s' => \$gaps1,
'g2:s' => \$gaps2,
'i1:s' => \$in1,
'i2:s' => \$in2,
'o1:s' => \$out1,
'o2:s' => \$out2,
'pref:s' => \$pref,
'num:s' => \$numsupers
);
my %gaps_h1 = ();
my %gaps_h2 = ();

open GAP1, "<$gaps1" or die "Couldn't open $gaps1\n";
while(<GAP1>){
    my @f = split;
    my $super = $f[0];
    $super=~s/Chr/SUPER_/;
    $gaps_h1{$super}++;
}
close GAP1;
open GAP2, "<$gaps2" or die "Couldn't open $gaps2\n";
while(<GAP2>){
    my @f = split;
    my $super = $f[0];
    $super=~s/Chr/SUPER_/;
    $gaps_h2{$super}++;
}
close GAP2;
if ($pref eq '$hap1'){
    for (my $i = 1; $i<=23; $i++){
        print STDERR $i,"\n";
        if (! exists($gaps_h1{"SUPER_".$i})){
            `cat $in1/SUPER_$i.fa >> $out1`;
            my $unlocs1 = "$in1/SUPER_$i"."_unloc*.fa";
            `cat $unlocs1 >> $out1`;

            `cat $in2/SUPER_$i.fa >> $out2`;
            my $unlocs2 = "$in2/SUPER_$i"."_unloc*.fa";
            `cat $unlocs2 >> $out2`;
        }elsif(! exists($gaps_h2{"SUPER_".$i})){
            `cat $in2/SUPER_$i.fa >> $out1`;
            my $unlocs2 = "$in2/SUPER_$i"."_unloc*.fa";
            `cat $unlocs2 >> $out1`;
            
            `cat $in1/SUPER_$i.fa >> $out2`;
            my $unlocs1 = "$in1/SUPER_$i"."_unloc*.fa";
            `cat $unlocs1 >> $out2`;
        }elsif($gaps_h2{"SUPER_".$i} >=  $gaps_h1{"SUPER_".$i}){ #hap1 is better super
            
            `cat $in1/SUPER_$i.fa >> $out1`;
            my $unlocs1 = "$in1/SUPER_$i"."_unloc*.fa";
            `cat $unlocs1 >> $out1`;

            `cat $in2/SUPER_$i.fa >> $out2`;
            my $unlocs2 = "$in2/SUPER_$i"."_unloc*.fa";
            `cat $unlocs2 >> $out2`;

        }else{ #hap2 is better super
            `cat $in2/SUPER_$i.fa >> $out1`;
            my $unlocs2 = "$in2/SUPER_$i"."_unloc*.fa";
            `cat $unlocs2 >> $out1`;
            
            `cat $in1/SUPER_$i.fa >> $out2`;
            my $unlocs1 = "$in1/SUPER_$i"."_unloc*.fa";
            `cat $unlocs1 >> $out2`;

        }
    }
    
    `cat $in1/*scaffold* >> $out1`;
    `cat $in2/*scaffold* >> $out2`;
}else{
    
    for (my $i = 1; $i<=$numsupers; $i++){
        print STDERR $i,"\n";
        if (! exists($gaps_h2{"SUPER_".$i})){
            `cat $in2/SUPER_$i.fa >> $out1`;
            my $unlocs2 = "$in2/SUPER_$i"."_unloc*.fa";
            `cat $unlocs2 >> $out1`;
            
            `cat $in1/SUPER_$i.fa >> $out2`;
            my $unlocs1 = "$in1/SUPER_$i"."_unloc*.fa";
            `cat $unlocs1 >> $out2`;
        }elsif(! exists($gaps_h1{"SUPER_".$i})){
            `cat $in1/SUPER_$i.fa >> $out1`;
            my $unlocs1 = "$in1/SUPER_$i"."_unloc*.fa";
            `cat $unlocs1 >> $out1`;

            `cat $in2/SUPER_$i.fa >> $out2`;
            my $unlocs2 = "$in2/SUPER_$i"."_unloc*.fa";
            `cat $unlocs2 >> $out2`;
        }elsif($gaps_h1{"SUPER_".$i} >=  $gaps_h2{"SUPER_".$i}){ #hap2 is better super
            `cat $in2/SUPER_$i.fa >> $out1`;
            my $unlocs2 = "$in2/SUPER_$i"."_unloc*.fa";
            `cat $unlocs2 >> $out1`;
            
            `cat $in1/SUPER_$i.fa >> $out2`;
            my $unlocs1 = "$in1/SUPER_$i"."_unloc*.fa";
            `cat $unlocs1>> $out2`;

        }else{ #hap1 is better super
            
            `cat $in1/SUPER_$i.fa >> $out1`;
            my $unlocs1 = "$in1/SUPER_$i"."_unloc*.fa";
            `cat $unlocs1 >> $out1`;

            `cat $in2/SUPER_$i.fa >> $out2`;
            my $unlocs2 = "$in2/SUPER_$i"."_unloc*.fa";
            `cat $unlocs2 >> $out2`;

        }
    }

    `cat $in2/*scaffold* >> $out1`;
    `cat $in1/*scaffold* >> $out2`;
}
