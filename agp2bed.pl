#!/usr/bin/env perl

while(<>){
    next if m/^#/;
    my @F=split;
    if($F[4] eq "W"){
        print join("\t",($F[5],$F[6]-1,$F[7])),"\n";
    }
}
