#!/usr/bin/env perl
# kj2 08.10.2020
# jmdw added in hap flag and changed output suffix for files

use strict;
use Bio::SeqIO;
use Getopt::Long;
no warnings 'uninitialized'; # turn that off for debugging

my $tpf;
my $fa;
my $out;
my $csvin;
my $hap;
my $help;

GetOptions (
    "fa:s"  => \$fa,
    "tpf:s" => \$tpf,
    "out:s" => \$out,
    # "csv:s" => \$csvin,
    "hap"   => \$hap,
    "h"     => \$help,
    "help"  => \$help,
);

if (($help) || (!$fa) || (!$tpf)) {
    print "This script takes an original assembly fasta file and a TPF file generated with rapid_split.pl and creates the assembly from the TPF file\n";
    print "Usage:\n";
    print "perl rapid_join_mod.pl -fa <fasta>\n";
    print "             -tpf <tpf>\n";
    print "             -out <outfile_fasta_prefix> \# unless specified the output will be written to <fasta>.curated.fasta\n";
    print "             -hap \# optional use only if generating haplotigs fasta\n";
    print "             -h/help # this message\n";   
    exit(0);
}

my $newout;
if ($out) {
    $newout = $out;
}
else {
    ($newout) = ($fa =~ /(\S+)\.fa/);
    $out = ${newout}.".intermediate.fa" unless ($out);
}

my $seqout;
if ($hap) {
	$seqout = Bio::SeqIO->new(-format => 'fasta',
				-file => "> ${newout}.additional_haplotigs.unscrubbed.fa");
}
else {
	$seqout = Bio::SeqIO->new(-format => 'fasta',
				-file => "> ${newout}.out.fa");
}

# storing sequences
my %seqhash;
my $seqin  = Bio::SeqIO->new('-format' => 'fasta',
                             '-file'   => $fa);
while (my $seqobj = $seqin->next_seq()) {
    die("ERROR: you are trying to run rapid_join.pl on a full curation file\n") if $seqobj->display_id =~/_ctg1/;
    $seqhash{$seqobj->display_id} = $seqobj;
}

open(TPF,$tpf);
my $seqobj;
my $lastscaff;
my $seqstring;
my %seen;
my $no;
my $lastgap;
my %created;
my %assembly;
open(CSVOUU,">${newout}.chromosome.list.csv");
while (<TPF>) {
    $no++;
    my $line = $_;
    next if $line =~/^\n/;
    chomp;
    
    # gaps
    if (/^GAP/) {
        # no gaps at start of component
        die("Sequence can't begin with gap in line $no\n") unless ($lastscaff);

        my $length;
        if (/^GAP\s+\S+\s+(\d+)/) {
            $length = $1;
        }
        else {
            $length = 200;
        }
        $seqstring .= "N"x$length;
        $lastgap++;
    }
    # sequence
    else {
        my ($undef,$ctg,$scaff,$ori) = split;
        my ($ctgname,$start,$end) = ($ctg =~ /(\S+)\:(\d+)-(\d+)/);
        
        # check whether scaffold name was used previously
        die("ERROR: attempt to reuse scaffold name $scaff, line $no\n") if (exists $created{$scaff});
        
        # check whether sequence is being reused
        push @{$seen{$ctgname}}, [$start,$end];
        
        # continuation of scaffold
        if ($lastscaff &&($scaff eq $lastscaff)) {
            
            undef $lastgap;
            my $oldobj = $seqhash{$ctgname};
            my $subseq;
            eval {
                $subseq = $oldobj->subseq($start,$end);
            };   
            die("ERROR: $scaff:$start-$end does not exist\n") if ($@);
            
                 
            # revcom
            if ($ori eq "MINUS") {
                $subseq = reverse($subseq);
                $subseq =~ tr/atcgATCG/tagcTAGC/;
            }
            
            $seqstring .= $subseq;
        }
        
        # new scaff
        else {
            
            # check whether there's a gap at start of last scaffold
            die("ERROR: Previous component ended in gap, line ",$no-1,"\n") if ($lastgap);
            
            # out with last
            if ($lastscaff) {
                $seqobj->seq($seqstring);
                $seqobj->display_id($lastscaff);
                $seqout->write_seq($seqobj);
                ### WRITE to CHR LIST
                my $seqname = $lastscaff;
                my $is_chrom = "yes";
                if ($seqname=~/unloc/){
                    $is_chrom = "no";
                }
                if ($seqname =~m/SUPER_/){
                    $seqname =~s/SUPER_//;
                    if ($seqname =~  m/^([^_]+)/){
                        my $chrnum = $1;
                        print CSVOUU "$lastscaff,$chrnum,$is_chrom\n";
                    }
                }
                
                $assembly{$seqobj->display_id($lastscaff)} = $seqobj;
                $created{$lastscaff}= length($seqstring);
                #print "Created $lastscaff with ", length($seqstring),"\n";

                undef $seqobj;
                undef $seqstring;
            }
            
            $seqobj = Bio::Seq->new();
            my $oldobj = $seqhash{$ctgname};
            $seqobj->display_id($scaff);
            
            my $subseq = $oldobj->subseq($start,$end);
            
            # revcom
            if ($ori eq "MINUS") {
                $subseq = reverse($subseq);
                $subseq =~ tr/atcgATCG/tagcTAGC/;
            }
            
            $seqstring .= $subseq;
            $lastscaff = $scaff;
            
            undef $lastgap;
        
        }
    }
}
# out with last
$seqobj->seq($seqstring);
$seqobj->display_id($lastscaff);
$seqout->write_seq($seqobj);
### WRITE to CHR LIST
my $seqname = $lastscaff;
my $is_chrom = "yes";
if ($seqname=~/unloc/){
    $is_chrom = "no";
}
if ($seqname =~  m/SUPER_/){
    $seqname =~  s/SUPER_//;
    if ($seqname =~  m/^([^_]+)/){
        my $chrnum = $1;
        print CSVOUU "$lastscaff,$chrnum,$is_chrom\n";
    }
}



$assembly{$seqobj->display_id($lastscaff)} = $seqobj;
$created{$lastscaff}= length($seqstring);
#print "Created $lastscaff with ", length($seqstring),"\n";

# check for overlaps
foreach my $ctg (keys %seen) {
    my @sorted = sort {$a->[0] <=> $b->[0]} @{$seen{$ctg}};
    my $max;
    foreach my $aref (@sorted) {
        my $start = $aref->[0];
        my $end = $aref->[1];
        my $tick = 0;
        if ($max) {
            if ($start < $max) {
                $tick++;    
            }
            elsif (($start < $max) && ($end > $max)) {
                $tick++;
            }
            elsif ($end < $max) {
                $tick++;
            }
        }
        if ($tick > 0) {
            die("ERROR: Sequence from ctg $ctg was used more than once, e.g. in region $start-$end\n");
        }
        $max = $end unless ($max && ($max > $end));
    }
}

