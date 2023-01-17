#! /usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

############################################################
##
## script to convert coordinates from a gff3 with an agp or lookup table.
##
## Author: Jèssica Gómez
## 
#############################################################


my ($annot, $remove_list,$cut_bed, $new_scaffolds_bed, $lookup, $join_scaffolds, $agp);

GetOptions(
           'gff:s'           => \$annot, 
           'l:s'                => \$remove_list,  ## list with the scaffols that need to be completely removed from the annotation
           'k:s'             => \$lookup,   # Lookup table with the new names of the scaffolds that are kept untouched but with new names
           'cut:s'              => \$cut_bed,   #Bed file with the coordinates of the scaffolds to cut, but keeping the same names. in column 4 you must specify the offset to change the coordinates of the remaining genes
           'new:s'               => \$new_scaffolds_bed, #Bed file with the coordinates of the regions that must be kept from the old scaffolds that we want to divide in parts. In column 4 you must specify the new name to give to the scaffold. 
           'join:s'              => \$join_scaffolds, ## output from "get_offsets_from_agouti.pl" script. 1-old scaffold name; 2-new offset; 3-end in the new scaffold; 5-new scaffold name; 6-"+" or "-" according to the orientation of the old scaffold in the new scaffold
           'agp:s'               => \$agp, ## agp file obtained to place the scaffols into pseudomolecules
           );

my %remove;
if ($remove_list) {
    open RM, "<", "$remove_list";
    while (<RM>){
        chomp;
        $remove{$_}++;
    }
    close RM;
}

my %lookup;
if ($lookup) {
    open IDS, "<", "$lookup";
    while (<IDS>){
        chomp;
        my @line = split /\t/, $_;
        $lookup{$line[1]} = $line[0];
    }
    close IDS;
}

my %cut; 
if ($cut_bed){
    open CUT, "<", "$cut_bed"; 
    while (<CUT>){
        chomp;
        my @line = split /\s/, $_; 
        $cut{$line[0]}->{start} = $line[1];
        $cut{$line[0]}->{end} = $line[2];
        $cut{$line[0]}->{offset} = $line[3];
    }
    close CUT;
}

my %new;
if ($new_scaffolds_bed){
    open NEW, "<", "$new_scaffolds_bed";
    #print "$new_scaffolds_bed\n";
    while (<NEW>){
        chomp;
        my @line = split /\t/, $_;
        $new{$line[0]}{$line[1]}{$line[2]} = $line[3];
    }
    close NEW;   
}

my %join;
if ($join_scaffolds){
    open JOIN, "<", "$join_scaffolds";
    while (<JOIN>){
        chomp;
        my @line = split /\t/, $_;
        $join{$line[0]}->{offset} = $line[1];
        $join{$line[0]}->{end} = $line[2];
        $join{$line[0]}->{name} = $line[3];
        $join{$line[0]}->{orientation} = $line[4];
    }
    close JOIN;
}

my %chrs;
if ($agp){
    open AGP, "<", "$agp";
    while (<AGP>){
        chomp;
        next if /^\#/o;
        my @line = split /\t/, $_;
        next if ($line[4] ne 'W');
        my $start = $line[6] - 1;
        my $offset = $line[1]-1;
        $line[8] = "+" if ($line[8] eq "?"); 
        $chrs{$line[5]}{$start}->{chromosome} = $line[0];
        $chrs{$line[5]}{$start}->{offset} = $offset;
        $chrs{$line[5]}{$start}->{new_end} = $line[2];
        $chrs{$line[5]}{$start}->{end} = $line[7];
        $chrs{$line[5]}{$start}->{orientation} = $line[8];
        #print "$line[5]\t$start\t$chrs{$line[5]}{$start}->{end}\n";
    }
    close AGP;
}

print "##gff-version 3\n";
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
printf "# date: %4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;
my %genes;
my %transcripts;
my $name;
open ANNOT, "<", "$annot";
while (<ANNOT>){
    chomp;
    next if /^\#/o;
    next if /^$/o;
    my @line = split /\t/, $_;
    next if (exists $remove{$line[0]});
    if (exists $cut{$line[0]}){
        if ($cut{$line[0]}->{offset} == 0){
            next if ($line[3] > $cut{$line[0]}->{start});
            print "$_\n";    
        }
        else {
            next if ($line[4]< $cut{$line[0]}->{end});
            $line[3] = $line[3] - $cut{$line[0]}->{offset};
            $line[4] = $line[4] - $cut{$line[0]}->{offset};
            my $line = join "\t", @line;
            print "$line\n";
        } 
    }
    elsif (exists $lookup{$line[0]}) {
        my $prev = $line[0];
        $line[0] = $lookup{$prev};
        my $line = join "\t", @line;
        print "$line\n";
    }
    elsif (exists $new{$line[0]}){
        my $i = 0;
        my $j = 0;
        my %key;
        my $id; 
        foreach (sort {$a<=>$b} keys %{$new{$line[0]}}){
            $key{$line[0]}{"part" . $i} = $_;
            $i++;
        }
        while ($j < $i){
            my $start = $key{$line[0]}{"part" . $j};
            my $end;
            foreach (keys %{$new{$line[0]}{$start}}) {$end = $_}
            if ($line[3] >= $start && $line[4] <= $end){
                if ($line[2] eq "gene"){
                    $name = "";
                    if ($line[8] =~ m/ID=([^;]+)/){
                        $genes{$1}++;
                        $id = $1;
                        
                    }
                    elsif ($line[8] =~ m/ID=/){
                        $genes{$'}++;
                        $id = $';
                       # print STDERR "Gene $id\n";
                    }
                }
                elsif ($line[2] eq "transcript" || $line[2] eq "mRNA" || $line[2] eq "Protein_Alignment"){
                    if ($line[8] =~ m/ID=([^;]+)/) {
                        $id = $1;
                        #print STDERR "Transcript $id\n";
                    }
                    elsif ($line[8] =~ m/ID=/){
                        $id = $';
                    }
                    if ($line[8] =~ m/Parent=([^;]+)/){
                        if (exists $genes{$1}){
                            $transcripts{$id}++;
                        }
                    }
                    elsif ($line[8] =~ m/Parent=/){
                        if (exists $genes{$'}){
                            $transcripts{$id}++;
                        }
                    }
                }
                else {
                    if ($line[8] =~ m/Parent=([^;]+)/){
                        $id = $1;
                        #print STDERR "$id\n";
                    }
                    elsif ($line[8] =~ m/Parent=/){
                        $id = $';
                    }
                }
                $line[8] =~ s/$line[0]/$new{$line[0]}{$start}{$end}/g;
                my $old = $line[0];
                my $old_start = $line[3];
                my $old_end = $line[4];
                $line[0] = $new{$line[0]}{$start}{$end};
                $line[3] = $old_start - $start;
                $line[4] = $old_end - $start;
                my $l = join "\t", @line;
                print "$l\n" if (exists $genes{$id} || exists $transcripts{$id});
                $j = $i;
            }
            else {
                $j++;
            }
        }
    }
    elsif (exists $chrs{$line[0]}){
     # print "$_\n";
      my $old = $line[0];
      my $old_start = $line[3];
      my $old_end = $line[4];
      foreach my $start (sort {$a<=>$b} keys %{$chrs{$old}}){
        my $end = $chrs{$old}{$start}->{end};
       # print "$old\t$old_start\t$old_end\n";
        my $id = "";
        if ($old_start >= $start && $old_end <= $chrs{$old}{$start}->{end}){
         # print "$old\t$old_start\t$old_end\n";
          if ($line[2] eq "gene" ){
            if ($line[8] =~ m/ID=([^;]+)/){
              $genes{$1}++;
              $id = $1;          
            }
            elsif ($line[8] =~ m/ID=/){
              $genes{$'}++;
              $id = $';
            }
          }
          elsif ($line[2] eq "transcript" || $line[2] eq "mRNA" || $line[2] eq "Protein_Alignment"|| $line[2] eq "ncRNA"){
            if ($line[8] =~ m/ID=([^;]+)/) {
              $id = $1;
              $transcripts{$id}++ if ($line[8] !~ m/Parent/);
            }
            elsif ($line[8] =~ m/ID=/){
              $id = $';
              $transcripts{$id}++ if ($line[8] !~ m/Parent/);
            }
            if ($line[8] =~ m/Parent=([^;]+)/){
              if (exists $genes{$1}){
                $transcripts{$id}++;
              }
            }
            elsif ($line[8] =~ m/Parent=/){
              if (exists $genes{$'}){
                $transcripts{$id}++;
              }
            }
          }
          else {
            if ($line[8] =~ m/Parent=([^;]+)/){
              $id = $1;
            }
            elsif ($line[8] =~ m/Parent=/){
              $id = $1;
            }
          }
         # $line[8] =~ s/$line[0]/$new{$line[0]}{$start}{$end}/g;
          my @new = @line;
          $new[0] = $chrs{$old}{$start}->{chromosome};
          my $offset_end = $old_end - $start;
          my $offset_start = $old_start - $start;
          if ($chrs{$old}{$start}->{orientation} eq "+" ){
            $new[3] = $chrs{$old}{$start}->{offset} + $offset_start;
            $new[4] = $chrs{$old}{$start}->{offset} + $offset_end;
          }
          elsif ($chrs{$old}{$start}->{orientation} eq "-" ){
            #print "$chrs{$old}{$start}->{new_end}\t$chrs{$old}{$start}->{offset}\t$old_start\t$old_end\t$offset_end\n";
            $new[3] = $chrs{$old}{$start}->{new_end} - $offset_end + 1;
            $new[4] = $chrs{$old}{$start}->{new_end} - $offset_start + 1;
            my $strand = $line[6];
            if ($strand eq "+"){
                $new[6] = "-";
            }
            else {
                $new[6] = "+";
            }
          }
          else {
            die "9th column in the join option input must be \"+\" or \"-\" or \"?\". In $old, $line[0]\n";
          }
          my $l = join "\t", @new;
          print "$l\n" if (exists $genes{$id} || exists $transcripts{$id} || $_ !~  m/Parent/);
        }
        elsif ($_ =~ m/Repeat/){
          my @new = @line;
          $new[0] = $chrs{$old}{$start}->{chromosome};
          if ($old_start < $start && $old_end > $start ){
            my $offset_end = $old_end - $start;
            if ($chrs{$old}{$start}->{orientation} eq "+" ){
              $new[3] = $chrs{$old}{$start}->{offset} + 1;
              $new[4] = $chrs{$old}{$start}->{offset} + $offset_end;
            }
            else{
              $new[3] = $chrs{$old}{$start}->{new_end} - $offset_end + 1;
              $new[4] = $chrs{$old}{$start}->{new_end};
              if ($new[6] eq "+"){
                $new[6] = "-";
              }
              else {
                $new[6] = "+";
              }
            }
            my $l = join "\t", @new;
            print "$l\n";
          #  print STDERR "!!@line\n";
          }
          elsif ($old_start < $end && $old_end > $end){
          #  print STDERR "!!@line\n";
            my $offset_start = $end - $old_start; 
            if ($chrs{$old}{$start}->{orientation} eq "+" ){
              $new[3] = $chrs{$old}{$start}->{new_end} - $offset_start;
              $new[4] = $chrs{$old}{$start}->{new_end};
            }
            else {
              $new[3] = $chrs{$old}{$start}->{offset} + 1;
              $new[4] = $chrs{$old}{$start}->{offset} + $offset_start + 1;
              if ($new[6] eq "+"){
                $new[6] = "-";
              }
              else {
                $new[6] = "+";
              }
            }
            my $l = join "\t", @new;
            print "$l\n";
          }
        }
      }
      
    }
    elsif (exists $join{$line[0]}){
        #print "!!!!$_\n";
        my $prev = $line[0];
        my $start = $line[3];
        my $end = $line[4];
        my $strand = $line[6];
        $line[0] = $join{$prev}->{name};
        if ($join{$prev}->{orientation} eq "+" ){
            $line[3] = $start + $join{$prev}->{offset};
            $line[4] = $end + $join{$prev}->{offset};
        }
        elsif ($join{$prev}->{orientation} eq "-"){
            $line[3] = $join{$prev}->{end} - $end + 1;
            $line[4] = $join{$prev}->{end} - $start + 1;
            if ($strand eq "+"){
                $line[6] = "-";
            }
            else {
                $line[6] = "+";
            }
        }
        else {
            die "5th column in the join option input must be \"+\" or \"-\". In $prev, $line[0]\n";
        }
        my $l = join "\t", @line;
        print "$l\n";
    }
    else {
        print "$_\n";
    }
}
close ANNOT;
