#! /usr/bin/perl

# Program to print out total yield (base pairs) from a fastq
# Fernando Cruz, 2019-07-23

use strict;
use Getopt::Long;
use File::Basename;

# define variables
my $min_len=0;
my $count=0;
my $read_id;
my $seq;
my $sign;
my $qual;
my $seq_len;

my $yield=0;
my $total_reads=0;
# reading STDIN
while(<STDIN>)
{
 chomp;
 $count++;
 

 if ($count == 1){
     $read_id=$_;
     $total_reads+=1;
 }

 if ($count == 2) {
     $seq=$_;
     $seq_len=(length($seq));
     #print "seq_len $seq_len DNA $seq\n";
 }

 if ($count == 3) {
     $sign=$_;
 }

 if ($count == 4) {
     $qual=$_;

    # print"$read_id\n$seq\n$sign\n$qual\n";
    # print "$seq_len\n";
     $yield+=$seq_len;
     $count=0;
 }
    

}# end STDIN


print "Total yield\: $yield\n";
print "Total reads\: $total_reads\n";
#print "$yield\n";
exit;
