#! /usr/bin/perl

# Program to look into a FASTA and flip the scaffolds in a list doing the reverse-complement of its sequence
 
# Fernando Cruz, 05-12-2022

use strict;
use Getopt::Long;
use File::Basename;

# define variables
my $file;
my $list;

my %fliplist=(); 
my @data=();


my $target_scaffold="";
my $id="";
my $seq="";
my $r="";

GetOptions(

           'fasta|f:s'   => \$file,
           'target_scaffolds|l:s'   => \$list  #
);


my ($basename,$path,$ext) = fileparse($file,qw(\.scaffolds\.fa));

#MAIN FUNCTIONS
#open the list file
target_list();
revcomp_scaffolds();

#SUB-ROUTINES
sub target_list {
    open (LIST, "$list") or die "print I cannot open this list\n";
    
    while(<LIST>){
     chomp;
     $target_scaffold=$_;
     $fliplist{$target_scaffold}+=1;    
   }

    close (LIST)  or die "print I cannot close this list\n";
#do nothin'
}# end subroutine


sub revcomp_scaffolds() {
   
    
    open (IN,"FastaToTbl $file |");
    open OUT, "|TblToFasta > $basename.revcomp.fa";
    while(<IN>){

	chomp;
        @data=split/\s/, $_;# the split recognize the space as separator. later i should use tab. there is a difference between perl and bash
	$id=$data[0];
        $seq=$data[1];
   
        
        if (exists ($fliplist{$id})){
   
	    print OUT "$id","_rc\t";
	    $seq=~tr/acgtACGT/tgcaTGCA/;
	     $r=reverse($seq);
	    print OUT "$r\n";
        }
	else{
	    print OUT "$id\t$seq\n";
   
	}

    }
    close IN;
    close OUT;


}# end of revcomp_routine
