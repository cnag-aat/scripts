#! /usr/bin/perl

# VERSION 2.0
# Program to look into a FASTA and look for contaminated scaffolds in a list. These will be output into a contaminant fasta and the rest in assembly.decont.fa
# This version expects a 3-column list with scaffold id, bestsum of a high taxonomic rank and the bestsum of species name (spaces replaced by underscore)

# Fernando Cruz, 06-02-2023

use strict;
use Getopt::Long;
use File::Basename;

# define variables
my $file;
my $list;

my %contamlist=(); 

my @data=();

my $species="";
my $target_scaffold="";

my $contam="yes";
my $line;
my $scaf_id;
my $contam_id;

GetOptions(

           'fasta|f:s'   => \$file,
           'contaminants_list|l:s'   => \$list  #
);

# USAGE: perl -f assembly.fa -l contaminant_list.txt

my ($basename,$path,$ext) = fileparse($file,qw(\.fa));

#MAIN FUNCTIONS
#open the list file
target_list();
select_fasta();

#SUB-ROUTINES
sub target_list {
    open (LIST, "$list") or die "print I cannot open this list\n";
    
    while(<LIST>){
     chomp;
     @data=split/\t/, $_; 
     $target_scaffold=$data[0];

     $contam_id=$data[0]."\|".$data[1]."\|".$data[2];
   
     $contamlist{$target_scaffold}=$contam_id;    
   }

    close (LIST)  or die "print I cannot close this list\n";
#do nothin'
}# end subroutine

sub select_fasta {
  # cleaning Input FASTA
open (CONT, ">contaminants.fa") or die "print I cannot write this list\n";
open (DECONT, ">$basename.decont.fa") or die "print I cannot write this list\n";

open (IN, "$file") or die print "I cannot open file $file \n";
{

  while(<IN>)
  {
      chomp;
      $line=$_;
      
      if ($line=~/^>/) # header
      {   
	  $scaf_id=$line;
          $scaf_id =~ s/\>//g;
          
	  if (exists $contamlist{$scaf_id}) {
	      #print "$scaf_id\n";
	               $contam_id=$contamlist{$scaf_id};
                       print CONT "\>$contam_id\n";
	  }
          else {
	               print DECONT "\>$scaf_id\n";
          }

      }# FASTA header

      else{# sequence
            if (exists $contamlist{$scaf_id}){ 
               print CONT "$line\n";    #print sequence
            }
            else {
		print DECONT "$line\n";
	    }

      }#FASTA sequence

  }
    



  close (IN) or die print "cannot close the file $file\n";
}#FASTA closed

close (CONT) or die "print I cannot close this list\n";
close (DECONT) or die "print I cannot close this list\n";

}#clean_fasta sub-routine
