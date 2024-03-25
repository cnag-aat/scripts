#! /usr/bin/perl
use strict;

# NOTE: This program estimates several statistics from the assembly.scaffolds.gaps.bed (e.g. sspaced.k10.scaffolds.gaps.bed)
# Date: 07-09-2017

#defining general variables:
my $tot=0;
my $sum_gap_length=0;
my $mean=0;
my $length;
#arrays & hashes
my @data;
# counter_intervals
my $gaps_1bp=0;
my $gaps_5bp=0;
my $gaps_10bp=0;
my $gaps_50bp=0;
my $gaps_100bp=0;
my $gaps_200bp=0;
my $gaps_300bp=0;
my $gaps_400bp=0;
my $gaps_500bp=0;
my $gaps_600bp=0;
my $gaps_700bp=0;
my $gaps_800bp=0;
my $gaps_900bp=0;
my $gaps_1000bp=0;

my $gaps_2k=0;
my $gaps_3k=0;
my $gaps_4k=0;
my $gaps_5k=0;
my $gaps_6k=0;
my $gaps_7k=0;
my $gaps_8k=0;
my $gaps_9k=0;
my $gaps_10k=0;

my $gaps_30k=0;
my $gaps_40k=0;
my $gaps_50k=0;

my $gaps_longer_than_50k=0;

my $syntax="\tUSAGE\: gaps-analysis.pl assembly.scaffolds.gaps.bed\n";



my $i=0;
my $f=0;

my $sd=0;
my $sum_sq_dev=0; # sum of square deviations

#command line input files and on the fly variable definitions using my
# @ARGV is a special variable that holds the arguments passed in from the command line

next unless ((exists $ARGV[0]) && ($ARGV[0]=~/\.gaps\.bed/)) or die print "*\.gaps\.bed file is required\!\!\!$syntax\n";
my $gap_file=$ARGV[0];


#workflow

mean_gap_length(); # estimates mean gap_lengthance between pairs
standard_dev();# estimates standard deviation, prints mean, sd, total gaps, and total length in gaps 
print_range_table();#print a tbale with a range of gap lengths
exit;

#subroutines:

sub mean_gap_length{

#it just loads true vars and their confidence


    $tot=0;
    open (IN, "$gap_file") or die print "can't find $gap_file\n";

    while (<IN>){

	chomp;
        @data=split/\t/,$_;
	$tot++;
       # print "$data[1]\t$data[2]\n";  
        $length= $data[2]-$data[1];
	$sum_gap_length+=$length;
        # Table covering PE range
        if ($length==1) {
	    $gaps_1bp++;
        }
        
        elsif (($length > 1) && ($length <= 5)){
            $gaps_5bp++;
        }
        
       
	elsif (($length > 5) && ($length <= 10)){
            $gaps_10bp++;
        }
       
	elsif (($length > 10) && ($length <= 50)){
            $gaps_50bp++;
        }

	elsif (($length > 50) && ($length <= 100)){
            $gaps_100bp++;
        }

	elsif (($length > 100) && ($length <= 200)){
            $gaps_200bp++;
        }

	elsif (($length > 200) && ($length <= 300)){
            $gaps_300bp++;
        }
       
	elsif (($length > 300) && ($length <= 400)){
            $gaps_400bp++;
        }

	elsif (($length > 400) && ($length <= 500)){
            $gaps_500bp++;
        }

        elsif (($length > 500) && ($length <= 600)){
            $gaps_600bp++;
	}

	elsif (($length > 600) && ($length <= 700)){
            $gaps_700bp++;
	}
      
        elsif (($length > 700) && ($length <= 800)){
            $gaps_800bp++;
	}
        elsif (($length > 800) && ($length <= 900)){
            $gaps_900bp++;
        }
	elsif (($length > 900) && ($length <= 1000)){
            $gaps_1000bp++;
        }


       # Table covering MP range
        elsif (($length > 1000) && ($length <= 2000)){
            $gaps_2k++;
        }
      	elsif (($length > 2000) && ($length <= 3000)){
            $gaps_3k++;
	}
	elsif (($length > 3000) && ($length <= 4000)){
            $gaps_4k++;
	}
	elsif (($length > 4000) && ($length <= 5000)){
            $gaps_5k++;
	}
	elsif (($length > 5000) && ($length <= 6000)){
            $gaps_6k++;
	}
	elsif (($length > 6000) && ($length <= 7000)){
            $gaps_7k++;
	}
	elsif (($length > 7000) && ($length <= 8000)){
            $gaps_8k++;
	}
	elsif (($length > 8000) && ($length <= 9000)){
            $gaps_9k++;
	}
	elsif (($length > 9000) && ($length <= 10000)){
            $gaps_10k++;
	}
	
       # Table covering FE range or longer...
	elsif (($length > 20000) && ($length <= 30000)){
            $gaps_30k++;
	}

	elsif (($length > 30000) && ($length <= 40000)){
            $gaps_40k++;
	}

	elsif (($length > 40000) && ($length <= 50000)){
            $gaps_50k++;
	}

      # Longer than 50Kb (out of range for MP and FE (such a long insert is very unlikely), outliers)
       elsif ($length > 50000){
       
             $gaps_longer_than_50k++;
       
       }

    }# while ENDS

    close (IN) or die print "unable to close $gap_file\n";

    $mean=$sum_gap_length/$tot;
#    print "mean\: $mean\n";


}# end calculate_gap_length


#######
sub  standard_dev {

    open (IN, "$gap_file") or die print "can't find $gap_file\n";

    while (<IN>){

	chomp;
	@data=split/\t/,$_;# split missing for some projects mean=sd
        
        $length= $data[2]-$data[1];
	
        $sum_sq_dev += (($length - $mean)**2);# square of deviations from the mean
      #  print "mean $mean value $_ sum_of_squared_deviation $sum_sq_dev \n";
       
    }#end while IN
   
    close (IN) or die print "cannot close the $gap_file\n";


       # unfold the loop and print                                                                                                                                                      #print "\n";
	   
	

	$sd=(($sum_sq_dev)/($tot-1))**(1/2); #square root is powered to 1/2...

        print "$gap_file\tmean\:\t$mean\tsd\:\t$sd\ttotal_gaps\:$tot\ttotal_bp\:$sum_gap_length\n";

} # end standard_dev

sub print_range_table {
 

    #open (OUT, ">".$vcf.".maf+allele_count.Rtable.txt") or die print "cannot print the outfile".$vcf.".maf.Rtable.txt\n";

    my $spacer = "\=" x 10;
    print "$spacer\n";
    #PE libs  
    print "1_bp\t$gaps_1bp\n";
    print "2\-5_bp\t$gaps_5bp\n";
    print "6\-10_bp\t$gaps_10bp\n";
    print "11\-50_bp\t$gaps_50bp\n";
    print "51\-100_bp\t$gaps_100bp\n";
    print "101\-200_bp\t$gaps_200bp\n";
    print "201\-300_bp\t$gaps_300bp\n";
    print "301\-400_bp\t$gaps_400bp\n";
    print "401\-500_bp\t$gaps_500bp\n";
    print "501\-600_bp\t$gaps_600bp\n";
    print "601\-700_bp\t$gaps_700bp\n";
    print "701\-800_bp\t$gaps_800bp\n";
    print "801\-900_bp\t$gaps_900bp\n";
    print "901\-1000_bp\t$gaps_1000bp\n";
    # MP
    print "1\-2kb\t$gaps_2k\n";
    print "2\-3kb\t$gaps_3k\n";
    print "3\-4kb\t$gaps_4k\n";
    print "4\-5kb\t$gaps_5k\n";
    print "5\-6kb\t$gaps_6k\n";
    print "6\-7kb\t$gaps_7k\n";
    print "7\-8kb\t$gaps_8k\n";
    print "8\-9kb\t$gaps_9k\n";
    print "9\-10kb\t$gaps_10k\n";
    # FE
    print "20\-30kb\t$gaps_30k\n";
    print "30\-40kb\t$gaps_40k\n";
    print "40\-50kb\t$gaps_50k\n";
    
    # LONG GAPS/OUTLIERS FOR ILLUMINA
    print "\>50kb\t$gaps_longer_than_50k\n";

}# end print_range_table

#

sub numerically
{

    $a<=>$b

}
