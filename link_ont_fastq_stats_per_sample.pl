#! /usr/bin/perl


# Note that this version was designed for the ONT reads of the CRE project.
# It names links to fastq files using the sample barcode ($sample) field as prefix. 
# Keeps the subproject name in the file name
# All the links are done inside the sample subdirectory
# Updated on 2020-03-11 to recongnize index "NA" in limsq table and look for the appropriate location. 
# Updated on 2020-06-04 to point to /scratch/project/devel/aateam
#                       to look for the ONT reads at both /scratch/production/shared/ONT/fastq and /project/production/shared/ONT/fastq
#                       NOTE that subproject CRE_13 points to both locations (project crashed and then everything started to be linked into /scratch/production)
#            2020-06-04 beside keeping just the passed libraries (libraryPassFail not "pass") it will discard failed sequencing runs (LanePassFail equal to "failed")
#            2020-06-17 contemplates the fact that some flowcells might have experienced pasue and reload at the gridion.
#                       this implies LIMS comment: "Flowcell with several non-indexed FLIs. The stats in the LIMS reflect the total evenly distributed amongst them. Only one fastq and fast5 will be found for the flowcell"
#                       In practice these FLIs are distinguished because they hold different library names, but both limsq lines point to the same fastq. So is IMPORTANT to avoid duplication of the data (having same fastq linked twice with different names)
#                       Solution: avoid including subproject and library in the soflink name. 
#                                 When duplicated then the system stderr will show up 
# Updated 2020-07-17 The program requires passing a subproject but the output is sent to a different subfolder per sample (sample barcode) inside the current working directory

use warnings;
use strict;

my $sp = $ARGV[0];
my $project_folder = "/project/production/shared/ONT/fastq"; #FAK60589/1/FAK60589_1_NB17.fastq.gz"
my $scratch_folder = "/scratch/production/shared/ONT/fastq"; #FAK60589/1/FAK60589_1_NB17.fastq.gz"
my $stats_folder = "/scratch/production/shared/ONT/stats";
#system "mkdir -p $sp";# output folder with subproject name
#system "module purge; module load jip-2.7/0.6; module load lims/1.2";
#print "$sp\n";

my $tbl = `/scratch/production/DAT/apps/LIMSQ/limsq -sp $sp`;
#print "$tbl\n";

my @f = split /\n/, $tbl;
my $total = scalar @f;
my $i = 1;
while ($i < $total) {
  #print "$f[$i]\n";
  my @l = split /;/, $f[$i];
  $i++;
  my ($sample,$sampleName, $subproject, $lib, $flowcell, $lane, $ind);
  $sampleName = $l[3];
  $sample = $l[4];
  $subproject = $l[1];# include subproject right before the library name
  $lib = $l[5];
  $flowcell = $l[9];
  $lane = $l[10];
  $ind = $l[11];#index or barcode to trim
  if ($ind eq "NA"){
      $ind=0;# allows locating the file, NA is the new convention in the limsq table
  }
  # Considers libraryPassFail status (pass) and discards failed sequencing runs (LanePassFail)   
  if (($l[12] ne "pass") || ($l[13] eq "fail")){ 
    print "$sampleName $lib $flowcell $lane $ind status is failed\n";
    next;
  }
  
  # Ensure the file link by looking into the new (/scratch) and old (/project) location:
  my $scratchFileToLink="$scratch_folder/$flowcell/$lane/$flowcell" ."_$lane"."_$ind" . ".fastq.gz";
  my $projectFileToLink = "$project_folder/$flowcell/$lane/$flowcell" ."_$lane"."_$ind" . ".fastq.gz";
  my $scratchStatDirToLink = "$stats_folder/$flowcell/$lane";

  #print STDERR "$fileToLink\n";
  if (-e $scratchFileToLink){
  	print STDERR "$scratchFileToLink exists\n";
  	system "mkdir -p $sample";
  	system "ln -s $scratchFileToLink $sample/$sample.$ind.$sampleName.$subproject.$flowcell.$lane.fastq.gz";  	
  }
  else{# link from old location
    
      if (-e $projectFileToLink){
	  print STDERR "$projectFileToLink exists\n";
	  system "mkdir -p $sample";
	  system "ln -s $projectFileToLink $sample/$sample.$ind.$sampleName.$subproject.$flowcell.$lane.fastq.gz";
      }
  }
if (-e $scratchStatDirToLink){
        print STDERR "$scratchStatDirToLink exists\n";
        system "mkdir -p $sample";
        system "ln -s $scratchStatDirToLink $sample/$sample.$ind.$sampleName.$subproject.$flowcell.$lane.stats";
  }


}#END while

print "WARNING:\nif symbolic links already exist, then, either you created some before or there is some Flowcell with several non-indexed FLIs.\n"
