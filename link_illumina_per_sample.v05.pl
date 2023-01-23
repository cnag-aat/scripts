#! /usr/bin/perl

# JUST LINK THE CORRESPONDING FLI FATSQ IF STATUS IS NOT FAILED "pass" or "waiting" are tolerated
# Note that this version was designed for the illumina reads of the CRE project.
# It names links to fastq files using the sample barcode ($sample) field as prefix. 
# Keeps the subproject name in the file name
# All the links are done inside the sample subdirectory (created in the current directory!)
# Requires subproject and sample name

# Last update  2021-06-09

use warnings;
use strict;

my $sp = $ARGV[0]; # subproject
my $target=$ARGV[1]; # subproject

#system "mkdir -p reads; mkdir -p reads/illumina";

#system "module purge; module load jip-2.7/0.6; module load lims/1.2";
#system "module load Python/2.7.15-foss-2018b";
print "linking reads for each sample belonging to $sp subproject\n"; # we need to specify the subproject (need to know if it's illumina or not)

my $tbl = `module load Python/2.7.15-foss-2018b; /home/groups/assembly/talioto/repos/scripts/limsq -sp $sp`;
#print "$tbl\n";
my @f = split /\n/, $tbl;
my $total = scalar @f;
my $i = 1;
while ($i < $total){
  #print "$f[$i]\n";
  my @l = split /;/, $f[$i];
  $i++;
  my ($sample, $sampleName, $subproject, $lib, $flowcell, $lane, $ind);
  
  $sample = $l[4];
  $sampleName = $l[3];

  #next unless ($sampleName eq $target); 

  $subproject = $l[1];
  $lib = $l[5];
  $flowcell = $l[9];
  $lane = $l[10];
  $ind = $l[11];
  # JUST LINK REGARDLESS OF THE FLI PassFail status
  # Lane PassFail is actually FLI PassFail 
  if ($l[13] eq "fail"){
    print "$sample $sampleName $lib $flowcell $lane $ind status is failed\n";
    next;
  }
  if ($l[21] eq "join"){
    #print "@l\n";
    $ind = 0;
  } 
  if (-e "/scratch/project/production/fastq/$flowcell/$lane/fastq/$flowcell" . "_$lane" . "_$ind" . "_1.fastq.gz"){
        system "mkdir -p $sample";# One subdirectory per sample, is good for LUSTRE
  	system "ln -s  /scratch/project/production/fastq/$flowcell/$lane/fastq/$flowcell" . "_$lane" . "_$ind" . "_1.fastq.gz $sample/$sample.$sampleName.$subproject.$lib.$flowcell.$lane.$ind.1.fastq.gz";
  	print "/scratch/project/production/fastq/$flowcell/$lane/fastq/$flowcell" . "_$lane" . "_$ind" . "_1.fastq.gz -> $sample/$sample.$sampleName.$subproject.$lib.$flowcell.$lane.$ind.1.fastq.gz\n";  
  }
  if  (-e "/scratch/project/production/fastq/$flowcell/$lane/fastq/$flowcell" . "_$lane" . "_$ind" . "_2.fastq.gz"){
        system "mkdir -p $sample";# One subdirectory per sample, is good for LUSTRE
	system "ln -s  /scratch/project/production/fastq/$flowcell/$lane/fastq/$flowcell" . "_$lane" . "_$ind" . "_2.fastq.gz $sample/$sample.$sampleName.$subproject.$lib.$flowcell.$lane.$ind.2.fastq.gz";
  	print " /scratch/project/production/fastq/$flowcell/$lane/fastq/$flowcell" . "_$lane" . "_$ind" . "_2.fastq.gz -> $sample/$sample.$sampleName.$subproject.$lib.$flowcell.$lane.$ind.2.fastq.gz\n";
  }
}


