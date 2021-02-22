#!/usr/bin/env perl
# Author: Tyler Alioto
# 18 July 2019
# Script to run Racon and Pilon to polish assemblies with ONT and Illumina reads, respectively
# Start by polishing with 2 rounds of racon and three rounds of pilon.
use Getopt::Long;
use File::Basename;
use warnings;
use strict;
use Getopt::Long 'HelpMessage';

=head1 NAME

Polish - polish with Racon and Pilon!

=head1 SYNOPSIS

  --input,-i      Input assembly (required)
  --ont,-o        ONT reads file (required)
  --pe1           Illumina read 1 (required)
  --pe2           Illumina read 2 (required)
  --rr            Number of rounds of racon (default: 2)
  --pr            Number of rounds of pilon (default: 3)
  --mr            Number of rounds of medaka (default: 1)
  -t              Number of threads
=head1 VERSION

1.0

=cut

my $PILON = "/apps/PILON/1.21/pilon";

my $usage = "usage:\n$0 -in <assembly.fasta> -ont <ONT_reads.fasta> -pe1 <illumina.read1.fastq> -pe2 <illumina.read2.fastq>\nReads may be in FASTA or FASTQ format, gzipped or not.\n";
GetOptions(
	   'i|input:s'  => \ my $in,
	   'ont|o:s'    => \ my $ont,
	   'pe1:s'      => \ my $pe1,
	   'pe2:s'      => \ my $pe2,
	   'rr:i'		=> \(my $rr = 0),
	   'pr:i'        => \(my $pr = 0),
	   't:i'        => \(my $threads = 1),
  	   'help'       => sub { HelpMessage(0) },
	  )or HelpMessage(1);
# die unless we got the mandatory argument
HelpMessage(1) unless ($in );
my ($b,$path,$ext) = fileparse($in,qw(\.fasta \.fa));

my $seq = $in;
if ($rr){
	for (my $i = 1; $i<=$rr; $i++){
		$seq = runRacon($seq,$i);
	}
}
if($pr){
	for (my $i = 1; $i<=$pr; $i++){
		$seq = runPilon($seq,$i);
	}
}
print "\nDone!\n";
sub runRacon {
	my $assembly = shift;
	my $round = shift;
	my ($base,$path,$ext) = fileparse($assembly,qw(\.fasta \.fa));
	print STDERR "\n###############################\n### Racon polishing round $round ###\n###############################\n\n";
	die "Input file $assembly not found! Aborting.\n" if !-e $assembly;
	if (! -e "$base.racon$round.fasta"){
		system("minimap2 -x map-ont -t 16 $assembly $ont | gzip > $base.racon$round.paf.gz");
		system("racon -t $threads $ont $base.racon$round.paf.gz $assembly > $base.racon$round.fasta");
		unlink "$base.racon$round.paf.gz";
	}else{
		print "Output file $base.racon$round.fasta already exists. Skipping this round.\n" 
	}
	return "$base.racon$round.fasta";
}
sub runPilon {
	my $assembly = shift;
	my $round = shift;
	my ($base,$path,$ext) = fileparse($assembly,qw(\.fasta \.fa));
	print STDERR "\n###############################\n### Pilon polishing round $round ###\n###############################\n\n";
	die "Input file $assembly not found! Aborting.\n" if !-e $assembly;
	if (! -e "$base.pilon$round.fasta"){
		system("bwa index $assembly");
		system("bwa mem -Y -t 16 $assembly $pe1 $pe2 | samtools view -Sb - | samtools sort -@ 8 -o $base.pilon$round.bam -");
		system("samtools index $base.pilon$round.bam");
		system("java -jar $PILON --genome $assembly --frags $base.pilon$round.bam --fix bases --changes --threads 16 --output $base.pilon$round");
		unlink ("$base.pilon$round.bam","$base.pilon$round.bam.bai");
		unlink glob "$assembly.*"
	}else{
		print "Output file $base.pilon$round.fasta already exists. Skipping this round.\n" 
	}
	return "$base.pilon$round.fasta";
}
