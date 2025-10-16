#!/usr/bin/env perl

use Data::Dumper;
use Getopt::Long;
use File::Basename qw( fileparse );
my $printhelp = 0;
my %manifest = ();
#################
#### EXAMPLE ####
#################
# STUDY   PRJEB77742
# SAMPLE  SAMEA114562287
# RUN_REF ERR13374505
# ASSEMBLYNAME    ibLobSubt1.1
# ASSEMBLY_TYPE   isolate
# COVERAGE        46
# PROGRAM NEXTDENOVO, HYPO, PURGE_DUPS, YAHS
# PLATFORM        ONT, Illumina, omniC
# MINGAPLENGTH    100
# MOLECULETYPE    genomic DNA
# FASTA   ibLobSubt1_v8.curated_primary.mt.scrubbed.fa.gz
# CHROMOSOME_LIST ibLobSubt1.1.chromosome.mt.list.tsv.gz
# UNLOCALISED_LIST        ibLobSubt1.1.unlocalized.tsv.gz 
$manifest{ASSEMBLY_TYPE}='isolate';
$manifest{PLATFORM}='ONT, Illumina, Omni-C';
$manifest{PROGRAM}='HIFIASM, YAHS';
$manifest{MINGAPLENGTH}='100';
$manifest{MOLECULETYPE}='genomic DNA';
my $assembler = 'HIFIASM, YAHS';
my $assembly_fasta = 0;
my $chrlist = 0;
my $unlocs = 0;
my $coverage = 0;
my $hic_protocol = "Omni-C";
my $species = 0;
my $name = 0;
my $mito_fasta = 0;
my $biosample = 0;
my $given_assembly_project = 0;
my $given_runs = 0;
my $assembly_project = 0;
my $runstring = '';
GetOptions(
    'name:s' => \$name,
	'assembly|fasta|fa|f:s' => \$assembly_fasta,
    'MT|mito|mt|m:s' => \$mito_fasta,
	'chr:s' => \$chrlist,
#	'unlocs|unloc|u:s' => \$unlocs,
    'cov:s' => \$coverage,
    'assembler:s' => \$assembler,
    'species:s' => \$species,
    'hic:s' => \$hic_protocol,
    'sample:s' => \$biosample,
    'project:s' => \$given_assembly_project,
    'runs:s' => \$given_runs,
    'h|help'   => \$printhelp
);
$manifest{ASSEMBLYNAME}=$name;
$manifest{COVERAGE}=$coverage;
# get sample and runs from https://genomes.cnag.cat/erga-stream/ena_runs/
# get project from 
die "first run:\n conda activate /software/assembly/conda/gfastats-1.3.6-3/\n" if system("gfastats");
die "Mandatory arguments: -name -fasta -chr -cov -species -sample\n" if !($name && $assembly_fasta && $chrlist && $coverage &&  $species && $biosample);
die "$mito_fasta doesn't exist" if $mito_fasta and ! -e $mito_fasta;
die "$assembly_fasta doesn't exist" if $assembly_fasta and ! -e $assembly_fasta;
die "$chrlist doesn't exist" if $chrlist and ! -e $chrlist;
$manifest{ASSEMBLYNAME}=$name;

print STDERR "Finding study, runs and samples...\n";
if ($given_assembly_project){
    $assembly_project = $given_assembly_project;
}else{
    $assembly_project = `curl -X 'GET'   'https://www.ebi.ac.uk/ena/portal/api/search?result=study&query=study_tree%28PRJEB61747%29&limit=0&includeMetagenomes=true'   -H 'accept: */*' | grep assembly | grep -v alternate | grep "$species" | cut -f 1 `;
    chomp $assembly_project;
}
print "$assembly_project\n";
$manifest{STUDY}=$assembly_project;
if ($given_runs){
    $runstring = $given_runs;
}else{
    my $data_project = `curl -X 'GET'   'https://www.ebi.ac.uk/ena/portal/api/search?result=study&query=study_tree%28PRJEB61747%29&limit=0&includeMetagenomes=true'   -H 'accept: */*' | grep data | grep "$species" | cut -f 1 `;
    chomp $data_project;
    print "$data_project\n";

    my $runfetchcmd = "curl -X 'GET' 'https://www.ebi.ac.uk/ena/portal/api/search?query=study_tree%28$data_project%29&result=read_run&fields=instrument_platform,run_accession,sample_accession&limit=0' | grep OXFORD_NANOPORE | ";
    print "$runfetchcmd\n";
    open RUNS, $runfetchcmd or die $!;
    my @runs = ();
    my %samples = ();
    my $linecount = 0;
    while (my $line = <RUNS>) {
        $linecount++;
        chomp $line;
        my @fields = split ' ',$line;
        if ($fields[0] =~/^ERR/){
            print STDERR $fields[0],"\n";
            push @runs, $fields[0];
            $samples{$fields[2]}++;
        }
    }
    if ($linecount < 1 || not scalar @runs){
        $runstring = '';
        print STDERR "WARNING: no runs found\n";
    }else{
        $runstring = join(",",@runs);
    }
}

if($biosample){
    $manifest{SAMPLE} = $biosample;
}else{
    print STDERR "WARNING: more than one sample found; should use specimen-level biosample or create a virtual sample first\n" if keys %samples > 1;
    $manifest{SAMPLE} = join(",",keys %samples);
}

$manifest{RUN_REF} = $runstring;
$manifest{PROGRAM}= uc($assembler);
$manifest{PLATFORM}= 'ONT, Illumina, '.$hic_protocol;
my $mingap = 1;
if ($assembly_fasta !~/fasta/){
    my $rename = $assembly_fasta;
    $rename =~ s/fa/fasta/;
    `ln -s $assembly_fasta $rename`;
    $mingap = `gfastats $rename | grep 'Smallest gap' | sed 's/.*\s: //'`;
    unlink($rename);
}else{
    $mingap = `gfastats $assembly_fasta | grep 'Smallest gap' | sed 's/.*\s: //'`;
}
chomp $mingap;
$manifest{MINGAPLENGTH}=$mingap;

# 2. Create the qqGluDors1.3.chromosome.list.txt.gz: 
# grep -v unloc qqGluDors1.3.chromosome.list.csv | sed 's/,/\t/g' | sed 's/yes/chromosome/' > qqGluDors1.3.chromosome.list.txt 
# echo -e  "qqGluDors_MT\tMT\tcircular-chromosome\tMitochondrion" >> qqGluDors1.3.chromosome.list.txt 
# gzip qqGluDors1.3.chromosome.list.txt 
`grep -v unloc $chrlist | sed 's/,/\t/g' | sed 's/yes/chromosome/' > $name.chromosome.list.txt`;
die "$assembly_fasta cannot be named $name.asm.fa.gz\n" if $assembly_fasta eq "$name.asm.fa.gz";
die "$chrlist cannot be named $name.asm.fa.gz\n" if $chrlist eq "$name.chromosome.list.txt";
if ($mito_fasta){
    `gzip -cdf $assembly_fasta $mito_fasta | gzip -c --fast  > $name.asm.fa.gz`;# if ! -e "$name.asm.fa.gz";
    my $mt_name = `head -n 1 $mito_fasta | sed s'/>//'`;
    chomp $mt_name; 
    print STDERR "Found mitogenome: $mt_name\n";
    `echo -e  "$mt_name\tMT\tcircular-chromosome\tMitochondrion" >> $name.chromosome.list.txt`;
}else{
    `gzip -cdf $assembly_fasta | gzip -c --fast > $name.asm.fa.gz`;# if ! -e "$name.asm.fa.gz";
}
$manifest{FASTA} = "$name.asm.fa.gz";

`gzip -f $name.chromosome.list.txt`;
$manifest{CHROMOSOME_LIST} = "$name.chromosome.list.txt.gz";

my $unlocs = `grep -c unloc $chrlist`;  #| gzip > qqGluDors1.3.unlocalized.txt `
chomp $unlocs;
print STDERR "$unlocs unlocs found\n";
if (scalar $unlocs){
    `grep unloc $chrlist | sed 's/,/\t/g' | cut -f 1,2 | gzip -c > $name.unlocalized.txt.gz `;
    $manifest{UNLOCALISED_LIST} = "$name.unlocalized.txt.gz";
}
# 3. Create the qqGluDors1.3.unlocalized.txt.gz: 
# grep unloc qqGluDors1.3.chromosome.list.csv | sed 's/,/\t/g' | cut -f 1,2 > qqGluDors1.3.unlocalized.txt 
# gzip qqGluDors1.3.unlocalized.txt 

print STDERR $manifest{ASSEMBLYNAME},"\n";
print_manifest(\%manifest);



sub print_manifest{
    my $mani = shift;
    open MANIFEST, ">".$mani->{ASSEMBLYNAME}.".assembly_manifest.txt" or die("couldn't write manifest\n");
    foreach my $k (qw(STUDY SAMPLE RUN_REF ASSEMBLYNAME ASSEMBLY_TYPE COVERAGE PROGRAM PLATFORM MINGAPLENGTH MOLECULETYPE FASTA CHROMOSOME_LIST UNLOCALISED_LIST)){
        print MANIFEST "$k\t".$mani->{$k}."\n";
        print STDERR "$k\t".$mani->{$k}."\n";
    }
    close MANIFEST;
    `mkdir -p submit`;
    my $manifestfile = $mani->{ASSEMBLYNAME}.".assembly_manifest.txt";
    print "Now check your manifest for errors, then run:\nconda activate /software/assembly/conda/JAVA/\n"; 
    print "java -jar /software/assembly/src/webin-cli-validator/webin-cli-8.1.0.jar -context genome -userName WEBIN-1543 -password BGWSEi8y -manifest $manifestfile -centerName 'CNAG'   -outputDir submit -submit \n";

}