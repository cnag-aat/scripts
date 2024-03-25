#!/usr/bin/env perl
use strict;
use Cwd;
use Getopt::Long;
use File::Spec;
use File::Basename qw( fileparse );
use File::Copy;

my $cwdir = getcwd();
my $pri = 'normal';
my $argfile = 0;
my $name = 0;
my $cpt = 1;
my $batchsize = 0;
my $partition = 'genD';
my $limit = "6:00:00";
my @ARGS = @ARGV;
GetOptions(
	   'n|name:s'        => \$name,
	   'a|args:s'        => \$argfile,
	   't|threads|cpus_per_task:s' => \$cpt,
	   'priority|qos:s'    => \$pri,
	   'l|limit:s'      => \$limit,
	   'b|batchsize:i' => \$batchsize,
	   'partition:s'   => \$partition
	  );


die "jobname is required: -n " if !$name;
my $class = getClass($pri);
my $tlimit = getLimit($limit);

my $numtasks=0;
if ($argfile) {
  open ARGS,"<$argfile";
  while (<ARGS>){
    $numtasks++ if m/\S/  && $_!~/^#/;
  }
  close ARGS;
  chmod '0444',$argfile;
}else{
  die "No command file given with -a switch\n";
}
my $batch='';
if($batchsize){$batch = '%'.$batchsize;}
my $timestring = "".localtime(time);
$timestring=~s/\s+/_/g;
$timestring=~s/:/-/g;
my $script = $name.'.'.$$.'.'.$timestring.".sh";
my $logdir = $name.'.'.$$.'.'.$timestring;
`mkdir -p $logdir`;
copy($argfile,"$logdir/$argfile");
open (SCR,">$logdir/$script") or die "couln't write to $logdir/$script\n";
my $outbase = "$logdir/$name".'_%A_%a';
my $out = "$logdir/$name".'_%A_%a.out';
my $err = "$logdir/$name".'_%A_%a.err';
print SCR qq(#!/bin/bash
#SBATCH --job-name=$name
#SBATCH --output=$out
#SBATCH --error=$err
#SBATCH --ntasks=1
#SBATCH --array=1-$numtasks$batch
#SBATCH --time=$tlimit
#SBATCH --partition $partition 
#SBATCH --qos=$class
#SBATCH --cpus-per-task $cpt
#SBATCH -D $cwdir 
);
my $jobfile = $name.'.$SLURM_ARRAY_TASK_ID.sh';
print SCR q/CMD=$(gawk "FNR==$SLURM_ARRAY_TASK_ID {print}" /,"$argfile)"; 
print SCR "\n",q(eval "$CMD"),"\n";
close SCR;

print STDERR "Running sbatch $logdir/$script\n";
#exit;
my ($commandname,$commandpath,$commandsuffix) = fileparse($0,(".pl"));
open CMDLINE, ">>$logdir/$commandname.log";
my @ar;
foreach my $a (@ARGS) {
  $b = '';
  if ($a!~/^-/) {
    $b = "'$a'";
  } else {
    $b=$a;
  }
  push @ar,$b;
}
print CMDLINE scalar(localtime),"\n$0 ",join(" ", @ar),"\n\n";
print STDERR scalar(localtime),"\n$0 ",join(" ", @ar),"\n";
close CMDLINE;



system("sbatch","$logdir/$script");

sub getLimit{
#Acceptable  time  formats  include  "minutes",  "minutes:seconds", "hours:minutes:seconds",  "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds".
  my $limit = shift;
  my $class = shift;
  my @lims;
  my $days = '0';
  my $hours = '00';
  my $minutes = '00';
  my $seconds = '00';
  my $time = 0;
#print STDERR "$limit\n";
  if ($limit){
    if($limit=~/(\d+)-(\S+)/){
      $days = $1;
      my $rest = $2;
      ($hours,$minutes,$seconds) = split ":",$rest;
      $hours+=($days * 24);
      $days=0;
$time = (defined($seconds)?$seconds:0) + (60 * (defined($minutes)?$minutes:0)) + (60 * 60 * (defined($hours)?$hours:0)) + (60 * 60 * 24 * $days);
    }elsif($limit=~/:/){
      my ($sec,$min,$hou) = reverse(split ":",$limit);
      $seconds = $sec if defined $sec;
      $minutes = $min if defined $min;
      $hours = $hou if defined $hou;
    }elsif($limit=~/^\d+$/){
      $minutes = $limit;
    }
    #convert to seconds
    $time = (defined($seconds)?$seconds:0) + (60 * defined($minutes)?$minutes:0) + (60 * 60 * defined($hours)?$hours:0) + (60 * 60 * 24 * $days);
  }
  #$limit= sprintf("%d:%02d:%02d",$hours,$minutes,$seconds);
  if ($class eq 'normal'){
    if (!$limit || $time > (60 * 60 * 24)) { 
      $limit = '24:00:00';
    }
  } elsif ($class eq 'lowprio') {
    if (!$limit || $time > (60 * 60 * 24 * 24)) { 
      $limit = '336:00:00';
    }
  } elsif ($class eq 'highprio') {
    if (!$limit || $time > (60 * 60 * 6)) { 
      $limit = '6:00:00';
    }
  } elsif ($class eq 'xlong') {
    if (!$limit || $time > (60 * 60 * 24 * 14)) {
      $limit = '336:00:00';
    }
  }elsif ($class eq 'debug'){
    if (!$limit || $time > (60 * 60 * 6)) { 
      $limit = '6:00:00';
    }
  } elsif ($class eq 'assembly') {
    if (!$limit || $time > (60 * 60 * 24 * 60)) {
      $limit = '1440:00:00';
    }
  }
  return $limit;
}

sub getClass{
  my $pri= shift;
  my $priority = "normal";
  if ($pri =~ /^n/i) {
    $priority = 'normal';
  } elsif ($pri =~ /^l/i) {
    $priority = 'lowprio';
  } elsif ($pri =~ /^h/i) {
    $priority = 'highprio';
  } elsif ($pri =~/^a/){
	$priority = 'assembly';
  } elsif ($pri =~/^x/){
	$priority = 'xlong';
  } 
  #print STDERR "$priority\n";
  return $priority;
}

exit;





























