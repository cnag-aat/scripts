#!/usr/bin/perl
use strict;
use warnings;

use Cwd;
use IPC::Open2;
use Time::HiRes qw(usleep nanosleep);
use File::Basename;
use File::Spec;
use Getopt::Long;
my $debug = 0;
my $usage = <<END;
usage: $0 [-debug] commands.pipe

   PIPE FORMAT is tab delimited with the following fields:
       job_name job_id dependencies partition class time tasks cpus_per_task initial_dir command
   
   job_name:      job name will show up in mnq
   job_id:        unique numeric id, usually sequential
   dependencies:  ':'-separated list of job_id's (field 2 of the pipe file), and/or the real id(s) of running job(s) prefixed by "j", for example,  "j16260898"
   partition:     No need to specify limited or projects (use a '.'), only himem, debug and development
                     For development, use this format: "development:XX" or "development:ExXX" where XX 
                     is your group (AA,AD,FB,StruG or StatG)
   class:         normal,highprio,lowprio
   time:          time in format minutes, hh:mm:ss, mm:ss, d-hh, d-hh:mm, d-hh:mm:ss
   tasks:         normally 1, but more if it's an MPI job
   cpus_per_task: number of cores per task
   initial_dir:   where to cd to before running command
   command:       self-explanatory (no tabs please!)

END
GetOptions(
	   'debug'      => \$debug,
	   );

my $limit = '24:00:00';
print $usage and exit if ! @ARGV;
my $pipeline = shift @ARGV;
my ($base,$path,$ext) = fileparse($pipeline,qr/\.[^.]*/);
my $name = $base;
my $time = scalar(localtime);
$time=~s/[\s:]/_/g;

my $class="normal";
my $mnsubmit="/opt/perf/bin/mnsubmit";

my $dir=cwd();
# File Format:
##jobname        id [integer]    dependencies [1:2:3]    class    time [minutes, hh:mm:ss, mm:ss, d-hh, d-hh:mm, d-hh:mm:ss] tasks cpu_per_task command  
open(CMDS,$pipeline) or die "couldn't open $pipeline: $!";
my @jobs;
while(<CMDS>){
  chomp;
  next if m/^#/;
  next if $_!~m/\S/;
  if(m/^>/){
    chomp;
    s/^>//;
    print STDERR $_,"\n";
    if (m/cd (\S+)/){
      chdir $1;
    }elsif (m/chdir (\S+)/){
      chdir $1;
    }else{
      system($_);
    }
    next;
  }
  my @f = split "\t",$_;
  die "Wrong number of fields\n$usage" if (scalar  @f) != 10;
  my($jobname,$id,$dep,$partition,$class,$time,$tasks,$cpt,$cwd,$command) = split "\t",$_;
  $command=~s/;$//;
  if ($class!~/[hldn]/i){$class='normal';}
  $jobname=~s/[^a-zA-Z0-9_-]+/_/g;
  push @jobs, {jobname=>$jobname,id=>$id,dep=>$dep,partition=>$partition,class=>$class,time=>$time,tasks=>$tasks,cpt=>$cpt,cwd=>$cwd,command=>$command};
}
close CMDS;


my %id;
foreach my $jobref (@jobs) {
  my %job = %$jobref;
  my $dependencies = '';
  if ($job{dep}=~/\d/){
    my @deps = split ":",$job{dep};
    my @jdeps;
    foreach my $d (@deps){
      if($d=~/j/i){
	$d=~s/[^\d]//;
	push @jdeps, $d;
      }else{
	push @jdeps, $id{$d} if exists $id{$d};
      }
    }
    $dependencies = '--dep afterok:'.join(":",@jdeps) if  scalar @jdeps;
  }
  if ($job{cwd} ne '.' && $job{cwd} ne '' && $job{cwd} ne '-' && $job{cwd} !~ /\$/){
    $dir=File::Spec->rel2abs($job{cwd});
  }
  my $logdir = "$dir/$base"."_log_".$time;
  `mkdir -p $logdir`;
  `cp $pipeline $logdir/`;
  my $logfile = "$logdir/$base.log";
  my $out =  "$logdir/$job{jobname}".".".$job{id}.'_%j.out';
  my $err =  "$logdir/$job{jobname}".".".$job{id}.'_%j.err';

  my $account = 0;
  my $constraint = 0;
  my $partition = $job{partition};
  if ($partition=~/development:(\S+)/){
    $partition = 'development';
    $constraint = $1;
    if($constraint=~/Ex(\S+)/){
      $account = $1;
    }else{
      $account = $constraint;
    }
  }
  my $partition_opt = '';
  if ($partition=~m/\w/){
    $partition_opt = " -p $partition ";
  }
  my $class = getClass($job{class});
  my $limit = getLimit($job{time},$class);
  my $jn = $job{jobname}.'_'.$job{id};
  $jn =~ s/[^A-Za-z0-9-_]//g;
  my $j = submit(
		 "echo $job{jobname}.$job{id} started: `date` >> $logfile; ".$job{command}."; echo $job{jobname}.$job{id} ended: `date` >> $logfile; ",
		 "$mnsubmit -job_name $jn $dependencies -cpus_per_task $job{cpt} -total_tasks $job{tasks} $partition_opt -c $class ",
		 $account,
		 $constraint,
		 $limit,
		 $dir,
		 $err,
		 $out,
		 $jn
	      );
  $id{$job{id}}=$j;
  print STDERR "$jn ($j) submitted\n";
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
    if (!$limit || $time > (60 * 60 * 24 * 14)) { 
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
  } elsif ($class eq 'himem') {
   if (!$limit || $time > (60 * 60 * 24 * 7)) { 
      $limit = '168:00:00';
    }
  } elsif ($class eq 'assembly') {
    if (!$limit || $time > (60 * 60 * 24 * 14)) {
      $limit = '336:00:00';
    }
  }
  return $limit;
}

sub submit {
  my ($scr,$cmd,$account,$constraint,$limit,$idir,$err,$out,$job)=@_;
  #print join("\t",($scr,$cmd,$account,$constraint,$limit,$err,$out,$job)),"\n" and exit;
  my $infile = $job.'.'.$$.'.'.time().".cmd";
  open (IN,">$infile") or die "couln't write to $infile\n";
  print IN <<EOF;
#!/bin/bash
# @ output           = $out
# @ error            = $err
# @ wall_clock_limit = $limit
EOF
  if($idir){
   print IN  '# @ initialdir = '."$idir\n";
  }
  if($account){
    print IN '# @ account = '."$account\n";
    if($constraint){
      print IN '# @ features = '."$constraint\n";
    }
  }
  # if ($class=~/himem|debug|limited|project/){
  #   print IN '# @ partition = '."$class\n";
  # }else{
  #   print IN '# @ class = '."$class\n";
  # }
  print IN "$scr\n";
  close IN;
  print STDERR "$cmd $infile\n";
  my @r = split " ", `$cmd $infile`;
  die "mnsubmit failed\n" if ! defined $r[3]; 
  my $jobid =$r[3];
  chomp $jobid;
  usleep(250000);
  if(!$debug){
	unlink $infile;
  }	
  return $jobid;
}
