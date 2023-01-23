#!/usr/bin/env perl
use strict;
use Cwd;
use Getopt::Long;
use File::Basename qw( fileparse );

my $cpt = 1;
my $name = 0;
my @com;
my $cmd = '';
my $initdir = getcwd();
#my $initdir = '.';
my $limit = 0;
my $pri = 'n';
my $priority = 'normal';
my $ext = 'out';
my $grp = 0;
my $cmdfile = 0;
my $mlimit = 5700000;
my $tpn = 8;
my $total_tasks = 1;
my $wait = 0;
my $tmpdir = 0;
my $partition_arg = "";
my $memarg = 0;
my $partition = 0;
my $mem = 0;
my ($commandname,$commandpath,$commandsuffix) = fileparse($0,(".pl"));
open CMDLINE, ">>$commandname.log";
my @ar;
foreach my $a (@ARGV) {
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
GetOptions(
	   'n|name:s'        => \$name,
	   'c|cmd|command:s' => \@com,
	   'd|dir:s' => \$initdir,
	   't|threads|cpus_per_task:s' => \$cpt,
	   'prio|priority:s'    => \$pri,
	   'l|limit:s'      => \$limit,
	   'oe:s'           => \$ext,
	   'group|g:s'      => \$grp,
	   'a:s'            => \$cmdfile,
	   'tpn:s'          => \$tpn,
	   'tasks:s'        => \$total_tasks,
	   'wait'           => \$wait,
	   'tmp:s'          => \$tmpdir,
	   'partition:s'   => \$partition,
           'mem:i'          => \$mem
	  );
if($partition){
  $partition_arg = " -partition $partition ";
}
if($mem){
  $memarg = " -memory $mem ";
}else{
  $memarg = " -memory 8G ";
}
$mlimit = $mlimit * $cpt;
my $most_tpn = int(8/$cpt);
if ($tpn > $most_tpn) {
  $tpn = $most_tpn;
}
### Let's handle maximum time limits for each queue
# lowprio 7 days
# normal 24 hours
# highprio 6 hours
my @lims;
if ($limit) {
  @lims = split ":",$limit;
  die "wrong limit format -- should be hh:mm:ss" if scalar(@lims)<3;
}
if ($pri =~ /^n/i) {
  $priority = 'normal';
  if (!$limit || $lims[0]>24) { 
    $limit = '24:00:00';
  }
}elsif ($pri =~ /^a/i) {
  $priority = 'assembly';
  if (!$limit || $lims[0]>336) {
    $limit = '336:00:00';
  }
} elsif ($pri =~ /^l/i) {
  $priority = 'lowprio';
  if (!$limit || $lims[0]>168) { 
    $limit = '168:00:00';
  }
} elsif ($pri =~ /^h/i) {
  $priority = 'highprio';
  if (!$limit || $lims[0]>6) { 
    $limit = '6:00:00';
  }
} elsif ($pri =~ /^d/i) {
  $priority = 'debug';
  if (!$limit || $lims[0]>6) { 
    $limit = '0:30:00';
  }
}
if ((scalar @com < 2)&& ($cmdfile)) {
  $cmd = $com[0] if scalar @com;
  @com = ();
  open (CMD,"<$cmdfile") or die "$!\n";
  while (<CMD>) {
    chomp;
    push @com, "$cmd $_";
  }
  close CMD;
}
$cmd = $com[0];
die "must provide a command" if !$cmd;

if (!$name) {
  my @args = split ' ', $cmd;
  $name = $args[0];
  #print STDERR "$name\n" and exit;
  $name =~ s/\.[^.\s]+$//;
}
$name = "$name.$$";
my $out =  $name.'_%j.'.$ext;
my $err =  $name.'_%j.err';
my $done =  $name.'.done';

# @ tasks_per_node   = $tpn
# ulimit -v $mlimit
open SCRIPT, ">$name.cmd";
print SCRIPT <<EOF;
#!/bin/bash
# @ job_name         = $name
# @ initialdir       = $initdir
# @ output           = $out
# @ error            = $err
# @ total_tasks      = $total_tasks
# @ wall_clock_limit = $limit
# @ cpus_per_task    = $cpt
date 1>&2
EOF
if ($tmpdir){
print SCRIPT "export TMPDIR=$tmpdir\n";
}
if ($grp) {
  print SCRIPT "newgrp $grp\n";
}
foreach my $cm (@com) {
  #print SCRIPT "echo '$cm' 1>&2\n";
  print SCRIPT "time $cm\n";
}
print SCRIPT <<EOF;
date 1>&2
echo '$name' > $done
EOF
close SCRIPT;
print STDERR "mnsubmit $memarg $partition_arg -c $priority $name.cmd\n";
my @r = split " ",`mnsubmit $memarg $partition_arg -c $priority $name.cmd`;
print STDERR $r[3],"\n";
print $r[3],"\n";
if ($wait) {
  my $jobid = $r[3];
  my $notdone = 1;
  while ($notdone) {
    $notdone = 0;
    sleep 30;
    open SQ,"squeue |" or die "couldn't open in pipe squeue!\n";
    while (<SQ>) {
      my @result = split;
      if ($result[0] eq $jobid) {
	$notdone++;
      }
    }
    close SQ;
    
  }
}

