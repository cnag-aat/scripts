#!/usr/bin/env perl
use File::Basename;
use Getopt::Long;
use strict;
#my $base = undef;
my $col = 1;
my $tag = 0;
my $file = undef;
GetOptions( 
	   'f|g:s'      => \$file,
	   'c|col|column:s' => \$col,
	   't|tag|tagvalue:s' => \$tag
	  );
### GFF DATA STRUCTURE ###

my $seqname  = 0;
my $source   = 1;
my $feature  = 2;
my $start    = 3;
my $end      = 4;
my $score    = 5;
my $strand   = 6;
my $frame    = 7;
my $group    = 8;

die "$!" if !defined $file;
print STDERR "splitting gff file on tag \"$tag\"\n" if $tag;
#$base = basename($file);
my ($base,$path,$ext) = fileparse($file,qw(\.gff \.gtf \.gff3 \.gff2 \.bed));
my %newfile;
my %fh;
if ($tag) {
  die "tag splitting only supported for gtf and gff3 files\n" if !($ext eq ".gff3" || $ext eq ".gtf" );
  open (IN, "<$file");
  while (<IN>) {
    next if  m/^track name/;
    if ($_ !~ m/^[# ]/) {
      chomp;
      my @l = split "\t",$_;
      if ($ext eq ".gff3") {
	#print STDERR "GFF3\n";
	my %tv = split(/[=;]/,$l[8]);
	#my %tv = (@tval);
	my $name = $tv{$tag};
	if (exists $fh{$name}) {
	} else {
	  open my $newfh, ">","$base.".$name."$ext";
	  $fh{$name}=$newfh;
	  #print $newfh "something\n";
	}
	my $handle = $fh{$name};
	print $handle $_;
      } elsif ($ext eq "gff2") {

      } elsif ($ext eq "gtf") {
	my %tv = split(/[ ;]/,$l[8]);
	my $name = $tv{$tag};
	if (exists $fh{$name}) {
	} else {
	  open my $newfh, ">","$base.".$name."$ext";
	  $fh{$name}=$newfh;
	  #print $newfh "something\n";
	}
	my $handle = $fh{$name};
	print $handle $_;
      }
    }
  }
  close IN;
} else {
  #my %fh;
  open (IN, "<$file");
  while (<IN>) {
    if ($_ !~ m/^[# ]/) {
      my @l = split "\t",$_;
      my $name = $l[$col-1];
      if (exists $fh{$name}) {
      } else {
	open my $newfh, ">","$base.".$name."$ext";
	$fh{$name}=$newfh;
	#print $newfh "something\n";
      }
      my $handle = $fh{$name};
      print $handle $_;
      #push @{$newfile{$l[$col-1]}},$_;
      #print STDERR $l[0],"\n";
    }
  }
  close IN;
}
