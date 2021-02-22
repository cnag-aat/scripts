#!/usr/bin/env perl
 
#use Pod::Usage;
use Getopt::Long;
use File::Basename;
 
$this_program = basename($0);
my $minmatch = 200;
my $minqual = 12;
my $output=0;
my $input=0;
my $calmd=0;
my $ref=0;
my $result = GetOptions(
			'i|input=s' => \$input,
			'q=i'   => \$minqual,
			'm|min=i'   => \$minmatch,
			'calmd' => \$calmd,
			'ref=s' => \$ref
			
		       );
my $usage = "usage: $this_program -i in.sam -o outname [-q minqual] [-m min_aligned_bases] [-calmd] [-ref]\n";
print $usage if (!$input || $input!~/\.sam$/) and exit;
print($usage,"Must specify ref if using calmd option\n") if ($calmd and !$ref) and exit;
my ($base,$path,$ext) = fileparse($input,qw(\.sam));
my $totalreadlength=0;
my $totalreadalnlength=0;
my $total_ievents=0;
my $total_devents=0;
my $total_alignedbases=0;
my $total_soft=0;
my $total_hard=0;
my $total_insert=0;
my $total_del=0;
my $total_subs=0;
my $totalreflength=0;
open (TABLE, ">$base.error_table.txt");
print TABLE join("\t",qw(read_id read_length accuracy substitution_error deletion_error deletion_rate insertion_error insertion_rate aligned softclip hardclip)),"\n";
if($calmd){
  open (SAM,"samtools calmd -S $input $ref | grep MD:Z | sam2pairwise |");
}else{
  open (SAM," grep MD:Z $input| sam2pairwise |");
}
while (my $samstring=<SAM>) {
  my $aln_ref = <SAM>;
  my $matches = <SAM>;
  my $aln_read= <SAM>;
  $t++;
  print STDERR " $t reads parsed.\r" unless ($t % 50000);
  my @samfields = split /\t/,$samstring;
  my $strand = strand($samfields[1]);
  my $contig = $samfields[2];
  next if $contig eq "*";
  my $read_id   = $samfields[0];
  my $from   = $samfields[3];
  my $mapq = $samfields[4];
  next if $mapq < $minqual;
  my $refalnlen    = 0;
  my $matches = 0;
  my $insert = 0;
  my $del = 0;
  my $soft = 0;
  my $hard = 0;
  my $intron = 0;
  my $ievent = 0;
  my $devent = 0;
  my $sub = 0;
  if ($samfields[4] < $minq) {
    $lowqual++;
  }
  while ($samfields[5]=~/(\d+)M/g) {
    $matches += $1;
  }
  while ($samfields[5]=~/(\d+)I/g) {
    $insert += $1;
    $ievent++;
  }
  while ($samfields[5]=~/(\d+)D/g) {
    $del += $1;
    $devent++;
  }
  while ($samfields[5]=~/(\d+)S/g) {
    $soft += $1;
  }
  while ($samfields[5]=~/(\d+)H/g) {
    $hard += $1;
  }
  while ($samfields[5]=~/(\d+)N/g) {
    $intron += $1;
  }
  $refalnlen = $matches + $del + $intron;
  $to = $from + $refalnlen -1;
  my $qalnlen = $matches + $insert;
  my $readlen = $matches + $insert + $soft +$hard;
  # M+I+S=ALNQLEN
  # M
  my $lefthard = 0;
  my $righthard = 0;
  my $leftsoft = 0;
  my $rightsoft = 0;
  my $qstart = 0;
  my $qend = $readlen;
  if ($samfields[5]=~s/^(\d+)H//) {
    $lefthard=$1;
  }
  if ($samfields[5]=~/^(\d+)S/) {
    $leftsoft=$1;
  }
  if ($samfields[5]=~s/(\d+)H$//) {
    $righthard=$1;
  }
  if ($samfields[5]=~/(\d+)S$/) {
    $rightsoft=$1;
  }
  if ($strand eq "+") {
    $qstart = 1 + $lefthard +$leftsoft;
    $qend = $readlen - $righthard -$rightsoft;
  } else {
    $qstart = 1 + $righthard + $rightsoft;
    $qend = $readlen - $lefthard - $leftsoft;
  }
  next if $refalnlen < $minmatch;
    
  my $read=$aln_read;
  $read=~s/-//g;
  my $readlength=length($read);
  $clipped_read=$read;
  $clipped_read=~s/N//g;
  my $readalnlen=length($clipped_read);
  my $ref=$aln_ref;
  $ref=~s/-//g;
  my $reflength=length($ref);
  

  my $ed = '0';
  my $same = 0;
  for (my $i =0; $i < length($aln_ref); ++$i) {
    my $q = substr($aln_read,$i,1);
    my $r = substr($aln_ref,$i,1);
    if ($q ne $r) {
      ++$ed;
      #annotated as a deletion into the reference sequence
      if ($q eq '-') {
      }
      #annotated as an insertion into the reference sequence
      elsif ($r eq '-') {
      }
      #annotated as soft clip
      elsif ($q eq 'N') {
      }
      #mismatch
      else {
	if ($q eq 'A') {
	  $sub++;
	} elsif ($q eq 'C') {
	  $sub++;
	} elsif ($q eq 'G') {
	  $sub++;
	} elsif ($q eq 'T') {
	  $sub++;
	} else {
	  die "That was unexpected q = $q que = $query ref = $reference\n";
	}
      }
    }else{
      $same++;
    }
  }
  my $alnlen = length($aln_read);
  $totalalnlen += $alnlen;
  $totalreadlength+=$readlength;
  $totalreflength+=$reflength;
  $totalreadalnlength+=$readalnlen;
  my $softpct = $soft/$readlength;
  $total_soft+=$soft;
  my $hardpct = $hard/($readlength+$hard); #hard clipping removes read sequence
  $total_hard+=$hard;
  my $insert_error = $insert/$refalnlen;
  $total_insert+=$insert;
  my $del_error = $del/$refalnlen;
  $total_del+=$del;
  my $subs_error = $sub/$refalnlen;
  $total_subs+=$sub;
  my $insrate = $ievent/$refalnlen;
  $total_ievents+=$ievent;
  my $delrate = $devent/$refalnlen;
  $total_devents+=$ievent;
  $total_alignedbases+=$matches;
  #print("$read\t$qseqlen{$read}\t$qstart\t$qend\t+\t$contig\t$refseqlen{$contig}\t$from\t$to\t$strand\t$mapq\t");printf("%0.2f\n",(1-(($insert+$del)/(($qalnlen+$refalnlen)/2))));
  my $accuracy=1-($subs_error + $del_error + $insert_error);
  print TABLE join("\t",($read_id,$readlen,$accuracy,$subs_error,$del_error,$delrate,$insert_error,$insrate, $readalnlen/$readlength,$softpct,$hardpct,)),"\n";
}

close SAM;
close TABLE;
sub strand{
  my $flag = shift;
  my $strand = '+';
  if ($flag & 0x10) {
    $strand = '-';
  }
  return $strand;		
}
