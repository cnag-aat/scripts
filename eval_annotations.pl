#!/usr/bin/env perl

#pseudo bed format output -> 1-based coords -need to change

use lib "/old_home/devel/talioto/myperlmods";
use Bio::Tools::GFF;
use strict;
use Getopt::Long;
use File::Path;
use Data::Dumper;
use FileHandle;
use File::Basename;
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

my $splice_window = 10;
my $exon_overlap = 10;
my $version = 3;
my $tx_start_end_near = 0;

my $ulevel = 'both';		#cds,exon
my $new_transcript = 1;
my $removestop = 0;
my $label=$$;
#my @files = grep /\.gff3/, @ARGV;
#print STDERR scalar @files,"\n";

my $aversion = 3;
my $pversion = 3;
my $source = 0;
my $annotation = "";
my $prediction = "";
my $expressed_annotation = 0;
my $evaluate=1;
my $eval_ss = 0;
my $column_label = '';
my $usage = $0;
$usage .=<<'END';
 [-ea expressed annotation] [-l label] [-ss] -p prediction -aa annotation
    -aa <annotation_file>
    -ea <expressed_annotation_file>
     -p <prediction_file>
    -ss evaluate splice sites
     -l <label>

METRIC DEFINITIONS
    TP = true positive
    FP = false positive
    FN = false negative

    Sn = sensitivity = TP/(TP+FN)
    Sp = specificity (TP + FP)/TP
    SnSp = (Sn + Sp)/2

    CTx = complete coding region (all CDS)
    ITx = "internal" transcript (all internal exons)
    Tx  = complete transcript (all exons)
    CDS = coding exon (First, Internal, or Terminal)
    Int = intron 
    Exo = exon
    IEx = internal exon
    Don = donor site
    Acc = acceptor site
    Nt  = nucleotide
    Cnt = coding nucleotide
    IeNt= nucleotides from internal exons

END
my $help = 0;
GetOptions(
	   'aa:s'      => \$annotation,
	   'ea:s'      => \$expressed_annotation,
	   'p:s'      => \$prediction,
	   'l:s'      => \$column_label,
	   'h'      => \$help,
	   #'av|version:s'      => \$aversion,
	   #'pv|version:s'      => \$pversion,
	   'eval!'    => \$evaluate,
	   'ss'       => \$eval_ss
	  );
if (!$annotation || !$prediction || $help){
    print $usage and exit;
}
my %aetranscripts; #annotated exon transcripts
my %actranscripts; #annotated cds transcripts
my %aietranscripts; #annotated interal exon transcripts
my %eaetranscripts; #expressed annotated exon transcripts
my %eactranscripts; #expressed annotated cds transcripts
my %eaietranscripts; #expressed annotated interal exon transcripts
my %petranscripts; #predicted exon transcripts
my %pctranscripts; #predicted cds transcripts
my %pietranscripts; #predicted internal exon transcripts

my $num_expressed_annotated_transcripts = undef;
my $num_annotated_transcripts = parseFile($annotation,'annotation',\%aetranscripts,\%actranscripts,\%aietranscripts,$eval_ss);
if ($expressed_annotation){
  $num_expressed_annotated_transcripts = parseFile($expressed_annotation,'expressed_annotation',\%eaetranscripts,\%eactranscripts,\%eaietranscripts,$eval_ss);
}else{
  $expressed_annotation = $annotation;
  $num_expressed_annotated_transcripts = $num_annotated_transcripts;
  %eaetranscripts = %aetranscripts;
  %eactranscripts = %actranscripts;
  %eaietranscripts = %aietranscripts;
}
exit if (!$evaluate);
my $num_predicted_transcripts = parseFile($prediction,'prediction',\%petranscripts,\%pctranscripts,\%pietranscripts,$eval_ss);
my $matching_cds = 0;
my $matching_transcripts = 0;
my $matching_itranscripts = 0;
my $matching_cds_expressed = 0;
my $matching_transcripts_expressed = 0;
my $matching_itranscripts_expressed = 0;
### Exact CDS (protein coding sequence)
foreach my $key (keys %actranscripts) {
    #print "cds:$key\n";
    $matching_cds++ if exists $pctranscripts{$key};
}
foreach my $key (keys %eactranscripts) {
    #print "cds:$key\n";
    $matching_cds_expressed++ if exists $pctranscripts{$key};
}
### All exons exact
foreach my $key (keys %aetranscripts) {
    #print "exon:$key\n";
    $matching_transcripts++ if exists $petranscripts{$key};
}
foreach my $key (keys %eaetranscripts) {
    #print "exon:$key\n";
    $matching_transcripts_expressed++ if exists $petranscripts{$key};
}
### All internal exons exact
foreach my $key (keys %aietranscripts) {
    #print "exon:$key\n";
    $matching_itranscripts++ if exists $pietranscripts{$key};
}
foreach my $key (keys %eaietranscripts) {
    #print "exon:$key\n";
    $matching_itranscripts_expressed++ if exists $pietranscripts{$key};
}
my ($abase,$apath,$aext) = fileparse($annotation,qw(.gff .gtf .gff2 .gff3));
my ($eabase,$eapath,$eaext) = fileparse($expressed_annotation,qw(.gff .gtf .gff2 .gff3));
my ($pbase,$ppath,$pext) = fileparse($prediction,qw(.gff .gtf .gff2 .gff3));
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);

printf "date\t%4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec;
print "ann\t$annotation\n";
print "exprsd_ann\t$expressed_annotation\n";
print "pred\t$prediction\n";
print "label\t$column_label\n";

printf ("CTx_Sn\t%4.1f\n",100*($matching_cds_expressed/($num_expressed_annotated_transcripts->{cds}))) if ($num_expressed_annotated_transcripts->{cds});
printf ("CTx_Sp\t%4.1f\n",100*($matching_cds/($num_predicted_transcripts->{cds}))) if ($num_predicted_transcripts->{cds});
printf ("CTx_SnSp\t%4.1f\n",100*(($matching_cds/($num_predicted_transcripts->{cds}))+($matching_cds_expressed/($num_expressed_annotated_transcripts->{cds})))/2) if ($num_expressed_annotated_transcripts->{cds} && $num_predicted_transcripts->{exon});

printf ("ITx_Sn\t%4.1f\n",100*($matching_itranscripts_expressed/($num_expressed_annotated_transcripts->{internalexon}))) if ($num_expressed_annotated_transcripts->{internalexon});
printf ("ITx_Sp\t%4.1f\n",100*($matching_itranscripts/($num_predicted_transcripts->{internalexon})))if($num_predicted_transcripts->{internalexon});
printf ("ITx_SnSp\t%4.1f\n",100*(($matching_itranscripts/($num_predicted_transcripts->{internalexon}))+($matching_itranscripts_expressed/($num_expressed_annotated_transcripts->{internalexon})))/2) if ($num_expressed_annotated_transcripts->{internalexon});

printf (" Tx_Sn\t%4.1f\n",100*($matching_transcripts_expressed/($num_expressed_annotated_transcripts->{exon}))) if ($num_expressed_annotated_transcripts->{exon});
printf (" Tx_Sp\t%4.1f\n",100*($matching_transcripts/($num_predicted_transcripts->{exon}))) if ($num_predicted_transcripts->{exon});
printf (" Tx_SnSp\t%4.1f\n",100*(($matching_transcripts/($num_predicted_transcripts->{exon}))+($matching_transcripts_expressed/($num_expressed_annotated_transcripts->{exon})))/2) if ($num_predicted_transcripts->{exon} && $num_expressed_annotated_transcripts->{exon});

my ($CDS_SPTP,@rest) = split " ",`intersectBed -a $pbase.cds4eval.uniq.bed -b $abase.cds4eval.uniq.bed -r -s -f 1 -u |wc -l`;
my ($CDS_SNTP,@rest2) = split " ",`intersectBed -a $pbase.cds4eval.uniq.bed -b $eabase.cds4eval.uniq.bed -r -s -f 1 -u |wc -l`;
my ($CDS_annotation,@resta) = split " ",`wc -l $eabase.cds4eval.uniq.bed`;
my ($CDS_prediction,@restb) = split " ",`wc -l $pbase.cds4eval.uniq.bed`;
printf ("CDS_Sn\t%4.1f\n",100*($CDS_SNTP/$CDS_annotation)) if $CDS_annotation;
printf ("CDS_Sp\t%4.1f\n",100*($CDS_SPTP/$CDS_prediction)) if $CDS_prediction;
printf ("CDS_SnSp\t%4.1f\n",100*(($CDS_SPTP/$CDS_prediction)+($CDS_SNTP/$CDS_annotation))/2) if ($CDS_annotation && $CDS_prediction);
#unlink "$pbase.cds4eval.uniq.bed";

my ($intron_SPTP,@rest) = split " ",`intersectBed -a $pbase.intron4eval.uniq.bed -b $abase.intron4eval.uniq.bed -r -s -f 1 -u |wc -l`;
my ($intron_SNTP,@rest) = split " ",`intersectBed -a $pbase.intron4eval.uniq.bed -b $eabase.intron4eval.uniq.bed -r -s -f 1 -u |wc -l`;
my ($intron_annotation,@resta) = split " ",`wc -l $eabase.intron4eval.uniq.bed`;
my ($intron_prediction,@restb) = split " ",`wc -l $pbase.intron4eval.uniq.bed`;
printf ("Int_Sn\t%4.1f\n",100*($intron_SNTP/$intron_annotation)) if $intron_annotation;
printf ("Int_Sp\t%4.1f\n",100*($intron_SPTP/$intron_prediction)) if $intron_prediction;
printf ("Int_SnSp\t%4.1f\n",100*(($intron_SPTP/$intron_prediction)+($intron_SNTP/$intron_annotation))/2) if ($intron_prediction && $intron_annotation);
#unlink "$pbase.intron4eval.uniq.bed";

my ($exon_SPTP,@restc) = split " ",`intersectBed -a $pbase.exon4eval.uniq.bed -b $abase.exon4eval.uniq.bed -r -s -f 1 -u |wc -l`;
my ($exon_SNTP,@restc) = split " ",`intersectBed -a $pbase.exon4eval.uniq.bed -b $eabase.exon4eval.uniq.bed -r -s -f 1 -u |wc -l`;
my ($exon_annotation,@restd) = split " ",`wc -l $eabase.exon4eval.uniq.bed`;
my ($exon_prediction,@reste) = split " ",`wc -l $pbase.exon4eval.uniq.bed`;
printf ("Exo_Sn\t%4.1f\n",100*($exon_SNTP/$exon_annotation));
printf ("Exo_Sp\t%4.1f\n",100*($exon_SPTP/$exon_prediction));
printf ("Exo_SnSp\t%4.1f\n",100*(($exon_SNTP/$exon_annotation)+($exon_SPTP/$exon_prediction))/2);
#unlink "$pbase.exon4eval.uniq.bed";

my ($internalexon_SPTP,@restf) = split " ",`intersectBed -a $pbase.internalexon4eval.uniq.bed -b $abase.internalexon4eval.uniq.bed -r -s -f 1 -u |wc -l`;
my ($internalexon_SNTP,@restf) = split " ",`intersectBed -a $pbase.internalexon4eval.uniq.bed -b $eabase.internalexon4eval.uniq.bed -r -s -f 1 -u |wc -l`;
my ($internalexon_annotation,@restg) = split " ",`wc -l $eabase.internalexon4eval.uniq.bed`;
my ($internalexon_prediction,@resth) = split " ",`wc -l $pbase.internalexon4eval.uniq.bed`;
printf ("IEx_Sn\t%4.1f\n",100*($internalexon_SNTP/$internalexon_annotation)) if $internalexon_annotation;
printf ("IEx_Sp\t%4.1f\n",100*($internalexon_SPTP/$internalexon_prediction)) if $internalexon_prediction;
printf ("IEx_SnSp\t%4.1f\n",100*(($internalexon_SNTP/$internalexon_annotation)+($internalexon_SPTP/$internalexon_prediction))/2) if ($internalexon_prediction && $internalexon_prediction);
#unlink "$pbase.internalexon4eval.uniq.bed";

if ($eval_ss) {
    my ($donors_SPTP,@resti) = split " ",`intersectBed -a $pbase.donors4eval.uniq.bed -b $abase.donors4eval.uniq.bed -r -s -f 1 -u |wc -l`;
    my ($donors_SNTP,@resti) = split " ",`intersectBed -a $pbase.donors4eval.uniq.bed -b $eabase.donors4eval.uniq.bed -r -s -f 1 -u |wc -l`;
    my ($donors_annotation,@restj) = split " ",`wc -l $eabase.donors4eval.uniq.bed`;
    my ($donors_prediction,@restk) = split " ",`wc -l $pbase.donors4eval.uniq.bed`;
    printf ("Don_Sn\t%4.1f\n",100*($donors_SNTP/$donors_annotation)) if $donors_annotation;
    printf ("Don_Sp\t%4.1f\n",100*($donors_SPTP/$donors_prediction)) if $donors_prediction;
    printf ("Don_SnSp\t%4.1f\n",100*(($donors_SPTP/$donors_prediction)+($donors_SNTP/$donors_annotation))/2) if ($donors_annotation && $donors_prediction);
    #unlink "$pbase.donors4eval.uniq.bed";

    my ($acceptors_SPTP,@restl) = split " ",`intersectBed -a $pbase.acceptors4eval.uniq.bed -b $abase.acceptors4eval.uniq.bed -r -s -f 1 -u |wc -l`;
    my ($acceptors_SNTP,@restl) = split " ",`intersectBed -a $pbase.acceptors4eval.uniq.bed -b $eabase.acceptors4eval.uniq.bed -r -s -f 1 -u |wc -l`;
    my ($acceptors_annotation,@restm) = split " ",`wc -l $eabase.acceptors4eval.uniq.bed`;
    my ($acceptors_prediction,@restn) = split " ",`wc -l $pbase.acceptors4eval.uniq.bed`;
    printf ("Acc_Sn\t%4.1f\n",100*($acceptors_SNTP/$acceptors_annotation)) if $acceptors_annotation;
    printf ("Acc_Sp\t%4.1f\n",100*($acceptors_SPTP/$acceptors_prediction)) if $acceptors_prediction;
    printf ("Acc_SnSp\t%4.1f\n",100*(($acceptors_SPTP/$acceptors_prediction)+($acceptors_SNTP/$acceptors_annotation))/2) if ( $acceptors_prediction && $acceptors_annotation);
    #unlink "$pbase.acceptors4eval.uniq.bed";
}
my ($ntsp,$nt_SNTP,@resto) = split " ",`coverageBed -a $pbase.exon4eval.uniq.projected.bed -b $eabase.exon4eval.uniq.projected.bed | gawk '{total+=\$6;cov+=\$5;}END{print cov/total "\t" cov;}'`;
my ($ntsp_all,$nt_SPTP,@resto2) = split " ",`coverageBed -a $pbase.exon4eval.uniq.projected.bed -b $abase.exon4eval.uniq.projected.bed | gawk '{total+=\$6;cov+=\$5;}END{print cov/total "\t" cov;}'`;
my $nt_pred_total = 0;
open(BED,"<$pbase.exon4eval.uniq.projected.bed") or die $!;while(<BED>){chomp;my @f=split "\t",$_;$nt_pred_total+=($f[2]-$f[1]);}close BED;
my $nt_aa_total = 0;
open(BED,"<$abase.exon4eval.uniq.projected.bed") or die $!;while(<BED>){chomp;my @f=split "\t",$_;$nt_aa_total+=($f[2]-$f[1]);}close BED;
my $ntsn = $nt_SNTP/$nt_aa_total;
printf ("Nt_Sn\t%4.1f\n",100*$ntsn);
printf ("Nt_Sp\t%4.1f\n",100*$ntsp_all);
printf ("Nt_SnSp\t%4.1f\n",100*($ntsn+$ntsp_all)/2);
#unlink "$pbase.exon4eval.uniq.projected.bed";

#print STDERR "coverageBed -a $pbase.cds4eval.uniq.projected.bed -b $eabase.cds4eval.uniq.projected.bed | gawk '{total+=\$6;cov+=\$5;}END{print cov/total ",'"\t"'," cov;}'\n";
my ($cdsntsp,$cdsnt_SNTP,@restp) = split " ",`coverageBed -a $pbase.cds4eval.uniq.projected.bed -b $eabase.cds4eval.uniq.projected.bed | gawk '{total+=\$6;cov+=\$5;}END{print cov/total "\t" cov;}'`;
my ($cdsntsp_all,$cdsnt_SPTP,@restp2) = split " ",`coverageBed -a $pbase.cds4eval.uniq.projected.bed -b $abase.cds4eval.uniq.projected.bed | gawk '{total+=\$6;cov+=\$5;}END{print cov/total "\t" cov;}'`;
my $cdsnt_pred_total = 0;
open(BED,"<$pbase.cds4eval.uniq.projected.bed") or die $!;while(<BED>){chomp;my @f=split "\t",$_;$cdsnt_pred_total+=($f[2]-$f[1]);}close BED;
my $cdsnt_aa_total = 0;
open(BED,"<$abase.cds4eval.uniq.projected.bed") or die $!;while(<BED>){chomp;my @f=split "\t",$_;$cdsnt_aa_total+=($f[2]-$f[1]);}close BED;
my $cdsntsn = $cdsnt_SNTP/$cdsnt_aa_total;
printf ("CNt_Sn\t%4.1f\n",100*$cdsntsn);
printf ("CNt_Sp\t%4.1f\n",100*$cdsntsp_all);
printf ("CNt_SnSp\t%4.1f\n",100*($cdsntsn+$cdsntsp_all)/2);
#unlink "$pbase.cds4eval.uniq.projected.bed";

my ($internalexonntsp,$internalexonnt_SPTP,@restq) = split " ",`coverageBed -a $pbase.internalexon4eval.uniq.projected.bed -b $eabase.internalexon4eval.uniq.projected.bed | gawk '{total+=\$6;cov+=\$5;}END{print cov/total "\t" cov;}'`;
my ($internalexonntsp_all,$internalexonnt_SPTP,@restq2) = split " ",`coverageBed -a $pbase.internalexon4eval.uniq.projected.bed -b $abase.internalexon4eval.uniq.projected.bed | gawk '{total+=\$6;cov+=\$5;}END{print cov/total "\t" cov;}'`;
my $internalexonnt_pred_total = 0;
open(BED,"<$pbase.internalexon4eval.uniq.projected.bed") or die $!;while(<BED>){chomp;my @f=split "\t",$_;$internalexonnt_pred_total+=($f[2]-$f[1]);}close BED;
my $internalexonnt_aa_total = 0;
open(BED,"<$abase.internalexon4eval.uniq.projected.bed") or die $!;while(<BED>){chomp;my @f=split "\t",$_;$internalexonnt_aa_total+=($f[2]-$f[1]);}close BED;
my $internalexonntsn = $nt_SNTP/$internalexonnt_aa_total;
printf ("IeNt_Sn\t%4.1f\n",100*$internalexonntsn);
printf ("IeNt_Sp\t%4.1f\n",100*$internalexonntsp_all);
printf ("IeNt_SnSp\t%4.1f\n",100*($internalexonntsn+$internalexonntsp_all)/2);
#unlink "$pbase.internalexon4eval.uniq.projected.bed";
#`rm $pbase.*.bed`;

sub parseFile{
    my $file = shift;
    my $category = shift;
    my $exontranscripts = shift;
    my $cdstranscripts = shift;
    my $internalexontranscripts = shift;
    my $evaluate_splicesites = shift;
    my ($base,$path,$ext) = fileparse($file,qw(.gff .gtf .gff2 .gff3));
    #my ($predbase,$predpath,$predext) = fileparse($prediction,qw(.gff .gtf .gff2 .gff3));
    my $cdsFileHandle;
    my $exonFileHandle;
    my $internalExonFileHandle;
    my $intronFileHandle;
    my $donorFileHandle;
    my $acceptorFileHandle;
    $cdsFileHandle = FileHandle->new(">$base.cds4eval.bed") if !-e "$base.cds4eval.bed";
    $exonFileHandle = FileHandle->new(">$base.exon4eval.bed") if !-e "$base.exon4eval.bed";
    $internalExonFileHandle = FileHandle->new(">$base.internalexon4eval.bed") if !-e "$base.internalexon4eval.bed";
    $intronFileHandle = FileHandle->new(">$base.intron4eval.bed") if !-e "$base.intron4eval.bed";
    if ($evaluate_splicesites) {
		$donorFileHandle = FileHandle->new(">$base.donors4eval.bed") if !-e "$base.donors4eval.bed";
		$acceptorFileHandle = FileHandle->new(">$base.acceptors4eval.bed") if !-e "$base.acceptors4eval.bed";
    }
    my $idcount = 1;
    my $totalnumber = 0;
    my $cdsunumber = 0;
    my $exonunumber = 0;
    my $internalexonunumber = 0;
    my %current;
	  my %all_transcripts = ();

    open FILE,"<$file" or die "couldn't open $category file $file:$!\n";
    #print STDERR "$file\n";
    while (<FILE>) {
		### skip sequence block if it exists
		if (m/^##FASTA/) {
			while (my $line = <>) {
				last if $line =~ m/^##gff-version/;
			}
		}elsif (m/^#/) {
			#skip comments and directives, including transcript/gene separater ###
		}elsif ($_ !~ m/\S/) {
			#skip comments and directives, including transcript/gene separater ###
		}else{
			my %record;
			makeFeature($_,\%record);
			my $transcript_id = '';
			if ($record{feature} eq 'transcript' || $record{feature} eq 'mRNA'){
				$transcript_id = $record{attributes}->{'ID'};
				$all_transcripts{$transcript_id}->{transcript}=\%record;
			}else{
				if ($record{feature} eq 'CDS' || $record{feature} eq 'exon' || $record{feature} =~ /utr/){
					$transcript_id = $record{attributes}->{'Parent'};
				}
				if ($record{feature} ne 'gene') {
					push @{$all_transcripts{$transcript_id}->{gff}},[split "\t",$_];
				}
	    		if ($record{feature} eq 'CDS') {
					push(@{$all_transcripts{$transcript_id}->{cds}},\%record);
	    		}
	    		if ($record{feature} eq 'exon') {
					push(@{$all_transcripts{$transcript_id}->{exons}},\%record);
					$all_transcripts{$transcript_id}->{exon_starts}->{$record{start}}++;
					$all_transcripts{$transcript_id}->{exon_ends}->{$record{end}}++;
	    		}
	    		if ($record{feature} eq 'UTR') {
					push(@{$all_transcripts{$transcript_id}->{utrs}},\%record);
					$all_transcripts{$transcript_id}->{utr_starts}->{$record{start}}++;
					$all_transcripts{$transcript_id}->{utr_ends}->{$record{end}}++;
	    		}
			}
		}
	}
	foreach my $t (sort keys %all_transcripts){
		my %current = %{$all_transcripts{$t}};
		$totalnumber++;
	    my $cdsustring = makeUniqueString(\%current,'cds');
	    
	    my $exonustring = makeUniqueString(\%current,'exon');
	    my $internalexonustring = makeUniqueString(\%current,'internalexon');
	    if (defined $cdsustring){
			#print STDERR "$cdsustring\n";
			if (!exists $cdstranscripts->{$cdsustring}) {
				$cdstranscripts->{$cdsustring}++;
				$cdsunumber++;
			}
			if (defined $cdsFileHandle){printCDS(\%current,$cdsFileHandle);}
	    }
	    if (defined $exonustring){
			if (!exists $exontranscripts->{$exonustring}) {
				$exontranscripts->{$exonustring}++;
				$exonunumber++;
			}
			if (defined $exonFileHandle){ printExons(\%current,$exonFileHandle);}
			if (defined $internalExonFileHandle){printInternalExons(\%current,$internalExonFileHandle);}
			if (defined $intronFileHandle){printIntrons(\%current,$intronFileHandle,$donorFileHandle,$acceptorFileHandle);}
	    }
	    if (defined $internalexonustring){
			#print STDERR "$cdsustring\n";
			if (!exists $internalexontranscripts->{$internalexonustring}) {
				$internalexontranscripts->{$internalexonustring}++;
				$internalexonunumber++;
			}
	    }



	}


#     print STDERR "total $category: $totalnumber\nunique $category CDS: $cdsunumber\n";
#     print STDERR "unique $category transcripts: $exonunumber\n";
#     print STDERR "unique $category itranscripts: $internalexonunumber\n";
    `sort -u -k1,1 -k2,3n -k6,6  $base.cds4eval.bed > $base.cds4eval.uniq.bed`;# if !-e "$base.cds4eval.uniq.bed";
    `sort -u -k1,1 -k2,3n -k6,6  $base.exon4eval.bed > $base.exon4eval.uniq.bed`;# if !-e "$base.exon4eval.uniq.bed";
    `sort -u -k1,1 -k2,3n -k6,6  $base.intron4eval.bed > $base.intron4eval.uniq.bed`;# if !-e "$base.intron4eval.uniq.bed";
    if($evaluate_splicesites){
	`sort -u -k1,1 -k2,3n -k6,6  $base.donors4eval.bed > $base.donors4eval.uniq.bed`;# if !-e "$base.donors4eval.uniq.bed";
	`sort -u -k1,1 -k2,3n -k6,6 $base.acceptors4eval.bed > $base.acceptors4eval.uniq.bed`;# if !-e "$base.acceptors4eval.uniq.bed";
    }
    `sort -u -k1,1 -k2,3n -k6,6 $base.internalexon4eval.bed > $base.internalexon4eval.uniq.bed`;# if !-e "$base.internalexon4eval.uniq.bed";
    `mergeBed -i $base.internalexon4eval.uniq.bed > $base.internalexon4eval.uniq.projected.bed`;# if !-e "$base.internalexon4eval.uniq.projected.bed";
    `mergeBed -i $base.cds4eval.uniq.bed > $base.cds4eval.uniq.projected.bed`;#  if !-e "$base.cds4eval.uniq.projected.bed";
    `mergeBed -i $base.exon4eval.uniq.bed > $base.exon4eval.uniq.projected.bed`;#  if !-e "$base.exon4eval.uniq.projected.bed";

    #unlink "$base.cds4eval.bed";
    #unlink "$base.exon4eval.bed";
    #unlink "$base.intron4eval.bed";
    #unlink "$base.donors4eval.bed";
    #unlink "$base.acceptors4eval.bed";
    #unlink "$base.internalexon4eval.bed";
    print STDERR "$file cdsunumber $cdsunumber\n";
    return {'cds'=>$cdsunumber,'exon'=>$exonunumber,'internalexon'=>$internalexonunumber};
    #print STDERR "Sorting transcripts\n";
}


sub transcriptsort
  {
      $a->{transcript}->{seqname} cmp $b->{transcript}->{seqname}
	||
	  $a->{transcript}->{strand} <=> $b->{transcript}->{strand}
	    ||
	      $a->{transcript}->{start} <=> $b->{transcript}->{start}
		||
		  $b->{transcript}->{end} <=> $a->{transcript}->{end}
	      }
sub makeUniqueString{
    my $t = shift;
    my $u = shift;
    my $ustring = undef;
    #print STDERR $t->{cds}->[0],"\n";;
    if ($u =~ m/(cds|both)/i && defined $t->{cds}) {
		my @CDS = sort {$a->{start} <=> $b->{start}} @{$t->{cds}};
		foreach my $c (@CDS) {
			if (defined $c){
				if (defined $c->{seqname}){
					$ustring .= join(':',($c->{seqname},$c->{start},$c->{end},$c->{strand})).':';
					#print STDERR "cds:$ustring\n";
				}
			}
		}
    }
    if ($u =~ m/^(exon|both)/i && defined $t->{exons}) {
		my @exons = sort {$a->{start} <=> $b->{start}} @{$t->{exons}};
		foreach my $e (@exons) {
			if (defined $e){
				if (defined $e->{seqname}){
					$ustring .= join(':',($e->{seqname},$e->{start},$e->{end},$e->{strand})).':';
					#print STDERR "exon:$ustring\n";
				}
			}
		}
    }
    if ($u =~ m/^internalexon/i && defined $t->{exons}) {
		my @exons = sort {$a->{start} <=> $b->{start}} @{$t->{exons}};
		if (scalar @exons > 2){
			shift @exons;
			pop @exons;
			foreach my $e (@exons) {
				if (defined $e) {
					if (defined $e->{seqname}) {
					$ustring .= join(':',($e->{seqname},$e->{start},$e->{end},$e->{strand})).':';
					#print STDERR "exon:$ustring\n";
					}
				}
			}
		}
    }
    #print STDERR "final:$ustring\n";
    return $ustring;
}
sub makeFeature{
    my $line = shift;
    my $r = shift;
    #print STDERR "$line";
    #print STDERR ;
    chomp $line;
    my @gff = split("\t",$line);
    #foreach my $f(@gff){print "$f\n";}
    #print STDERR $gff[$feature],"\t";
    $r->{seqname}=$gff[$seqname];
    $r->{source}=$gff[$source];
    $r->{feature}=$gff[$feature];
    $r->{start}=$gff[$start];
    $r->{end}=$gff[$end];
    $r->{score}=$gff[$score];
    $r->{strand}=$gff[$strand];
    $r->{frame}=$gff[$frame];
    $r->{group}=$gff[$group];
    my %att = split(/[=;]/,$gff[$group]);
    $r->{attributes}=\%att;
    $r->{gffstring}=$line;
}



sub gffsort
  {
      $a->{CDS}->[0]->seq_id cmp $b->{CDS}->[0]->seq_id
	||
	  $a->{CDS}->[0]->start <=> $b->{CDS}->[0]->start
	    ||
	      $a->{CDS}->[0]->end <=> $b->{CDS}->[0]->end
		||
		  $a->{CDS}->[0]->strand <=> $b->{CDS}->[0]->strand
	      }
sub merge_clusters
  {
      my $clusterhashref = shift;
      my $clusterin = shift;
      my @inclust = @$clusterin;
      my $res_clust = $inclust[0];
      #print STDERR join("\t",@inclust),"\n";
      if (scalar (@inclust < 2)) {
		#print STDERR "RESULT CLUSTER $res_clust\n";
		return $res_clust;
      } else {
		for (my $i = 1;$i< scalar @inclust;$i++) {
			foreach my $transcript (@{$clusterhashref->{$inclust[$i]}->{transcripts}}) {
				push (@{$clusterhashref->{$res_clust}->{transcripts}},$transcript);
				#push (@{$clusterhashref->{$res_clust}->{transcripts}},$transcript);
				if (defined $transcript->{cds}) {
					push (@{$clusterhashref->{$res_clust}->{cds}},@{$transcript->{cds}});
				}
				if (defined $transcript->{exons}) {
					push (@{$clusterhashref->{$res_clust}->{exons}},@{$transcript->{exons}});
				}
				foreach my $start (sort keys %{$clusterhashref->{$inclust[$i]}->{starts}}) {
					$clusterhashref->{$res_clust}->{starts}->{$start}++;
				}
				foreach my $end (sort keys %{$clusterhashref->{$inclust[$i]}->{ends}}) {
					$clusterhashref->{$res_clust}->{ends}->{$end}++;
				}
				delete($clusterhashref->{$inclust[$i]});
			}
		}
		#print STDERR "RESULT CLUSTER $res_clust\n";
		return $res_clust;	
      }
  }
sub printCDS{
    my $t = shift;
    my $cdsFh = shift;
    return if !defined $t->{cds};
    my @CDS = sort {$a->{start} <=> $b->{start}} @{$t->{cds}};
    foreach my $c (@CDS) {
	print $cdsFh join("\t",$c->{seqname},$c->{start}-1,$c->{end},$t->{id},$c->{score},$c->{strand}),"\n";
    }
}
sub printExons{
    my $t = shift;
    my $exonsFh = shift;
    return if !defined $t->{exons};
    my @exons = sort {$a->{start} <=> $b->{start}} @{$t->{exons}};
    foreach my $e (@exons) {
	print $exonsFh join("\t",$e->{seqname},$e->{start}-1,$e->{end},$t->{id},$e->{score},$e->{strand}),"\n";
    }
}
sub printInternalExons{
    my $t = shift;
    my $internalExonsFh = shift;
    return if !defined $t->{exons};
    my @exons = sort {$a->{start} <=> $b->{start}} @{$t->{exons}};
    if (scalar @exons > 2) {
		shift @exons;
		pop @exons;
		foreach my $e (@exons) {
			print $internalExonsFh join("\t",$e->{seqname},$e->{start}-1,$e->{end},$t->{id},$e->{score},$e->{strand}),"\n";
		}
    }
}
sub printIntrons{
    my $t = shift;
    my $intronFh = shift;
    my $donorFh = shift;
    my $acceptorFh = shift;
    return if !defined $t->{exons};
    my @exons = sort {$a->{start} <=> $b->{start}} @{$t->{exons}};
    my $n = 0;
    my $END = 0;
    if (scalar @exons > 1) {
	foreach my $e (@exons) {
	    $n++;
	    if ($n != 1) {
			die "previous exon end $END is greater than current exon start ".$e->{start}."\n" if $END > $e->{start}-1; 
			print  $intronFh join("\t",$e->{seqname},$END +1 -1,$e->{start}-1,$t->{id},'.',$e->{strand}),"\n";
			if (defined $donorFh && defined $acceptorFh){
				print  $donorFh join("\t",$e->{seqname},$END -1,$END+1,$t->{id},'.',$e->{strand}),"\n";
				print  $acceptorFh join("\t",$e->{seqname},$e->{start}-1,$e->{start},$t->{id},'.',$e->{strand}),"\n";
			}
	    }
	    $END = $e->{end};
	}
    }
}
