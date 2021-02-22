#!/usr/bin/env perl
my $max = shift @ARGV;
my @lines = ();
my %seen;
my $count = 0;
while (<STDIN>) {
	$count++;
   	push @lines,$_;
}	

my $tot = scalar @lines;

my $c = 0;
print $lines[0];
$seen{0}++;
while ($c < $max){
	my $r = int rand $tot; 
	if (!(exists $seen{$r})){
		print $lines[$r];
		$c++
	}
	$seen{$r}++;
}
