#!/usr/bin/env perl
use Getopt::Long;
use File::Copy;
my $dir = 0;
my $op = 'link';
my $num = 4000;
my $usage = "usage: $0 -d <source dir> -o <operation:cp,ln,mv> -n <number of files per subdir: default=4000\n";
GetOptions(
	   'd|dir:s'      => \$dir,
	   'o|operation:s'      => \$op,
	   'n:s'            => \$num
	   );
print $usage and exit if !$dir;
opendir(my $dh, $dir) || die "Can't opendir $dir: $!";
my @nondots = grep /^[^\.]/,readdir($dh);
closedir $dh;
my $subdir=0;
mkdir($subdir,0755);
my $count = 1;
print STDERR "Status: $op files to $subdir ...\n";
foreach my $file (@nondots){
  #print STDERR "$dir/$file\n";
  if($op =~ /l/){
    symlink("$dir/$file","$subdir/$file");
  }elsif($op =~ /m/){
    move("$dir/$file","$subdir/$file") or die "Move failed: $!";
  }elsif($op =~ /c/){
    copy("$dir/$file","$subdir/$file") or die "Copy failed: $!";
  }
  #++$count;
  if($count==$num){
    $count = 1;
    $subdir++;
    print STDERR "Status: $op files to $subdir ...\n";
    mkdir($subdir,0755);
  }else{
    $count++;
  }	
}
