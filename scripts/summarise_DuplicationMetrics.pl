#!/usr/bin/perl

use warnings;
use strict;

my $version = 0.9;

usage() if (scalar(@ARGV) != 2);

my $dir = $ARGV[0];
my $filter = $ARGV[1];
opendir(my $dh, $dir) or doe $!;

my %ht = ();
while (my $file = readdir($dh)) {
    if ($file =~ /$filter/) {
	my $id = $file;
	$file = join("/", $dir, $file);
	open(my $fh, $file) or die $!;
	while (my $line = <$fh>) {
	    if ($line =~ /PERCENT_DUPLICATION/) {
		my @arr = split("\t", <$fh>);
		$ht{$id} = [@arr];
	    }
	}
	close($fh);
    }
}

my @cn = ("Sample", "Read pairs", "Duplicates", "Percent dupes");
printf("%s\n", join("\t", @cn));
foreach my $key (sort {$a cmp $b} keys %ht) {
    my @arr = @{$ht{$key}};
    printf("%s\n", join("\t", $key, $arr[2], $arr[6], $arr[8]));
    
}


sub usage {
    printf("summarise_DuplicationMetrics version %s by Maurits Evers (maurits.evers\@anu.edu.au)\n", $version);
    printf("Summarise picard-tools MarkDuplicates files.\n");
    printf("Usage:\n");
    printf("  summaris_DuplicationMetrics.pl <dir> <filter>\n");
    printf("\n");
    printf("  <dir>     Directory that contains picard-tools files.\n");
    printf("  <filter>  Regular expression to filter files.\n"); 
   exit; 
}
