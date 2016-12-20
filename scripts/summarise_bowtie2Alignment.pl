#!/usr/bin/perl

# Author: Maurits Evers (maurits.evers@anu.edu.au)

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
	my $nTotal = 0;
	my $nMapped = 0;
	my $nStep = 0;
	while (my $line = <$fh>) {
	    if ($line =~ /(\d+) reads; of these:/) {
		$nTotal = $1;
		$nStep++;
	    }
	    if ($line =~ /(\d+)\s.*aligned concordantly exactly 1 time/) {
		$nMapped += $1;
		$nStep++;
	    }
	    if ($line =~ /(\d+)\s.*aligned concordantly >1 time/) {
		$nMapped += $1;
		$nStep++;
	    }
	    my @arr = ($nTotal, $nMapped);
	    $ht{$id} = [@arr];
	    last if ($nStep == 3);
	}
	close($fh);
    }
}

my @cn = ("Sample", "Number of reads", "Number of mapped reads", "Percentage mapped");
printf("%s\n", join("\t", @cn));
foreach my $key (sort {$a cmp $b} keys %ht) {
    my @arr = @{$ht{$key}};
    printf("%s\n", join("\t", $key, $arr[0] * 2, $arr[1] * 2, sprintf("%3.2f", $arr[1] / $arr[0] * 100.0)));
}


sub usage {
    printf("summarise_bowtie2Alignment version %s by Maurits Evers (maurits.evers\@anu.edu.au)\n", $version);
    printf("Summarise bowtie2 alignment stats.\n");
    printf("Usage:\n");
    printf("  summaris_bowtie2Alignment.pl <dir> <filter>\n");
    printf("\n");
    printf("  <dir>     Directory that contains picard-tools files.\n");
    printf("  <filter>  Regular expression to filter files.\n"); 
   exit; 
}
