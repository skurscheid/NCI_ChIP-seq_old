#!/usr/bin/perl

use warnings;
use strict;

my $version = 0.9;

usage() if (scalar(@ARGV) != 1);

my $fn = $ARGV[0];
open(my $fh, $fn) or die(sprintf("[ERROR] Could not open %s.", $fn));

while (<$fh>) {
    chomp($_);
    my @arr_bed = split("\t", $_);
    my @arr_gtf = ($arr_bed[0], "macs2", "peak", $arr_bed[1] + 1, $arr_bed[2]);
    push(@arr_gtf, ($arr_bed[4], $arr_bed[5], "."));
    my @id = sprintf("peak_id \"%s\"", $arr_bed[3]);
    push(@id, sprintf("fold_enrichment \"%s\"", $arr_bed[6]));
    push(@id, sprintf("log10_qval \"%4.3f\"", -$arr_bed[8]));
    push(@id, sprintf("rel_summit_pos \"%i\"", $arr_bed[9]));
    push(@arr_gtf, join("; ", @id));
    printf("%s\n", join("\t", @arr_gtf));
}

sub usage {
    printf("bed2gtf.pl version %s by Maurits Evers (maurits.evers\@anu.edu.au)\n", $version);
    printf("Convert BED-formatted annotation file to GTF-formatted file.\n");
    printf("Usage:\n");
    printf("  bed2gtf.pl <BED>\n");
    printf("\n");
    printf("  <BED>     BED-formatted input file.\n");
   exit; 
}
