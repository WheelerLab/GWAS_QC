#!/usr/bin/perl
use warnings;
use strict;

###usage: replaceXYMT.pl HM3_ASN_CEU_YRI_Unrelated_hg18.coords

open(A, "$ARGV[0]");

my $outfile = $ARGV[0] . ".XYM";

open(OUT, ">$outfile");
my %fail;
while(<A>){
    chomp;
    my ($c, $p1, $p2) = split(/\s+/);
    if($c eq "chr23"){
	$c =~ s/$c/chrX/;
    }elsif($c eq "chr24"){
	$c =~ s/$c/chrY/;
    }elsif($c eq "chr25"){
	$c =~ s/$c/chrXY/;
    }elsif($c eq "chr26"){
	$c =~ s/$c/chrM/;
    }
    print OUT "$c $p1 $p2\n";
}
