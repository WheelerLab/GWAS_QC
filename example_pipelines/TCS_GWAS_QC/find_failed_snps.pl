#!/usr/bin/perl
use warnings;
use strict;

###usage: liftOver_sucesses.pl *.failures *.coords.bim outfile

open(A, "$ARGV[0]");
my %fail;
while(<A>){
    chomp;
    my ($chr, $pos) = split(/\s+/);
    if($chr ne "#Deleted"){
	my $info = $chr . $pos;   
	$fail{$info} = $pos;
    }
}

open(B, "$ARGV[1]");
open(OUT, ">$ARGV[2]");

while(<B>){
    chomp;
    my ($chr, $pos1, $pos2, $c, $rs) = split(/\s+/);
    my $in = $chr . $pos1;
    if(defined($fail{$in})){
	print OUT "$rs\n";
    }
}
