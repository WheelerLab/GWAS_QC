#!/usr/bin/perl
use warnings;
use strict;

while(<>){
    chomp;
    my ($c, $pos, $p2, $chr, $rs, $cm, $old, $a1, $a2) = split(/\s+/);
    print "$chr\t$rs\t$cm\t$pos\t$a1\t$a2\n";
}
