use warnings;
use strict;

while(<>){
    chomp;
    my ($chr, $rs, $cm, $bp, $a1, $a2) = split(/\s+/);
    if($chr >= 1 && $chr <= 24){
	print "$rs\n";
    }
}
