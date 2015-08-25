use warnings;
use strict;

open(A, "T13b.HapMap3.euro.keep.list");
my %keep;
while(<A>){
    chomp;
    my ($id) = split(/\s+/);
    $keep{$id} = 1;
}

open(B, "T13b.HapMap3.forPCA.fam.placeholder");
open(OUT, ">T13b.fullid.HapMap3.euro.keep.list");

while(<B>){
    chomp;
    my ($num, $fid, $iid) = split(/\s+/);
    if(defined($keep{$num})){
	print OUT "$fid $iid\n";
    }
}
