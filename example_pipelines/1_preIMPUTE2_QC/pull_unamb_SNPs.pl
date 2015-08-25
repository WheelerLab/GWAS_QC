use warnings;
use strict;

###retrieves list of non-A/T or C/G (unambiguous strand) SNPs from a *bim file###
###use: perl 8_pull_unamb_SNPs.pl GEUVADIS.SNP.bim > GEUVADIS.unamb.SNPlist


while(<>){
    chomp;
    my ($chr, $snp, $cm, $bp, $a1, $a2) = split(/\s+/);	
    if($a1 eq "A" && $a2 eq "T"){
	next;
    }elsif($a1 eq "T" && $a2 eq "A"){
	next;
    }elsif($a1 eq "C" && $a2 eq "G"){
	next;
    }elsif($a1 eq "G" && $a2 eq "C"){
	next;
    }else{
#	print "$chr\t$snp\t$cm\t$bp\t$a1\t$a2\n";
	print "$snp\n";
    }
}
