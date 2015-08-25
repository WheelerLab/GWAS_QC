use warnings;
use strict;

###flags sample only if sex mis-match, not if missing X-chr data###

while(<>){
    chomp;
    my ($sp, $FID, $IID, $PEDSEX, $SNPSEX, $STATUS, $F) = split(/\s+/);
    if($STATUS eq "PROBLEM" && $SNPSEX != 0){
	print "$FID\t$IID\t$PEDSEX\t$SNPSEX\t$STATUS\t$F\n";
    }
}
