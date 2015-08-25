use warnings;
use strict;

my $bfile = $ARGV[0];
my $num = $ARGV[1];

print "genotypename: $bfile.ped\n";
print "snpname: $bfile.map\n";
print "indivname: $bfile.fam\n";
print "evecoutname: $bfile.rmout.evec\n";
print "evaloutname: $bfile.rmout.eval\n";
print "outliername: $bfile.rmout.outlier\n";
print "numoutevec: 10\n";
print "numoutlieriter: $num\n";
print "numoutlierevec: 2\n";
print "outliersigmathresh: 6\n";

