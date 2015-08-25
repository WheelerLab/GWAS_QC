use warnings;
use strict;

my $bfile = $ARGV[0];
my $num = $ARGV[1];

print "genotypename: $bfile.ped\n";
print "snpname: $bfile.map\n";
print "indivname: $bfile.fam\n";
print "evecoutname: $bfile.evec\n";
print "evaloutname: $bfile.eval\n";
print "outliername: $bfile.outlier\n";
print "numoutevec: 10\n";
print "numoutlieriter: $num\n";
print "numoutlierevec: 2\n";
print "outliersigmathresh: 6\n";

