#!/usr/bin/perl
use warnings;
use strict;

##convert genotypes from RIKEN to tped/tfam files for plink

##get SNP position info
open(A, "/group/dolan-lab/hwheeler/ThePlatinumStudy/RIKEN_data/HumanOmniExpressExome-8v1-2_A.csv");
my %pos;
while(<A>){
    chomp;
    my (@data) = split(/,/);
    my $name = $data[1];
    my $chr = $data[9];
    my $pos = $data[10];
    my $info = $chr . "," . $pos;
    $pos{$name} = $info;
}

##get Sample ID info
open(B, "/group/dolan-lab/hwheeler/ThePlatinumStudy/RIKEN_genotyping/RIKEN_layout_20150624_FINAL.txt");
my %sample;
while(<B>){
    chomp;
    my ($PATNO, $pid, $info, $plateid, $ptid, $pat, $mat, $sex) = split(/\t/);
    my $platepos = "N88_" . $plateid;
    $sample{$platepos} = $PATNO . " 0 0 " . $sex; ##don't include hapmap relative information since FID are not the same
}

#for testing: open(C, "/group/dolan-lab/hwheeler/ThePlatinumStudy/RIKEN_data/head_N88_Recluster_TOP_20150911_FinalReport.csv"); 
open(C, "/group/dolan-lab/hwheeler/ThePlatinumStudy/RIKEN_data/FinalReport/N88_Recluster_TOP_20150911_FinalReport.csv");
open(TPED, ">genotypes/N88_Recluster_TOP_20150911_FinalReport.tped");
open(TFAM, ">genotypes/N88_Recluster_TOP_20150911_FinalReport.tfam");

while(<C>){
    chomp;
    my ($snp, @rest) = split(/,/);
    if($snp eq ""){
        foreach my $ind (@rest){
	    my $patno = $sample{$ind};
	    print TFAM "$ind $patno -9\n";
        }
	next;
    }
    if(defined($pos{$snp})){
        my ($c,$bp) = split(/,/,$pos{$snp});
        print TPED "$c $snp 0 $bp"; 
	my $a1;
	my $a2;
	foreach my $id (@rest){
	    if($id =~ m/(\w)(\w)/){
		$a1 = $1;
		$a2 = $2;
	    }else{
		$a1 = 0;
		$a2 = 0;
	    }
	    print TPED " $a1 $a2";
	}
	print TPED "\n";
    }
}    
