#!/bin/bash

###################################################################
#    THE PLATINUM STUDY RIKEN SET 2 GWAS QC WORKFLOW              #
#    Heather E. Wheeler 2017-01-21                                #
###################################################################

##load software on tarbell
module load R
module load plink/1.90 ##plink2
module load vcftools
module load tabix/0.2.6

$DIR = /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/
cd $DIR

#1. Convert tped/tfam made in 01_make_tped.pl to PLINK binary, set any heterozygous haploid genotypes (male X,Y,MT) to missing
plink --tfile F13_Recluster_TOP_20161226_FinalReport --set-hh-missing --make-bed --out F13_Recluster_TOP_20161226_FinalReport
### 474 people (471 males, 3 females) loaded from .fam.
### Warning: 37109 het. haploid genotypes present (see
### F13_Recluster_TOP_20161226_FinalReport.hh ).
### Warning: Nonmissing nonmale Y chromosome genotype(s) present.
### Total genotyping rate is 0.998724.

#2. Check sex and calculate call rates for flagging poorly called SNPs and individuals
plink --bfile F13_Recluster_TOP_20161226_FinalReport --check-sex --missing --out F13_Recluster_TOP_20161226_FinalReport
### --missing: Sample missing data report written to
### F13_Recluster_TOP_20161226_FinalReport.imiss, and variant-based missing data
### report written to F13_Recluster_TOP_20161226_FinalReport.lmiss.
### 958497 variants and 474 people pass filters and QC.
### --check-sex: 17864 Xchr and 0 Ychr variant(s) scanned, no problems detected.
### Report written to F13_Recluster_TOP_20161226_FinalReport.sexcheck .

#3. Recalculate individual call rates after removing SNPs with call rates <99%
plink --bfile F13_Recluster_TOP_20161226_FinalReport --geno 0.01 --make-bed --out F13_Recluster_TOP_20161226_FinalReport.geno0.01
### 8725 variants removed due to missing genotype data (--geno).
### 949772 variants and 474 people pass filters and QC.
plink --bfile N88_Recluster_TOP_20150911_FinalReport.geno0.01 --missing --out N88_Recluster_TOP_20150911_FinalReport.geno0.01
### Total genotyping rate is 0.999721.
### --missing: Sample missing data report written to
### F13_Recluster_TOP_20161226_FinalReport.geno0.01.imiss, and variant-based
### missing data report written to
### F13_Recluster_TOP_20161226_FinalReport.geno0.01.lmiss.

### looks great, all individuals now have >99.4% call rates (see 03_GWAS_QC_plots.Rmd output)

#4. Calculate HWE statistics to flag SNPs later
plink --bfile F13_Recluster_TOP_20161226_FinalReport.geno0.01 --hardy --out F13_Recluster_TOP_20161226_FinalReport.geno0.01
### remove SNPs with HWE P<1e-06
plink --bfile F13_Recluster_TOP_20161226_FinalReport.geno0.01 --hwe 1e-06 --make-bed --out F13_Recluster_TOP_20161226_FinalReport.geno0.01.hwe1e-06
### --hwe: 566 variants removed due to Hardy-Weinberg exact test.
### 949206 variants and 474 people pass filters and QC.

#5. LD prune (rm 1 SNP if r2>0.3 in 50 SNP window) for relationship check and heterozygosity calculation
plink --bfile F13_Recluster_TOP_20161226_FinalReport.geno0.01.hwe1e-06 --indep-pairwise 50 5 0.3 --out F13_Recluster_TOP_20161226_FinalReport.geno0.01.hwe1e-06 --maf 0.01
### 268268 variants removed due to minor allele threshold(s)
### Pruning complete.  494303 of 680720 variants removed.
### Marker lists written to
### F13_Recluster_TOP_20161226_FinalReport.geno0.01.hwe1e-06.prune.in and
### F13_Recluster_TOP_20161226_FinalReport.geno0.01.hwe1e-06.prune.out .

#6. Relationship check
plink --bfile F13_Recluster_TOP_20161226_FinalReport.geno0.01.hwe1e-06 --extract F13_Recluster_TOP_20161226_FinalReport.geno0.01.hwe1e-06.prune.in --genome --min 0.05 --out F13_Recluster_TOP_20161226_FinalReport.geno0.01.hwe1e-06.LD0.3

### checked for expected duplicates and known hapmap relationships and found 5 pairs of unexpected duplicates in 40_GWAS_riken2_QC_plots.Rmd, made list of one of a known duplicate pair (17), hapmap samples (6), all unexpected duplicates (10), and one of pair with pi-hat > 0.125 (0). Rerun relationship check with these excluded. 
plink --bfile F13_Recluster_TOP_20161226_FinalReport.geno0.01.hwe1e-06 --extract F13_Recluster_TOP_20161226_FinalReport.geno0.01.hwe1e-06.prune.in --remove RIKEN2_hapmapDuplicateRelated0.125_List.txt --genome --out F13_Recluster_TOP_20161226_FinalReport.geno0.01.hwe1e-06.LD0.3.rmDupsHapmap.rmPIHAT0.125
### 949206 variants loaded from .bim file.
### 474 people (471 males, 3 females) loaded from .fam.
### --extract: 186417 variants remaining.
### --remove: 441 people remaining.

#7. Check heterozygosity (across all autosomal SNPs) -- look at that distribution across individuals to check for and rm outliers (F: mean +/-6 sd), see 40_GWAS_riken2_QC_plots.Rmd.
plink --bfile F13_Recluster_TOP_20161226_FinalReport.geno0.01.hwe1e-06 --extract F13_Recluster_TOP_20161226_FinalReport.geno0.01.hwe1e-06.prune.in --het --out F13_Recluster_TOP_20161226_FinalReport.geno0.01.hwe1e-06.LD0.3.rmDupsHapmap.rmPIHAT0.125 --remove RIKEN2_hapmapDuplicateRelated0.125_List.txt
### --extract: 186417 variants remaining.
### --remove: 441 people remaining.
### Total genotyping rate in remaining samples is 0.999718.
### 186417 variants and 441 people pass filters and QC.
### --het: 181942 variants scanned, report written to
### F13_Recluster_TOP_20161226_FinalReport.geno0.01.hwe1e-06.LD0.3.rmDupsHapmap.rmPIHAT0.125.het

### add 3 to remove list (see 40_GWAS_riken2_QC_plots.Rmd): RIKEN2_hapmapDuplicate_pihat0.125_hetOutlier_exclusion_list.txt

#8. Prepare genotype files for PCA
#make list of chr0 SNPs to remove
plink --bfile F13_Recluster_TOP_20161226_FinalReport.geno0.01.hwe1e-06 --chr 0 --make-bed --out F13_chr0
cut -f 2 F13_chr0.bim > chr0.SNPlist

plink --bfile F13_Recluster_TOP_20161226_FinalReport.geno0.01.hwe1e-06 --remove RIKEN2_hapmapDuplicate_pihat0.125_hetOutlier_exclusion_list.txt --exclude chr0.SNPlist --make-bed --out F13_Recluster_TOP_20161226_FinalReport.forPCA
###Total genotyping rate in remaining samples is 0.999722.
###948628 variants and 438 people pass filters and QC.

#merge with hg19 HAPMAP3
plink --bfile F13_Recluster_TOP_20161226_FinalReport.forPCA --bmerge /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/HAPMAP3/HM3_ASN_CEU_YRI_Unrelated_hg19_noAmbig.bed /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/HAPMAP3/HM3_ASN_CEU_YRI_Unrelated_hg19_noAmbig.bim /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/HAPMAP3/HM3_ASN_CEU_YRI_Unrelated_hg19_noAmbig.fam --make-bed --out F13_HM3
#add warning SNPs to *missnp list and remove from HM3
grep -o 'rs[0-9][0-9]*' F13_HM3.log >> F13_HM3-merge.missnp 
plink --bfile /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/HAPMAP3/HM3_ASN_CEU_YRI_Unrelated_hg19_noAmbig --exclude F13_HM3-merge.missnp --make-bed --out /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/HAPMAP3/HM3_hg19_for_F13_PCA

##try merge again
plink --bfile F13_Recluster_TOP_20161226_FinalReport.forPCA --bmerge /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/HAPMAP3/HM3_hg19_for_F13_PCA.bed /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/HAPMAP3/HM3_hg19_for_F13_PCA.bim /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/HAPMAP3/HM3_hg19_for_F13_PCA.fam --out F13_HM3_try2
#####see F13_HM3_try2.log for lots of variants in F13 with the same position, should get filtered out in LD-prune prior to PCA
# filter merged file to SNPs with >99% genotypes and MAF > 0.05 and LD prune
plink --bfile F13_HM3_try2 --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.3 --out F13_HM3_LDpruned
# make ped and map file for smartpca 
plink --bfile F13_HM3_try2 --geno 0.01 --maf 0.05 --extract F13_HM3_LDpruned.prune.in --recode --out F13_HM3_LDpruned
# make fam file (no -9s) for smartpca  
awk '{print $1,$2,$3,$4,$5,1}' F13_HM3_try2.fam > F13_HM3_LDpruned.fam
# make parfile for smartpca 
perl /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/PCA/make_par_file.pl /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/F13_HM3_LDpruned 0 > /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/PCA/F13_HM3_LDpruned.par

#9. Run PCA
module load eigensoft/5.0.1
cd /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/PCA/
qsub run_smartpca.sh

#plot PCs and keep samples w/in 12 sd of CEU cluster, see 40_GWAS_QC_plots.Rmd, less stringent prior to combining with RIKEN1
#rerun PCA with Euro samples only
cd /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/
plink --bfile F13_Recluster_TOP_20161226_FinalReport.forPCA --keep F13.euro.GWAS.PCAs --make-bed --out F13_Recluster_TOP_20161226_FinalReport.postPCA.euro
plink --file F13_HM3_LDpruned --keep F13.euro.GWAS.PCAs --recode --out F13_euro_LDpruned
awk '{print $1,$2,$3,$4,$5,1}' F13_Recluster_TOP_20161226_FinalReport.postPCA.euro.fam > F13_euro_LDpruned.fam
perl /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/PCA/make_par_file.pl /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/F13_euro_LDpruned 0 > /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/PCA/F13_euro_LDpruned.par
cd /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/PCA/
qsub run_smartpca.sh
#plot PCs, PtStudy and Norway are interspersed, see 03_GWAS_QC_plots.Rmd

#10. Plate effects analysis
##make pseudo-phenotype file for plink
./41_make_plate_phenofile_riken2.py
##run GWAS's with each plate as "cases", rm maf<0.05 SNPs to improve the stability of test statistics
plink --bfile F13_Recluster_TOP_20161226_FinalReport.postPCA.euro --maf 0.05 --assoc --pheno plate_check_phenotypes_riken2 --all-pheno --adjust --out F13_plate_effects_maf0.05
##view QQ plots, make list of SNPs with p<1e-4 in plate GWAS's in 42_plate_effects_riken2.Rmd
##QQ plots look good, only 374 SNPs with P<1e-4, remove from dataset

#11. Make post-qc bed/bim/fam 
#rm plate effect SNPs and the 6504 SNPs in exm-N88_HM3.list to simplify merging with N88
cat F13_plate_effects_excludeSNPs_p1e-4.txt exm-N88_HM3.list > F13_plate_exm_exclude.SNPlist
plink --bfile F13_Recluster_TOP_20161226_FinalReport.postPCA.euro --exclude F13_plate_exm_exclude.SNPlist --make-bed --out F13_Recluster_TOP_20161226_FinalReport.postQC

#12. prep for merge with N88 (riken1)
##rm ambiguous-strand SNPs (A/T and C/G)
perl /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/08_pull_unamb_SNPs.pl F13_Recluster_TOP_20161226_FinalReport.postQC.bim > F13.unamb.snplist
plink --bfile F13_Recluster_TOP_20161226_FinalReport.postQC --extract F13.unamb.snplist --maf 0.01 --make-bed --out F13.forImputation