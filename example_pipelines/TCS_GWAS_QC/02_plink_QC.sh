#!/bin/bash

###################################################################
#    THE PLATINUM STUDY RIKEN GWAS QC WORKFLOW                    #
#    Heather E. Wheeler 2015-10-09                                #
###################################################################

##load software on tarbell
module load R
module load plink/1.90 ##plink2
module load vcftools
module load tabix/0.2.6

$DIR = /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/
cd $DIR

#1. Convert tped/tfam made in 01_make_tped.pl to PLINK binary, set any heterozygous haploid genotypes (male X,Y,MT) to missing
plink --tfile N88_Recluster_TOP_20150911_FinalReport --set-hh-missing --make-bed --out N88_Recluster_TOP_20150911_FinalReport
### 1145 people (1142 males, 3 females) loaded from .fam.
### Warning: 323965 het. haploid genotypes present (see
### N88_Recluster_TOP_20150911_FinalReport.hh ).
### Warning: Nonmissing nonmale Y chromosome genotype(s) present.
### Total genotyping rate is 0.996293.

#2. Check sex and calculate call rates for flagging poorly called SNPs and individuals
plink --bfile N88_Recluster_TOP_20150911_FinalReport --check-sex --missing --out N88_Recluster_TOP_20150911_FinalReport
### --missing: Sample missing data report written to
### N88_Recluster_TOP_20150911_FinalReport.imiss, and variant-based missing data
### report written to N88_Recluster_TOP_20150911_FinalReport.lmiss.
### 964193 variants and 1145 people pass filters and QC.
### --check-sex: 18460 Xchr and 0 Ychr variant(s) scanned, no problems detected.
### Report written to N88_Recluster_TOP_20150911_FinalReport.sexcheck .

#3. Recalculate individual call rates after removing SNPs with call rates <99%
plink --bfile N88_Recluster_TOP_20150911_FinalReport --geno 0.01 --make-bed --out N88_Recluster_TOP_20150911_FinalReport.geno0.01
### 30486 variants removed due to missing genotype data (--geno).
### 933707 variants and 1145 people pass filters and QC.
plink --bfile N88_Recluster_TOP_20150911_FinalReport.geno0.01 --missing --out N88_Recluster_TOP_20150911_FinalReport.geno0.01
### Total genotyping rate is 0.999497.
### --missing: Sample missing data report written to
### N88_Recluster_TOP_20150911_FinalReport.geno0.01.imiss, and variant-based
### missing data report written to
### N88_Recluster_TOP_20150911_FinalReport.geno0.01.lmiss.

### looks great, all individuals now have >99.2% call rates (see 03_GWAS_QC_plots.Rmd output)

#4. Calculate HWE statistics to flag SNPs later
plink --bfile N88_Recluster_TOP_20150911_FinalReport.geno0.01 --hardy --out N88_Recluster_TOP_20150911_FinalReport.geno0.01
### 1842/932292 SNPs (0.2%) have P<1e-06, remove next

#5. LD prune (rm 1 SNP if r2>0.3 in 50 SNP window) for relationship check and heterozygosity calculation
plink --bfile N88_Recluster_TOP_20150911_FinalReport.geno0.01 --indep-pairwise 50 5 0.3 --out N88_Recluster_TOP_20150911_FinalReport.geno0.01

#5. Relationship check
plink --bfile N88_Recluster_TOP_20150911_FinalReport.geno0.01 --extract N88_Recluster_TOP_20150911_FinalReport.geno0.01.prune.in --genome --min 0.05 --out N88_Recluster_TOP_20150911_FinalReport.geno0.01.LD0.3
### checked for expected duplicates and known hapmap relationships and found 8 pairs of unexpected duplicates in 03_GWAS_QC_plots.Rmd, made list of one of a known duplicate pair, hapmap samples, and all 16 unexpected duplicates. Rerun relationship check with these excluded.
plink --bfile N88_Recluster_TOP_20150911_FinalReport.geno0.01 --extract N88_Recluster_TOP_20150911_FinalReport.geno0.01.prune.in --remove hapmapDuplicateList.txt --genome --min 0.05 --out N88_Recluster_TOP_20150911_FinalReport.geno0.01.LD0.3.rmKnownDupsHapmap 
### remove one of pair w/pi-hat > 0.05 (keep sample with audiometry if possible), run full IBD analysis:
plink --bfile N88_Recluster_TOP_20150911_FinalReport.geno0.01 --extract N88_Recluster_TOP_20150911_FinalReport.geno0.01.prune.in --remove hapmapDuplicate_pihat0.05_exclusion_list.txt --genome --out N88_Recluster_TOP_20150911_FinalReport.geno0.01.LD0.3.rmKnownDupsHapmap.rmPIHAT0.05

#5. Check heterozygosity (across all autosomal SNPs) -- look at that distribution across individuals to check for and rm outliers (F: mean +/-3 sd), see 03_GWAS_QC_plots.Rmd.
plink --bfile N88_Recluster_TOP_20150911_FinalReport.geno0.01 --het --extract N88_Recluster_TOP_20150911_FinalReport.geno0.01.prune.in --remove hapmapDuplicate_pihat0.05_exclusion_list.txt --out N88_Recluster_TOP_20150911_FinalReport.geno0.01.LD0.3.rmKnownDupsHapmap.rmPIHAT0.05
### --extract: 294247 variants remaining.
### --remove: 1036 people remaining.
### Total genotyping rate in remaining samples is 0.999575.
### 294247 variants and 1036 people pass filters and QC.
### --het: 282186 variants scanned, report written to
### N88_Recluster_TOP_20150911_FinalReport.geno0.01.LD0.3.rmKnownDupsHapmap.rmPIHAT0.05.het

#6. Prepare genotype files for PCA
#make list of chr0 SNPs to remove
plink --bfile N88_Recluster_TOP_20150911_FinalReport.forPCA --chr 0 --make-bed --out N88_chr0
cut -f 2 N88_chr0.bim > chr0.SNPlist

plink --bfile N88_Recluster_TOP_20150911_FinalReport.geno0.01 --remove hapmapDuplicate_pihat0.05_hetOutlier_exclusion_list.txt --exclude chr0.SNPlist --make-bed --out N88_Recluster_TOP_20150911_FinalReport.forPCA
### Total genotyping rate in remaining samples is 0.99952.
### 933707 variants and 1014 people pass filters and QC.

#merge dataset with HapMap3 to run smartpca
module load eigensoft/5.0.1 

###liftOver hg18 to hg19 coordinated in HapMap3 data
cd /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/HAPMAP3
module load liftOver
awk '{print "chr"$1,$4,$4+1}' /group/dolan-lab/hwheeler/nas40t0/HAPMAP3/HM3_ASN_CEU_YRI_Unrelated.bim > HM3_ASN_CEU_YRI_Unrelated_hg18.coords
#change chr23 -> chrX, chr24 -> chrY, chr26 -> chrM.
perl replaceXYM.pl HM3_ASN_CEU_YRI_Unrelated_hg18.coords
liftOver HM3_ASN_CEU_YRI_Unrelated_hg18.coords.XYM /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/hg18ToHg19.over.chain.gz HM3.hg18-to-hg19.successes HM3.hg18-to-hg19.failures
paste HM3_ASN_CEU_YRI_Unrelated_hg18.coords.XYM /group/dolan-lab/hwheeler/nas40t0/HAPMAP3/HM3_ASN_CEU_YRI_Unrelated.bim > HM3_ASN_CEU_YRI_Unrelated.coords.bim
perl find_failed_snps.pl HM3.hg18-to-hg19.failures HM3_ASN_CEU_YRI_Unrelated.coords.bim HM3.failures
plink --bfile /group/dolan-lab/hwheeler/nas40t0/HAPMAP3/HM3_ASN_CEU_YRI_Unrelated --exclude HM3.failures --make-bed --out HM3_ASN_CEU_YRI_Unrelated_hg19
paste HM3.hg18-to-hg19.successes HM3_ASN_CEU_YRI_Unrelated_hg19.bim > prebim
perl update_bim.pl prebim > HM3_ASN_CEU_YRI_Unrelated_hg19.bim
#remove A/T and C/G SNPs (ambiguous strand) from file
perl pull_unamb_SNPs.pl HM3_ASN_CEU_YRI_Unrelated_hg19.bim > HM3.unamb.SNPlist
plink --bfile HM3_ASN_CEU_YRI_Unrelated_hg19 --extract HM3.unamb.SNPlist --make-bed --out HM3_ASN_CEU_YRI_Unrelated_hg19_noAmbig

#merge with hg19 HAPMAP3
cd /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/
plink --bfile N88_Recluster_TOP_20150911_FinalReport.forPCA --bmerge ../HAPMAP3/HM3_ASN_CEU_YRI_Unrelated_hg19_noAmbig.bed ../HAPMAP3/HM3_ASN_CEU_YRI_Unrelated_hg19_noAmbig.bim ../HAPMAP3/HM3_ASN_CEU_YRI_Unrelated_hg19_noAmbig.fam --make-bed --out N88_HM3
#add warning SNPs to *missnp list and remove from HM3
#Warning: Multiple positions seen for variant 'rs1886730'.
#Warning: Multiple positions seen for variant 'rs11573991'.
#Warning: Multiple positions seen for variant 'rs2234167'.
#Warning: Multiple chromosomes seen for variant 'rs12043679'.
#Warning: Multiple chromosomes seen for variant 'rs10182914'.
#Warning: Multiple chromosomes seen for variant 'rs12496398'.
#Warning: Multiple chromosomes seen for variant 'rs11942835'.
#Warning: Multiple chromosomes seen for variant 'rs4291004'.
#Warning: Multiple chromosomes seen for variant 'rs2569201'.
#Warning: Multiple chromosomes seen for variant 'rs4714901'.
#Warning: Multiple chromosomes seen for variant 'rs16942804'.
#Warning: Multiple chromosomes seen for variant 'rs10106770'.
#Warning: Multiple chromosomes seen for variant 'rs10957824'.
#Warning: Multiple chromosomes seen for variant 'rs1578263'.
#Warning: Multiple chromosomes seen for variant 'rs12804886'.
#Warning: Multiple chromosomes seen for variant 'rs12900938'.
#Warning: Multiple chromosomes seen for variant 'rs11857958'.
#Warning: Multiple chromosomes seen for variant 'rs17728665'.
cat warninglist N88_HM3-merge.missnp > warning_N88_HM3-merge.missnp
cd /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/HAPMAP3
plink --bfile HM3_ASN_CEU_YRI_Unrelated_hg19_noAmbig --extract ../genotypes/warning_N88_HM3-merge.missnp --make-bed --out HM3_hg19_forPCA

##try merge again
cd /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/
plink --bfile N88_Recluster_TOP_20150911_FinalReport.forPCA --bmerge ../HAPMAP3/HM3_hg19_forPCA.bed ../HAPMAP3/HM3_hg19_forPCA.bim ../HAPMAP3/HM3_hg19_forPCA.fam --out N88_HM3
##see N88_HM3.log for lots of variants in N88 with the same position, exclude one of each pair (exm-*) before GWAS
# filter merged file to SNPs with >90% genotypes
plink --bfile N88_HM3 --geno 0.2 --maf 0.05 --make-bed --out N88_HM3_geno0.2_maf0.05
# (no exm variants remain b/c none in hapmap)
# make ped and map file and fam file (no -9s) for smartpca
plink --bfile N88_HM3_geno0.2_maf0.05 --indep-pairwise 50 5 0.2 --recode --out N88_HM3_LDpruned
awk '{print $1,$2,$3,$4,$5,1}' N88_HM3_geno0.2_maf0.05.fam > N88_HM3_LDpruned.fam
# make parfile for smartpca
cd /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/PCA
perl make_par_file.pl /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/N88_HM3_LDpruned 0 > N88_HM3_LDpruned.par

#6. Run PCA
module load eigensoft/5.0.1
qsub run_smartpca.sh

#plot PCs and keep samples w/in 10 sd of CEU cluster, see 03_GWAS_QC_plots.Rmd
#rerun PCA with Euro samples only
cd /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/
plink --bfile N88_Recluster_TOP_20150911_FinalReport.forPCA --keep N88.euro.GWAS.PCAs --make-bed --out N88_Recluster_TOP_20150911_FinalReport.postPCA.euro
plink --file N88_HM3_LDpruned --keep N88.euro.GWAS.PCAs --recode --out N88_euro_LDpruned
awk '{print $1,$2,$3,$4,$5,1}' N88_Recluster_TOP_20150911_FinalReport.postPCA.euro.fam > N88_euro_LDpruned.fam
cd /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/PCA
perl make_par_file.pl /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/N88_euro_LDpruned 0 > N88_euro_LDpruned.par
qsub run_smartpca.sh
#plot PCs, PtStudy and Norway are interspersed, see 03_GWAS_QC_plots.Rmd

#7. Plate effects analysis
##make pseudo-phenotype file for plink
./04_make_plate_phenofile.py
##run GWAS's with each plate as "cases", rm maf<0.05 SNPs to improve the stability of test statistics
plink --bfile N88_Recluster_TOP_20150911_FinalReport.postPCA.euro --maf 0.05 --assoc --pheno plate_check_phenotypes --all-pheno --adjust --out N88_plate_effects_maf0.05
##view QQ plots, make list of SNPs with p<1e-4 in plate GWAS's in 05_plate_effects.Rmd
##QQ plots look good, only 950 SNPs with P<1e-4, remove from dataset

#8. Make post-qc bed/bim/fam 
plink --bfile N88_Recluster_TOP_20150911_FinalReport.postPCA.euro --exclude N88_plate_effects_excludeSNPs_p1e-4.txt --make-bed --out N88_Recluster_TOP_20150911_FinalReport.postQC

