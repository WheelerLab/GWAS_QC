#############################################################
# NEURP IMPUTATION WORKFLOW
# Heather Wheeler 3/25/2014
# 	  updated: 10/14/2014
# GWAS QC - PHASING with SHAPEIT - IMPUTATION with IMPUTE2
#
##############################################################

#1. Check call rates and rm poorly called SNPs.	See plink.sh in 1_preIMPUTE2_QC/
	plink --noweb --bfile ${DATA_DIR}$bfile --missing --out ${DATA_DIR}$bfile
        plink --noweb --bfile ${DATA_DIR}$bfile --geno 0.1	--set-hh-missing --make-bed --out ${DATA_DIR}$bfile.geno
        plink --noweb --bfile ${DATA_DIR}$bfile.geno --mind 0.1 --make-bed	--out ${DATA_DIR}$bfile.geno.mind

#2. Check sex and rm discrepencies. See plink.sh in 1_preIMPUTE2_QC/
	plink --noweb --bfile ${DATA_DIR}$bfile.geno.mind --check-sex --out ${DATA_DIR}$bfile.geno.mind
        perl parse_sexcheck.pl ${DATA_DIR}$bfile.geno.mind.sexcheck > ${DATA_DIR}$bfile.wrongsex.samples
        plink --noweb --bfile ${DATA_DIR}$bfile.geno.mind --remove ${DATA_DIR}$bfile.wrongsex.samples --make-bed --out ${DATA_DIR}$bfile.geno.mind.sex

#3. Calculate HWE statistics and rm SNPs with P < e-8. See plink.sh in 1_preIMPUTE2_QC/
	plink --noweb --bfile ${DATA_DIR}$bfile.geno.mind.sex --hardy --out ${DATA_DIR}$bfile.geno.mind.sex
	plink --noweb --bfile ${DATA_DIR}$bfile.geno.mind.sex --hwe 0.00000001 --make-bed --out ${DATA_DIR}$bfile.geno.mind.sex.hwe

#4. Check heterozygosity (across all autosomal SNPs) -- look at that distribution across individuals to check for and rm outliers (F > 0.1 or F < -0.05). See plink.sh in 1_preIMPUTE2_QC/
	perl pull_autosome_snplist.pl ${DATA_DIR}$bfile.geno.mind.sex.hwe.bim > ${DATA_DIR}$bfile.autosome.snplist
	plink --noweb --bfile ${DATA_DIR}$bfile.geno.mind.sex.hwe --extract ${DATA_DIR}$bfile.autosome.snplist --make-bed --out ${DATA_DIR}$bfile.geno.mind.sex.hwe.auto
	plink --noweb --bfile ${DATA_DIR}$bfile.geno.mind.sex.hwe.auto --het --out ${DATA_DIR}$bfile.geno.mind.sex.hwe.auto

	R --vanilla < pull_het_outliers.r --args ${DATA_DIR}$bfile.geno.mind.sex.hwe.auto.het
	#keep X and Y chrs
	perl ${DATA_DIR}pull_chr1-24_snplist.pl ${DATA_DIR}$bfile.geno.mind.sex.hwe.bim > ${DATA_DIR}$bfile.chr1-24.snplist                                                                                                                           
        plink --noweb --bfile ${DATA_DIR}$bfile.geno.mind.sex.hwe --extract ${DATA_DIR}$bfile.chr1-24.snplist --remove ${DATA_DIR}$bfile.geno.mind.sex.hwe.auto.het.outlier.list --make-bed --out ${DATA_DIR}$bfile.geno.mind.sex.hwe.chr1-24.het

#5. Relationship check -- rm one of pair with pi_hat>0.05, rm the one with no phenotype if possible. See plink.sh in 1_preIMPUTE2_QC/
	plink --noweb --bfile ${DATA_DIR}$bfile.geno.mind.sex.hwe.chr1-24.het --genome --min 0.05 --out ${DATA_DIR}$bfile.geno.mind.sex.hwe.chr1-24.het
	#make ${DATA_DIR}$bfile.related.to.remove by hand
	plink --noweb --bfile ${DATA_DIR}$bfile.geno.mind.sex.hwe.chr1-24.het --remove ${DATA_DIR}$bfile.related.to.remove --make-bed --out ${DATA_DIR}$bfile.geno.mind.sex.hwe.chr1-24.het

#5b. Remove A/T and C/G SNPs prior to merging with hapmap and for imputation
	perl ${DATA_DIR}pull_unamb_SNPs.pl ${DATA_DIR}$bfile.geno.mind.sex.hwe.chr1-24.het.bim > ${DATA_DIR}$bfile.unamb.SNPlist
	plink --noweb --bfile ${DATA_DIR}$bfile.geno.mind.sex.hwe.chr1-24.het --extract ${DATA_DIR}$bfile.unamb.SNPlist --make-bed --out ${DATA_DIR}$bfile.unamb

#6. Compute PCs and rm outliers With HapMap: See plink.sh in 1_preIMPUTE2_QC/

	#merge dataset with HapMap3 to run smartpca
	plink --noweb --bfile ${DATA_DIR}$bfile.unamb --bmerge /nas40t0/HAPMAP3/HM3_ASN_CEU_YRI_Unrelated.bed /nas40t0/HAPMAP3/HM3_ASN_CEU_YRI_Unrelated.bim /nas40t0/HAPMAP3/HM3_ASN_CEU_YRI_Unrelated.fam --make-bed --out ${DATA_DIR}$bfile.HapMap3
	plink --noweb --bfile ${DATA_DIR}$bfile.unamb --exclude ${DATA_DIR}$bfile.HapMap3.missnp --make-bed --out ${DATA_DIR}$bfile.rm.missnp
	plink --noweb --bfile /nas40t0/HAPMAP3/HM3_ASN_CEU_YRI_Unrelated --exclude ${DATA_DIR}$bfile.HapMap3.missnp --make-bed --out HapMap3
	plink --noweb --bfile ${DATA_DIR}$bfile.rm.missnp --indep-pairwise 50 5 0.2 --out ${DATA_DIR}$bfile.rm.missnp 
	plink --noweb --bfile ${DATA_DIR}$bfile.rm.missnp --extract ${DATA_DIR}$bfile.rm.missnp.prune.in --make-bed --out ${DATA_DIR}$bfile.ldpruned
	awk '{print $2}' ${DATA_DIR}$bfile.ldpruned.bim > ${DATA_DIR}$bfile.merge.snplist
	plink --noweb --bfile HapMap3 --extract ${DATA_DIR}$bfile.merge.snplist --make-bed --out HapMap3
	plink --noweb --bfile ${DATA_DIR}$bfile.ldpruned --bmerge HapMap3.bed HapMap3.bim HapMap3.fam --recode --out ${DATA_DIR}$bfile.HapMap3.forPCA
	cut -d" " -f 1-5 ${DATA_DIR}$bfile.HapMap3.forPCA.ped > fam
	awk '{print $1,$2,$3,$4,$5,1}' fam > ${DATA_DIR}$bfile.HapMap3.forPCA.fam

	perl make_par_file.pl ${DATA_DIR}$bfile.HapMap3.forPCA 0 > ${DATA_DIR}$bfile.HapMap3.forPCA.par
	smartpca -p ${DATA_DIR}$bfile.HapMap3.forPCA.par

	#plot in R and choose genetic Europeans (PC1 & PC2 w/in 10 SDs of CEU means) for followup
	#make pop vector for legend, color
	tail -n +2 ${DATA_DIR}$bfile.HapMap3.forPCA.evec > ${DATA_DIR}$bfile.HapMap3.evec
	awk '{print "GWAS"}' ${DATA_DIR}$bfile.ldpruned.fam > ${DATA_DIR}$bfile.pop
	cat ${DATA_DIR}$bfile.pop HapMap3.pop > ${DATA_DIR}$bfile.HapMap3.pop
	R --vanilla < plot.pca.choose.euros.r --args ${DATA_DIR}$bfile.HapMap3
	perl split.colon.pl ${DATA_DIR}$bfile.HapMap3.euro.GWAS.PCAs > ${DATA_DIR}$bfile.HapMap3.euro.keep.list

	###upon examing PCA plots, I should keep all E5103 and GoKinD samples, manually update *keep.list files
	cp E5103_p1.ldpruned.fam E5103_p1.HapMap3.euro.keep.list
	cp E5103_p2.ldpruned.fam E5103_p2.HapMap3.euro.keep.list
	cp GoKinD.ldpruned.fam GoKinD.HapMap3.euro.keep.list

#7. Compute PCs and rm ouliers With Genetic Europeans: See plink.sh in 1_preIMPUTE2_QC/

	plink --noweb --bfile ${DATA_DIR}$bfile.ldpruned --keep ${DATA_DIR}$bfile.HapMap3.euro.keep.list --recode --out ${DATA_DIR}$bfile.ldpruned.euro.forPCA
	cut -d" " -f 1-5 ${DATA_DIR}$bfile.ldpruned.euro.forPCA.ped > fam
	awk '{print $1,$2,$3,$4,$5,1}' fam > ${DATA_DIR}$bfile.ldpruned.euro.forPCA.fam
	perl make_par_file.pl ${DATA_DIR}$bfile.ldpruned.euro.forPCA 0 > ${DATA_DIR}$bfile.ldpruned.euro.forPCA.par
	smartpca -p ${DATA_DIR}$bfile.ldpruned.euro.forPCA.par
	perl make_par_file_rmoutlier.pl ${DATA_DIR}$bfile.ldpruned.euro.forPCA 5 > ${DATA_DIR}$bfile.ldpruned.euro.forPCA.rmout.par
	smartpca -p ${DATA_DIR}$bfile.ldpruned.euro.forPCA.rmout.par 	

	#plot homogeneous PCA results in R
	tail -n +2 ${DATA_DIR}$bfile.ldpruned.euro.forPCA.evec > ${DATA_DIR}$bfile.ldpruned.euro.evec
	R --vanilla < plot.pca.r --args ${DATA_DIR}$bfile.ldpruned.euro
	tail -n	+2 ${DATA_DIR}$bfile.ldpruned.euro.forPCA.rmout.evec > ${DATA_DIR}$bfile.ldpruned.euro.rmout.evec
        R --vanilla < plot.pca.r --args ${DATA_DIR}$bfile.ldpruned.euro.rmout 

	###for consistency with Baldwin et al, don't remove flagged outliers,  make final QC'd bed/bim/bam for SHAPEIT
	plink --noweb --bfile ${DATA_DIR}$bfile.unamb --keep ${DATA_DIR}$bfile.HapMap3.euro.keep.list --make-bed --out /nas40t0/hwheeler/NEURP/IMPUTE2_NEURP/2_prePHASE_shapeit/${DATA_DIR}$bfile.QC

#8. liftOver SNP coordinates hg18/B36 to hg19/B37 See shapeit.sh in 2_prePHASE_shapeit/
	#only necessary for GoKinD, the rest are B37
	awk '{print "chr"$1,$4,$4+1}' GoKinD.QC.bim > GoKinD.QC.B36.coords 
	~/bin/liftOver GoKinD.QC.B36.coords ~/bin/hg18ToHg19.over.chain.gz GoKinD.B36toB37.successes GoKinD.B36toB37.failures
	paste GoKinD.QC.B36.coords GoKinD.QC.bim > GoKinD.QC.coords.bim.merged
	perl ~/bin/find_failed_snps.pl GoKinD.QC.coords.bim.merged GoKinD.B36toB37.failures > GoKinD.failures
	plink --noweb --bfile GoKinD.QC --exclude GoKinD.failures --make-bed --out GoKinD.QC
	paste GoKinD.B36toB37.successes GoKinD.QC.bim > prebim
	perl ~/bin/update_bim.pl prebim > GoKinD.QC.bim

#9. merge E5103_p1 and E5103_p2 See shapeit.sh in 2_prePHASE_shapeit/
	plink --noweb --bfile E5103_p1.QC --bmerge E5103_p2.QC.bed E5103_p2.QC.bim E5103_p2.QC.fam --make-bed --out E5103.QC

#10. split plink bfiles by chr
	#add to 2_prePHASE_shapeit/print22.pl:
	print OUT "plink --bfile ${DATA_DIR}$bfile.QC --chr $i --make-bed --out ${DATA_DIR}$bfile.chr$i\n";

	perl print22.pl
	source qsub.txt

#11. strand alignment with TGP
	#add to 2_prePHASE_shapeit/print22.pl, comment out previous: ###next time update to ref1k.201312.plink ##alignment should be same b/c both B37
	print OUT "plink --noweb --bfile /nas40t2/anuarStuff2/ref1kg/ref1k.201312.plink/ref1k.201312.$i --bmerge ${DATA_DIR}$bfile.chr$i.bed ${DATA_DIR}$bfile.chr$i.bim ${DATA_DIR}$bfile.chr$i.fam --make-bed --out TGP.${DATA_DIR}$bfile.merged.chr$i\n";
	print OUT "plink --noweb --bfile ${DATA_DIR}$bfile.chr$i --flip TGP.${DATA_DIR}$bfile.merged.chr$i.missnp --make-bed --out ${DATA_DIR}$bfile.B37.chr$i.unphased\n";
	perl print22.pl
        source qsub.txt

#12. prephase with shapeit
	#add to 2_prePHASE_shapeit/print22.pl, comment out previous:
	print OUT "~/bin/shapeit.v2.r644.linux.x86_64 -B ${DATA_DIR}$bfile.B37.chr$i.unphased -M genetic_map_chr$i\_combined_b37.txt -O ${DATA_DIR}$bfile.B37.chr$i.phased\n";
	perl print22.pl
        source qsub.txt

#13. run IMPUTE2
	#from cwd, replace called script (impute2.sh) with appropriate cohort name (${DATA_DIR}$bfile)
	###old way### source qsub.impute2.set1.txt #run sets 1-6, each contains 100 jobs (max # of jobs allowed to run per user on cluster)
	#updated to run 10 chr regions per run (reduces total #jobs to 58, so don't need to keep checking cluster)
	source qsub.all.txt

#14. concatenate IMPUTE2 results 
	#make bfile.list.subset if not everything is done running. From 3_IMPUTE2_results/
	qsub concat_impute2_results.sh

#15. pull out SNPs with info score > 0.3
	#within 3_IMPUTE2_results/cat_results/
	qsub pull_snps_info.sh		

#16. pull out SNPs from HapMap2 for PrediXcan (predictors currently only have HapMap2 SNPs, allows Haky's imputeProb2Dosage.py to work), see /nas40t2/hwheeler/predict_EXP
	#hapmap.snplist is from /userhome/hwheeler/Paclitaxel_GWAS/CEU_and_YRI1_merge.bim
	#within 3_IMPUTE2_results/cat_results/
	qsub pull_hapmap_snps.sh

#17. run SNPTEST on allele dosages to get GWAS results within 4_SNPTEST_GWAS
	#for resid phenotypes add to print22.pl:
	print OUT "~/bin/snptest_v2.4.1 -data /nas40t0/hwheeler/NEURP/IMPUTE2_NEURP/3_IMPUTE2_results/cat_results/${DATA_DIR}$bfile";
	print OUT "_chr$i.info_gt0.3_impute2.gz ${DATA_DIR}$bfile.sample -o output/${DATA_DIR}$bfile.chr$i.SNPTEST.resid.nocovar.out -pheno resid -frequentist 1 -method expected\n";
	
	source qsub.txt

	#for GoKinD add to print22.pl:
	snptest -data /nas40t0/hwheeler/NEURP/IMPUTE2_NEURP/3_IMPUTE2_results/GoKinD_chr22.pos15000001-20000000.impute2 GoKinD.noPCA.sample -o GK.covar.noPCA -pheno NEURP -frequentist 1 -method expected -cov_all
