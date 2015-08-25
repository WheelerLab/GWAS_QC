#!/bin/bash
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N plink
#$ -j y


#0 plink --file GoKinDAndEdic --make-bed --out GoKinD

for bfile in `cat bfile.list.subset`;

do

#1
	plink --noweb --bfile $bfile --missing --out $bfile
#1
	plink --noweb --bfile $bfile --geno 0.1 --set-hh-missing --make-bed --out $bfile.geno
#1
	plink --noweb --bfile $bfile.geno --mind 0.1 --make-bed --out $bfile.geno.mind

#2
	plink --noweb --bfile $bfile.geno.mind --check-sex --out $bfile.geno.mind
#2
	perl parse_sexcheck.pl $bfile.geno.mind.sexcheck > $bfile.wrongsex.samples
#2
	plink --noweb --bfile $bfile.geno.mind --remove $bfile.wrongsex.samples --make-bed --out $bfile.geno.mind.sex

#3
	plink --noweb --bfile $bfile.geno.mind.sex --hardy --out $bfile.geno.mind.sex
#3
      plink --noweb --bfile $bfile.geno.mind.sex --hwe 0.00000001 --make-bed --out $bfile.geno.mind.sex.hwe

#4
      perl pull_autosome_snplist.pl $bfile.geno.mind.sex.hwe.bim > $bfile.autosome.snplist
#4
      plink --noweb --bfile $bfile.geno.mind.sex.hwe --extract $bfile.autosome.snplist --make-bed --out $bfile.geno.mind.sex.hwe.auto
#4
      plink --noweb --bfile $bfile.geno.mind.sex.hwe.auto --het --out $bfile.geno.mind.sex.hwe.auto	
#4
	R --vanilla < pull_het_outliers.r --args $bfile.geno.mind.sex.hwe.auto.het
#4
      plink --noweb --bfile $bfile.geno.mind.sex.hwe.auto --remove $bfile.geno.mind.sex.hwe.auto.het.outlier.list --make-bed --out $bfile.geno.mind.sex.hwe.auto.het

#5	plink --noweb --bfile $bfile.geno.mind.sex.hwe.auto.het	--genome --min 0.1 --out $bfile.geno.mind.sex.hwe.auto.het
#5	plink --noweb --bfile $bfile.geno.mind.sex.hwe.auto.het --remove $bfile.related.to.remove --make-bed --out $bfile.geno.mind.sex.hwe.auto.het

#6	plink --noweb --bfile $bfile.geno.mind.sex.hwe.auto.het --bmerge /nas40t0/HAPMAP3/HM3_ASN_CEU_YRI_Unrelated.bed /nas40t0/HAPMAP3/HM3_ASN_CEU_YRI_Unrelated.bim /nas40t0/HAPMAP3/HM3_ASN_CEU_YRI_Unrelated.fam --make-bed --out $bfile.HapMap3
#6      plink --noweb --bfile $bfile.geno.mind.sex.hwe.auto.het --exclude $bfile.HapMap3.missnp --make-bed --out $bfile.rm.missnp
#6      plink --noweb --bfile /nas40t0/HAPMAP3/HM3_ASN_CEU_YRI_Unrelated --exclude $bfile.HapMap3.missnp --make-bed --out HapMap3
#6      plink --noweb --bfile $bfile.rm.missnp --indep-pairwise 50 5 0.2 --out $bfile.rm.missnp
#6      plink --noweb --bfile $bfile.rm.missnp --extract $bfile.rm.missnp.prune.in --make-bed --out $bfile.ldpruned
#6      awk '{print $2}' $bfile.ldpruned.bim > $bfile.merge.snplist
#6      plink --noweb --bfile HapMap3 --extract $bfile.merge.snplist --make-bed --out HapMap3
#6      plink --noweb --bfile $bfile.ldpruned --bmerge HapMap3.bed HapMap3.bim HapMap3.fam --recode --out $bfile.HapMap3.forPCA
#6      cut -d" " -f 1-5 $bfile.HapMap3.forPCA.ped > fam
#6      awk '{print $1,$2,$3,$4,$5,1}' fam > $bfile.HapMap3.forPCA.fam

#6      perl make_par_file.pl $bfile.HapMap3.forPCA 0 > $bfile.HapMap3.forPCA.par
#6      smartpca -p $bfile.HapMap3.forPCA.par

#6        #plot in R and choose genetic Europeans (PC1 & PC2 w/in 10 SDs of CEU means) for followup
#6        #make pop vector for legend, color
#6        tail -n +2 $bfile.HapMap3.forPCA.evec > $bfile.HapMap3.evec
#6        awk '{print "GWAS"}' $bfile.ldpruned.fam > $bfile.pop
#6        cat $bfile.pop HapMap3.pop > $bfile.HapMap3.pop
#6        R --vanilla < plot.pca.choose.euros.r --args $bfile.HapMap3
#6        perl split.colon.pl $bfile.HapMap3.euro.GWAS.PCAs > $bfile.HapMap3.euro.keep.list

#With Genetic Europeans:
#7        plink --noweb --bfile $bfile.ldpruned --keep $bfile.HapMap3.euro.keep.list --recode --out $bfile.ldpruned.euro.forPCA
#7        cut -d" " -f 1-5 $bfile.ldpruned.euro.forPCA.ped > fam
#7        awk '{print $1,$2,$3,$4,$5,1}' fam > $bfile.ldpruned.euro.forPCA.fam
#7        perl make_par_file.pl $bfile.ldpruned.euro.forPCA 0 > $bfile.ldpruned.euro.forPCA.par
#7        smartpca -p $bfile.ldpruned.euro.forPCA.par
#7        perl make_par_file_rmoutlier.pl $bfile.ldpruned.euro.forPCA 5 > $bfile.ldpruned.euro.forPCA.rmout.par
#7        smartpca -p $bfile.ldpruned.euro.forPCA.rmout.par

	tail -n +2 $bfile.ldpruned.euro.forPCA.evec > $bfile.ldpruned.euro.evec
        R --vanilla < plot.pca.r --args $bfile.ldpruned.euro
        tail -n +2 $bfile.ldpruned.euro.forPCA.rmout.evec > $bfile.ldpruned.euro.rmout.evec
        R --vanilla < plot.pca.r --args $bfile.ldpruned.euro.rmout

#	plink --noweb --bfile $bfile.geno.mind.sex.hwe.auto.het --keep $bfile.HapMap3.euro.keep.list --make-bed --out /nas40t0/hwheeler/NEURP/IMPUTE2_NEURP/2_prePHASE_shapeit/$bfile.QC



done

