#!/bin/sh
#$ -pe smp 1                # Request cores
#$ -l avx2                  # faster nodes, also some things like gcta64 only run on this instruction set
#$ -l h_vmem=1G
#$ -l h_rt=240:0:0          # longer runtime as now might also be doing the gender splits
#$ -m a 	             # Send email if job aborted
#$ -M d.vanheel@qmul.ac.uk  # The email address to notify
#$ -cwd 


# wash up and cat together each GWAS

# CHANGE THESE
daterun=2021_05_11
phenoarray=( 	ANA_A2_V2	ANA_B1_V2	ANA_B2_V2	ANA_C2_v2	)
filepath="/data/Blizard-VanHeelLab/GenesandHealth/analysis_hmw208/2021_04_covid19hgiGWAS_2"
prefixarray=( "." )

## prefixarray=( "." "_FemaleOnly." "_MaleOnly." )





cd "$filepath"

for pheno in "${phenoarray[@]}"
do

   for prefix in "${prefixarray[@]}"
   do

   head -n 1 saige2_"$pheno"_chr1_"$daterun""$prefix"gwas > saige3_"$pheno"_"$daterun""$prefix"completegwas.txt


#fix chr 23 output which doesnt have rsID when vcf DS is used not bgen

mv saige2_"$pheno"_chr23_"$daterun""$prefix"gwas saige2_"$pheno"_chr23_"$daterun""$prefix"gwas_columnsshifted

awk <saige2_"$pheno"_chr23_"$daterun""$prefix"gwas_columnsshifted '{print $1,$2,$3,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25}' >saige2_"$pheno"_chr23_"$daterun""$prefix"gwas


      for chr in {1..23}
      do

      tail -n+2 saige2_"$pheno"_chr"$chr"_"$daterun""$prefix"gwas >> saige3_"$pheno"_"$daterun""$prefix"completegwas.txt
      done

   gzip saige3_"$pheno"_"$daterun""$prefix"completegwas.txt
   done

done

# gsutil cp UKBB.Doe.ANA_A2_V2.5.ALL.ALL.EUR.154.1341.SAIGE.20201214.txt.gz gs://covid19-hg-upload-your_study/

cp saige3_ANA_A2_V2_2021_05_11.completegwas.txt.gz GNH.vanHeel.ANA_A2_V2.7.ALL.ALL.SAS.151.34013.SAIGE.20210512.txt.gz
cp saige3_ANA_B1_V2_2021_05_11.completegwas.txt.gz GNH.vanHeel.ANA_B1_V2.7.ALL.ALL.SAS.339.4812.SAIGE.20210512.txt.gz
cp saige3_ANA_C2_v2_2021_05_11.completegwas.txt.gz GNH.vanHeel.ANA_C2_V2.7.ALL.ALL.SAS.5151.29013.SAIGE.20210512.txt.gz
cp saige3_ANA_B2_V2_2021_05_11.completegwas.txt.gz GNH.vanHeel.ANA_B2_V2.7.ALL.ALL.SAS.339.33825.SAIGE.20210512.txt.gz
gsutil cp GNH.vanHeel.ANA_A2_V2.7.ALL.ALL.SAS.151.34013.SAIGE.20210512.txt.gz gs://covid19-hg-upload-gnh/
gsutil cp GNH.vanHeel.ANA_B1_V2.7.ALL.ALL.SAS.339.4812.SAIGE.20210512.txt.gz gs://covid19-hg-upload-gnh/
gsutil cp GNH.vanHeel.ANA_C2_V2.7.ALL.ALL.SAS.5151.29013.SAIGE.20210512.txt.gz gs://covid19-hg-upload-gnh/
gsutil cp GNH.vanHeel.ANA_B2_V2.7.ALL.ALL.SAS.339.33825.SAIGE.20210512.txt.gz gs://covid19-hg-upload-gnh/



