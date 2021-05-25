#!/bin/sh

daterun=2021_05_11
phenoarray=( 	ANA_A2_V2	ANA_B1_V2	ANA_B2_V2	ANA_C2_v2	)
prefixarray=( "." )
filepath="/data/Blizard-VanHeelLab/GenesandHealth/analysis_hmw208/2021_04_covid19hgiGWAS_2"


for pheno in "${phenoarray[@]}"
do

for prefix in "${prefixarray[@]}"
do

for chr in {1..22}
do
filetoqsub=saige2_"$pheno"_chr"$chr"_"$daterun""$prefix"qsub

echo >"$filetoqsub" "

#!/bin/sh
#$ -pe smp 1               # Request cores. 1 core on nxv was only 25% slower than 8 cores
#$ -l h_rt=240:0:0           # Runtime, alas takes just below/above 1 hour.
#$ -l avx2		     # faster nodes
#$ -l h_vmem=1G
#$ -m a 	             # Send email if aborted
#$ -M d.vanheel@qmul.ac.uk  # The email address to notify
#$ -cwd 
cd "$filepath"

module load singularity/3.5.3

SINGULARITYENV_OMP_NUM_THREADS=1 singularity exec /data/home/hmw208/saige:0.41.sif   \
step2_SPAtests.R   \
--bgenFile=/data/Blizard-VanHeelLab/GenesandHealth/GSAv3EAMD/Dec2020_interim_36k_TOPMED-r2_Imputation_b38/Dec2020_TOPMED-r2_bgen_pgen/chr"$chr".dose.merged_bangla_pak_INFO0.3_MAF0.0001.bgen \
--sampleFile=/data/Blizard-VanHeelLab/GenesandHealth/GSAv3EAMD/Dec2020_interim_36k_TOPMED-r2_Imputation_b38/Dec2020_TOPMED-r2_bgen_pgen/Dec2020_34164_plink2_bgen_samplefileSAIGE.txt     \
--minMAC=5     \
--GMMATmodelFile=saige1_"$pheno"_"$daterun""$prefix"rda    \
--varianceRatioFile=saige1_"$pheno"_"$daterun""$prefix"varianceRatio.txt           \
--numLinesOutput=1000             \
--IsOutputNinCaseCtrl=TRUE \
--IsOutputHetHomCountsinCaseCtrl=TRUE \
--IsOutputAFinCaseCtrl=TRUE       \
--LOCO=FALSE     \
--SAIGEOutputFile=saige2_"$pheno"_chr"$chr"_"$daterun""$prefix"gwas 

"

qsub "$filetoqsub"

qstat | wc -l

# sleep 1s

done


chr=23

filetoqsub=saige2_"$pheno"_chr"$chr"_"$daterun""$prefix"qsub

echo >"$filetoqsub" "

#!/bin/sh
#$ -pe smp 1               # Request cores. 1 core on nxv was only 25% slower than 8 cores
#$ -l h_rt=240:0:0           # Runtime, alas takes just below/above 1 hour.
#$ -l avx2		     # faster nodes
#$ -l h_vmem=1G
#$ -m a 	             # Send email if aborted
#$ -M d.vanheel@qmul.ac.uk  # The email address to notify
#$ -cwd 
cd "$filepath"

module load singularity/3.5.3

SINGULARITYENV_OMP_NUM_THREADS=1 singularity exec /data/home/hmw208/saige:0.41.sif   \
step2_SPAtests.R   \
--vcfFile=/data/Blizard-VanHeelLab/GenesandHealth/GSAv3EAMD/Dec2020_interim_36k_TOPMED-r2_Imputation_b38/Dec2020_TOPMED_Rsq0.1_imputed_vcfs/merged_bangla_pak/chr23.dose.merged_bangla_pak_INFO0.3_MAF0.0001.vcf.gz_HAPLOID/chr23.dose.merged_bangla_pak_INFO0.3_MAF0.0001.vcf.gz   \
--vcfFileIndex=/data/Blizard-VanHeelLab/GenesandHealth/GSAv3EAMD/Dec2020_interim_36k_TOPMED-r2_Imputation_b38/Dec2020_TOPMED_Rsq0.1_imputed_vcfs/merged_bangla_pak/chr23.dose.merged_bangla_pak_INFO0.3_MAF0.0001.vcf.gz_HAPLOID/chr23.dose.merged_bangla_pak_INFO0.3_MAF0.0001.vcf.gz.tbi \
--vcfField=DS \
--chrom=chrX \
--minMAC=5     \
--GMMATmodelFile=saige1_"$pheno"_"$daterun""$prefix"rda    \
--varianceRatioFile=saige1_"$pheno"_"$daterun""$prefix"varianceRatio.txt           \
--numLinesOutput=1000             \
--IsOutputNinCaseCtrl=TRUE \
--IsOutputHetHomCountsinCaseCtrl=TRUE \
--IsOutputAFinCaseCtrl=TRUE       \
--LOCO=FALSE     \
--SAIGEOutputFile=saige2_"$pheno"_chr"$chr"_"$daterun""$prefix"gwas 

"

qsub "$filetoqsub"

done
done

