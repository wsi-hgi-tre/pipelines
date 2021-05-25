## this script for BINARY traits 
daterun=2021_05_11

phenoarray=( 	ANA_A2_V2	ANA_B1_V2	ANA_B2_V2	ANA_C2_v2	)

phenofile="/data/Blizard-VanHeelLab/GenesandHealth/analysis_hmw208/2021_04_covid19hgiGWAS_2/2021_CovarsPhenos_34164_forcovid19hgi_v3_deident.txt"

#pruned grm file with 105754 variants 34344 samples
fileforgrm="/data/Blizard-VanHeelLab/GenesandHealth/GSAv3EAMD/Dec2020_interim_36k_TOPMED-r2_Imputation_b38/Dec2020_PlinkHardGenotypedGSAdata_CallRateQConly/GenesAndHealth_GSAv3EAMD_SubProj12345678_Clust4_IDupdated3_251120QCremoved_validNHS_NoDups_No1000Goutliers_Updates12_LDpruned.PrunedForSAIGEGRM.autosome"



## RUN THIS SCRIPT FROM WHERE YOU WANT THE SAIGE ANALYSIS OUTPUT


for pheno in "${phenoarray[@]}"
do

echo >saige1_"$pheno"_"$daterun".qsub    "
#!/bin/sh
#$ -pe smp 8                # Request cores
#$ -l avx2                  # faster nodes, also some things like gcta64 only run on this instruction set
#$ -l h_vmem=8G
#$ -l h_rt=1:0:0          # longer runtime as now also doing the gender splits
#$ -m a 	             # Send email if job aborted
#$ -M d.vanheel@qmul.ac.uk  # The email address to notify
#$ -cwd 

module load singularity/3.5.3 

## ALL genders
singularity exec /data/home/hmw208/saige:0.41.sif step1_fitNULLGLMM.R             \
--plinkFile=$fileforgrm    \
--phenoFile="$phenofile"   \
--phenoCol="$pheno"             \
--covarColList=EthnicTightPCAPakBang_Dec2020,age_at_recruit,age_at_recruit_sq,S1QST_Gender,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20  \
--sampleIDColinphenoFile=SAMPLE_ID      \
--nThreads=\$NSLOTS             \
--LOCO=FALSE      \
--traitType=binary         \
--outputPrefix=saige1_"$pheno"_"$daterun"  


"

qsub saige1_"$pheno"_"$daterun".qsub

done

