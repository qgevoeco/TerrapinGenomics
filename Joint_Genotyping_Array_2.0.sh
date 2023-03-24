#!/bin/bash
#SBATCH --job-name=J_G_Array_2.0
#SBATCH --ntasks=20
#SBATCH --mem=40GB
#SBATCH --partition=general          # name of partition to submit job
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL              # will send email for begin,end,fail
#SBATCH --mail-user=sbw0033@auburn.edu
#SBATCH --output=J_G_%A_%a._2.0output 		#Changes the output to correspond to each subjob
#SBATCH --error=J_G_%A_%a._2.0error 		#Changes the error to correspond to each subjob
#SBATCH --array=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25

##########################################################################################################
#This script performs a few main tasks
#1) Generates vcfs with only SNPs from a directory containing all Gvcfs
#2) Filters the VCFs to contain only high-confidence variants
#3) Sets up base recalibration
##########################################################################################################

#For testing on an interactive node
#srun -N1 -n4 --pty /bin/bash

#Load modules using specific version to ensure future compatibility
module load gatk/4.1.9.0
module load samtools/1.11
module load bcftools
module load vcftools

#Set Paths
GvcfInputDirectory="/scratch/sbw0033/TerrapinGenomics/Data/Gvcfs"
SampleMapFile="/scratch/sbw0033/TerrapinGenomics/Data/cohort.sample_map_0_WGS"
GenomicDataBaseDirectory="/scratch/sbw0033/TerrapinGenomics/Data/ReducedDatabases"
VCFDirectory="/scratch/sbw0033/TerrapinGenomics/Data/VCFs"
InputBamDirectory="/scratch/sbw0033/TerrapinAllAlignments/ReducedTerrapinFinalAlignments"
RecalDirectory="/scratch/sbw0033/TerrapinGenomics/Data/Recalibration"
NewBamDirectory="/scratch/sbw0033/TerrapinAllAlignments/Validation_Round1_bams"

#cd /scratch/sbw0033/TerrapinGenomics/Data/
#mkdir VCFs
cd $GvcfInputDirectory

#Gotta make a GATK database
gatk GenomicsDBImport --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
  --max-num-intervals-to-import-in-parallel 5 \
  --sample-name-map $SampleMapFile \
  --genomicsdb-workspace-path $GenomicDataBaseDirectory/Chr_${SLURM_ARRAY_TASK_ID} \
  --intervals SUPER_${SLURM_ARRAY_TASK_ID}

##########################################################################################################
#Get a VCF for each chromosome
##########################################################################################################
cd $VCFDirectory
# Get SNPs
gatk GenotypeGVCFs --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
  -R /scratch/sbw0033/TerrapinGenomics/Data/TerrapinReferenceOnlySupers.fasta \
  -V gendb://$GenomicDataBaseDirectory/Chr_${SLURM_ARRAY_TASK_ID} \
  -O $VCFDirectory/genotyped_Chr_${SLURM_ARRAY_TASK_ID}_0.vcf

##########################################################################################################
#Make chromosome level VCFs that only include variants
##########################################################################################################
gatk SelectVariants -R /scratch/sbw0033/TerrapinGenomics/Data/TerrapinReferenceOnlySupers.fasta \
  -V $VCFDirectory/genotyped_Chr_${SLURM_ARRAY_TASK_ID}_0.vcf \
  --select-type-to-include SNP \
  -O JustSNPs_Chr_${SLURM_ARRAY_TASK_ID}_0.vcf

##########################################################################################################
#Filter variants for different parameters
##########################################################################################################
#Worth noting that the output for this KEEPS sites which do not meet the criteria. Sites not up to snuff will only be tagged in the INFO column of the output VCF
gatk VariantFiltration \
  -R /scratch/sbw0033/TerrapinGenomics/Data/TerrapinReferenceOnlySupers.fasta \
	-V JustSNPs_Chr_${SLURM_ARRAY_TASK_ID}_0.vcf \
  -O filtered_Chr_${SLURM_ARRAY_TASK_ID}_0.vcf \
  --filter-name "SOR" \
  --filter-expression "SOR > 3.0" \
  --filter-name "QD" \
  --filter-expression "QD < 2.0" \
  --filter-name "MQ" \
  --filter-expression "MQ < 40.0" \
  --filter-name "MQRankSum" \
  --filter-expression "MQRankSum < -12.5" \
  --filter-name "FS" \
  --filter-expression "FS > 60.0" \
  --filter-name "ReadPosRankSum" \
  --filter-expression "ReadPosRankSum < -5.0"

##########################################################################################################
#Merge the VCFs together (has to be done after running the array- maybe start a new job with this?
##########################################################################################################
#Use bcftools, with a list of filenames (no clue if I can supply a text file, this was how they wrote it on the manual)
#bcftools concat -o /scratch/sbw0033/TerrapinGenomics/Data/total_chroms_0.vcf filtered_Chr_1_0.vcf filtered_Chr_2_0.vcf filtered_Chr_3_0.vcf filtered_Chr_4_0.vcf filtered_Chr_5_0.vcf filtered_Chr_6_0.vcf filtered_Chr_7_0.vcf filtered_Chr_8_0.vcf filtered_Chr_9_0.vcf filtered_Chr_10_0.vcf filtered_Chr_11_0.vcf filtered_Chr_12_0.vcf filtered_Chr_13_0.vcf filtered_Chr_14_0.vcf filtered_Chr_15_0.vcf filtered_Chr_16_0.vcf filtered_Chr_17_0.vcf filtered_Chr_18_0.vcf filtered_Chr_19_0.vcf filtered_Chr_20_0.vcf filtered_Chr_21_0.vcf filtered_Chr_22_0.vcf filtered_Chr_23_0.vcf filtered_Chr_24_0.vcf filtered_Chr_25_0.vcf

#Again, make a concatenated gile containing ALL SITES
#bcftools concat -o /scratch/sbw0033/TerrapinGenomics/Data/total_chroms_invariant_0.vcf genotyped_Chr_1_0.vcf genotyped_Chr_2_0.vcf genotyped_Chr_3_0.vcf genotyped_Chr_4_0.vcf genotyped_Chr_5_0.vcf genotyped_Chr_6_0.vcf genotyped_Chr_7_0.vcf genotyped_Chr_8_0.vcf genotyped_Chr_9_0.vcf genotyped_Chr_10_0.vcf genotyped_Chr_11_0.vcf genotyped_Chr_12_0.vcf genotyped_Chr_13_0.vcf genotyped_Chr_14_0.vcf genotyped_Chr_15_0.vcf genotyped_Chr_16_0.vcf genotyped_Chr_17_0.vcf genotyped_Chr_18_0.vcf genotyped_Chr_19_0.vcf genotyped_Chr_20_0.vcf genotyped_Chr_21_0.vcf genotyped_Chr_22_0.vcf genotyped_Chr_23_0.vcf genotyped_Chr_24_0.vcf genotyped_Chr_25_0.vcf

#Now that we've got our big old vcf, we can delete the old ones
#rm filtered_Chr_${SLURM_ARRAY_TASK_ID}_0.vcf
#rm JustSNPs_Chr_${SLURM_ARRAY_TASK_ID}_0.vcf
#rm genotyped_Chr_${SLURM_ARRAY_TASK_ID}_0.vcf

##########################################################################################################
#Use vcftools to explore coverage, missingness, allele frequencies, and heterozygosity
##########################################################################################################
#Make a table with all the variants to look at in R
#gatk VariantsToTable -V /scratch/sbw0033/TerrapinGenomics/Data/total_chroms.vcf -O /scratch/sbw0033/TerrapinGenomics/Data/Terrapin.sitetable -F CHROM -F POS -F QUAL -F DP -F FS -F QD -F SOR -F MQ -F ReadPosRankSum -F BaseQRankSum -F MQRankSum

#Put it on my computadora
#scp -r sbw0033@easley.auburn.edu:/scratch/sbw0033/Mitochondria_Project/Data/Mito_VCFs/filtered_0.vcf /Users/samweaver/Docs/TerrapinRProject/Data
#scp -r sbw0033@easley.auburn.edu:/scratch/sbw0033/Mitochondria_Project/Data/Mito_VCF_summaries /Users/samweaver/Docs/TerrapinRProject/Data
#scp -r sbw0033@easley.auburn.edu:/scratch/sbw0033/TerrapinGenomics/Data/Terrapin.sitetable /Users/samweaver/Docs/TerrapinRProject/Data
