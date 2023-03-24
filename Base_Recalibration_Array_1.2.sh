#!/bin/bash
#SBATCH --job-name=Mito_Recal_1.2
#SBATCH --ntasks=20
#SBATCH --mem=40GB
#SBATCH --partition=general          # name of partition to submit job
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL              # will send email for begin,end,fail
#SBATCH --mail-user=sbw0033@auburn.edu
#SBATCH --output=Mito_Recal_1.2.output 		#Changes the output to correspond to each subjob
#SBATCH --error=Mito_Recal_1.2.error 		#Changes the error to correspond to each subjob

##########################################################################################################
#This script performs a few main tasks on our new dataset
#1) Generates vcfs with only SNPs from a directory containing all Gvcfs
#2) Filters the VCFs to contain only high-confidence variants
#3) Sets up second round base recalibration
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
GenomicDataBaseDirectory="/scratch/sbw0033/Mitochondria_Project/Data/Mito_Database"
VCFDirectory="/scratch/sbw0033/Mitochondria_Project/Data/Mito_VCFs"
ReferenceGenome="/scratch/sbw0033/Mitochondria_Project/Data/Terrapin_Mitochondrial_Reference.fasta"

##########################################################################################################
#Make the GATK database
##########################################################################################################
cd $GvcfInputDirectory
gatk GenomicsDBImport --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
  --max-num-intervals-to-import-in-parallel 20 \
  --sample-name-map $SampleMapFile \
  --genomicsdb-workspace-path $GenomicDataBaseDirectory/Mito_DB_test \
  --intervals KX774423.1 #This is the name of the only contig

##########################################################################################################
#Call genotypes
##########################################################################################################
cd $VCFDirectory
# Get SNPs
gatk GenotypeGVCFs --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
  -R $ReferenceGenome \
  -V gendb://$GenomicDataBaseDirectory/Mito_DB_test \
  -O $VCFDirectory/AllVariants_genotyped_1.vcf

##########################################################################################################
#Make a VCF with only variants
##########################################################################################################
gatk SelectVariants -R $ReferenceGenome \
  -V $VCFDirectory/AllVariants_genotyped_1.vcf \
  --set-filtered-gt-to-nocall \
  -O AllVariants_1.vcf

##########################################################################################################
#Filter variants for different parameters
##########################################################################################################
gatk VariantFiltration \
  -R $ReferenceGenome -V AllVariants_1.vcf -O AllVariants_Filtered_1.vcf \
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
