#!/bin/bash
#SBATCH --job-name=GATKtestArray
#SBATCH --mem=40GB
#SBATCH --ntasks=20
#SBATCH --partition=general          # name of partition to submit job
#SBATCH --time=36:00:00
#SBATCH --mail-type=ALL              # will send email for begin,end,fail
#SBATCH --mail-user=sbw0033@auburn.edu
#SBATCH --output=GATK_test_%A_%a.output 		#Changes the output to correspond to each subjob
#SBATCH --error=GATK_test_%A_%a.error 		#Changes the error to correspond to each subjob
#SBATCH --array=1,14,141,153,17,18,2,22,23,24,44,47,6,60,64,7,76,9,93,96

ReferenceGenome="/scratch/sbw0033/TerrapinGenomics/Data/TerrapinReferenceOnlySupers.fasta"
ReferenceIndexPrefix="/scratch/sbw0033/TerrapinGenomics/Data/TerrapinGenomeBWAIndex/Terrapin_Reduced_Reference_index"
InputBamDirectory="/scratch/sbw0033/TerrapinAllAlignments/ReducedTerrapinFinalAlignments"
OutputGvcfDirectory="/scratch/sbw0033/TerrapinGenomics/Data/Gvcfs"

#Load modules using specific version to ensure future compatibility
module load gatk/4.1.9.0
module load samtools/1.11
module load picard/2.23.9

#Make the fasta index for GATK
#cd /scratch/sbw0033/TerrapinGenomics/Data/
#gzip -d TerrapinReferenceHaplotype1.fasta.gz
#java -Xms2g -Xmx16g -jar /tools/picard-2.23.9/libs/picard.jar CreateSequenceDictionary -R TerrapinReferenceOnlySupers.fasta
#samtools faidx TerrapinReferenceOnlySupers.fasta

########################
#Get into the directory with all the Cleaned reads
cd "$InputBamDirectory"

#Run haplotypeCaller on every single sample individually
gatk HaplotypeCaller --input "$InputBamDirectory"/"${SLURM_ARRAY_TASK_ID}"_0.bam --output "$OutputGvcfDirectory"/"${SLURM_ARRAY_TASK_ID}"_0.g.vcf.gz --reference "$ReferenceGenome" --native-pair-hmm-threads 20 -ERC GVCF

cd "$OutputGvcfDirectory"

#Make a cohort map with samples associated with files
echo "Sample"_"${SLURM_ARRAY_TASK_ID}"'\t'"${SLURM_ARRAY_TASK_ID}"_0.g.vcf.gz >> /home/sbw0033/cohort.sample_map_New_2
