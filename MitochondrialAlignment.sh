#!/bin/bash
#SBATCH --job-name=Mito_Alignment
#SBATCH --ntasks=14
#SBATCH --partition=general          # name of partition to submit job
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL              # will send email for begin,end,fail
#SBATCH --mail-user=sbw0033@auburn.edu
#SBATCH --output=Mito_Alignment_%A_%a.output 	#Changes the output to correspond to each subjob
#SBATCH --error=Mito_Alignment_%A_%a.error 	#Changes the error to correspond to each subjob
#SBATCH --array=1,14,141,153,17,18,2,22,23,24,44,47,6,60,64,7,76,9,93

RawReadDirectory="/scratch/sbw0033/TerrapinGenomics/Data/FastqCopiesFinal/clean"
ReferenceIndexPrefix="/scratch/sbw0033/Mitochondria_Project/Data/Mitochondiral_BWA_Index/Terrapin_Mitochondrial_Reference_index"
AlignmentDirectory="/scratch/sbw0033/Mitochondria_Project/Data/Alignments/RawAlignments"
SortedBamDirectory="/scratch/sbw0033/Mitochondria_Project/Data/Alignments/SortedAlignments"
CleanedBamDirectory="/scratch/sbw0033/Mitochondria_Project/Data/Alignments/CleanedAlignments"
FinalAlignmentDirectory="/scratch/sbw0033/Mitochondria_Project/Data/Alignments/FinalAlignments"

#Get that lil genome onto easley
#scp -r /Users/samweaver/Docs/TerrapinRProject/Data/Terrapin_Mitochondrial_Reference.fasta sbw0033@easley.auburn.edu:/scratch/sbw0033/Mitochondria_Project/Data #Transfer the annotation

#/scratch/sbw0033/Mitochondria_Project/Data

#Load modules using specific version to ensure future compatibility
module load bwa/0.7.17
module load samtools/1.11
module load picard/2.23.9
module load java/15.0.1
module load bcftools

#index the reference genome
#cd /scratch/sbw0033/Mitochondria_Project/Data/Terrapin_Mitochondrial_Reference.fasta
#bwa index Terrapin_Mitochondrial_Reference.fasta -p Terrapin_Mitochondrial_Reference_index
#mv Terrapin_Mitochondrial_Reference_index.* /scratch/sbw0033/Mitochondria_Project/Data/Mitochondiral_BWA_Index/

############################################################################################################
#Run BWA
############################################################################################################
#Get into the directory with all the raw reads
cd "$RawReadDirectory"
bwa mem -M -t 14 \
	"${ReferenceIndexPrefix}" \
	"${RawReadDirectory}"/"${SLURM_ARRAY_TASK_ID}"/"${SLURM_ARRAY_TASK_ID}"_FW_reads.fq.gz "${RawReadDirectory}"/"${SLURM_ARRAY_TASK_ID}"/"${SLURM_ARRAY_TASK_ID}"_RV_reads.fq.gz \
	2> /scratch/sbw0033/bwaSample_"${SLURM_ARRAY_TASK_ID}".err \
	> "${AlignmentDirectory}"/"${SLURM_ARRAY_TASK_ID}".sam

############################################################################################################
#Convert sams, ditch unmapped reads, sort alignments
############################################################################################################
cd "${AlignmentDirectory}"
samtools view -bS -@ 14 "${SLURM_ARRAY_TASK_ID}".sam > "${SLURM_ARRAY_TASK_ID}".bam #Converts to .bam
#samtools view -bS -@ 14 1.sam > 1.bam #Converts to .bam
samtools view -b -F 4 "${SLURM_ARRAY_TASK_ID}".bam > Mapped_only_"${SLURM_ARRAY_TASK_ID}".bam #Discards unmapped reads
samtools sort "${AlignmentDirectory}"/Mapped_only_${SLURM_ARRAY_TASK_ID}.bam --threads 14 -o "${SortedBamDirectory}"/"${SLURM_ARRAY_TASK_ID}"_sorted.bam #Sorts bam
#Delete intermediate files
rm Mapped_only_"${SLURM_ARRAY_TASK_ID}".bam
rm "${SLURM_ARRAY_TASK_ID}".bam

############################################################################################################
#Now, we want to get rid of any bad alignments from our sorted bam files
############################################################################################################
cd "${SortedBamDirectory}"
#Now, keep alignments with a map quality score greater than or equal to 30
samtools view -@ 14 -q 30 -b "${SLURM_ARRAY_TASK_ID}"_sorted.bam > ${SLURM_ARRAY_TASK_ID}_sorted.q30.bam
#Remove secondary alignments
samtools view -@ 14 -h -F 0x900 "${SLURM_ARRAY_TASK_ID}"_sorted.q30.bam > "${CleanedBamDirectory}"/"${SLURM_ARRAY_TASK_ID}"_sorted.q30.primary_only.bam
#Delete intermediate files
rm ${SLURM_ARRAY_TASK_ID}_sorted.q30.bam

############################################################################################################
#Use picard to add read group info and mark duplicates
############################################################################################################
cd "${CleanedBamDirectory}"
java -Xms2g -Xmx16g -jar /tools/picard-2.23.9/libs/picard.jar AddOrReplaceReadGroups \
	I="${SLURM_ARRAY_TASK_ID}"_sorted.q30.primary_only.bam \
	O="${SLURM_ARRAY_TASK_ID}"_IDed.bam \
	RGID=ID_"${SLURM_ARRAY_TASK_ID}" \
	RGLB="lib1" \
	RGPL="ILLUMINA" \
	RGPU=Barcode_"${SLURM_ARRAY_TASK_ID}" \
	RGSM=RGSample_"${SLURM_ARRAY_TASK_ID}"
#Use picard to mark duplicates
java -Xms2g -Xmx16g -jar /tools/picard-2.23.9/libs/picard.jar MarkDuplicates I="${SLURM_ARRAY_TASK_ID}"_IDed.bam O="${FinalAlignmentDirectory}"/"${SLURM_ARRAY_TASK_ID}"_0.bam M="${FinalAlignmentDirectory}"/"${SLURM_ARRAY_TASK_ID}"_marked_dup_metrics.txt CREATE_INDEX=true REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE

############################################################################################################
#Add indexes to all bam files
############################################################################################################
cd "${FinalAlignmentDirectory}"
#Index the bam file
samtools index -@ 14 "${SLURM_ARRAY_TASK_ID}"_0.bam
#Delete intermediate files
rm ${CleanedBamDirectory}/"${SLURM_ARRAY_TASK_ID}"_IDed.bam

############################################################################################################
#Generate concensus sequence from our samples
############################################################################################################
mkdir /scratch/sbw0033/TerrapinGenomics/Data/Concensus_SQs/
# Get consensus fastq file (vcfutils.pl is part of bcftools)
#samtools mpileup -uf /scratch/sbw0033/Mitochondria_Project/Data/Terrapin_Mitochondrial_Reference.fasta "${SLURM_ARRAY_TASK_ID}"_0.bam | bcftools call -c | vcfutils.pl vcf2fq > /scratch/sbw0033/TerrapinGenomics/Data/Concensus_SQs/${SLURM_ARRAY_TASK_ID}_cns.fastq
# Convert .fastq to .fasta and set bases of quality lower than 20 to N
#cd /scratch/sbw0033/TerrapinGenomics/Data/Concensus_SQs/
#seqtk seq -aQ64 -q20 -n N ${SLURM_ARRAY_TASK_ID}_cns.fastq > ${SLURM_ARRAY_TASK_ID}_cns.fasta
