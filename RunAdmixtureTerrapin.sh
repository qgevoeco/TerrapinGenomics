#!/bin/bash
#SBATCH --job-name=J_G_Array_2.0
#SBATCH --ntasks=14
#SBATCH --mem=20GB
#SBATCH --partition=general          # name of partition to submit job
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL              # will send email for begin,end,fail
#SBATCH --mail-user=sbw0033@auburn.edu
#SBATCH --output=TerrapinAdmixture_%A_%a.0output 		#Changes the output to correspond to each subjob
#SBATCH --error=TerrapinAdmixture_%A_%a.error 		#Changes the error to correspond to each subjob
#SBATCH --array=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25

# Run an ADMIXTURE analysis for various values of K
# From the PCA, it looks like we have three or four major clusters, so we will
# not run all the way up to the 27 basin number.
# Set the limits of the run here
KMIN=2
KMAX=10
THREADS="14"

# Set paths to program, output directory, and input file
ADMIXTURE="/home/mcgaughs/shared/Software/admixture_linux-1.3.0/admixture"
OUT_DIR="/home/mcgaughs/weave271/Reanalyses/ADMIXTURE/Output"
IN_DIR="/home/mcgaughs/weave271/Reanalyses/ADMIXTURE/"

# Call the cd command because ADMIXTURE always writes output to the CWD
mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

# Copy the input files into the CWD. This is only a convenience because we have
# to make the chromosome numbers numeric again
cp "${IN_DIR}/New_pruned_Terrapin.bed" .
cp "${IN_DIR}/New_pruned_Terrapin.fam" .
# Remove the "StacksLocus_" prefix from the chromosome names
sed -e 's/scaffold_//g' "${IN_DIR}/New_pruned_Terrapin.bim" > "New_pruned_Terrapin.bim"

# For each K, run ADMIXTURE
for ((i=${KMIN};i<=${KMAX};i++))
do
    "${ADMIXTURE}" \
        --cv \
        "New_pruned_Terrapin.bed" \
        "${i}" \
        -j"${THREADS}" \
        | tee "ADMIXTURE_K${i}.log"
done
