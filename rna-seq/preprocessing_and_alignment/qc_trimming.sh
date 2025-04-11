#!/bin/bash

#SBATCH -J qc_trim # A job-name
#SBATCH -n 32 # Number of cores
#SBATCH --mem 400000 # Memory request (80GB)
#SBATCH -o ./logs.out/slurm.%A_%a.out # Output file, indexed by array ID
#SBATCH -e ./logs.out/slurm.%A_%a.err # Error file, indexed by array ID
#SBATCH --array=0-1%1 # Array job: 16'

module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-Miniconda3-23.9.0-0 
module load Trim_Galore/0.6.7-GCCcore-11.3.0

# Input parameters
DATASET2ANALYZE=$1

# Paths 
PATH2RAW='data/fastqs/raw' # Path to raw files
PATH2TRIMMED='data/fastqs/trimmed' # Path to trimmed files
PATH2FASTQC='data/fastqs/fastqc' # Path to FastQC reports

mkdir -p $PATH2FASTQC
mkdir -p $PATH2TRIMMED

# Parameters
TRIM_Q=20 # Read quality threshold
TRIM_TYPE='--illumina'

cat << \
.
            QC & TRIMMING
=========================================

# PIPELINE GLOBAL OPTIONS
- Dataset Name: $DATASET2ANALYZE

# DATA FOLDER
Raw fastq folder: $PATH2RAW

# OUTPUT FOLDERS
Trimmed reads output folder: $PATH2TRIMMED
FastQC and MultiQC reports folder: $PATH2FASTQC

# TRIM GALORE PARAMETERS
Quality threshold: $TRIM_Q
Adapter type: $TRIM_TYPE

=========================================
.
\

if [ $# -ne 1 ]
then
    cat << \
.
    echo "Please, give:"
    echo "1) Naming for the dataset to analyze"   
.
\

fi

# Get list of files
files_R1=(${PATH2RAW}/F0004*R1_001.fastq.gz)
files_R2=(${PATH2RAW}/F0004*R2_001.fastq.gz)

# file per index array
file_R1=${files_R1[$SLURM_ARRAY_TASK_ID]} 
file_R2=${files_R2[$SLURM_ARRAY_TASK_ID]}

sample=$(basename $file | sed 's/_R[12]_001.*//')

# FastQC Raw
echo "Performing FastQC on ${sample}"
fastqc $file_R1 -o $PATH2FASTQC
fastqc $file_R2 -o $PATH2FASTQC
echo "FastQC done for ${sample}"


# Trimming and FastQC

echo "Trimming ${file} using Trim Galore!"
trim_galore --quality "$TRIM_Q" "$TRIM_TYPE" --paired "$file_R1" "$file_R2" -o "$PATH2TRIMMED" --fastqc_args "-o $PATH2FASTQC"
echo "Trimming done for ${sample}"



# MultiQC
if [ $SLURM_ARRAY_TASK_ID -eq $SLURM_ARRAY_TASK_MAX ] && [ "$DO_FASTQC_RAW" = true -o "$DO_FASTQC_TRIMMED" = true ]
then
echo "Merging FastQC reports with MultiQC"
source activate multiqc_1.9
multiqc $PATH2FASTQC/* -n ${DATASET2ANALYZE}_all -f -o $PATH2FASTQC 
conda deactivate
echo "MultiQC analysis finished!"
fi

echo "Merging FastQC reports with MultiQC"
source activate multiqc_1.9
multiqc $PATH2TRIMMED/*.txt -n ${DATASET2ANALYZE}_trimming -f -o $PATH2FASTQC 
conda deactivate
echo "MultiQC analysis finished!"