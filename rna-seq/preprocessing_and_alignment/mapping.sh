#!/bin/bash

#SBATCH -J star # A job-name
#SBATCH -n 16 # Number of cores 
#SBATCH --mem 400000 #(GB per node)
#SBATCH -o ./jobs.out/slurm.%j.out
#SBATCH -e ./jobs.out/slurm.%j.err
#SBATCH --array=0-1%1 # Array job'

: 'Mapping transcriptomic reads to genome with STAR'

#load modules in hpc
module load gzip/1.10-GCCcore-11.2.0 
module load STAR/2.7.10b-GCC-11.3.0
module load MultiQC/1.9-Miniconda3-23.9.0-0 

# global options
DATASET_NAME=$1

# paths
PATH2READS='data/fastqs/trimmed'
PATH2FASTA='data/references/fasta/GRCm38.p6.genome.fa'
PATH2ANNOTATION='data/references/annotation/gencode.vM25.primary_assembly.annotation.gtf'
PATH2GENOME_IDX='data/references/fasta/index'

MAPPING_OUT='star_out'
MULTIQC_OUT='star_out/reports'

# STAR parameters
THREADS=16
# index generation
GENOME_SA_INDEXBASES=12
GENOME_SA_SPARSE=3
#mapping
LOAD_MODE='NoSharedMemory'
READ_FORMAT='zcat'
GENE_ID='gene_name'
FILTER_MULTIMAP=10
OUTSAM_FORMAT='None'
QUANTMODE='GeneCounts' #Transcriptome SAM

cat << \
.
            MAPPING
=========================================

# DATA FOLDER
Raw fastq folder: $PATH2READS

# REFERENCES
Fasta: $(echo $PATH2FASTA | sed 's:.*/::')
Anntotation: $(echo $PATH2ANNOTATION | sed 's:.*/::')

# OUTPUT FOLDERS
Mapping output: $MAPPING_OUT
MultiQC report: $MULTIQC_OUT

# MAPPING PARAMETERS
# Memory
Threads: $THREADS
# memory for index generation
SA pre-indexingn length: $GENOME_SA_INDEXBASES
Indexing distance : $GENOME_SA_SPARSE
#mapping
Genome loading mode: $LOAD_MODE
Read format= $READ_FORMAT
Gene id in annotation: $GENE_ID
Multimappers filter: $FILTER_MULTIMAP
SAM/BAM Output format: $OUTSAM_FORMAT 
STAR quantMode: $QUANTMODE  


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


# 1. GENOME INDEX GENERATION
if [[ ! -f ${PATH2GENOME_IDX}/genomeParameters.txt ]] 
then
    mkdir -p ${PATH2GENOME_IDX}
    echo '
    Starting to index the genome...
    '

    STAR \
    --runMode genomeGenerate \
    --runThreadN $THREADS \
    --genomeDir $PATH2GENOME_IDX \
    --genomeFastaFiles $PATH2FASTA \
    --sjdbGTFfile $PATH2ANNOTATION \
    --genomeSAindexNbases $GENOME_SA_INDEXBASES \
    --genomeSAsparseD $GENOME_SA_SPARSE
else
    echo 'Genome index already computed.'
fi 


# 2.ALIGNMENT
files_R1=(${PATH2READS}/*_val_1.fq.gz)
files_R2=(${PATH2READS}/*_val_2.fq.gz)

# file per index array
file_R1=${files_R1[$SLURM_ARRAY_TASK_ID]} 
file_R2=${files_R2[$SLURM_ARRAY_TASK_ID]}

condition
condition=$(echo $file_R2 | sed 's:.*/::' | cut -d '.' -f 1 | cut -d '_' -f 1)

# if [[ ! -f $MAPPING_OUT/$condition/${condition}_Log.final.out ]]; then
#     echo '
#     Starting to map '${condition} 'reads
#     '
#     STAR \
#     --genomeDir $PATH2GENOME_IDX \
#     --genomeLoad $LOAD_MODE \
#     --runThreadN $THREADS \
#     --readFilesCommand $READ_FORMAT \
#     --readFilesIn $file_R1 $file_R2 \
#     --sjdbGTFfile $PATH2ANNOTATION \
#     --sjdbGTFtagExonParentGene $GENE_ID \
#     --outFilterMultimapNmax $FILTER_MULTIMAP \
#     --outFileNamePrefix $MAPPING_OUT/$condition/${condition}_isoform_ \
#     --outSAMtype $OUTSAM_FORMAT \
#     --quantMode $QUANTMODE

# fi

if [[ ! -f $MAPPING_OUT/$condition/${condition}_isoform_Log.final.out ]]; then
    echo '
    Starting to map '${condition} 'reads
    '
    STAR \
    --genomeDir $PATH2GENOME_IDX \
    --genomeLoad $LOAD_MODE \
    --runThreadN $THREADS \
    --readFilesCommand $READ_FORMAT \
    --readFilesIn $file_R1 $file_R2 \
    --sjdbGTFfile $PATH2ANNOTATION \
    --sjdbGTFtagExonParentGene $GENE_ID \
    --outFilterMultimapNmax $FILTER_MULTIMAP \
    --outFileNamePrefix $MAPPING_OUT/$condition/${condition}_isoform_ \
    --outSAMtype $OUTSAM_FORMAT \
    --quantMode $QUANTMODE

fi
echo 'Mapping done for '${condition}

# store mapping report in MultiQC folder
mkdir -p ${MULTIQC_OUT}

# for res in $MAPPING_OUT/*/*_Log.final.out; do
#     condition=$(basename $(dirname "$res"))  # Extracts the subdirectory name as the condition
#     cp "${res}" "$MULTIQC_OUT/${condition}_Log.final.out"
# done

# mkdir -p ${MAPPING_OUT}/counts
# for res in $MAPPING_OUT/*/*_ReadsPerGene.out.tab; do
#     condition=$(basename $(dirname "$res"))  # Extracts the subdirectory name as the condition
#     cp "${res}" "$MAPPING_OUT/counts/${condition}_counts.tab"
# done

# run MultiQC on alignment reports
source activate multiqc_1.9
multiqc $MULTIQC_OUT/* -n $DATASET_NAME -f -s -o $MULTIQC_OUT/

echo '
==> STAR ALIGNMENT FINISHED <==
--------------------------------

Number of files analyzed: '$(ls ${PATH2READS}/*R1* -1 | wc -l)'
Alignment output stored in: '${MAPPING_OUT}'
MultiQC and STAR reports stored in: '${MULTIQC_OUT}'
'


