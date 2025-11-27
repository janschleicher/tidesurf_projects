#!/bin/bash

WD=$(pwd)
NUM_CORES=12
MAX_MEM_CR=120  # max memory for cellranger in GB
DATA_DIR="../data/stewart"
FASTQ_DIR="${DATA_DIR}/fastq"
SAMPLES=("Naive" "Transitional" "IgM-Memory" "Classical-Memory" "DN")

# Download the data if its not already there
FTP_FOLDERS=(
  "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR975/ERR9751750"
  "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR975/ERR9751751"
  "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR975/ERR9751752"
  "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR975/ERR9751753"
  "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR975/ERR9751754"
)
if [ ! -d "${FASTQ_DIR}" ]; then
    echo "Creating ${FASTQ_DIR}"
    mkdir -p ${FASTQ_DIR}
fi
for FOLDER in "${FTP_FOLDERS[@]}"; do
    wget -N -P ${FASTQ_DIR} "${FOLDER}/*"
done

# Create splici reference
REFERENCE="refdata-gex-GRCh38-2020-A"
SPLICI_REF="../data/reference_genomes/GRCh38_splici_rl100_ftl5"
SPLICI_IDX="${SPLICI_REF}_idx"
if [ ! -d ${SPLICI_IDX} ]; then
    if [ ! -d ${SPLICI_REF} ]; then
        echo "Creating ${SPLICI_REF}"
        conda run -n alevin --no-capture-output pyroe make-splici ../data/reference_genomes/${REFERENCE}/fasta/genome.fa \
          ../data/reference_genomes/${REFERENCE}/genes/genes.gtf 100 $SPLICI_REF \
          --flank-trim-length 5 --filename-prefix splici
    fi
    conda run -n alevin --no-capture-output salmon index -t ${SPLICI_REF}/splici_fl95.fa -i ${SPLICI_IDX} -p ${NUM_CORES}
fi

# Run alevin-fry
for i in {1..5}; do
    SAMPLE=${SAMPLES[$i-1]}
    OUT_DIR=${DATA_DIR}/alevin/${SAMPLE}
    conda run -n alevin --no-capture-output python3 salmon_script.py \
      --fastq-prefix=${FASTQ_DIR}/${SAMPLE}_S${i}_L00${i} \
      --idx-path=${SPLICI_IDX} \
      --tgm-path=${SPLICI_REF}/splici_fl95_t2g_3col.tsv \
      --output-dir=${OUT_DIR} \
      --five-prime
    conda run -n alevin --no-capture-output python3 make_adata.py \
      -q ${OUT_DIR}/quant_res \
      -e2g ${SPLICI_REF}/gene_id_to_name.tsv \
      -of raw_all
done

# Run Cell Ranger
mkdir ${DATA_DIR}/cellranger
cd ${DATA_DIR}/cellranger || exit

for SAMPLE in "${SAMPLES[@]}"; do
    cellranger count --id=$SAMPLE \
        --transcriptome="../../reference_genomes/${REFERENCE}" \
        --fastqs=../fastq \
        --sample=$SAMPLE \
        --localcores=${NUM_CORES} \
        --localmem=${MAX_MEM_CR}
done

# Run tidesurf
cd ${WD} || exit
conda run -n tidesurf --no-capture-output python3 run_tidesurf_multisample.py \
  -o ${DATA_DIR}/tidesurf \
  --orientation antisense \
  --filter_cells \
  --whitelist cellranger \
  ${DATA_DIR}/cellranger \
  ../data/reference_genomes/${REFERENCE}/genes/genes.gtf

# Run velocyto
for SAMPLE in "${SAMPLES[@]}"; do
    { time conda run -n velocyto --no-capture-output velocyto run10x \
      -@ ${NUM_CORES} \
      ${DATA_DIR}/cellranger/${SAMPLE} \
      ../data/reference_genomes/${REFERENCE}/genes/genes.gtf ; } 2> ${DATA_DIR}/cellranger/${SAMPLE}_time.txt
done

# Run STARsolo
mkdir ${DATA_DIR}/starsolo
cd ${DATA_DIR}/starsolo || exit

for i in {1..5}; do
    SAMPLE=${SAMPLES[$i-1]}
    OUT_DIR="${SAMPLE}/"
    mkdir $OUT_DIR
    cd $OUT_DIR || exit

    { time 
    STAR --runThreadN ${NUM_CORES} \
      --genomeDir ../../../reference_genomes/GRCh38_STAR \
      --readFilesIn ../../fastq/${SAMPLE}_S${i}_L00${i}_R1_001.fastq.gz ../../fastq/${SAMPLE}_S${i}_L00${i}_R2_001.fastq.gz \
      --soloType CB_UMI_Simple \
      --soloCBwhitelist ~/cellranger-7.1.0/lib/python/cellranger/barcodes/737K-august-2016.txt \
      --outFilterScoreMin 30 \
      --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
      --soloUMIfiltering MultiGeneUMI_CR \
      --soloUMIdedup 1MM_CR \
      --soloCellFilter EmptyDrops_CR \
      --soloBarcodeMate 1 \
      --clip5pNbases 39 0 \
      --soloCBstart 1 \
      --soloCBlen 16 \
      --soloUMIstart 17 \
      --soloUMIlen 10 \
      --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
      --outSAMtype BAM SortedByCoordinate \
      --soloFeatures Gene Velocyto \
      --readFilesCommand zcat ; } 2> time.txt
    cd ..
done
