#!/bin/bash

WD=$(pwd)
NUM_CORES=50
MAX_MEM_CR=700  # max memory for cellranger in GB
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