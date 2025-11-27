#!/bin/bash

WD=$(pwd)
NUM_CORES=10
MAX_MEM_CR=120  # max memory for cellranger in GB
DATA_DIR="../data/fu"
FASTQ_DIR="${DATA_DIR}/fastq"

# Samples, IDs and read lengths
SAMPLES=(
  "SRR23050737"
  "SRR23050726"
  "SRR23050721"
  "SRR23050720"
  "SRR23050719"
  "SRR23050715"
  "SRR23050718"
  "SRR23050738"
  "SRR23050736"
  "SRR23050735"
  "SRR23050717"
  "SRR23050716"
)
IDS=(
  "MJ001"
  "MJ002"
  "MJ003"
  "MJ016"
  "MJ017"
  "MJ005"
  "MJ006"
  "MJ007"
  "MJ008"
  "MJ009"
  "MJ018"
  "MJ019"
)
READ_LENGTHS=(91 91 91 91 91 101 101 101 101 101 101 101)

mkdir -p ${FASTQ_DIR}
mkdir -p ${DATA_DIR}/cellranger
mkdir -p ${DATA_DIR}/starsolo

# Create splici references
REFERENCE="refdata-gex-GRCh38-2020-A"
SPLICI_REF_PREFIX="../data/reference_genomes/GRCh38_splici_rl"
# SPLICI_IDX="${SPLICI_REF}_idx"
for RL in 91 101; do
    FL=$((RL-5))
    SPLICI_REF="${SPLICI_REF_PREFIX}${RL}_ftl5"
    SPLICI_IDX="${SPLICI_REF}_idx"
    if [ ! -d ${SPLICI_IDX} ]; then
        if [ ! -d ${SPLICI_REF} ]; then
            echo "Creating ${SPLICI_REF}"
            conda run -n alevin --no-capture-output pyroe make-splici ../data/reference_genomes/${REFERENCE}/fasta/genome.fa \
              ../data/reference_genomes/${REFERENCE}/genes/genes.gtf $RL $SPLICI_REF \
              --flank-trim-length 5 --filename-prefix splici
        fi
        conda run -n alevin --no-capture-output salmon index -t ${SPLICI_REF}/splici_fl${FL}.fa -i ${SPLICI_IDX} -p ${NUM_CORES}
    fi
done

# Iterate over samples
for i in {0..11}; do
	cd ${FASTQ_DIR} || exit

	SAMPLE=${SAMPLES[i]}
	ID=${IDS[i]}
	READ_LENGTH=${READ_LENGTHS[i]}
	FLANK_LENGTH=$((READ_LENGTH-5))

	# Download bam-file
	wget -N https://sra-pub-src-1.s3.amazonaws.com/${SAMPLE}/${ID}_GEX.bam.1
	mv "${ID}_GEX.bam.1" "${ID}_GEX.bam"

	# Convert to fastq files
	cellranger bamtofastq --nthreads=${NUM_CORES} --reads-per-fastq=9999999999999999999 "${ID}_GEX.bam" "${ID}"
	rm "${ID}_GEX.bam"

	cd "${WD}" || exit
	mv ${FASTQ_DIR}/${ID}/*/*.fastq.gz ${FASTQ_DIR}/${ID}

	# Run alevin-fry
	FASTQ_FILES=( ${FASTQ_DIR}/${ID}/*_R1_001.fastq.gz )
	echo "${FASTQ_FILES[@]}"
	FASTQ_PREFIX=()
	for f in "${FASTQ_FILES[@]}"; do
    FASTQ_PREFIX+=( "${f::-16}" )
  done
	# FASTQ_PREFIX=${FASTQ_FILES[0]::-16}
	echo "${FASTQ_PREFIX[@]}"

	SPLICI_REF="${SPLICI_REF_PREFIX}${READ_LENGTH}_ftl5"
	SPLICI_IDX="${SPLICI_REF}_idx"
	OUT_DIR="${DATA_DIR}/alevin/${ID}"
	conda run -n alevin --no-capture-output python3 salmon_script.py \
      --fastq-prefix $(echo "${FASTQ_PREFIX[@]}") \
      --idx-path "${SPLICI_IDX}" \
      --tgm-path "${SPLICI_REF}/splici_fl${FLANK_LENGTH}_t2g_3col.tsv" \
      --output-dir "${OUT_DIR}" \
      --five-prime
  conda run -n alevin --no-capture-output python3 make_adata.py \
      -q "${OUT_DIR}/quant_res" \
      -e2g "${SPLICI_REF}/gene_id_to_name.tsv" \
      -of raw_all

	# Run cellranger
	cd ${DATA_DIR}/cellranger || exit
	cellranger count --id="${ID}" \
	    --fastqs="../fastq/${ID}" \
	    --sample=bamtofastq \
	    --transcriptome="../../reference_genomes/${REFERENCE}" \
	    --localcores=${NUM_CORES} \
	    --localmem=${MAX_MEM_CR}

	# Run tidesurf
  conda run -n tidesurf --no-capture-output tidesurf \
    -o "${DATA_DIR}/tidesurf/${ID}" \
    --orientation antisense \
  --filter_cells \
  --whitelist cellranger \
    "${DATA_DIR}/cellranger/${ID}" \
    ../data/reference_genomes/${REFERENCE}/genes/genes.gtf

  # Run velocyto
  { time conda run -n velocyto --no-capture-output velocyto run10x \
      -@ ${NUM_CORES} \
      "${DATA_DIR}/cellranger/${ID}" \
      ../data/reference_genomes/${REFERENCE}/genes/genes.gtf ; } 2> "${DATA_DIR}/cellranger/${ID}_time.txt"

  # Remove velocyto sorted BAM file
  rm "${DATA_DIR}/cellranger/${ID}/outs/cellsorted_possorted_genome_bam.bam"

  # Run STARsolo
  cd $WD || exit
  cd ${DATA_DIR}/starsolo || exit
  OUT_DIR="${ID}/"
  mkdir $OUT_DIR
  cd $OUT_DIR || exit

  R1_FILES=( ../../fastq/${ID}/*_R1_001.fastq.gz )
  printf -v R1_STRING '%s,' "${R1_FILES[@]}"
  R1_STRING="${R1_STRING%,}"
  R2_FILES=( ../../fastq/${ID}/*_R2_001.fastq.gz )
  printf -v R2_STRING '%s,' "${R2_FILES[@]}"
  R2_STRING="${R2_STRING%,}"

  if [ $READ_LENGTH -eq 91 ]; then
    { time 
      STAR --runThreadN ${NUM_CORES} \
        --genomeDir ../../../reference_genomes/GRCh38_STAR \
        --readFilesIn $R1_STRING $R2_STRING \
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
        --outSAMtype None \
        --soloFeatures Gene Velocyto \
        --readFilesCommand zcat ; } 2> time.txt

        # To create sorted BAM files, add the following two arguments:
        # --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
        # --outSAMtype BAM SortedByCoordinate \
        # Omitted because it crashes due to insufficient RAM (when using FASTQ files from bamtofastq)
        # For unsorted BAM output:
        # --outSAMtype BAM Unsorted \
  else
    # Paired-end alignment fails with Segmentation fault for samples with read length 101
    # -> single-end alignment
    { time 
      STAR --runThreadN ${NUM_CORES} \
        --genomeDir ../../../reference_genomes/GRCh38_STAR \
        --readFilesIn $R2_STRING $R1_STRING \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist ~/cellranger-7.1.0/lib/python/cellranger/barcodes/737K-august-2016.txt \
        --outFilterScoreMin 30 \
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
        --soloUMIfiltering MultiGeneUMI_CR \
        --soloUMIdedup 1MM_CR \
        --soloCellFilter EmptyDrops_CR \
        --soloStrand Reverse \
        --soloBarcodeReadLength 0 \
        --soloCBstart 1 \
        --soloCBlen 16 \
        --soloUMIstart 17 \
        --soloUMIlen 10 \
        --outSAMtype None \
        --soloFeatures Gene Velocyto \
        --readFilesCommand zcat ; } 2> time.txt
  fi

	# Remove fastqs
  cd "${WD}" || exit
	rm -r "${FASTQ_DIR:?}/${ID:?}"
done
