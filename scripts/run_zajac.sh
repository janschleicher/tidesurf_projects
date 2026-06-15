#!/bin/bash

WD=$(pwd)
NUM_CORES=10
MAX_MEM_CR=120  # max memory for cellranger in GB
DATA_DIR="../data/zajac"
FASTQ_DIR_CR="${DATA_DIR}/fastq"
PACBIO_DIR="${DATA_DIR}/pacbio"

# Samples
SAMPLES=(
  "G16_315_ccRCC"
  "G16_667_ccRCC"
  "G18_109_ccRCC"
  "G18_626_Normal"
  "G18_626_ccRCC"
)

mkdir -p ${FASTQ_DIR_CR}
mkdir -p ${DATA_DIR}
CELLRANGER_DIR=${DATA_DIR}/cellranger
mkdir -p ${CELLRANGER_DIR}
mkdir -p ${CELLRANGER_DIR}/original
mkdir -p ${CELLRANGER_DIR}/filtered
mkdir -p ${DATA_DIR}/starsolo

# Create splici references
REFERENCE="GRCh38_GENCODEv39"
SPLICI_REF="../data/reference_genomes/GRCh38_GENCODEv39_splici_rl91_ftl5"
SPLICI_IDX="${SPLICI_REF}_idx"
if [ ! -d ${SPLICI_IDX} ]; then
    if [ ! -d ${SPLICI_REF} ]; then
        echo "Creating ${SPLICI_REF}"
        conda run -n alevin --no-capture-output pyroe make-splici ../data/reference_genomes/${REFERENCE}/fasta/genome.fa \
          ../data/reference_genomes/${REFERENCE}/genes/genes.gtf 91 $SPLICI_REF \
          --flank-trim-length 5 --filename-prefix splici
    fi
    conda run -n alevin --no-capture-output salmon index -t ${SPLICI_REF}/splici_fl86.fa -i ${SPLICI_IDX} -p ${NUM_CORES}
fi

# Iterate over samples
for i in {0..4}; do
	SAMPLE=${SAMPLES[i]}
    echo $SAMPLE

    # Run cellranger (first time, all reads)
    cd ${CELLRANGER_DIR}/original || exit
    cellranger count --id="${SAMPLE}" \
        --fastqs="../../fastq/original" \
        --sample="${SAMPLE}" \
        --transcriptome="../../../reference_genomes/${REFERENCE}" \
        --localcores=${NUM_CORES} \
        --localmem=${MAX_MEM_CR}
  
    # Match and filter cell barcodes and UMIs
    cd ${WD} || exit
    conda run -n tidesurf --no-capture-output python3 umi_matching.py \
        -cr "${CELLRANGER_DIR}/original/${SAMPLE}/outs/possorted_genome_bam.bam" \
        -pb "${PACBIO_DIR}/original/${SAMPLE}/scisoseq.mapped.bam" \
        -o "${PACBIO_DIR}/original/${SAMPLE}/" || exit
    conda run -n tidesurf --no-capture-output python3 filter_bam.py \
        -f "${PACBIO_DIR}/original/${SAMPLE}/common_tags.csv" \
        -i "${PACBIO_DIR}/original/${SAMPLE}/scisoseq.mapped.bam" \
        -o "${PACBIO_DIR}/original/${SAMPLE}/scisoseq.mapped_filtered.bam" \
        -pb || exit
    conda run -n tidesurf --no-capture-output python3 filter_bam.py \
        -f "${PACBIO_DIR}/original/${SAMPLE}/common_tags.csv" \
        -i "${CELLRANGER_DIR}/original/${SAMPLE}/outs/possorted_genome_bam.bam" \
        -o "${CELLRANGER_DIR}/original/${SAMPLE}/outs/possorted_genome_bam_filtered.bam" || exit

    # Run IsoSeq pipeline
    cd ${PACBIO_DIR} || exit
    mkdir filtered/${SAMPLE}
    isoseq collapse --log-level INFO \
        --log-file filtered/${SAMPLE}/isoseq_collapse.log \
        -j 7 \
        --keep-non-real-cells \
        --max-5p-diff 1000 \
        --max-3p-diff 100 \
        --alarms filtered/${SAMPLE}/alarms.json \
        original/${SAMPLE}/scisoseq.mapped_filtered.bam \
        filtered/${SAMPLE}/scisoseq.mapped_transcripts.collapse.gff || exit
    pigeon sort --log-level INFO \
        --log-file filtered/${SAMPLE}/pigeon_sort.log \
        -o filtered/${SAMPLE}/scisoseq.mapped_transcripts.sorted.gff \
        filtered/${SAMPLE}/scisoseq.mapped_transcripts.collapse.gff || exit
    pigeon classify --log-level INFO \
        --log-file filtered/${SAMPLE}/pigeon_classify.log \
        -j 7 \
        --out-dir filtered/${SAMPLE}/ \
        --out-prefix scisoseq \
        --flnc filtered/${SAMPLE}/scisoseq.mapped_transcripts.collapse.abundance.txt \
        --ref absolutized.referenceset.xml \
        filtered/${SAMPLE}/scisoseq.mapped_transcripts.sorted.gff || exit
    pigeon filter --log-level INFO \
        --log-file filtered/${SAMPLE}/pigeon_filter.log \
        filtered/${SAMPLE}/scisoseq_classification.txt
    pigeon report --log-level INFO \
        --log-file filtered/${SAMPLE}/pigeon_report.log \
        filtered/${SAMPLE}/scisoseq_classification.txt \
        filtered/${SAMPLE}/scisoseq_saturation.txt || exit
    pigeon make-seurat --keep-ribo-mito-genes \
        --log-level INFO \
        --log-file filtered/${SAMPLE}/pigeon_make_seurat.log \
        --annotations absolutized.referenceset.xml \
        --dedup original/${SAMPLE}/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta \
        --group filtered/${SAMPLE}/scisoseq.mapped_transcripts.collapse.group.txt \
        --out-dir filtered/${SAMPLE}/matrices/ \
        --out-prefix scisoseq \
        filtered/${SAMPLE}/scisoseq_classification.txt || exit
    
    # Generate FASTQ files from filtered BAM files
    cd ${WD} || exit
    cellranger bamtofastq "${CELLRANGER_DIR}/original/${SAMPLE}/outs/possorted_genome_bam_filtered.bam" \
        "${FASTQ_DIR_CR}/filtered/${SAMPLE}" || exit
    mv ${FASTQ_DIR_CR}/filtered/${SAMPLE}/*/* ${FASTQ_DIR_CR}/filtered/${SAMPLE}/
    rm -r ${FASTQ_DIR_CR}/filtered/${SAMPLE}/${SAMPLE}*/


    # Run alevin-fry
    FASTQ_FILES=( ${FASTQ_DIR_CR}/filtered/${SAMPLE}/*_R1_001.fastq.gz )
    echo "${FASTQ_FILES[@]}"
    FASTQ_PREFIX=()
    for f in "${FASTQ_FILES[@]}"; do
        FASTQ_PREFIX+=( "${f::-16}" )
    done
    
    echo "${FASTQ_PREFIX[@]}"

    OUT_DIR="${DATA_DIR}/alevin/${SAMPLE}"
    conda run -n alevin --no-capture-output python3 salmon_script.py \
        --fastq-prefix $(echo "${FASTQ_PREFIX[@]}") \
        --idx-path "${SPLICI_IDX}" \
        --tgm-path "${SPLICI_REF}/splici_fl86_t2g_3col.tsv" \
        --output-dir "${OUT_DIR}" || exit
    conda run -n alevin --no-capture-output python3 make_adata.py \
        -q "${OUT_DIR}/quant_res" \
        -e2g "${SPLICI_REF}/gene_id_to_name.tsv" \
        -of raw_all || exit

    # Run cellranger (second time, filtered reads)
    cd ${CELLRANGER_DIR}/filtered || exit
    cellranger count --id="${SAMPLE}" \
        --fastqs="../../fastq/filtered/${SAMPLE}" \
        --sample=bamtofastq \
        --transcriptome="../../../reference_genomes/${REFERENCE}" \
        --localcores=${NUM_CORES} \
        --localmem=${MAX_MEM_CR} || exit

    # Run tidesurf
    cd ${WD} || exit
    conda run -n tidesurf --no-capture-output tidesurf \
        -o "${DATA_DIR}/tidesurf/${SAMPLE}" \
        --orientation sense \
        --no_filter_cells \
        --export_umi_tables \
        "${DATA_DIR}/cellranger/filtered/${SAMPLE}" \
        ../data/reference_genomes/${REFERENCE}/genes/genes.gtf || exit

    # Run velocyto
    { time conda run -n velocyto --no-capture-output velocyto run10x \
        -@ ${NUM_CORES} \
        "${DATA_DIR}/cellranger/filtered/${SAMPLE}" \
        ../data/reference_genomes/${REFERENCE}/genes/genes_mt.gtf || exit ; } 2> "${DATA_DIR}/cellranger/filtered/${SAMPLE}_time.txt"

    # Remove velocyto sorted BAM file
    rm "${DATA_DIR}/cellranger/filtered/${SAMPLE}/outs/cellsorted_possorted_genome_bam.bam"

    # Run STARsolo
    cd $WD || exit
    cd ${DATA_DIR}/starsolo || exit
    OUT_DIR="${SAMPLE}/"
    mkdir $OUT_DIR
    cd $OUT_DIR || exit

    R1_FILES=( ../../fastq/filtered/${SAMPLE}/*_R1_001.fastq.gz )
    printf -v R1_STRING '%s,' "${R1_FILES[@]}"
    R1_STRING="${R1_STRING%,}"
    R2_FILES=( ../../fastq/filtered/${SAMPLE}/*_R2_001.fastq.gz )
    printf -v R2_STRING '%s,' "${R2_FILES[@]}"
    R2_STRING="${R2_STRING%,}"

    { time 
    STAR --runThreadN ${NUM_CORES} \
        --genomeDir ../../../reference_genomes/GRCh38_STAR_GENCODEv39 \
        --readFilesIn $R2_STRING $R1_STRING \
        --soloType CB_UMI_Simple \
        --clipAdapterType CellRanger4 \
        --soloCBwhitelist ~/cellranger-7.1.0/lib/python/cellranger/barcodes/737K-august-2016.txt \
        --outFilterScoreMin 30 \
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
        --soloUMIfiltering MultiGeneUMI_CR \
        --soloUMIdedup 1MM_CR \
        --soloCellFilter EmptyDrops_CR \
        --soloStrand Forward \
        --soloBarcodeReadLength 0 \
        --soloCBstart 1 \
        --soloCBlen 16 \
        --soloUMIstart 17 \
        --soloUMIlen 12 \
        --outSAMtype None \
        --soloFeatures Gene Velocyto \
        --readFilesCommand zcat  || exit ; } 2> time.txt

        # To create sorted BAM files, add the following two arguments:
        # --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
        # --outSAMtype BAM SortedByCoordinate \
        # Omitted because it crashes due to insufficient RAM (when using FASTQ files from bamtofastq)
        # For unsorted BAM output:
        # --outSAMtype BAM Unsorted \
    cd ${WD} || exit
done
