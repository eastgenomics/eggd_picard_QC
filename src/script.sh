#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

# err() - Output error messages to STDERR
#
# Arguments:
#   $* - The error message to display
err() {
    # Write error message to STDERR
    echo "$*" >&2
}

# create_interval_file() - Converts a BED file to a Picard Interval List. 
# See https://gatk.broadinstitute.org/hc/en-us/articles/360037593251-BedToIntervalList-Picard for details
#
# Arguments:
#   $1 - Path to BED file
#   $2 - Path to sorted BAM file
#   $3 - Output filename
#   $4 - Java maximum heap size
function create_interval_file() {
    local BEDFILE_PATH=$1
    local SORTED_BAM=$2
    local OUTPUT_TARGETS=$3
    local MAXHEAP=$4

	docker exec picard_image java -Xmx"${MAXHEAP}" -jar /usr/picard/picard.jar BedToIntervalList \
	    -I="${BEDFILE_PATH}" \
	    -O="${OUTPUT_TARGETS}" \
	    -SD="${SORTED_BAM}"
}

# collect_targeted_pcr_metrics() - Collects targeted PCR metrics
# See https://gatk.broadinstitute.org/hc/en-us/articles/360037438131-CollectTargetedPcrMetrics-Picard for details 
#
# Arguments:
#   $1 - Path to sorted BAM file
#   $2 - Path to reference genome
#   $3 - Path to picard interval list
#   $4 - Path to output directory
#   $5 - Java maximum heap size
collect_targeted_pcr_metrics() {
    local SORTED_BAM=$1
    local REF_GENOME=$2
    local TARGETS_FILE=$3
    local OUTPUT_DIR=$4
    local MAXHEAP=$5

    local SORTED_BAM_PREFIX
    SORTED_BAM_PREFIX=$(basename "${SORTED_BAM}" .bam)

    docker exec picard_image java -Xmx"${MAXHEAP}" -jar /usr/picard/picard.jar CollectTargetedPcrMetrics  \
        -I="${SORTED_BAM}" \
        -R="${REF_GENOME}" \
	    -O="${OUTPUT_DIR}/${SORTED_BAM_PREFIX}.targetPCRmetrics.txt" \
        -AI="${TARGETS_FILE}" \
        -TI="${TARGETS_FILE}" \
        --PER_TARGET_COVERAGE="${OUTPUT_DIR}/${SORTED_BAM_PREFIX}.perTargetCov.txt"
}

# collect_multiple_metrics() - Collected multiple metrics. Note that not all outputs are relevant for all type of sequencing.
# See https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectMultipleMetrics-Picard for details
#
# Arguments:
#   $1 - Path to sorted BAM file
#   $2 - Path to reference genome
#   $3 - Path to output directory
#   $4 - Java maximum heap size
#
# TODO: Investigate CollectSequencingArtifactMetrics - it errors out with TSO500 BAMs due to
# "Record contains library that is missing from header", and so is not used (fix unclear)
collect_multiple_metrics() {
    local SORTED_BAM=$1
    local REF_GENOME=$2
    local OUTPUT_DIR=$3
    local MAXHEAP=$4

    local SORTED_BAM_PREFIX
    SORTED_BAM_PREFIX=$(basename "${SORTED_BAM}" .bam)

    docker exec picard_image java -Xmx"${MAXHEAP}" -jar /usr/picard/picard.jar CollectMultipleMetrics \
        -I="${SORTED_BAM}" \
        -R="${REF_GENOME}" \
	    --PROGRAM=null \
	    --PROGRAM=CollectAlignmentSummaryMetrics \
	    --PROGRAM=CollectInsertSizeMetrics \
	    --PROGRAM=QualityScoreDistribution \
	    --PROGRAM=MeanQualityByCycle \
	    --PROGRAM=CollectBaseDistributionByCycle \
	    --PROGRAM=CollectGcBiasMetrics \
	    --PROGRAM=CollectQualityYieldMetrics \
	    -O="${OUTPUT_DIR}/${SORTED_BAM_PREFIX}"
}

# collect_hs_metrics() - Collect hybrid-selection (HS) metrics.
# See https://gatk.broadinstitute.org/hc/en-us/articles/360036856051-CollectHsMetrics-Picard for details
#
# Arguments:
#   $1 - Path to sorted BAM file
#   $2 - Path to targets file
#   $3 - Path to reference genome
#   $4 - Path to output directory
#   $5 - Java maximum heap size
collect_hs_metrics() {
    local SORTED_BAM=$1
    local TARGETS_FILE=$2
    local REF_GENOME=$3
    local OUTPUT_DIR=$4
    local MAXHEAP=$5

    local SORTED_BAM_PREFIX
    SORTED_BAM_PREFIX=$(basename "${SORTED_BAM}" .bam)

    docker exec picard_image java -Xmx"${MAXHEAP}" -jar /usr/picard/picard.jar CollectHsMetrics \
        --BI="${TARGETS_FILE}" \
        --TI="${TARGETS_FILE}" \
        --I="${SORTED_BAM}" \
        --O="${OUTPUT_DIR}/${SORTED_BAM_PREFIX}.hsmetrics.tsv" \
        --R="${REF_GENOME}" \
        --PER_TARGET_COVERAGE="${OUTPUT_DIR}/${SORTED_BAM_PREFIX}.pertarget_coverage.tsv"\
        --COVERAGE_CAP=100000
}

# collect_rnaseq_metrics() - Collect RNA-seq metrics
# See https://gatk.broadinstitute.org/hc/en-us/articles/360037057492-CollectRnaSeqMetrics-Picard for details
#
# Arguments:
#   $1 - Path to sorted BAM file
#   $2 - Path to refFlat file
#   $3 - Path to output directory
#   $4 - Java maximum heap size
collect_rnaseq_metrics() {
    local SORTED_BAM=$1
    local REF_FLAT=$2
    local OUTPUT_DIR=$3
    local MAXHEAP=$4

    local SORTED_BAM_PREFIX
    SORTED_BAM_PREFIX=$(basename "${SORTED_BAM}" .bam)

    docker exec picard_image java -Xmx"${MAXHEAP}" -jar /usr/picard/picard.jar CollectRnaSeqMetrics \
        -I="${SORTED_BAM}" \
        -O="${OUTPUT_DIR}/${SORTED_BAM_PREFIX}.RNAmetrics.tsv" \
        --REF_FLAT="${REF_FLAT}" \
        -STRAND=SECOND_READ_TRANSCRIPTION_STRAND
}

# collect_variant_calling_metrics() - Collect variant calling metrics
# See https://gatk.broadinstitute.org/hc/en-us/articles/360037057132-CollectVariantCallingMetrics-Picard for details
#
# Arguments:
#   $1 - Path to VCF
#   $2 - Path to dbSNP VCF
#   $3 - Path to sequence dictionary file
#   $4 - Path to output directory
#   $5 - Java maximum heap size
collect_variant_calling_metrics() {
    local VCF=$1
    local DBSNP_VCF=$2
    local SEQ_DICT=$3
    local OUTPUT_DIR=$4
    local MAXHEAP=$5

    local VCF_PREFIX
    VCF_PREFIX=$(basename "${VCF}" .vcf.gz)

    docker exec picard_image java -Xmx"${MAXHEAP}" -jar /usr/picard/picard.jar CollectVariantCallingMetrics \
        --DBSNP="${DBSNP_VCF}" \
        --INPUT="${VCF}" \
        --OUTPUT="${OUTPUT_DIR}/${VCF_PREFIX}.variantcallingmetrics" \
        --SEQUENCE_DICTIONARY="${SEQ_DICT}" \
        --GVCF_INPUT true
}

# make_refflat() - generates refFlat file 
# See https://biopet.github.io/gtftorefflat/develop/ for details
#
# Arguments:
#   $1 - path to input GTF
#   $2 - path to output refFlat file
#   $3 - Java maximum heap size
make_refflat() {
    local GTF=$1
    local OUTPUT_REFFLAT=$2
    local MAXHEAP=$3

    docker exec picard_image java -Xmx"${MAXHEAP}" -jar /usr/bin/GtftoRefflat-assembly-0.1.jar \
        -g "${GTF}" \
        -r "${OUTPUT_REFFLAT}"
}

main() {
    ## Sanity checks
    # Exit if no picard functions selected
    if [[ "$run_CollectTargetedPcrMetrics" == "false" ]] && \
        [[ "$run_CollectHsMetrics" == "false" ]] && \
        [[ "$run_CollectMultipleMetrics" == "false" ]] && \
        [[ "$run_CollectRnaSeqMetrics" == "false" ]] && \
        [[ "$run_CollectVariantCallingMetrics" == "false" ]]; then
        err "No picard functions selected!"
        exit 1
    fi

    # Exit if input args don't align with selected functions
    if [[ ( "$run_CollectTargetedPcrMetrics" == "true" || \
            "$run_CollectHsMetrics" == "true" || \
            "$run_CollectMultipleMetrics" == "true" ) && \
            ( -z "$sorted_bam" || -z "$fasta_index" || -z "$bedfile" ) ]] ; then
        err "One of run_CollectTargetedPcrMetrics, run_CollectHsMetrics or run_CollectMultipleMetrics was requested, but one or more of sorted_bam, fasta_index or bedfile are missing. Exiting..."
        exit 1
    fi

    if [[ "$run_CollectRnaSeqMetrics" == "true" && -z "$sorted_bam" ]]; then
        err "run_CollectRnaSeqMetrics was requested, but sorted_bam is missing. Exiting..."
        exit 1
    fi

    if [[ "$run_CollectVariantCallingMetrics" == "true" && \
        ( -z "$vcf" || -z "$dbsnp_vcf" ) ]]; then
        err "run_CollectVariantCallingMetrics was requested, but one or more of vcf or dbsnp_vcf are missing. Exiting..."
        exit 1
    fi

    ## Setup 
    ### puts inputs in /home/dnanexus/in/
    dx-download-all-inputs

    # Calculate 90% of memory size for java
    MEM=$(head -n1 /proc/meminfo | awk '{print int($2*0.9)}')
    MEM_IN_MB="$(("${MEM}"/1024))m"

    tar zxvf "$fasta_index_path" -C in/ 
    OUTPUT_DIR="${HOME}/out/eggd_picard_stats/QC"
    mkdir -p "$OUTPUT_DIR"

    # Start container process in background, and use docker exec to run multiple commands.
    # This is because we need the docker container to persist between steps
    ## note that mount point for GTF needs to be made first
    GENOME_BUILD_NAME=$(basename $fasta_index_name .tar.gz)
    CTAT_LIBDIR="/home/dnanexus/${GENOME_BUILD_NAME}/ctat_genome_lib_build_dir"

    DOCKER_IMAGENAME=$(find /image -name "*.tar.gz")
    sudo docker load -i "${DOCKER_IMAGENAME}"
    DOCKER_IMAGE=$(docker image ls -q)
    docker run \
        --name picard_image \
        --mount type=bind,source=/home/dnanexus/in/,target=/in \
        --mount type=bind,source="${CTAT_LIBDIR}",target=/ctat_libdir \
        --mount type=bind,source="${OUTPUT_DIR}",target=/out \
        --entrypoint /bin/bash \
        -itd "${DOCKER_IMAGE}"

    # Create the interval file if required
    if [[ "$run_CollectMultipleMetrics" == true ]] || \
        [[ "$run_CollectHsMetrics" == true ]] || \
        [[ "$run_CollectTargetedPcrMetrics" == true ]]; then
        echo "Generating interval file"
        create_interval_file "/in/${bedfile_name}" "/in/${sorted_bam_name}" "/in/targets.picard" "${MEM_IN_MB}"
    fi

    ## Run picard commands
    if [[ "$run_CollectMultipleMetrics" == true ]]; then
        collect_multiple_metrics "/in/${sorted_bam_name}" "/in/genome.fa" "/out/" "${MEM_IN_MB}"
    fi

    if [[ "$run_CollectHsMetrics" == true ]]; then
        collect_hs_metrics "/in/${sorted_bam_name}" "/in/targets.picard" "/in/genome.fa" "/out/" "${MEM_IN_MB}"
    fi

    if [[ "$run_CollectTargetedPcrMetrics" == true ]]; then
        collect_targeted_pcr_metrics "/in/${sorted_bam_name}" "/in/genome.fa" "/in/targets.picard" "/out/" "${MEM_IN_MB}"
    fi

    if [[ "$run_CollectRnaSeqMetrics" == true ]]; then
        # Create refFlat file if not provided by user
        if [ -z "$ref_annot_refflat" ]; then
            echo "No refFlat file provided - creating GTF from refFlat file in CTAT bundle"
            REF_FLAT="${GENOME_BUILD_NAME}_ref_annot.refflat"
            OUTPUT_FILE="${OUTPUT_DIR}/${GENOME_BUILD_NAME}_ref_annot.refflat"
            make_refflat "/ctat_libdir/ref_annot.gtf" "/in/${REF_FLAT}" "${MEM_IN_MB}"
        else
            REF_FLAT="${ref_annot_refflat_name}"
        fi
        collect_rnaseq_metrics "/in/${sorted_bam_name}" "/in/${REF_FLAT}" "/out/" "${MEM_IN_MB}"
    fi

    if [[ "$run_CollectVariantCallingMetrics" == true ]]; then
        collect_variant_calling_metrics "/in/${vcf_name}" "/in/${dbsnp_vcf_name}" "/in/genome.dict" "/out/" "${MEM_IN_MB}"
    fi

    dx-upload-all-outputs --parallel
}
