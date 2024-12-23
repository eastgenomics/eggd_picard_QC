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
	    -I "${BEDFILE_PATH}" \
	    -O "${OUTPUT_TARGETS}" \
	    -SD "${SORTED_BAM}"
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
        -I "${SORTED_BAM}" \
        -R "${REF_GENOME}" \
	    -O "${OUTPUT_DIR}/${SORTED_BAM_PREFIX}.targetPCRmetrics.txt" \
        -AI "${TARGETS_FILE}" \
        -TI "${TARGETS_FILE}" \
        --PER_TARGET_COVERAGE "${OUTPUT_DIR}/${SORTED_BAM_PREFIX}.perTargetCov.txt"
}

# collect_multiple_metrics() - Collected multiple metrics. Note that not all outputs are relevant for all type of sequencing.
# See https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectMultipleMetrics-Picard for details
#
# Arguments:
#   $1 - Path to sorted BAM file
#   $2 - Path to reference genome
#   $3 - Path to output directory
#   $4 - Java maximum heap size
collect_multiple_metrics() {
    local SORTED_BAM=$1
    local REF_GENOME=$2
    local OUTPUT_DIR=$3
    local MAXHEAP=$4

    local SORTED_BAM_PREFIX
    SORTED_BAM_PREFIX=$(basename "${SORTED_BAM}" .bam)

    docker exec picard_image java -Xmx"${MAXHEAP}" -jar /usr/picard/picard.jar CollectMultipleMetrics \
        -I "${SORTED_BAM}" \
        -R "${REF_GENOME}" \
	    --PROGRAM null \
	    --PROGRAM CollectAlignmentSummaryMetrics \
	    --PROGRAM CollectInsertSizeMetrics \
	    --PROGRAM QualityScoreDistribution \
	    --PROGRAM MeanQualityByCycle \
	    --PROGRAM CollectBaseDistributionByCycle \
	    --PROGRAM CollectGcBiasMetrics \
	    --PROGRAM CollectQualityYieldMetrics \
	    -O "${OUTPUT_DIR}/${SORTED_BAM_PREFIX}"
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
        --BI "${TARGETS_FILE}" \
        --TI "${TARGETS_FILE}" \
        --I "${SORTED_BAM}" \
        --O "${OUTPUT_DIR}/${SORTED_BAM_PREFIX}.hsmetrics.tsv" \
        --R "${REF_GENOME}" \
        --PER_TARGET_COVERAGE "${OUTPUT_DIR}/${SORTED_BAM_PREFIX}.pertarget_coverage.tsv"\
        --COVERAGE_CAP 100000
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
        -I "${SORTED_BAM}" \
        -O "${OUTPUT_DIR}/${SORTED_BAM_PREFIX}.RNAmetrics.tsv" \
        --REF_FLAT "${REF_FLAT}" \
        -STRAND SECOND_READ_TRANSCRIPTION_STRAND
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
        --DBSNP "${DBSNP_VCF}" \
        --INPUT "${VCF}" \
        --OUTPUT "${OUTPUT_DIR}/${VCF_PREFIX}.variantcallingmetrics" \
        --SEQUENCE_DICTIONARY "${SEQ_DICT}" \
        --GVCF_INPUT true
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

    if [[ "$run_CollectRnaSeqMetrics" == "true" && ( -z "$sorted_bam" || -z "$ref_annot_refflat" ) ]]; then
        err "run_CollectRnaSeqMetrics was requested, but one or more of sorted_bam or ref_annot_refflat are missing. Exiting..."
        exit 1
    fi

    if [[ "$run_CollectVariantCallingMetrics" == "true" && \
        ( -z "$vcf" || -z "$vcf_index" || -z "$dbsnp_vcf" ) ]]; then
        err "run_CollectVariantCallingMetrics was requested, but one or more of vcf, vcf_index, or dbsnp_vcf are missing. Exiting..."
        exit 1
    fi

    ## Setup 
    ### puts inputs in /home/dnanexus/in/
    dx-download-all-inputs

    ### move all inputs to flat /home/dnanexus/input/ directory
    mkdir -p ~/input/
    find ~/in -type f -name "*" -print0 | xargs -0 -I {} mv {} ~/input/

    # Calculate 90% of memory size for java
    MEM=$(head -n1 /proc/meminfo | awk '{print int($2*0.9)}')
    MEM_IN_MB="$(("${MEM}"/1024))m"

    tar zxvf ~/input/${fasta_index_name} -C ~/input/ 
    OUTPUT_DIR="${HOME}/out/eggd_picard_stats/QC"
    mkdir -p "$OUTPUT_DIR"

    DOCKER_IMAGENAME=$(find /image -name "*.tar.gz")
    sudo docker load -i "${DOCKER_IMAGENAME}"
    DOCKER_IMAGE=$(docker image ls -q)
    docker run \
        --name picard_image \
        --mount type=bind,source=/home/dnanexus/input/,target=/input \
        --mount type=bind,source="${OUTPUT_DIR}",target=/out \
        --entrypoint /bin/bash \
        -itd "${DOCKER_IMAGE}"

    # Create the interval file if required
    if [[ "$run_CollectMultipleMetrics" == true ]] || \
        [[ "$run_CollectHsMetrics" == true ]] || \
        [[ "$run_CollectTargetedPcrMetrics" == true ]]; then
        echo "Generating interval file"
        create_interval_file "/input/${bedfile_name}" "/input/${sorted_bam_name}" "/input/targets.picard" "${MEM_IN_MB}"
    fi

    ## Run picard commands
    if [[ "$run_CollectMultipleMetrics" == true ]]; then
        collect_multiple_metrics "/input/${sorted_bam_name}" "/input/genome.fa" "/out/" "${MEM_IN_MB}"
    fi

    if [[ "$run_CollectHsMetrics" == true ]]; then
        collect_hs_metrics "/input/${sorted_bam_name}" "/input/targets.picard" "/input/genome.fa" "/out/" "${MEM_IN_MB}"
    fi

    if [[ "$run_CollectTargetedPcrMetrics" == true ]]; then
        collect_targeted_pcr_metrics "/input/${sorted_bam_name}" "/input/genome.fa" "/input/targets.picard" "/out/" "${MEM_IN_MB}"
    fi

    if [[ "$run_CollectRnaSeqMetrics" == true ]]; then
        collect_rnaseq_metrics "/input/${sorted_bam_name}" "/input/${ref_annot_refflat_name}" "/out/" "${MEM_IN_MB}"
    fi

    if [[ "$run_CollectVariantCallingMetrics" == true ]]; then
        collect_variant_calling_metrics "/input/${vcf_name}" "/input/${dbsnp_vcf_name}" "/input/genome.dict" "/out/" "${MEM_IN_MB}"
    fi

    dx-upload-all-outputs --parallel
}
