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

create_interval_file() {
	# Converts a BED file to a Picard Interval List
	# See https://gatk.broadinstitute.org/hc/en-us/articles/360037593251-BedToIntervalList-Picard-
	# Note that SD can take any of the follow. Here we are using the BAM.
	# - A file with .dict extension generated using Picard's CreateSequenceDictionaryTool
	# - A reference.fa or reference.fasta file with a reference.dict in the same directory
	# - Another IntervalList with @SQ lines in the header from which to generate a dictionary
	# - A VCF that contains #contig lines from which to generate a sequence dictionary
	# - A SAM or BAM file with @SQ lines in the header from which to generate a dictionary
	$java -jar /picard.jar BedToIntervalList \
	    -I="$bedfile_path" \
	    -O=targets.picard \
	    -SD="$sorted_bam_path"
}

collect_targeted_pcr_metrics() {
	# Call Picard CollectMultipleMetrics. Requires the co-ordinate sorted BAM file given to the app
	# as input. The file is referenced in this command with the option 'I=<input_file>'. Here, the
	# downloaded BAM file path is accessed using the DNAnexus helper variable $sorted_bam_path.
	# All outputs are saved to $output_dir (defined in main()) for upload to DNAnexus.
	$java -jar /picard.jar CollectTargetedPcrMetrics  \
        -I="$sorted_bam_path" \
        -R=genome.fa \
	    -O="$output_dir/$sorted_bam_prefix.targetPCRmetrics.txt" \
        -AI=targets.picard \
        -TI=targets.picard \
        --PER_TARGET_COVERAGE="$output_dir/$sorted_bam_prefix.perTargetCov.txt"
}

collect_multiple_metrics() {
	# Call Picard CollectMultipleMetrics. Requires the co-ordinate sorted BAM file given to the app
	# as input. The file is referenced in this command with the option 'I=<input_file>'. Here, the
	# downloaded BAM file path is accessed using the DNAnexus helper variable $sorted_bam_path.
	# All outputs are saved to $output_dir (defined in main()) for upload to DNAnexus.
	# Note that not all outputs are relevent for all types of sequencing
	# e.g. some aren't applicable for amplicon NGC
    # TODO: Investigate CollectSequencingArtifactMetrics
	# Note that CollectSequencingArtifactMetrics errors out with TSO500 BAMs due to 
	# "Record contains library that is missing from header" and so not used (fix unclear)
	$java -jar /picard.jar CollectMultipleMetrics \
        -I="$sorted_bam_path" \
        -R=genome.fa \
	    --PROGRAM=null \
	    --PROGRAM=CollectAlignmentSummaryMetrics \
	    --PROGRAM=CollectInsertSizeMetrics \
	    --PROGRAM=QualityScoreDistribution \
	    --PROGRAM=MeanQualityByCycle \
	    --PROGRAM=CollectBaseDistributionByCycle \
	    --PROGRAM=CollectGcBiasMetrics \
	    --PROGRAM=CollectQualityYieldMetrics \
	    -O="$output_dir/$sorted_bam_prefix"
}

collect_hs_metrics() {
	echo "collect_hs_metrics"
	# Call Picard CollectHsMetrics. Requires the co-ordinate sorted BAM file given to the app as
	# input (I=). Outputs the hsmetrics.tsv and pertarget_coverage.tsv files to $output_dir
	# (defined in main()) for upload to DNAnexus. Note that coverage cap is set to 100000 (default=200).
	$java -jar /picard.jar CollectHsMetrics \
        --BI=targets.picard \
        --TI=targets.picard \
        --I="$sorted_bam_path" \
        --O="$output_dir/${sorted_bam_prefix}.hsmetrics.tsv" \
        --R=genome.fa \
        --PER_TARGET_COVERAGE="$output_dir/${sorted_bam_prefix}.pertarget_coverage.tsv"\
        --COVERAGE_CAP=100000
}

collect_rnaseq_metrics() {
	echo "collect_rnaseq_metrics"
	# Call Picard CollectRnaSeqMetrics. akes a SAM/BAM file containing
	# the aligned reads from an RNAseq experiment and produces metrics
	# describing the distribution of the bases within the transcripts
	$java -jar /picard.jar CollectRnaSeqMetrics \
        -I="$sorted_bam_path" \
        -O="$output_dir/${sorted_bam_prefix}.RNAmetrics.tsv" \
        --REF_FLAT="$ref_flat" \
        -STRAND=SECOND_READ_TRANSCRIPTION_STRAND
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

    ## Setup 
    dx-download-all-inputs

    # Calculate 90% of memory size for java
    mem_in_mb=$(head -n1 /proc/meminfo | awk '{print int($2*0.9/1024)}')

    # Set java command with the calculated maximum memory usage
    java="java -Xmx${mem_in_mb}m"

    tar zxvf "$fasta_index_path"
    output_dir=$HOME/out/eggd_picard_stats/QC
    mkdir -p "$output_dir"

    # Create the interval file if required
    if [[ "$run_CollectMultipleMetrics" == true ]] || \
        [[ "$run_CollectHsMetrics" == true ]] || \
        [[ "$run_CollectTargetedPcrMetrics" == true ]]; then
        echo "Generating interval file"
        create_interval_file
    fi

    ## Run picard commands
    if [[ "$run_CollectMultipleMetrics" == true ]]; then
        collect_multiple_metrics
    fi

    if [[ "$run_CollectHsMetrics" == true ]]; then
        collect_hs_metrics
    fi

    if [[ "$run_CollectTargetedPcrMetrics" == true ]]; then
        collect_targeted_pcr_metrics
    fi

    if [[ "$run_CollectRnaSeqMetrics" == true ]]; then
        # Create refFlat file if not provided by user
        if [ -z "$ref_annot_refflat" ]; then
            echo "No refFlat file provided - creating GTF from refFlat file in CTAT bundle"
            lib_dir=$(echo $fasta_index_name | cut -d "." -f 1,2)
            ref_annot_gtf="/home/dnanexus/${lib_dir}/ctat_genome_lib_build_dir/ref_annot.gtf"
            $java -jar /GtftoRefflat-assembly-0.1.jar \
                -g "$ref_annot_gtf" \
                -r "${lib_dir}_ref_annot.refflat"
            ref_flat="${lib_dir}_ref_annot.refflat"
        else
            ref_flat="$ref_annot_refflat_path"
        fi
        collect_rnaseq_metrics
    fi

    dx-upload-all-outputs --parallel
}
