#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

create_interval_file() {
	echo "create_interval_file"
	# Converts a BED file to a Picard Interval List
	# List https://gatk.broadinstitute.org/hc/en-us/articles/360037593251-BedToIntervalList-Picard-
	$java -jar /picard.jar BedToIntervalList \
	I="$bedfile_path" \
	O=targets.picard \
	SD="$sorted_bam_path"
}

collect_targeted_pcr_metrics() {
	echo "collect_targeted_pcr_metrics"
	# Call Picard CollectMultipleMetrics. Requires the co-ordinate sorted BAM file given to the app
	# as input. The file is referenced in this command with the option 'I=<input_file>'. Here, the
	# downloaded BAM file path is accessed using the DNAnexus helper variable $sorted_bam_path.
	# All outputs are saved to $output_dir (defined in main()) for upload to DNAnexus.
	$java -jar /picard.jar CollectTargetedPcrMetrics  I="$sorted_bam_path" R=genome.fa \
	O="$output_dir/$sorted_bam_prefix.targetPCRmetrics.txt" AI=targets.picard TI=targets.picard \
	PER_TARGET_COVERAGE="$output_dir/$sorted_bam_prefix.perTargetCov.txt"
}

collect_multiple_metrics() {
	echo "collect_multiple_metrics"
	# Call Picard CollectMultipleMetrics. Requires the co-ordinate sorted BAM file given to the app
	# as input. The file is referenced in this command with the option 'I=<input_file>'. Here, the
	# downloaded BAM file path is accessed using the DNAnexus helper variable $sorted_bam_path.
	# All outputs are saved to $output_dir (defined in main()) for upload to DNAnexus.
	# Note that not all outputs are relevent for all types of sequencing
	# e.g. some aren't applicable for amplicon NGC
	# Note that CollectSequencingArtifactMetrics errors out with TSO500 stitched BAMs due to 
	# "Record contains library that is missing from header" and so not used
	$java -jar /picard.jar CollectMultipleMetrics I="$sorted_bam_path" R=genome.fa \
	PROGRAM=null \
	PROGRAM=CollectAlignmentSummaryMetrics \
	PROGRAM=CollectInsertSizeMetrics \
	PROGRAM=QualityScoreDistribution \
	PROGRAM=MeanQualityByCycle \
	PROGRAM=CollectBaseDistributionByCycle \
	PROGRAM=CollectGcBiasMetrics \
	PROGRAM=CollectQualityYieldMetrics \
	O="$output_dir/$sorted_bam_prefix"
	# PROGRAM=CollectSequencingArtifactMetrics \
}

calculate_hs_metrics() {
	echo "collect_hs_metrics"
	# Call Picard CollectHsMetrics. Requires the co-ordinate sorted BAM file given to the app as
	# input (I=). Outputs the hsmetrics.tsv and pertarget_coverage.tsv files to $output_dir
	# (defined in main()) for upload to DNAnexus. Note that coverage cap is set to 100000 (default=200).
	$java -jar /picard.jar CollectHsMetrics BI=targets.picard TI=targets.picard I="$sorted_bam_path" \
	O="$output_dir/${sorted_bam_prefix}.hsmetrics.tsv" R=genome.fa \
	PER_TARGET_COVERAGE="$output_dir/${sorted_bam_prefix}.pertarget_coverage.tsv" \
	COVERAGE_CAP=100000
}

main() {

##### SETUP #####

# Download input files from inputSpec to ~/in/. Allows the use of DNA Nexus bash helper variables.
dx-download-all-inputs

# Calculate 90% of memory size for java
mem_in_mb=`head -n1 /proc/meminfo | awk '{print int($2*0.9/1024)}'`
# Set java command with the calculated maximum memory usage
java="java -Xmx${mem_in_mb}m"

# Unpack the reference genome for Picard. Produces genome.fa, genome.fa.fai, and genome.dict files.
tar zxf $fasta_index_path

# Create directory for Picard stats files to be uploaded from the worker
output_dir=$HOME/out/eggd_picard_stats/QC
mkdir -p $output_dir

##### MAIN #####

# Create the interval file
create_interval_file

# if it's a capture panel call the relevant modules
if [[ "$enrichment_method" == "Hybridisation" ]]; then
# Call Picard CollectMultipleMetrics
collect_multiple_metrics
# Call Picard CalculateHSMetrics
calculate_hs_metrics

# if it's a amplicon panel call the relevant modules
elif [[ "$enrichment_method" == "Amplicon" ]]; then
collect_targeted_pcr_metrics
else
echo "unknown capture type"
fi

##### CLEAN UP #####

# Upload all results files and directories in $HOME/out/eggd_picard_stats/
dx-upload-all-outputs --parallel
}