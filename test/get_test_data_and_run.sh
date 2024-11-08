
if [ ! -f "test/input_sorted.bam" ]; then
  wget -O input_sorted.bam "https://nf-core-awsmegatests.s3-eu-west-1.amazonaws.com/viralrecon/results-3731dd3a32a67a2648ea22c2bd980c224abdaee2/platform_illumina/variants/bowtie2/SAMPLE_02.sorted.bam"
fi

./scripts/run_bam_to_bed_to_graph.sh \
  -i test/input_sorted.bam \
  -o test/results/ \
  -d false \
  -p "SAMPLE: TEST" \
  -c mapping_status_test

