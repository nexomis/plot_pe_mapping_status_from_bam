
if [ ! -f "test/input_sorted.bam" ]; then
  wget -O test/input_sorted.bam "https://nf-core-awsmegatests.s3-eu-west-1.amazonaws.com/viralrecon/results-3731dd3a32a67a2648ea22c2bd980c224abdaee2/platform_illumina/variants/bowtie2/SAMPLE_02.sorted.bam"
  wget -O test/input.gff.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz
  gunzip test/input.gff.gz
  sed -i "s/NC_045512.2/MN908947.3/g" test/input.gff
fi

./scripts/run_bam_to_bed_to_graph.sh \
  -scripts_dir_path scripts/ \
  -input_bam test/input_sorted.bam \
  -output_dir test/results/ \
  -delete_bam_files false \
  -plot_title "SAMPLE: TEST" \
  -output_bn_depth_pdf mapping_status_test \
  -annot_gff_file test/input.gff \
  -annot_feat_id_regex ".*gene=([^; ]+).*" \
  -annot_feat_id_catch "\\1" \
  -annot_interest_type "gene,CDS" \
  -max_pb_by_A4_width 10000 \
  -plotly_out_dir test/results/plotly/
