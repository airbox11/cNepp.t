
rm -rf 1.vcf
cat /omics/groups/OE0422/internal/yanhong/all_in_one_pipeline/vcf_header.txt >> 1.vcf
cat *org | grep -v '^#' >> 1.vcf