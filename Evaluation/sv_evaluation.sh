truvari bench -c $1 \
-b benchmark_data/HG002_SVs_Tier1_v0.6.vcf.gz\
 -o SV_eval --giabreport \
 -r 1000 -f benchmark_data/hs37d5.fa \
 --passonly --includebed benchmark_data/HG002_SVs_Tier1_v0.6.bed