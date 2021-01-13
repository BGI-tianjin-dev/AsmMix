for type in snp indel
  do
    bcftools view -O z --type ${type}s $1 > $type.vcf.gz
    tabix -p vcf -f $type.vcf.gz
    
    rtg RTG_MEM=8G vcfeval \
      -c $type.vcf.gz \
      -b benchmark_data/$type.vcf.gz \
      -e benchmark_data/HG002_NA24385_highconf.bed \
      -t benchmark_data/hs37d5.sdf \
      -o ${1}_${type}
    rtg RTG_MEM=8G vcfeval \
      -c $type.vcf.gz \
      -b benchmark_data/$type.vcf.gz \
      -e benchmark_data/HG002_NA24385_highconf.bed \
      -t benchmark_data/hs37d5.sdf \
      -o ${1}_${type}_squash --squash-ploidy
done
