#mkdir -p /nascngb/gccnt/wangou/shichang/LFR_WGS_data/HG002_1.5ng/output/09.stat/variant/HG002_1.5ng/evaluation
for type in snp indel
  do
    /ldfssz1/MGI_BIT/RUO/qiuweijing/Project/01.stLFR/v1/Publish_v1/latest_wat/lib/bcftools view -O z --type ${type}s $1 > $type.vcf.gz
    /ldfssz1/MGI_BIT/RUO/qiuweijing/Project/01.stLFR/v1/Publish_v1/latest_wat/lib/tabix -p vcf -f $type.vcf.gz
    #rm -fr /nascngb/gccnt/wangou/shichang/LFR_WGS_data/HG002_1.5ng/output/09.stat/variant/HG002_1.5ng/evaluation/$type

    /ldfssz1/MGI_BIT/RUO/qiuweijing/Project/01.stLFR/v1/Publish_v1/latest_wat/tools/rtg-tools-3.7.1/rtg RTG_MEM=8G vcfeval \
      -c $type.vcf.gz \
      -b /nascngb/gccnt/wangou/shichang/LFR_WGS_data/HG002_1.5ng/output/09.stat/variant/HG002_1.5ng/evaluation/$type.vcf.gz \
      -e /ldfssz1/MGI_BIT/RUO/lizhanqing/01.Pipeline/01.WGS/Database/AshkenazimTrio/HG002_NA24385_son/latest/GRCh37/HG002_NA24385_highconf.bed \
      -t /ldfssz1/MGI_BIT/RUO/qiuweijing/Project/01.stLFR/v1/Publish_v1/latest_wat/database/hs37d5/hs37d5.fa.SDF \
      -o ${1}_${type}
      #-o ${1}_${type}_squash --squash-ploidy
      #2> /nascngb/gccnt/wangou/shichang/LFR_WGS_data/HG002_1.5ng/output/09.stat/variant/HG002_1.5ng/evaluation/$type.log
    #for type2 in tp fp fn
    #do
     # gzip -cd /nascngb/gccnt/wangou/shichang/LFR_WGS_data/HG002_1.5ng/output/09.stat/variant/HG002_1.5ng/evaluation/$type/$type2.vcf.gz | grep -v ^# | awk 'BEGIN{OFS="\t"}{a=$2-1;print $1,a,$2}' > /nascngb/gccnt/wangou/shichang/LFR_WGS_data/HG002_1.5ng/output/09.stat/variant/HG002_1.5ng/evaluation/$type/$type2.bed
    #done
done
