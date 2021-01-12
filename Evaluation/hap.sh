REF=~/wupei1/ref_genome/hs37d5.fa
HAP1=$1
HAP2=$2
#Bcftools=/ldfssz1/MGI_BIT/RUO/qiuweijing/Project/01.stLFR/v1/Publish_v1/latest_wat/lib/bcftools
Bcftools=/home/wupei1/wupei1/software/bcftools/bcftools-1.11/bcftools
Tabix=/ldfssz1/MGI_BIT/RUO/qiuweijing/Project/01.stLFR/v1/Publish_v1/latest_wat/lib/tabix
Bgzip=/ldfssz1/MGI_BIT/RUO/qiuweijing/Project/01.stLFR/v1/Publish_v1/latest_wat/lib/bgzip

#~/wupei1/software/minimap2 -c --cs -x asm5 $REF $HAP1 -t 30 | sort -k6,6 -k8,8n > hap1.paf
#~/wupei1/software/minimap2 -c --cs -x asm5 $REF $HAP2 -t 30 | sort -k6,6 -k8,8n > hap2.paf
#~/wupei1/software/minimap2-master/misc/paftools.js call -f $REF hap1.paf > hap1.vcf
#~/wupei1/software/minimap2-master/misc/paftools.js call -f $REF hap2.paf > hap2.vcf
$Bcftools sort hap1.vcf > hap1_sort.vcf
$Bcftools sort hap2.vcf > hap2_sort.vcf
$Bcftools norm -f $REF hap1_sort.vcf > hap1_norm.vcf
$Bcftools norm -f $REF hap2_sort.vcf > hap2_norm.vcf
$Bcftools sort hap1_norm.vcf > hap1_sort.vcf
$Bcftools sort hap2_norm.vcf > hap2_sort.vcf
$Bgzip hap1.vcf
$Bgzip hap2.vcf
$Tabix hap1.vcf.gz
$Tabix hap2.vcf.gz
$Bcftools merge hap1.vcf.gz hap2.vcf.gz -o merge.vcf --force-samples
perl ~/wupei1/ASM_mix/replace_raw/SV_eval/vcf_reformat.pl merge.vcf > hap.vcf
