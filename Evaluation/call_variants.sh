HAP1=$1
HAP2=$2
#here to config, reference genome and other tools
REF=hs37d5.fa 
Bcftools=bcftools
Tabix=tabix
Bgzip=bgzip
Minimap=minimap2
Paftools=paftools.js #see https://github.com/lh3/minimap2/tree/master/misc on how to config paftools

$Minimap2 -c --cs -x asm5 $REF $HAP1 -t 30 | sort -k6,6 -k8,8n > hap1.paf
$Minimap2 -c --cs -x asm5 $REF $HAP2 -t 30 | sort -k6,6 -k8,8n > hap2.paf
$Paftools call -f $REF hap1.paf > hap1.vcf
$Paftools call -f $REF hap2.paf > hap2.vcf
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
perl vcf_reformat.pl merge.vcf > hap.vcf
$Bgzip hap.vcf
$Tabix hap.vcf.gz