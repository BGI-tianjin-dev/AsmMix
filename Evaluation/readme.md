# dependencies

minimap2

rtgtools

truvari

bcftools

# usage

run `bash build_data.sh` to download all nessesary data 

to generate diploid vcf, run `bash call_variants.sh hap1.fa hap2.fa`, produce a vcf file hap.vcf.gz

evaluation of SNVs: `bash snp_evaluation.sh hap.vcf.gz`

evaluation of SVs: `bash sv_evaluation.sh hap.vcf.gz`

evaluation of phasing: `bash phase_evaluation.sh hap.vcf.gz`
