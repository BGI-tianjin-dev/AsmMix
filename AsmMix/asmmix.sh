#!/bin/bash
#####################################################################
#
#   brief :  main pipeline script of asmmix.
#   usage :  use -h/--help for details.
#   author :  xumengyang@genomics.cn && liuchao3@genomics.cn
#
############################################################3########

#####################################################################
#
#   version 
#
############################################################3########
VERSION="1.0.0"
RELEASE_TIME="2024-02-20"
HOOM_DIR=`dirname $0`
echo "INFO  :   Run asmmix. from $HOOM_DIR ;"
echo "          Version : $VERSION ;"
echo "          Release time : $RELEASE_TIME ."
echo ""
#####################################################################
#
#   function
#
############################################################3########
function print_info_line()
{
    echo "              -   $1"
}
function print_info()
{
    echo "";
    echo "INFO  :   $1"
}
function print_fatal()
{
    echo "";
    echo "FATAL :   $1 ! exit ...";
    echo " try -h/--help for usage."
    exit 1 ;
}
function print_help()
{
    echo "Usage:"
    echo "      asmmix  --tgsasm SCAFF_FILE --slrasm TGS_READS_FILE --output OUT_PREFIX --quast quast/quast.py [options...]"
    echo "      required :"
    echo "          --tgsasm     <tgs_assmbly_file>      input TGS assembly."
    echo "          --slrasm     <slr_assembly_file>     input SLR assembly."
    echo "          --output     <output_prefix>         output prefix."
    echo "          --quast      <quast_prefix>          installed quast path."
    echo "      optional:"
    echo "          --max-length-to-replace  <maxlen args>      for example, --max-length-to-replace 50 "
    echo "          --num-thread             <thread args>      for example, --num-thread 32 "

}



#command: example.sh stLFR_asm.fa SMS.fa
#output will be mix.fa
#unzip fasta files before using, file prefix should not contain "."

mkdir AsmMix_run
quast.py $1 -r $2 -o AsmMix_run/quast -s -m 10000 -x7000 -t 32
fa_name="${1%.*}"

create_dummy_scaff_from_quast.py AsmMix_run/quast/contigs_report_${fa_name}.stdout AsmMix_run/mix.scaff AsmMix_run/mix.blat
evaluate_scaff_struct.py AsmMix_run/mix.scaff AsmMix_run/mix.blat > AsmMix_run/mix.eval
replace_sequence_parallel.py --max-length-to-replace 50 --num-thread 32 \
  $1 $2 AsmMix_run/mix.eval mix.fa





  
