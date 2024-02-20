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
    echo "      required:"
    echo "          --tgsasm     <tgs_assmbly_file>      input TGS assembly."
    echo "          --slrasm     <slr_assembly_file>     input SLR assembly."
    echo "          --output     <output_prefix>         output prefix."
    echo "          --quast      <quast_prefix>          installed quast path."
    echo "      optional:"
    echo "          --max-length-to-replace  <maxlen args>      for example, --max-length-to-replace 50 "
    echo "          --num-thread             <thread args>      for example, --num-thread 32 "

}
function check_arg_null() {
    if [[ -z $2 ]] ; then 
        print_fatal "Arg $1 is not assigned !!!";
    fi
}

function check_file()
{
    if [[ ! -e $1 ]] ; then
        print_fatal "File $1 does not exist !!!"
    fi
}

function check_file_exe()
{
    check_file $1
    if [[ ! -x $1 ]] ; then
        print_fatal "File $1 is not executable !!!"
    fi
}
function check_arg_exist()
{
    check_arg_null $1 $2
    check_file $2
}

function check_arg_exe()
{
    check_arg_null $1 $2
    check_file_exe $2
}

#####################################################################
#
#   check args
#
############################################################3########
INPUT_TGS_ASM=""
INPUT_SLR_ASM=""
OUT_PREFIX="asmmix"
QUAST=""

MAXLEN=50
THREAD=32


print_info "Parsing args starting ..."
if [[ $# -lt 1 ]] ; then 
    print_help
    exit 0 ;
fi
ARGS=`getopt -o h  --long tgsasm:,slrasm:,output:,quast:,max-length-to-replace:,num-thread  -- "$@"`
eval set -- "$ARGS"
while true; do
    case "$1" in
        -h|--help)
            shift;
            print_help;
            exit 0;
        ;;
        --tgsasm)
            shift;
            INPUT_TGS_ASM=$1
            shift;
            echo  "             --tgsasm $INPUT_TGS_ASM"
        ;;
        --slrasm)
            shift;
            INPUT_SLR_ASM=$1
            shift;
            echo  "             --slrasm $INPUT_SLR_ASM"
        ;;
        --output)
            shift;
            OUT_PREFIX=$1
            shift;
            echo  "             --output $OUT_PREFIX"
        ;;
        --quast)
            shift;
            QUAST=$1
            shift;
            echo  "             --quast $QUAST"
        ;;
        --max-length-to-replace)
            shift;
            MAXLEN=$1
            shift;
            echo  "             --max-length-to-replace $MAXLEN"
        ;;
        --num-thread)
            shift;
            THREAD=$1
            shift;
            echo  "             --num-thread $THREAD"
        ;;
        --)
            shift;
            break;
        ;;
    esac
done

print_info "Parsing args end ."
#####################################################################
#
#   check env
#
############################################################3########

# check binary
Candidate=$HOOM_DIR"/AsmMix/tgsgapcandidate"
GapCloser=$HOOM_DIR"/AsmMix/tgsgapcloser"
SeqGen=$HOOM_DIR"/AsmMix/tgsseqgen"
SeqSplit=$HOOM_DIR"/AsmMix/tgsseqsplit"

print_info "Checking basic args & env ..."

check_file_exe $Candidate
check_file_exe $GapCloser
check_file_exe $SeqGen
check_file_exe $SeqSplit
#check_file_exe $MiniMap2
# check input args.
check_arg_exist "input_scaff" $INPUT_SCAFF
check_arg_exist "reads"  $TGS_READS
check_arg_null "output" $OUT_PREFIX

if [[ $MIN_NREAD -lt 1 ]] ; then 
    echo "Error : min_nread < 1. exit ..."
    exit 1
fi

if [[ $NE == "no" ]] ; then
    if [[ $NGS_READS != "" ]] ; then 
        check_arg_exist "ngs" $NGS_READS
        check_arg_exist "pilon" $PILON
        check_arg_exe "samtools" $SAMTOOL
        check_arg_exe "java" $JAVA
        print_info_line "Doing error correction by pilon with ngs reads. "
        USE_RACON="no"
    else
        check_arg_exe "racon" $RACON
        print_info_line "Doing error correction by racon with tgs reads."
        USE_RACON="yes"
    fi
else 
    print_info_line "No error correction by --ne option"
fi
# pacbio special default values.
if [[ $TGS_TYPE == "pb" ]] ; then 
    if [[ $USE_NEW_MINMAP2_ARG != "yes" ]] ; then 
        MINIMAP2_PARAM=" -x ava-pb "
    fi
    if [[ $MIN_IDY_USER == "no" ]] ; then
        MIN_IDY="0.2"
    fi
    if [[ $MIN_MATCH_USER == "no" ]] ; then
        MIN_MATCH="200"
    fi
fi
print_info_line "TGS reads type is $TGS_TYPE . MINIMAP2_PARAM is $MINIMAP2_PARAM  MIN_IDY is $MIN_IDY . MIN_MATCH is $MIN_MATCH ."

print_info "Checking basic args & env end."
#####################################################################
#
#   step 1 , split input scaffold
#
############################################################3########


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





  
