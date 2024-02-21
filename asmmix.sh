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
chmod 777 $HOOM_DIR/asmmix.sh
echo "INFO  :   Run asmmix.sh from $HOOM_DIR ;"
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
    echo "      sh asmmix.sh  --tgsasm INPUT_TGS_ASM --slrasm INPUT_SLR_ASM --output OUT_PREFIX --quast quast/quast.py [options...]"
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

INPUT_TGS_ASM=$(realpath $INPUT_TGS_ASM)
INPUT_SLR_ASM=$(realpath $INPUT_SLR_ASM)

print_info "Parsing args end ."
#####################################################################
#
#   check env
#
############################################################3########

# check python scripts
CREATDUMMY=$HOOM_DIR"/AsmMix/create_dummy_scaff_from_quast.py"
EVASCAFF=$HOOM_DIR"/AsmMix/evaluate_scaff_struct.py"
REPSEQ=$HOOM_DIR"/AsmMix/replace_sequence_parallel.py"
chmod 777 $CREATDUMMY $EVASCAFF $REPSEQ

print_info "Checking basic args & env ..."

check_file_exe $CREATDUMMY
check_file_exe $EVASCAFF
check_file_exe $REPSEQ

# check input args.
check_arg_exist "tgsasm" $INPUT_TGS_ASM
check_arg_exist "slrasm"  $INPUT_SLR_ASM
check_arg_exist "quast"  $QUAST
check_arg_null "output" $OUT_PREFIX


print_info "Checking basic args & env end."
#####################################################################
#
#   step 1, obtain mapping relationship of assemblies using quast
#
############################################################3########

print_info "Step 1 , run QUAST to obtain the mapping relationship of assemblies. "
if [[ ! -e 'done_step1_tag' ]] ; then
    mkdir AsmMix_run
    ln -s $INPUT_TGS_ASM AsmMix_run/tgsasm.fasta
    ln -s $INPUT_SLR_ASM AsmMix_run/slrasm.fasta
    $QUAST AsmMix_run/slrasm.fasta -r AsmMix_run/tgsasm.fasta -o AsmMix_run/quast -s -m 10000 -x7000 -t $THREAD

    check_file AsmMix_run/quast/contigs_reports/contigs_report_slrasm.stdout
    date >>'done_step1_tag'
else 
    echo 'skip step1 since done_step1_tag exists'
fi
print_info "Step 1, done ."


#####################################################################
#
#   step 2, create dummy scaff from quast
#
############################################################3########

if [[ ! -e 'done_step2_tag' ]] ; then
    $CREATDUMMY AsmMix_run/quast/contigs_reports/contigs_report_slrasm.stdout AsmMix_run/mix.scaff AsmMix_run/mix.blat

    check_file AsmMix_run/mix.blat
    date >>'done_step2_tag'
else 
    echo 'skip step2 since done_step2_tag exists'
fi
print_info "Step 2, done ."


#####################################################################
#
#   step 3, split input scaffold
#
############################################################3########

if [[ ! -e 'done_step3_tag' ]] ; then
    $EVASCAFF AsmMix_run/mix.scaff AsmMix_run/mix.blat > AsmMix_run/mix.eval

    check_file AsmMix_run/mix.eval
    date >>'done_step3_tag'
else 
    echo 'skip step3 since done_step3_tag exists'
fi
print_info "Step 3, done ."


#####################################################################
#
#   step 4, split input scaffold
#
############################################################3########

if [[ ! -e 'done_step4_tag' ]] ; then
    $REPSEQ --max-length-to-replace $MAXLEN --num-thread $THREAD \
  AsmMix_run/slrasm.fasta AsmMix_run/tgsasm.fasta AsmMix_run/mix.eval ${OUT_PREFIX}.fasta

    check_file ${OUT_PREFIX}.fasta
    date >>'done_step4_tag'
else 
    echo 'skip step4 since done_step4_tag exists'
fi
print_info "Step 4, done ."

#####################################################################
#
# ALL DONE
#
#####################################################################
print_info "ALL DONE !!! "


  
