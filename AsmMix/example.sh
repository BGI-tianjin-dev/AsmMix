mkdir AsmMix_run
quast.py $1 -r $2 -o AsmMix_run/quast -s -m 10000 -x7000 -t 32
fa_name="${1%.*}"
create_dummy_scaff_from_quast.py AsmMix_run/quast/${fa_name}.stdout AsmMix_run/mix.scaff AsmMix_run/mix.blat
evaluate_scaff_struct.py  AsmMix_run/mix.scaff AsmMix_run/mix.blat > AsmMix_run/mix.eval
replace_sequence_parallel.py --max-length-to-replace 50 --num-thread 32 \
$1 $2 AsmMix_run/mix.eval mix.fa