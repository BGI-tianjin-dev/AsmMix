open IN1, "$ARGV[0]";

my (@chr,@pos,@phase);
my $long_penalty=-5;
my $short_penalty=-1;

my ($total,$short_cnt,$long_cnt);
my (@phase_block_len);
while(<IN1>)
{
	chomp;
	my @array = split(/\t/);
	push @chr, $array[0];
	push @pos, $array[1];
	push @phase, $array[2];
}

for (my $start=0; $start <= $#chr; $start++) 
{
	my $stop;
	for($stop=$start; $stop<= $#chr; $stop++)
	{
		if($chr[$start] ne $chr[$stop]){last;}
	}
	my $score0=($phase[$start]==0)?0:$short_penalty;
	my $score1=($phase[$start]==0)?$short_penalty:0;
	my @trace_back0;
	my @trace_back1;
	for(my $loop=$start+1;$loop<$stop;$loop++)
	{
		my $score00=$score0+(($phase[$loop]==0)?0:$short_penalty);
		my $score01=$score0+(($phase[$loop]==1)?0:$short_penalty)+$long_penalty;
		my $score10=$score1+(($phase[$loop]==0)?0:$short_penalty)+$long_penalty;
		my $score11=$score1+(($phase[$loop]==1)?0:$short_penalty);

		if($score00>$score10){$score0=$score00; push @trace_back0,0;}
		else {$score0=$score10; push @trace_back0,1;}

		if($score01>$score11){$score1=$score01;push @trace_back1,0;}
		else {$score1=$score11;push @trace_back1,1;}
	}

	my @path;
	if($score0>$score1){unshift(@path,0);}
	else{unshift(@path,1);}

	for(my $loop=$#trace_back0;$loop>=0;$loop--)
	{
		if($path[0]==0){unshift(@path,$trace_back0[$loop]);}
		elsif($path[0]==1){unshift(@path,$trace_back1[$loop]);}
	}
 	#print "$start\t$stop\t$#path\t$#trace_back0\t$#trace_back1\n";
 	$phase_block_start=0;
	for(my $loop=$start;$loop<$stop;$loop++)
	{
		if($phase[$loop] != $path[$loop-$start]){$short_cnt++;}
		if($loop>$start && $path[$loop-$start]!=$path[$loop-$start-1]) 
		{
			$long_cnt++;
			push @phase_block_len, $pos[$loop-1]-$phase_block_start;
			$phase_block_start=$pos[$loop];
		}
		$total++;
	}

			push @phase_block_len, $pos[$stop-1]-$phase_block_start;
	# for(my $loop=$start;$loop<$stop;$loop++)
	# {
	# 	print "$phase[$loop]\t$path[$loop-$start]\n";
	# }
	# print "\n";
	# print "$start -> $stop\n";
	$start=$stop-1;
	#last;
}
print "$short_cnt/$total\t$long_cnt/$total\t$total\n";


my @sorted = sort { $a <=> $b } @phase_block_len;
my $sum=0;
foreach ( @sorted) {
	$sum+=$_;
}

my $sum2=0;
foreach ( @sorted) {
	$sum2+=$_;
	if(2*$sum2>$sum){print "$_\n"; last;}
}
close IN1;
