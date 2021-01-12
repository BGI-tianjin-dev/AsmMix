my $vcf_file1 =  $ARGV[0];
my $vcf_file2 =  $ARGV[1];
my @array;
open IN1, "zcat $vcf_file1| sort -k1,1n -k2,2n |";
open IN2, "zcat $vcf_file2| sort -k1,1n -k2,2n |";

my (@chr1, @pos1,@ref1,@alt1,@phase1,@info1);
my (@chr2, @pos2,@ref2,@alt2,@phase2);
my $chr_num=0;
my %chr_hash;

while(<IN1>)
{
	chomp;
 	$line=$_;
    if($line=~/^#/){ next;}
    else
	{
		@array = split /\t/, $line;
		if($array[0]!~/\d+/){next;}
		if($array[9]!~/0\|1/ && $array[9]!~/1\|0/ ) {next;}
		push @chr1, int($array[0]);
		push @pos1, int($array[1]);
		push @ref1, $array[3];
		push @alt1, $array[4];
		push @info1, $array[7];
		if($array[9]!~/0\|1/) 
		{
			push @phase1, 0;
		}
		else{push @phase1, 1;}
	}
}

while(<IN2>)
{
chomp;
 	$line=$_;
    if($line=~/^#/){ next;}
    else
	{
		@array = split /\t/, $line;
		if($array[0]!~/\d+/){next;}
		if($array[9]!~/0\|1/ && $array[9]!~/1\|0/ ) {next;}
		push @chr2, int($array[0]);
		push @pos2, int($array[1]);
		push @ref2, $array[3];
		push @alt2, $array[4];
		if($array[9]!~/0\|1/) 
		{
			push @phase2, 0;
		}
		else{push @phase2, 1;}
	}
}

print STDERR "$#chr1\t$#pos1\t$#ref1\t$#alt1\t$#phase1\n";
print STDERR "$#chr2\t$#pos2\t$#ref2\t$#alt2\t$#phase2\n";

close IN1;
close IN2;
my ($i,$j);
$i=0;
$j=0;

while( $i <= $#chr1 && $j<=$#chr2)
{
    my $cmp1=($chr1[$i] <=> $chr2[$j] );
    my $cmp2=($pos1[$i] <=> $pos2[$j] );
    my $cmp3=($ref1[$i] cmp $ref2[$j] );
    my $cmp4=($alt1[$i] cmp $alt2[$j] );
    my $cmp= ($cmp1*9+$cmp2*3+$cmp3);
    if($cmp>0)
    {
      #  print "fp:$chr2[$j]\t$pos2[$j]\n";
        $j++;
    }
    elsif($cmp<0)
    {
      #  print "fn:$chr1[$i]\t$pos1[$i]\n";
    	$i++;
    }
    else
    {
    if($cmp4==0)
    {
        my $cmpx=($phase1[$i] == $phase2[$j])?0:1;
        $info1[$i]=~/QNAME=(.+);QSTART=(.+);/;
        print "$chr1[$i]\t$pos1[$i]\t$cmpx\n";
    }
    $i++;$j++;
    }
}
