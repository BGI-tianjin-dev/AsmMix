#!/usr/bin/perl
# Author: Wu Pei

my $vcf_file1 =  $ARGV[0];
my @array;
open IN1, "<$vcf_file1";
while(<IN1>)
{
        $line=$_;
        if(/^##/){print ;}
        elsif(/^#/){print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n";}
        else
        {
                @array = split /\t/, $line;
                my $n_cnt= () = $array[4] =~ /N/g;
                if(length($array[4])*4 < $n_cnt*5) {next;}
                print "$array[0]\t$array[1]\t$array[2]\t$array[3]\t$array[4]\t$array[5]\t$array[6]\t$array[7]\t$array[8]\t";
                if($array[9] eq "./.") {print "0|";}
                elsif($array[9] eq "1/1") {print "1|";}

                chomp($array[10]);
                if($array[10] eq "2/2") {print "2\n";}
                elsif($array[10] eq "1/1"){print "1\n";}
                elsif($array[10] eq "./."){print "0\n";}
        }
}
