#!/usr/bin/perl

# written by Miguel Brown 2010-Jan-21. This script takes a dataset with a header row and label column and calculates the relative frequency of each element against the sum of the column

use strict;
if(@ARGV !=2){
    print STDERR "Usage: $0 {file}{data start column - 0-based array style}\n";
    exit(1);
}
my $input= $ARGV[0];
my $s = $ARGV[1];
open(INPUT, "<", $input) or die ("Cannot open file $input\n");
my @sums;
my @data;
my $header=<INPUT>;
my @header=split /\t/, $header;
my $len=@header;
print $header;
my $i=0;

while(my $line=<INPUT>){
    chomp $line;
    my @line=split /\t/,$line;
    $data[$i]=$line;
    for(my $j=$s;$j<$len;$j++){
	$sums[$j]+=$line[$j];
    }
    $i++;
}
close(INPUT);

foreach my $cur(@data){
    my @cur=split /\t/, $cur;
    print $cur[0];
    for ($i=1;$i<$s;$i++){
	print "\t$cur[$i]";
    }

    for(my $j=$s;$j<$len;$j++){
	if($sums[$j]){
	    my $rf=$cur[$j]/$sums[$j];
	    print "\t$rf";
	}
	else{
	    print "\t0";
	}
    }
    print "\n";
}
