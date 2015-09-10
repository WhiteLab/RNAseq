#!/usr/bin/perl

# written by Miguel Brown 2010-Nov-30. Uses sample/mir tables to select either by top percentage, fixed value, or read count

use strict;

if(@ARGV != 4){
    print STDERR "Usage: $0 {table} {option: percent, fixed or count} {cutoff value} {output for mirs by samp}\n";
    exit(1);
}

my $input= $ARGV[0];
my $opt=$ARGV[1];
my $value=$ARGV[2];
my $out=$ARGV[3];

open(INPUT, "<", $input) or die ("Cannot open file $input\n");
my $head=<INPUT>;
chomp $head;
my @head=split /\t/, $head;
my %temp;
my %sum;

while(my $line=<INPUT>){
    chomp $line;
    my @line=split /\t/,$line;
    if($opt eq "count"){
	count(\@line,\$value);
    }
    else{
	sort_vals(\%temp,\@line,\@head,\%sum);
    } 
}
close(INPUT);

if($opt ne "count"){
    calc_top(\%temp,\%sum,\$value);
}

sub count{
    my ($info,$val)=@_;
    my $len=@$info;
    my $sum=0;
    for(my $i=1;$i<$len;$i++){
	$sum+=$$info[$i];
    }
    if($sum >= $$val){
	print "$$info[0]\n";
    }
}

sub sort_vals{
    my ($temp,$info,$head,$sum)=@_;
    my $len=@$info;
    for(my $i=1;$i<$len;$i++){
	if($opt eq "percent"){
	    $$sum{$head[$i]}+=$$info[$i];
	}
	$temp{$head[$i]}{$$info[$i]}{$$info[0]}=0;
    }
}

sub calc_top{
    my ($temp,$sum,$val)=@_;
    my %list;
    my $lim=0;
    if($opt eq "fixed"){
	$lim=1;
    }
    open(OUT, ">", $out);
    foreach my $samp(sort keys %$temp){
	print OUT "$samp";
	foreach my $ct(sort {$b <=> $a} keys %{$$temp{$samp}}){
	    if($lim>$$val){
		$lim=0;
		print OUT "\n";
		if($opt eq "fixed"){
		    $lim=1;
		}
		last;
		
	    }
	    
	    elsif($opt eq "fixed"){
		$lim+=1;
		foreach my $mir(sort keys %{$$temp{$samp}{$ct}}){
		    print OUT "\t$mir";
		    if(!exists $list{$mir}){
			$list{$mir}=1;
		    }
		}
	    }
	    
	    
	    else{
		foreach my $mir(sort keys %{$$temp{$samp}{$ct}}){
		    print OUT "\t$mir";
		    $lim+=$ct/$$sum{$samp};
		    if(!exists $list{$mir}){
			$list{$mir}=1;
		    }
		}
	    }
	    
	}
    }
    foreach my $mir(sort keys %list){
	print "$mir\n";
    }
    close OUT;
}
