#!/usr/bin/perl
#modified by Miguel Brown, 2009-Nov-19
# used to extract samples of interest organized by row
use strict;
if(@ARGV != 5){
    print STDERR "Usage: $0 {table} {selection list} {y or n to indicate whether the table has a header row} {Start position of labels, array style} {y or n to remove extra spaces}\n";
    exit(1);
}
my $table=$ARGV[0];
my $list=$ARGV[1];
my $h=$ARGV[2];
my $s=$ARGV[3];
my $t=$ARGV[4];
open (LIST,"<",$list) or die("Fail $list\n");

#my $discard=<LIST>;

my %list;
while(my $line=<LIST>){
    chomp $line;
    my @line=split /\t/,$line;
    $list{$line[0]}=undef;
}
close(LIST);

open (TABLE,"<",$table) or die ("Fail $table\n");
if($h =~ /Y/i){
    my $header=<TABLE>;
    print $header;
#    $header=<TABLE>;
#    print $header;
}
while(my $line=<TABLE>){
    chomp $line;
    if($t =~ /Y/i){
	$line =~ s/\s+/\t/g;
    }
    my @line=split /\t/,$line;
    if(exists($list{$line[$s]})){print $line."\n";}
#    else{
#	print STDERR "$line[$s] not found\n";
#    }
}
close(TABLE);
