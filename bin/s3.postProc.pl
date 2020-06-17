#!/usr/bin/perl -w

use Getopt::Long;
use strict;
use vars qw/$opt_help $opt_canon $opt_leader $opt_project $opt_sample $opt_range/;

my $MAX_VAR=100;

$opt_leader=74;
$opt_project="PRJXXX";
$opt_sample="SRXXXX";
$opt_range=8;
GetOptions("help|h", "canon|c", "leader|l=s", "project|p=s", "sample|s=s", "range|r=s");

if ($opt_help) {
    print "perl -.pl sgRNA.out [-c]\n";
    print "-c\tCanon only\n";
    print "-h\tHelp info\n";
    print "-p\tProject name\n";
    print "-s\tSample name\n";
    print "-l\tLeader pos";
    print "-r\tRange to cluster together\n";
    exit;

}

my %H;
open(IN,$ARGV[0]);
while(my $line = <IN>){
    next if ($line=~ /^@/ || $line=~/^#/ || $line=~/^$/);
    chomp($line);
#SRR1942990.10000101     R       66      28535   TTGATTTTAACGAA  ResultSM,24S77M,24S77M
    my ($name,$s,$l,$r,$o,$i) =split(/\t/, $line);
    if($opt_canon && abs($l-$opt_leader) > $opt_range){
	next;
    }
    $line="$s\t$l\t$r\t$o";
#    print "$line\n";
    if(defined($H{$line})){
	$H{$line}+=1;
    }else{
	$H{$line}=1;
    }
}
my @sorted = sort{$H{$b} <=> $H{$a}} keys %H;
my $n=$#sorted;
if($n>$MAX_VAR){
    $n=$MAX_VAR;
}
@sorted = @sorted[0..$n];
my @keep;
for my $i (0..$#sorted){
    $keep[$i]=1;
}
for my $i (0..$#sorted){
#    print "$i $sorted[$i] $H{$sorted[$i]}\n";
    for my $j (($i+1)..$#sorted){
	my $diff=compareBPS($sorted[$i],$sorted[$j]);
	if($diff==0){
	    $keep[$j]=0;
	    $H{$sorted[$i]}+=$H{$sorted[$j]};
	}
    }

}
my @final;
for my $i (0..$#keep){
    if($keep[$i]==1){
	push(@final, $sorted[$i]);
	print "$opt_project\t$opt_sample\t$sorted[$i]\t$H{$sorted[$i]}\n";
    }

}
#seek IN, 0, 0;
#while(my $line=<IN>){

#}
close(IN);

sub compareBPS{
    my ($b1, $b2)=@_;
    my ($s1, $l1, $r1, $o1)=split(/\t/,$b1);
    my ($s2, $l2, $r2, $o2)=split(/\t/,$b2);
    if(abs($l1-$l2) <= $opt_range  && abs($r1-$r2) <= $opt_range){
	return(0);
    }else{
	return(1);
    }
}

