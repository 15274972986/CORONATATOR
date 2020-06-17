#!/usr/bin/perl -w

use Getopt::Long;
use strict;
use vars qw/$opt_help $opt_anno $opt_range/;

$opt_range=8;
GetOptions("help|h", "anno|a=s", "range|r=s");
my $in=$ARGV[0];

if ($opt_help) {
    print "perl -.pl postProc.out [-a anno-file] [-r cluster-range]\n";
    print "-a\tAnno file\n";
    print "-h\tHelp info\n";
    print "-r\tCluster range [$opt_range]\n";
    exit;
}

open(ANNO, $opt_anno) || die "Bad annotation file\n";

my @Ag;
my %Hg;
while(my $line=<ANNO>){
    next if ($line=~ /^@/ || $line=~/^#/ || $line=~/^$/);
    chomp($line);
    my ($gname, $pos) =split(/\t/, $line);
    $Hg{$gname} = $pos;
    push(@Ag,$gname);
    
}
close(ANNO);

my %Hsample;
open(IN, $in) || die "Bad input\n";
while(my $line=<IN>){
    next if ($line=~ /^@/ || $line=~/^#/ || $line=~/^$/);
    chomp($line);
    my ($prj, $srx, $s, $l, $r, $o, $cnt) =split(/\t/, $line);
    my $sample=$prj."_".$srx;
    if(defined($Hsample{$sample})){
	$Hsample{$sample}+=$cnt;
    }else{
	$Hsample{$sample}=$cnt;
    }
}


my %Hmatrix;
seek IN, 0, 0;
while(my $line=<IN>){
    next if ($line=~ /^@/ || $line=~/^#/ || $line=~/^$/);
    chomp($line);
    my ($prj, $srx, $s, $l, $r, $o, $cnt) =split(/\t/, $line);
    my $x=annotateSG($r,$o);
    my $sample=$prj."_".$srx;
    my $ratio=$cnt/$Hsample{$sample};
    print "$line\t$ratio\t$x\n";
    if(defined($Hmatrix{$sample})){
	if(defined($Hmatrix{$sample}->{$x})){
	    $Hmatrix{$sample}->{$x}+=$cnt;
	}else{
	    $Hmatrix{$sample}->{$x}=$cnt;
	}
    }else{
	$Hmatrix{$sample}={};
	$Hmatrix{$sample}->{$x}=$cnt;
    }
    
}
close(IN);
my $out=$in.".matrix";
open(OUT,">$out") || die "Cannot open output file to write\n";
print OUT "SAMPLE\tORF1";
for my $i (@Ag){
    print OUT "\t$i";
}
print OUT "\tOTHER\n";
for my $k (sort(keys %Hmatrix)){
    print OUT $k;
    my $sum=0;
    if(defined($Hmatrix{$k}->{"ORF1"})){
	$sum=$Hmatrix{$k}->{"ORF1"};
    }
    print OUT "\t$sum";
    for my $i (@Ag){
	my $x=0;
	if(defined $Hmatrix{$k}->{$i}){
	    $x=$Hmatrix{$k}->{$i};
	}
	$sum+=$x;
	print OUT "\t$x";
    }
    my $other=$Hsample{$k}-$sum;
    print OUT "\t$other\n";
}

#for my $k (keys %Hsrx){
#    print "$k $Hsrx{$k}\n";
#}

sub annotateSG{ #annotation hash hard coded Hg
    my ($r,$o)=@_;
    if($o=~/LEADER/){
	return("ORF1");
    }
    for my $k (keys %Hg){
	if(abs($r-$Hg{$k})<= $opt_range){
	    return($k);
	}
    }
    if(length($o) >4) {
	return("Novel");
    }else{
	return("LowQual");
    }
}
