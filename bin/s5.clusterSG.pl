#!/usr/bin/perl -w

use Getopt::Long;
use strict;
use vars qw/$opt_help $opt_range/;

$opt_range=8;
GetOptions("help|h", "range|r=s");
my $in=$ARGV[0];

if ($opt_help) {
    print "perl -.pl postProc.out [-r cluster-range]\n";
    print "-h\tHelp info\n";
    print "-r\tCluster range [$opt_range]\n";
    exit;
}
#PRJNA615032     SRX8089273      F       77      33      AACCAACTTT      39      0.000110972948209779    Novel

my %H;
my $cnt_sum=0;
open(IN, $in) || die "Bad input\n";
while(my $line=<IN>){
    next if ($line=~ /^@/ || $line=~/^#/ || $line=~/^$/);
    chomp($line);
    my ($prj, $srx, $s, $l, $r, $o, $cnt, $ratio, $anno) =split(/\t/, $line);
    my $sfactor=sizeFactor($cnt/$ratio);
    $cnt_sum+=$cnt*$sfactor;
    if(defined($H{$r})){
	$H{$r}->{'cnt'}+=$cnt*$sfactor;
	my $q=getQual($anno,$o);
	if($H{$r}->{'q'} < $q){
	    $H{$r}->{'q'}=$q;
	    $H{$r}->{'anno'}=$anno;
	    $H{$r}->{'o'}=$o;
	}
	$H{$r}->{'s'}=addUniq($H{$r}->{'s'},$srx);
	$H{$r}->{'p'}=addUniq($H{$r}->{'p'},$prj);
	
    }else{
	my $href;
	$href->{'anno'}=$anno;
	$href->{'cnt'}=$cnt*$sfactor;
 	$href->{'q'}=getQual($anno,$o);
	$href->{'o'}=$o;
	$href->{'s'}=[];
	push(@{$href->{'s'}},$srx);
	$href->{'p'}=[];
	push(@{$href->{'p'}},$prj);
	$H{$r}=$href;
    }
}
close(IN);

my @Csort=sort {$H{$b}->{'cnt'} <=> $H{$a}->{'cnt'}} (keys %H);
for my $i (0..$#Csort){
    next unless(defined($H{$Csort[$i]}));
    for my $j (($i+1)..$#Csort){
	next unless(defined($H{$Csort[$j]}));
	if(abs($Csort[$i]-$Csort[$j]) <= $opt_range){
	    $H{$Csort[$i]}->{'cnt'} += $H{$Csort[$j]}->{'cnt'};
	    if($H{$Csort[$j]}->{'q'} > $H{$Csort[$i]}->{'q'}){
		$H{$Csort[$i]}->{'q'}=$H{$Csort[$j]}->{'q'};
	    }
	    for my $x (@{$H{$Csort[$j]}->{'p'}}){
		$H{$Csort[$i]}->{'p'} = addUniq($H{$Csort[$i]}->{'p'}, $x);
	    }
	    for my $x (@{$H{$Csort[$j]}->{'s'}}){
		$H{$Csort[$i]}->{'s'} = addUniq($H{$Csort[$i]}->{'s'}, $x);
	    }
	    delete($H{$Csort[$j]});
	}
    }
}

for my $k (sort {$a<=>$b} (keys %H)){
    my $anno=$H{$k}->{'anno'};
    my $cnt=$H{$k}->{'cnt'};
    my $q=$H{$k}->{'q'};
    my $o=$H{$k}->{'o'};
    my $s=join(",",@{$H{$k}->{'s'}});
    my $p=join(",",@{$H{$k}->{'p'}});
    my $ratio=$cnt/$cnt_sum;
    print "$k\t$cnt\t$ratio\t$anno\t$o\t$p\t$s\n";
    
    
}




sub getQual{
    my ($anno, $o)=@_;
    if($anno eq "LowQual" ){
	return(0);
    }elsif($anno eq "Novel"){
	return(1);
    }else{
	return(2);
    }
}

sub addUniq{
    my ($aref, $x)=@_;
    my @ret=@{$aref};
    for my $i (@ret){
	if($i eq $x){
	    return(\@ret);
	}
    }
    push(@ret, $x);
    return(\@ret);
}

sub sizeFactor{ # from 0 to 1, 1 means count everything.
    my ($tot)=@_;
    my $max=100;
    if($tot<$max){
	return(1);
    }else{
	return($max/$tot);
    }
}
