#!/usr/bin/perl -w

use Getopt::Long;
use strict;
use vars qw/$opt_help/;
#use Array::Utils qw(:all);
GetOptions("help|h");

my $flagReverse=16;
my $flagFirst=64;
my $flagSuppl=2048;


if ($opt_help) {
    print "Usage: samtools view in.bam | perl -.pl \n";
    print "-h\tHelp info\n";
    exit;

}
my $Fcount=0;
my $Rcount=0;
while(my $line = <STDIN>){
    next if ($line=~ /^@/ || $line=~/^#/ || $line=~/^$/);
    chomp($line);
    my ($name,$flag,$rname,$pos,$qmap,$cigar,$mate_ref,$mate_pos,$insert,$seq) =split(/\t/, $line);
    my $isRev=(($flag & $flagReverse) == $flagReverse)?"isRev":"notRev";
    my $isF=(($flag & $flagFirst) == $flagFirst)?"isF":"notF";
    my $strand;
    if($isF eq "isF"){
	if($isRev eq "isRev"){
	    $strand="R";
	}else{
	    $strand="F";
	}
    }else{
	if($isRev eq "isRev"){
	    $strand="F";
	}else{
	    $strand="R";
	}
    }
    if($strand eq "F"){
	$Fcount+=1;
    }else{
	$Rcount+=1;
    }
#    print "$strand\n";
}

my $x=-1;
if(($Fcount+$Rcount) >0 ){
    $x=$Fcount/($Fcount+$Rcount);
}
print $x;
