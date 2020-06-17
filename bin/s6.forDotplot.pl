#!/usr/bin/perl -w

use Getopt::Long;
use strict;
#1       2.67904045136522        0.00208810635336338     LowQual N       PRJNA605907     SRX7705831,SRX7705832

while(my $line=<STDIN>){
    next if ($line=~ /^@/ || $line=~/^#/ || $line=~/^$/);
    chomp($line);
    my ($r,$cnt,$ratio,$anno,$o,$prj,$srx) =split(/\t/, $line);
    my @cols=split(/,/, $srx);
    for my $i (@cols){
	print "$r\t$cnt\t$ratio\t$i\n";
    }
}
