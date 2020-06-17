#!/usr/bin/perl -w

use Getopt::Long;
use strict;
use vars qw/$opt_ref $opt_help $opt_klen $opt_leader $opt_minh $opt_debug/;
#use Array::Utils qw(:all);

my $flagReverse=16;
my $flagFirst=64;
my $flagSuppl=2048;

$opt_leader=-100;  #specify 3' position of leader, used to quantify gRNA in relation to sgRNA, use -100 as a hack to ignore
$opt_klen=8;    #kmer length used to do pseudo alignment
$opt_minh=8;     #minimal overhang
GetOptions("ref|r=s","help|h","klen|k=s", "leader|l=s", "minh|m=s", "debug|d");

if ($opt_help) {
    print "Usage: samtools view in.bam | perl -.pl -r consensus.fasta [-k kmer-length] [-l leader-3-position] [-m min-overhang]\n";
    print "-r\tConsensus reference fasta\n";
    print "-k\tKmer length [$opt_klen]\n";
    print "-l\tLeader 3' position, used to check for gRNA, unless =0 \n";
    print "-m\tMinimal overhang [$opt_minh]\n";
    print "-d\tDebug info\n";
    exit;

}


open(REF, $opt_ref) or die ("Bad ref file $opt_ref");
my $Ref=""; #to hold the consensus ref
my $Href; #only hash ref can be returned by function??
while(<REF>){
    next if(/^>/);
    chomp;
    $Ref.=uc($_);
}
$Href=buildHash($Ref, $opt_klen);
#my $x=substr($ref, 28255, 10);
#print "Debug $x\n";

#M04943:62:000000000-CLTCJ:1:1105:5373:9058      163     MT121215.1  15      60      151M    =       38      174     CCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACG EEDEEFFFFEFFGGGGGGGGGGGGHGHHHHGGHHHHHHHHHHHHHHHHHHHHHHHHHGGGGHHHHHHHHHHHHHHHGHHGHHHHHHGGGGGGHHHHHHHHHHHHHHHHHHHGGFHHHHHHHHHHHHHHHHHHHHHHHGGGHHHHHGHHHHF NM:i:0  MD:Z:151        MC:Z:151M       AS:i:151        XS:i:0
while(my $line = <STDIN>){
    next if ($line=~ /^@/ || $line=~/^#/ || $line=~/^$/);
    chomp($line);
    my ($name,$flag,$rname,$pos,$qmap,$cigar,$mate_ref,$mate_pos,$insert,$seq) =split(/\t/, $line);
    my $pos0=$pos-1;     #pos start from 1
    my $isRev=(($flag & $flagReverse) == $flagReverse)?"isRev":"notRev";
    my $isFirst=(($flag & $flagFirst) == $flagFirst)?"isF":"notF";
    my $isSuppl=(($flag & $flagSuppl) == $flagSuppl)?"isS":"notS";
    my $rlen=lenCIGAR($cigar);
    if($cigar=~/^(\d+)M$/){  #perfect match, some will support gRNA
	my $lmatch=$1;
#	my $min_match=30; #match to right of bp to call support
	if(($opt_leader-$pos)>$opt_minh+14 && ($lmatch+$pos-$opt_leader) > $opt_minh+2){
	    print "$name\t$isFirst\t$isRev\t$isSuppl\t$pos\t$cigar\tResult0\t$opt_leader,*\t$opt_leader,*\tLEADERTRS\n";
	}
    }else{
	my ($btype, $in);
	my $lhang=0;
	my $rhang=0;
	my $match=[];
	if($cigar =~ /^(\d+)S(\d+)M/){	#left 5' break point
	    $lhang = $1;
	    my $mlen =$2;
	    if($lhang > $opt_minh){
		$in=substr($seq, 0, $lhang);
		$match = findMatch($Href, $in);
		if(scalar(@$match) == 1){ ## should be more check, should be left of $pos ??
		    $match->[0] += $lhang-1;
		    my $rleft=$lhang-1;
		    my $rright=$lhang;
		    my $overlap=rightExtension(substr($seq,$lhang,$mlen),$match->[0] +1);
		    $match->[0] +=$overlap;
		    $rleft +=$overlap;
		    my $soverlap=substr($seq,$rright,$overlap);
		    print "$name\t$isFirst\t$isRev\t$isSuppl\t$pos\t$cigar\tResultSM\t@$match,$rleft\t$pos0,$rright\t$soverlap\n";
		    if($opt_debug){
			my $s1="Ref:".substrSafe($Ref,$match->[0]-4,5);
			$s1.="X".substrSafe($Ref,$match->[0] + 1, 5);
			my $s2="Read:".substrSafe($seq, $rleft-4, 5);
			$s2.="X".substrSafe($seq,$rleft+1,5);
			print "DebugSeq $s1 $s2\n";
			$s1="Ref:".substrSafe($Ref,$pos0-5,5);
			$s1.="X".substrSafe($Ref,$pos0, 5);
			$s2="Read:".substrSafe($seq, $rright-5, 5);
			$s2.="X".substrSafe($seq,$rright,5);
			print "DebugSeq $s1 $s2\n";
		    }
		}else{# when no match or more than one match, same as hard clip
		    my $soverlap=substr($seq,$lhang,20);
		    print "$name\t$isFirst\t$isRev\t$isSuppl\t$pos\t$cigar\tResultSM\t*,*\t$pos0,$lhang\t$soverlap\n";
		}
	    }
	}elsif($cigar =~/^(\d+)H\d+M/){
	    $lhang=$1;
	    my $soverlap=substr($seq,0,20);
	    print "$name\t$isFirst\t$isRev\t$isSuppl\t$pos\t$cigar\tResultHM\t*,*\t$pos0,$lhang\t$soverlap\n";
	    if($opt_debug){
		my $s1="Ref:".substrSafe($Ref,$pos0-5,5);
		$s1.="X".substrSafe($Ref,$pos0 , 5);
		my $s2="Read:NNNNN";
		$s2.="X".substrSafe($seq,0,5);
		print "DebugSeq $s1 $s2\n";
	    }
	}

	if($cigar =~ /(\d+)M(\d+)S$/){	#right 3' break point
	    my $mlen=$1;
	    $rhang = $2;
	    if($rhang > $opt_minh){
		$in=substr($seq, length($seq)-$rhang, $rhang);
		$match = findMatch($Href, $in);
		if(scalar(@$match) == 1){
		    my $rleft=$rlen-$rhang-1;
		    my $rright=$rlen-$rhang;
		    my $bleft=$pos0+$rlen-$rhang-$lhang-1;
		    my $overlap=leftExtension(substr($seq,$lhang,$mlen),$match->[0] -1);
		    $match->[0]-=$overlap;
		    $rright-=$overlap;
		    my $soverlap=substr($seq, $rleft-$overlap+1,$overlap);
		    print "$name\t$isFirst\t$isRev\t$isSuppl\t$pos\t$cigar\tResultMS\t$bleft,$rleft\t@$match,$rright\t$soverlap\n";
		    if($opt_debug){
			my $s1="Ref:".substrSafe($Ref,$bleft-4,5);
			$s1.="X".substrSafe($Ref,$bleft + 1, 5);
			my $s2="Read:".substrSafe($seq, $rleft-4, 5);
			$s2.="X".substrSafe($seq,$rleft+1,5);
			print "DebugSeq $s1 $s2\n";
			$s1="Ref:".substrSafe($Ref,$match->[0] -5,5);
			$s1.="X".substr($Ref,$match->[0], 5);
			$s2="Read:".substrSafe($seq, $rright-5, 5);
			$s2.="X".substrSafe($seq,$rright,5);
			print "DebugSeq $s1 $s2\n";
		    }
		}else{ # when no match or more than one match, same as hard clip
		    my $rleft=$rlen-$rhang-1;
		    my $bleft=$pos0+$rlen-$rhang-$lhang-1;
		    my $soverlap=substr(reverse($seq), $rhang, 20);
		    print "$name\t$isFirst\t$isRev\t$isSuppl\t$pos\t$cigar\tResultMS\t$bleft,$rleft\t*,*\t$soverlap\n";
		}
	    }
	}elsif($cigar =~ /\d+M(\d+)H$/){
	    $rhang=$1;
	    my $rleft=$rlen-$rhang-1;
	    my $bleft=$pos0+$rlen-$rhang-$lhang-1;
	    my $soverlap=substr(reverse($seq),0,20);
	    print "$name\t$isFirst\t$isRev\t$isSuppl\t$pos\t$cigar\tResultMH\t$bleft,$rleft\t*,*\t$soverlap\n";
	    if($opt_debug){
		my $s1="Ref:".substrSafe($Ref,$bleft-4,5);
		$s1.="X".substrSafe($Ref,$bleft + 1, 5);
		my $rleftseq=$rleft;
		if($cigar =~ /^\d+H/){ # hard clipped seq has to discount lhang
		    $rleftseq=$rleftseq-$lhang;
		}
		my $s2="Read:".substrSafe($seq, $rleftseq-4, 5);
		$s2.="X".substrSafe($seq,$rleft+1,5);
		print "DebugSeq $s1 $s2\n";
	    }
	}
    }

}

#for my $k (keys %{$Href}){
#    my $aref=$Href -> {$k};
#    my $n=scalar(@{$aref});
#    print "Debugx $k $n\n";
#}

#my $x=leftExtension("XCXXF",5,"ABCDEF");
#print "DebugExtension x$x\n";


sub lenCIGAR{
    my ($s)=@_;
    my $ret=0;
    while($s =~ /(\d+)([MSHID])/g){
	$ret+=$1 unless($2 eq "D");   # D I could be aligner specific?? use D for bwa??
    }
    return($ret);
}

sub buildHash{
    my ($r, $k)=@_;
    my %h;
    if(length($r) <= $k) {
	return();
    }
    for(my $i=0; $i<(length($r)-$k); $i++){
	my $mer= substr($r, $i, $k);
	if(defined($h{$mer})){
	    push(@{$h{$mer}}, $i);
	}else{
	    $h{$mer}=[$i];
	}
    }
    return(\%h);
}


sub findMatch{  # pseudo alignment
    my ($href, $in)=@_;
    my %h;
    my $max=(length($in)-$opt_klen)>10?10:(length($in)-$opt_klen);  
    for(my $i=0; $i <= $max; $i+=2){
	my $mer=substr($in, $i, $opt_klen);
	if(defined($href->{$mer})){
	    for my $x (@{$href->{$mer}}){
		my $y=$x-$i;
		if(defined($h{$y})){
		    $h{$y}+=1;
		}else{
		    $h{$y}=1;
		}
	    }
	}
    }
    my $ret=pickMatch(\%h);
    $ret=pickMatch2($ret,$in); # used hard coded Ref to save space??
    return($ret);
}
sub pickMatch{
    my ($href)=@_;
    my $ret =[];
    my $best =0;
#    my $tmp=printHash($href);
#    print "DebugMatch $tmp\n";
    for my $k (keys %{$href}){
	if($href->{$k} > $best){
	    $best = $href->{$k};
	}
    }
    for my $k (keys %{$href}){
	if($href->{$k} == $best){
	    push(@{$ret}, $k);
	}
    }
    return($ret);
}
sub pickMatch2{  # used hard coded ref
    my ($aref,$s)=@_;
    my $ret=[];
    my ($ss,$score,$x);
    my $best=0;
    my $slen=length($s);
    my %h;
    for my $i (0..$#{$aref}){
	$x=$aref->[$i];
	if($x<0){
	    $ss=substr($s,-$x,);
	    $score=length($ss)-hd($ss, substr($Ref,0,length($ss)));
	}elsif($x+$slen > length($Ref)){
	    $ss=substr($s,0,length($Ref)-$x);
	    $score=length($ss)-hd($ss, substr($Ref,$x));
	}else{
	    $ss=$s;
	    $score=$slen-hd($ss, substr($Ref, $x, $slen));
	}
#	my $tmp=substr($Ref, $x, $slen);
#	print "DebugMatch3 $ss $tmp\n";
	$score=$score/length($ss);
	$h{$x}=$score;
	if($score > $best){
	    $best=$score;
	}		
    }
#    my $tmp=printHash(\%h);
#    print "DebugMatch2 $tmp\n";

    if($best < 0.9){
	return($ret);
    }
    for my $k (keys %h){
	if($h{$k}==$best){
	    push(@$ret, $k);
	}
    }
    return($ret);
}
sub rightExtension{   # again, ref hard coded, allow one mismatch, return length, did not consider when go over ref
    my ($in, $p)=@_;
    my $bad=0;
    my $lastgood=0;
    for my $i (0..(length($in)-1)){
	my $xbase=substr($in, $i, 1); 
	my $rbase=substrSafe($Ref,$p+$i, 1); #potential bug
#	print "DebugrightExt $xbase $rbase\n";
	if($xbase eq $rbase){
	    $lastgood=$i+1;
	}else{
	    if($bad==1){
		last;
	    }else{
		$bad+=1;
	    }
	}
    }
    if($lastgood >=2 && substr($in,$lastgood-2,1) ne substr($Ref,$p+$lastgood-2,1)){  #hack to handle mistch one base from edge
	return($lastgood-2);
    }else{
	return($lastgood);
    }
}
sub leftExtension{   # again, ref hard coded, allow one mismatch, return length, did not consider when go over ref
    my ($in, $p)=@_;
    my $bad=0;
    my $lastgood=0;
    $in=reverse($in);
    for my $i (0..(length($in)-1)){
	my $xbase=substr($in, $i, 1);
	my $rbase=substrSafe($Ref,$p-$i, 1);
#	print "DebugleftExt $xbase $rbase\n";
	if($xbase eq $rbase){
	    $lastgood=$i+1;
	}else{
	    if($bad==1){
		last;
	    }else{
		$bad+=1;
	    }
	}
    }
    if($lastgood>=2 && substr($in, $lastgood-2,1) ne substr($Ref,$p-$lastgood+2,1)){
	return($lastgood-2);  # hack to handle mismatch one base to edge.
    }else{
	return($lastgood);
    }
}

#for(my $i=3; $i <= 3; $i++){
#    print "DDDD  $i\n";
#}
sub printHash{
    my ($href)=@_;
    my $ret="";
    for my $k (keys %$href){
	my $v=$href->{$k};
	$ret.=",$k:$v";
    }
    return($ret);
}

sub hd { #credit to perl monk
    my ($k, $l) = @_;
    my $diff = $k ^ $l;
    my $num_diff = $diff =~ tr/\0//c;
    return($num_diff);
}

sub substrSafe{
    my ($s,$p,$l)=@_;
    if($p<0){
	$s=("N"x(-$p)).$s;
	$p=0;
    }
    if($p+$l>length($s)){
	$s.=("N"x($p+$l-length($s)));
    }
    return(substr($s,$p,$l));
}
