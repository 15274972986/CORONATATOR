#!/usr/bin/perl -w

use Getopt::Long;
use strict;
use vars qw/$opt_help $opt_leader $opt_canon/;

GetOptions("help|h","leader|l=s", "canon|c");


$opt_leader=74;
if ($opt_help) {
    print "Usage: sort in.bp | perl -.pl [-l leader-position]\n";
    print "-h\tHelp\n";
    print "-l\tLeader position [$opt_leader]\n";
    print "-c\tCanonical only\n";
    exit;
}
#"$name\t$isFirst\t$isRev\t$isSuppl\t$pos\t$cigar\tRESULTSL\t@$match,$rleft\t$pos,$rright\n"
my @Aread;
my @Apair;
my $prename="";
my $preisF="";
my $prestrand="";
while(my $line=<STDIN>){
    chomp($line);
    my ($name,$isF,$isRev,$isS,$pos,$cigar,$type,$lstr,$rstr, $soverlap) = split(/\t/,$line);
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
    if($name ne $prename){
	if($prename ne ""){
	    my $d0= dumpRead(\@Aread);
	    push(@Apair,@$d0);
	    my $d = dumpPair(\@Apair);
	    for my $x (@$d){
		my $l=$x->{'l'};
		my $r=$x->{'r'};
		my $o=$x->{'o'};
		my $i=$x->{'i'};
		if($o eq ""){
		    $o="*";
		}
		if( ! $opt_canon  || abs($l-$opt_leader)<4 ){
		    print "$prename\t$prestrand\t$l\t$r\t$o\t$i\n";
		}
	    }
	}
	$prestrand=$strand;  # could be trouble
	$prename=$name;
	$preisF=$isF;
	@Aread=();
	@Apair=();
    }elsif($preisF ne $isF){
	my $d0=dumpRead(\@Aread);
	push(@Apair,@$d0);
	@Aread=();
	$preisF=$isF;
    }
    my ($lref, $lread)=split(/,/,$lstr);
    my ($rref, $rread)=split(/,/,$rstr);
    if($lref eq "*" && $lread eq "*"){ #add right half bp to @Aread
	my %h=('dir' => "R",'cigar'=> $cigar, 'ref'=> $rref, 'read'=> $rread, 'o'=> $soverlap);
	push(@Aread, \%h);
    }elsif($rref eq "*" && $rread eq "*"){ #add left half bp to @Aread
	my %h=('dir' => "L",'cigar'=> $cigar,'ref'=> $lref, 'read'=> $lread, 'o'=> $soverlap);
	push(@Aread, \%h);
    }else{ # full bp add to @Apair
	my %h=('l'=>$lref,'r'=>$rref,'o'=>$soverlap, 'i'=>$type.",".$cigar.",".$cigar);
	push(@Apair, \%h);	
    }
}

my $d0=dumpRead(\@Aread);
push(@Apair, @$d0);
my $d=dumpPair(\@Apair);
for my $x (@$d){
    my $l=$x->{'l'};
    my $r=$x->{'r'};
    my $o=$x->{'o'};
    my $i=$x->{'i'};
    if($o eq ""){
	$o="*";
    }
    if( ! $opt_canon  || abs($l-$opt_leader)<4 ){
	print "$prename\t$prestrand\t$l\t$r\t$o\t$i\n";
    }
}


#my @a=(10,20,28, 100,110,90);
#my $b=dumpRead(\@a);
#print "Debug a@a\nb@$b\n";

sub dumpRead{ #cluster number in array
    my ($aref)=@_;
    my @ret;
    for my $i (0..$#{$aref}){
	for my $j (($i+1)..$#{$aref}){
	    my $x=combineBP($aref->[$i], $aref->[$j]);
	    if($x){
		push(@ret, $x);
	    }
	}
    }
    return(\@ret);
}
sub combineBP{ # combine half bps
    my ($href1, $href2)=@_;
    if($href1->{'cigar'} eq $href2->{'cigar'} || $href1->{'dir'} eq $href2->{'dir'}){
	return();
    }
    my ($l, $r);
    if($href1->{'dir'} eq "R"){
	$r=$href1;
	$l=$href2;
    }else{
	$r=$href2;
	$l=$href1;
    }
    if(leftHang($l->{'cigar'}) < leftHang($r->{'cigar'}) && rightHang($l->{'cigar'}) > rightHang($r->{'cigar'}) ){
	my $loverlap=$l->{'read'} - $r->{'read'} +1;
	my $soverlap;
	if($loverlap>=0){
	    $soverlap=substrSafe($r->{'o'} , 0, $loverlap);
	}else{
	    $soverlap="N"x(-$loverlap);
	}
	my $i="Combined,".$l->{'cigar'}.",".$r->{'cigar'};
	my %h=('l'=> $l->{'ref'}, 'r'=> $r->{'ref'}, 'o'=>$soverlap, 'i'=>$i);
#	my $tmp=printHash(\%h);
#	print "Debug combine $tmp\n";
	return(\%h);
    }else{
	return();
    }
}
sub leftHang{
    my ($cigar)=@_;
    if($cigar =~ /^(\d+)[HS]/){
	return($1);
    }else{
	return(0);
    }
}
sub rightHang{
    my ($cigar)=@_;
    if($cigar =~ /(\d+)[HS]$/){
	return($1);
    }else{
	return(0);
    }
}
sub dumpPair{
    my ($aref)=@_;
    my @ret;
    my @keep;
#    for my $i (0..$#{$aref}){
#	$keep[$i]=1;
#    }
    for my $i (0..$#{$aref}){
	$keep[$i]=1;
	for my $j (($i+1)..$#{$aref}){
	    my $x=compareBP($aref->[$i],$aref->[$j]);
	    if($x==0){
	    }else{
		$keep[$i]=0;
		last;
	    }
	}
    }
    for my $i (0..$#keep){
	if($keep[$i]==1){
	    push(@ret, $aref->[$i]);
	}
    }
    return(\@ret);
}

sub compareBP{ #compare full bps
    my ($href1, $href2)=@_;
    if(abs($href1->{'l'} - $href2->{'l'}) > 5){
	return(0);
    }elsif(abs($href1->{'r'} - $href2->{'r'}) > 5){
	return(0);
    }elsif($href1->{'o'} eq $href2->{'o'}){
	return(1);
    }else{
	return(2);
    }
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
sub printHash{
    my ($href)=@_;
    my $ret="";
    for my $k (keys %$href){
	my $v=$href->{$k};
	$ret.=",$k:$v";
    }
    return($ret);
}
