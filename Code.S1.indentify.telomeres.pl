#!/usr/bin/perl
use File::Basename;

my ($infile)=@ARGV;
die "[USAGE] perl $0 input.fasta\n"  unless( -s $ARGV[0] );

print "#"x"39"."\n".
"#"." "x"37"."#\n".
"#        Telomere_identifier          #\n".
"#"." "x"37"."#\n".
"#"x"39"."\n";

# load Fasta file
print "[Loading] loading contigs from $infile\n";
my %ctgs;
load_fasta($infile,\%ctgs);
my $total_ctgs=keys %ctgs;
print "          loaded $total_ctgs contigs\n";
die "Error: Unable to load sequences from $infile\n!" unless((keys %ctgs)>0);

# check $output_list
print "[Output] the outputs will be saved at $infile.Telomeres.tsv\n";
open OUT,">","$infile.Telomeres.tsv" or die $!;

# search Telomeres
print "[Searching] ...\n\n";
my $TSVinfo="ctgID\ttelo_count\t5'telo-size\t3'telo-size\n";
my %tel_type_stat;
my %tel5Len_stat;
my %tel3Len_stat;
for my $ctgname(sort keys %ctgs){
	my $ctg=uc $ctgs{$ctgname};
	   $ctg=~s/-/N/g;
	my ($teltype,$tel5Len,$tel3Len)=($search_mode eq "p") ? detect_with_partial_mode($ctg,0.8,100) : detect_with_end2end_mode($ctg) ;
	$tel_count{$ctgname}=($teltype==8) ? 2 : ($teltype>0) ? 1 : 0;
	$tel_type_stat{$teltype}++ if($teltype>0);
	$tel5Len_stat{$tel5Len}++ if($tel5Len>0);
	$tel3Len_stat{$tel3Len}++ if($tel3Len>0);
	$TSVinfo.="$ctgname\t$tel_count{$ctgname}\t$tel5Len\t$tel3Len\n";
}


# output result to STDOUT

#my $total_ctgs=keys %ctgs;
my $two_tels=$tel_type_stat{8};
my $single_tels5=$tel_type_stat{5};
my $single_tels3=$tel_type_stat{3};
my $single_tels=$tel_type_stat{5}+$tel_type_stat{3};
my $none_tels=$total_ctgs-$two_tels-$single_tels;
my $two_per=sprintf("%.3f",$two_tels/$total_ctgs);
my $single_per=sprintf("%.3f",$single_tels/$total_ctgs);
my $single5_per=sprintf("%.3f",$single_tels5/$total_ctgs);
my $single3_per=sprintf("%.3f",$single_tels3/$total_ctgs);
my $none_per=sprintf("%.3f",$none_tels/$total_ctgs);
$two_per*=100;
$single_per*=100;
$none_per*=100;
$single5_per*=100;
$single3_per*=100;
print "[Summary]\n";
print "Contigs with both telomeres: $two_tels [$two_per\%]\n".
    "Contigs with single telomeres: $single_tels [$single_per\%]\n".
	"        with 5'telomere: $single_tels5 [$single5_per\%]\n".
	"        with 3'telomere: $single_tels3 [$single3_per\%]\n".
    "Contigs without any telomere: $none_tels [$none_per\%]\n";
print OUT "#both telomeres: $two_tels [$two_per\%];\n".
			"#single telomere: $single_tels [$single_per\%]\n".
			"#none: $none_tels [$none_per\%]\n".
			$TSVinfo;



sub load_fasta {
    my $file=shift;
    my $hash=shift;
	if($file=~/gz$/){
		open Genome,"zcat $file|" or die $!;
	}else{
		open Genome,"<","$file" or die $!;
	}
    {
		local $/=undef;
        my $a=<Genome>;
        $a=~s/\n/!/g;
        my @b=split />/,$a;
        for(@b){
			next unless(m/^(.*?)!(.*?)$/);
			my $name=$1;
			my $seq=$2;
			$name=$1 if($name=~/^(.*?)(\s|\t|\|)/);
			$seq=~s/!//g;
			$$hash{$name}=$seq;
		}
    }
	close Genome;

	if((keys %$hash)==0){
		open Genome,"<","$file" or die $!;
		my $name;
		while(<Genome>){
			chomp;
			if(m/>(.*)$/i){
				$name=$1;
				$name=$1 if($name=~/^(.*?)(\s|\t|\|)/);
			}else{
				$$hash{$name}.=$_;
			}
		}
		close Genome;
	}
	return $hash;
}


sub detect_with_partial_mode {
	my ($ctg,$cutoff,$fakeLengthLimit)=@_;
	my $ctgLen=length $ctg;
    my $effect_seqLen=0;
	my $teltype=0;
	my $tel5Len=0;
	my $tel3Len=0;
	my $FP_5Len=0;
	my $FP_3Len=0;
	my $ctgRC=get_reverse_complementry($ctg);
	
	if($ctg=~/^([ATCGN]*?)([AC]*CCCCAA[AC]*AACCCC)([ATCGN]+?)$/i and (length $1)<=$fakeLengthLimit and (length $3)>=$ctgLen*$cutoff ){
		$FP_5Len=length $1;
		$tel5Len=length $2;
		$teltype+=5;
    }
	if($ctgRC=~/^([ATCGN]*?)([AC]*CCCCAA[AC]*AACCCC)([ATCGN]+?)$/i and (length $1)<=$fakeLengthLimit and (length $3)>=$ctgLen*$cutoff ){
		$tel3Len=length $2;
		$FP_3Len=length $1;
		$teltype+=3;
    }
	$tel5Len+=$FP_5Len;
	$tel3Len+=$FP_3Len;
	return ($teltype,$tel5Len,$tel3Len);
}



sub detect_with_end2end_mode {
	my $ctg=shift;
	my $teltype=0;
	my $tel5Len=0;
	my $tel3Len=0;
	if($ctg=~/^([AC]*CCCCAA[AC]*AACCCC)([ATCGN]+?)(GGGGTT[GT]*TTGGGG[GT]*)$/){
		$tel5Len=length $1;
		$tel3Len=length $3;
		$teltype=8
	}elsif($ctg=~/^([AC]*CCCCAA[AC]*AACCCC)([ATCGN]+?)$/){
		$tel5Len=length $1;
		$teltype=5
	}elsif($ctg=~/^([ATCGN]+?)(GGGGTT[GT]*TTGGGG[GT]*)$/){
		$tel3Len=length $2;
		$teltype=3
	}
	return ($teltype,$tel5Len,$tel3Len);
}

sub get_reverse_complementry {
	my $seq=shift;
	$seq=~tr/ATCG/TAGC/;
	$seq=reverse $seq;
	return $seq;
}
