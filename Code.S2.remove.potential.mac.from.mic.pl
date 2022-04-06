#!/usr/bin/perl
=head1 NAME

 Code.S2.remove.potential.mac.from.mic.pl 
 
=head1 DESCRIPTION

 This script is used for remove potential MAC contigs from MDA assemblies
 BLAST+ need to be installed.

=head1 AUTHOR

 Liping Lyu [liping_lv@foxmail.com]

=head1 USAGE

 [USAGE]
  perl Assembly_Procedure_vx.x.pl  [options] -in in.fasta
 [options]
	General options:
	-in [file]: input.MDA.assembly.fasta
	-p [name]: the prefix of outputs. Default [assembly]
	-o [path]: the directory of outputs. Default is current directory
	-t [int]: threads. [2]
	-mac [file]: specify the MAC.fasta

=cut 

#===================================================================================================================================================#
use File::Basename;
use Getopt::Long; 
use strict;

my $infile="";
my $outPrefix="assembly";
my $outPath="./";
my $thread=2;
my $macFile;
my $help;
GetOptions(
  'in=s'    =>\$infile,
  'p=s' =>\$outPrefix,
  'o=s'    =>\$outPath,
  't=s'   =>\$thread,
  'mac=s'  =>\$macFile,
  'h' =>\$help,
);
chop $outPath if($outPath=~/\/$/);

die `pod2text $0` unless(-s $infile and -s $macFile and !$help );

# creat output Path
mkdir $outPath unless(-e $outPath);
# load genome fasta
print "\tLoading Genome.fasta\n";
my %FastaSeq;
load_fasta($infile,\%FastaSeq);

my $step_code="MIC";
print "Running Step $step_code : remove the MIC sequences with high identity with MAC-sequences\n";
my $step_out_dir="$outPath/Step.$step_code.Remove.high-identity.to.MAC";
mkdir $step_out_dir unless(-e $step_out_dir);
my $outFile="$step_out_dir/$outPrefix.no.MAC-similar.sequences.in.MIC.fasta";

my $blast_database_dir=$step_out_dir."/MAC_index";
my $blast_database=$blast_database_dir."/MAC";
my $map_result=$step_out_dir."/MIC_map2_MAC.tsv";
my $evalue_limit=1e-10;
my $Pi_limit=95;
my $qCov_limit=97;

mkdir $blast_database_dir unless(-e $blast_database_dir);
my $cmd="makeblastdb -in $macFile -dbtype nucl -parse_seqids -out $blast_database ";
`$cmd`;
my $cmd="blastn -task megablast -num_threads $thread -evalue $evalue_limit -max_target_seqs 5 -ungapped  "
        ."-db $blast_database "
        ."-query $infile "
        ."-out  $map_result "
        ."-outfmt '6 qseqid qlen qcovs sseqid slen qstart qend sstart send pident length bitscore evalue' ";
`$cmd`;
my $deleted_mic_size=filt_mac_out(\%FastaSeq,$map_result,$Pi_limit,$qCov_limit);
my $deleted_mic_size=int($deleted_mic_size/1000000);
print "\t$deleted_mic_size Mb in $infile were removed\n";
output_a_sequence_hash(\%FastaSeq,$outFile);

sub load_fasta {
    my $file=shift;
    my $hash=shift;
    open Genome,"<","$file" or die $!;
    {
        local $/=undef;
        my $a=<Genome>;
            $a=~s/\n/!/g;
        my @b=split />/,$a;
        for(@b){
            next unless(m/^(.*?)!(.*?)$/);
            my $name=$1;
            my $seq=$2;
            $seq=~s/!//g;
            $$hash{$name}=$seq;
        }
    }
    close Genome;
    return $hash;
}

sub filt_mac_out {
	my ($sequences,$map_result,$Pi_limit,$qCov_limit)=@_;
	my $deleted_mic_size=0;
	open BlastOut,"<","$map_result" or die $!;
	while(<BlastOut>){
		chomp;
		my ($qseqid,$qlen,$qcovs,$sseqid,$slen,$qstart,$qend,$sstart,$send,$pident,$length,$bitscore,$evalue)=split /\t/,;
		next if(  $pident<$Pi_limit );
		my $qcov=$length/$qlen*100;
		if($qcov>$qCov_limit){
			$deleted_mic_size+=$qlen if(exists $$sequences{$qseqid});
			delete $$sequences{$qseqid};
		}
	}
	return $deleted_mic_size;
}

sub output_a_sequence_hash {
	my $sequences=shift;
	my $outFile=shift;
	open OUT,">","$outFile" or die $!;
	for my $ctgname(sort keys %{$sequences}){
		print OUT ">$ctgname\n$$sequences{$ctgname}\n";
    	}
	close OUT;
}
