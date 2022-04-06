use Getopt::Long;

my $pointers_path;
my $GFFs_PATH;
my $Proteins_PATH;
my $Genome_PATH;
my $OrthoVenn_File;
my $mds_gff_PATH;
 GetOptions(
        'p=s'  =>\$pointers_path,
        'gff=s'  =>\$GFFs_PATH,
        'prot=s'	=>\$Proteins_PATH,
		'ort=s'	=>\$OrthoVenn_File,
        'g=s'   =>\$Genome_PATH,
        'mds=s' =>\$mds_gff_PATH,
 );


die "perl $0 
    -p /path/of/*.pointers.csv [contains information of pointers;format:ctgID,seq,size,start,stop]
    -gff /path/of/*.mac.genome.gff3
    -prot /path/of/*.proteins.fa
    -ort orthoClouster.txt [each line includes all the ID of homolugues genes]
    -g /path/of/*.genome.fa
    -mds /path/to/intergrated.gff3   \n" unless(-s $OrthoVenn_File);

our %Halteria_adjust;

# load gff3
chop $GFFs_PATH if($GFFs_PATH=~/\/$/);
chomp(my $GFFs_file=`ls $GFFs_PATH`);
my %gff_inf;
for my $gff (split /\n/,$GFFs_file){
    if($gff=~/^(.*?)\.mac.*?/){
        my $species=$1;
        my $type=($species =~/Hal/i) ? "locus_tag" : "Parent" ;
        load_GFF("$GFFs_PATH/$gff",\%{$gff_inf{$species}},$type);
    }else{
        die "Fail to identify $gff\n";
    }
}
print "Load GFF done\n";

# load mds.gff3
chop $mds_gff_PATH if($mds_gff_PATH=~/\/$/);
chomp(my $mds_gff=`ls $mds_gff_PATH`);
my %mds_inf;
for my $mdsFile (split /\n/,$mds_gff){
    if($mdsFile=~/^(.*?)\.Data.*?/){
        my $species=$1;
        load_mds_GFF("$mds_gff_PATH/$mdsFile",\%{$mds_inf{$species}});
    }else{
        die "Fail to identify $gff\n";
    }
}
print "Load GFF done\n";

# load orthovenn file
my %OrthoCluster;
load_orthoClusters($OrthoVenn_File,\%OrthoCluster);

# load proteins
chop $Proteins_PATH if($Proteins_PATH=~/\/$/);
chomp(my $Prot_files=`ls $Proteins_PATH`);
my %proteins;
for my $prot (split /\n/,$Prot_files){
    if($prot=~/^(.*?)\.mac.*?/){
        my $species=$1;
        load_fasta("$Proteins_PATH/$prot",\%{$proteins{$species}});
    }else{
        die "Fail to identify $prot\n";
    }
}
print "Load proteins done\n";

# load genomes
chop $Genome_PATH if($Genome_PATH=~/\/$/);
chomp(my $Geno_files=`ls $Genome_PATH`);
my %genome;
for my $genome_file (split /\n/,$Geno_files){
    if($genome_file=~/^(.*?)\.mac.*?/){
        my $species=$1;
        load_fasta("$Genome_PATH/$genome_file",\%{$genome{$species}});
    }else{
        die "Fail to identify $genome_file\n";
    }
}
print "Load genome done\n";

# load pointers
chop $pointers_path if($pointers_path=~/\/$/);
chomp(my $Pointer_files=`ls $pointers_path`);
my %pointers_pos;
for my $pointer (split /\n/,$Pointer_files){
    if($pointer=~/^(.*?)\.Data.*?/){
        my $species=$1;
        load_pointers("$pointers_path/$pointer",\%{$pointers_pos{$species}},\%{$gff_inf{$species}});
    }else{
        die "Fail to identify $pointer\n";
    }
}
print "Load pointers done\n";


#################################

open LOG,">","$0.dealing.log" or die $!;
open adjusted_site,">","$0.adjsuted.site.tsv";
print adjusted_site "clusterID\tspecies\tmacID\tgeneID\tmacStart\tmacStop\tsite.of.MSA\tsite.of.macContig\n";
for my $clusterID (sort keys %{$OrthoCluster{"cluster_content"}}){
    next unless( (keys %{$OrthoCluster{"cluster_content"}{$clusterID}})>=1);
    my %qualified;
    for my $species (sort keys %{$OrthoCluster{"cluster_content"}{$clusterID}}){
        for my $geneID (sort keys %{$OrthoCluster{"cluster_content"}{$clusterID}{$species}}){
            my $correspond_macID=$gff_inf{$species}{"geneID2macID"}{$geneID};
            $qualified{$species}++ if( (keys %{$pointers_pos{$species}{$correspond_macID}})>0);
        }
    }   
    next unless ((keys %qualified)>=1);
    print LOG "dealing $clusterID\n";
    my $cluster_member_path="./OUT.01.globalpair.cluster.members.seq";
    mkdir $cluster_member_path unless(-e $cluster_member_path);
    my $fasta_name="$cluster_member_path/$clusterID.fa";
    my $genome_name="$cluster_member_path/$clusterID.genome.fa";
    open OUT,">",$fasta_name or die $!;
    my %outputed_genomeID;
    open genome,">",$genome_name or die $!;
    for my $species (sort keys %{$OrthoCluster{"cluster_content"}{$clusterID}}){
        for my $geneID (sort keys %{$OrthoCluster{"cluster_content"}{$clusterID}{$species}}){
            print OUT ">$species\|$geneID\n$proteins{$species}{$geneID}\n";
            my $correspond_macID=$gff_inf{$species}{"geneID2macID"}{$geneID};
            print genome ">$species\|$correspond_macID\n$genome{$species}{$correspond_macID}\n" if($outputed_genomeID{"$species\|$correspond_macID"}==0);
            $outputed_genomeID{"$species\|$correspond_macID"}++;
        }
    }
    close OUT;
    close genome;

    my $mafft_cmd="mafft --maxiterate 1000 --thread 8 --globalpair $fasta_name > $fasta_name.MAFFT.fa";
    print "Runing $mafft_cmd\n";
    my $tmp=`$mafft_cmd`;
    #die "test over\n";

    my $blatcmd="blat $genome_name $fasta_name -t=dnax -q=prot -noHead $fasta_name.MAFFT.psl";
    my $tmp=`$blatcmd`;
    my %protPos_link2_nuclPos;
    load_PSL("$fasta_name.MAFFT.psl",\%protPos_link2_nuclPos);

    my %MSA_fasta;
    load_MSA_fasta("$fasta_name.MAFFT.fa",\%MSA_fasta);
    my %MSA_pointer_site;
    my %geneIDlist;
    open ADJUST,">","$cluster_member_path/$clusterID.adjsuted.txt";
    for my $species (sort keys %MSA_fasta){
        for my $geneID (sort keys %{$MSA_fasta{$species}}){
            $geneIDlist{$geneID}++;
            my $pointers_site_line;
            my $correspond_macID=$gff_inf{$species}{"geneID2macID"}{$geneID};
            my $MSA_seq=$MSA_fasta{$species}{$geneID};
            my $valid_position=0;
            for(my $i=0;$i<(length $MSA_seq);$i++){
                my $char=substr($MSA_seq,$i,1);
                if($char=~/-/){
                    $pointers_site_line.=" ";
                }elsif((keys %{$protPos_link2_nuclPos{$species}{$geneID}})>0){
                    $valid_position++;
                    #$protPos_link2_nuclPos{$species}{$geneID}{$i}>0
                    my ($contig_start,$contig_stop)=split /\./,$protPos_link2_nuclPos{$species}{$geneID}{$valid_position};
                    my $overlap=check_overlaps(\%{$pointers_pos{$species}{$correspond_macID}},$contig_start,$contig_stop);
                    my $mds_annotated_site=check_MDS_annotated(\%{$mds_inf{$species}{$correspond_macID}},$contig_start,$contig_stop);
                    if($overlap){
                        my $k=$i+1;
                        $MSA_pointer_site{$k}{"count1"}{$species}{$geneID}++;
                        $MSA_pointer_site{$k}{"count2"}{$geneID}++;
                        $pointers_site_line.="*";
						print adjusted_site "$clusterID\t$species\t$correspond_macID\t$geneID\t$contig_start\t$contig_stop\t$i\t$valid_position\n";
                    }elsif( $mds_annotated_site){
                        $pointers_site_line.="-";
                    }
                    else{
                        $pointers_site_line.="?";
                    }
                }
            }
            print ADJUST ">$species\|$geneID\|$correspond_macID\n$pointers_site_line\n$MSA_seq\n";
        }
    }
    
    close ADJSUT;
    for my $MSAsite (sort {$a<=>$b} keys %MSA_pointer_site){
        if((keys %{$MSA_pointer_site{$MSAsite}{"count2"}}) == (keys %geneIDlist)){
            print "Class.A $clusterID:$MSAsite\n";
        }elsif( (keys %{$MSA_pointer_site{$MSAsite}{"count1"}}) >=2){
            print "Class.B $clusterID:$MSAsite\n";
        }elsif( (keys %{$MSA_pointer_site{$MSAsite}{"count1"}}) >=1 and (keys %{$MSA_pointer_site{$MSAsite}{"count2"}})>=2){
            print "Class.C $clusterID:$MSAsite\n";
        }
        for my $species (sort keys %{$MSA_pointer_site{$MSAsite}{"count1"}}){
            for my $geneID (sort keys %{$MSA_pointer_site{$MSAsite}{"count1"}{$species}}){
                print "$species\t$geneID\t$MSAsite\n";
            }
        }
    }
}

sub load_PSL {
	my $file=shift;
	my $hash=shift;
    my %tmp_records;
    open INF,"<",$file or die $!;
    while(<INF>){
        my $Qspecies;
        my $Tspecies;
        my ($match,$mismatch,$repmatch,$N,$QgapCount,$QgapBases,$TgapCount,$TgapBases,$strand,$Qname,$Qsize,$Qstart,$Qend,$Tname,$Tsize,$Tstart,$Tend,$blockCount,$blockSizes,$qStarts,$tStarts)=split /\s+|\t/,;
        if($Qname=~/^(.*?)\|(.*?)$/){
            $Qspecies=$1;
            $Qname=$2;
        }
        if($Tname=~/^(.*?)\|(.*?)$/){
            $Tspecies=$1;
            $Tname=$2;
        }       
        next unless($Qspecies eq $Tspecies);
        $tmp_records{$Qspecies}{$Qname}{$_}=$match;
    }
    for my $Qspecies (sort keys %tmp_records){
        for my $Qname1 (sort keys %{$tmp_records{$Qspecies}}){
            my $best_one= (sort {$tmp_records{$Qspecies}{$Qname1}{$b} <=> $tmp_records{$Qspecies}{$Qname1}{$a}} keys %{$tmp_records{$Qspecies}{$Qname1}})[0];
            my ($match,$mismatch,$repmatch,$N,$QgapCount,$QgapBases,$TgapCount,$TgapBases,$strand,$Qname2,$Qsize,$Qstart,$Qend,$Tname,$Tsize,$Tstart,$Tend,$blockCount,$blockSizes,$qStarts,$tStarts)=split /\s+|\t/,$best_one;
            my $Qname2=$2 if($Qname2=~/^(.*?)\|(.*?)$/);
            my @BlockSizes=split /,/,$blockSizes;
            my @QStarts=split /,/,$qStarts;
            my @Tstarts=split /,/,$tStarts;
            for(my $i=0;$i<@BlockSizes;$i++){
                my $Qstart=1+$QStarts[$i];
                my $Tstart=1+$Tstarts[$i];
                my $BlockSize=$BlockSizes[$i];

                my $Qend=$Qstart+$BlockSize-1;

                for my $Qpos ($Qstart..$Qend){
                    my $distance=$Qpos-$Qstart;
                    my $Tstart_here=$Tstart+$distance*3;
                    my $Tend_here=$Tstart_here+2;
                    if($strand=~/.-/){
                        my $tmp=$Tsize-$Tstart_here+1;
                        $Tstart_here=$Tsize-$Tend_here+1;
                        $Tend_here=$tmp;
                    }
                    $$hash{$Qspecies}{$Qname2}{$Qpos}="$Tstart_here.$Tend_here";
                }
            }
        }
    } 
}

sub check_overlaps {
    my $refs=shift;
    my $query_start=shift;
    my $query_stop=shift;
    my $overlap=0;
    for my $refPos (sort {$a<=>$b} keys %$refs){
        my ($refStart,$refStop)=split /\./,$refPos;
        if( $query_stop<$refStart or $query_start>$refStop ){
        }else{
            $overlap++;
            last;
        }
    }
    return $overlap;
}

sub load_pointers {
	my $file=shift;
	my $hash=shift;
    my $geneIDs=shift;
    open INF,"<",$file or die $!;
    while(<INF>){
        chomp;
        my @inf=split /,/;
        my ($macID,$pointerSeq,$pointerSize,$pointerStart,$pointerStop)=@inf;
        $$hash{$macID}{"$pointerStart.$pointerStop"}=$pointerSeq;# if($pointerSize<=30);
    }
    close INF;
}

sub load_orthoClusters {
	my $file=shift;
	my $hash=shift;
	my $clusterID=1;
	open OrthoVenn,"<",$file or die $!;
	while(<OrthoVenn>){
		chomp;
		my @inf=split /\t/,;
		for my $geneID (@inf){
			if($geneID=~/^(.*?)\|.*?_(.*?)$/i){
				my $species=$1;
				my $gene_realID=$2;
                $gene_realID=$1 if($gene_realID=~/^(.*?)_fr/);
                $gene_realID=$Halteria_adjust{$gene_realID} if(exists $Halteria_adjust{$gene_realID});
                $species=$1 if($species=~/^(.*?)(_TSA)$/);
				my $clusterID_here="Clustr.$clusterID";
				$$hash{"cluster_content"}{$clusterID_here}{$species}{$gene_realID}++;
				$$hash{"gene_details"}{$species}{$gene_realID}=$clusterID_here;
			}
		}
		$clusterID++;
	}
	close OrthoVenn;
}

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
            $name=$1 if($name=~/^(.*?)(\t|\s+|\|)/);
			$name=$1 if($name=~/^(.*?)_fr/);
            $$hash{$name}=$seq;
        }
    }
    close Genome;
}

sub load_MSA_fasta {
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
            if($name=~/^(.*?)\|(.*?)$/){
                $$hash{$1}{$2}=$seq;
            }
        }
    }
    close Genome;
}


sub load_GFF {
    #my $type=($species =~/Hal/i) ? "locus_tag" : "Parent" ;
    #: ($species =~/Oxy|Hal/i) ? "gene" : ""; Hal=>locus_tag;else=>Parent
	my $file=shift;
	my $hash=shift;
	my $type=shift;
    my %exon;
	open GFF,"<",$file or die $!;
	while(<GFF>){
		chomp;
		my @inf=split /\t/,;
        if($inf[2]=~/cds/i and $inf[-1]=~/locus_tag=(FGO.*?);.*;protein_id=(TNV.*?);/){
            $Halteria_adjust{$2}=$1;
        }
        if($inf[2]=~/exon/i){
            $inf[-1].=";";
		    my $geneID=$1 if($inf[-1]=~/$type=(.*?)\;/);
            $geneID=$1 if($geneID=~/^(.*?).t\d+/);
		    $$hash{"exons_pos"}{$inf[0]}{$geneID}{"$inf[3].$inf[4]"}++;
            $$hash{"geneID2macID"}{$geneID}=$inf[0];
        }
	}
	close GFF;
}

sub load_mds_GFF {
    #my $type=($species =~/Hal/i) ? "locus_tag" : "Parent" ;
    #: ($species =~/Oxy|Hal/i) ? "gene" : ""; Hal=>locus_tag;else=>Parent
	my $file=shift;
	my $hash=shift;
	open GFF,"<",$file or die $!;
	while(<GFF>){
		chomp;
		my @inf=split /\t/,;
        if($inf[2]=~/mds/i){
            $$hash{$inf[0]}{"$inf[3].$inf[4]"}++;
        }
	}
	close GFF;
}

sub check_MDS_annotated {
    my $database=shift;
    my $annotated=0;
    while(my $check_site=shift @_){
        for my $databaseSection (sort {$a<=>$b} keys %{$database}){
            my ($databaseStart,$databaseStop)=split /\./,$databaseSection;
            if($databaseStart<=$check_site and $check_site<=$databaseStop){
                $annotated++;
                last;
            }
        }
        last if($annotated)
    } 
    return $annotated;
}
