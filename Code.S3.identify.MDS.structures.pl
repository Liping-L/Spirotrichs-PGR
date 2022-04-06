# Merge intergrated MDSs
# load options and prepare output path
  use Getopt::Long;
my $outPath="./Output";
my $MACfile;
my $MICfile;
my $telo_file;
my $input_mode;
my $inPut_File;
our $Pident_limit||=97;
our $Evalue_limit||=1E-10;
our $insertion_length_limit||=0.05; # percentage of whole MDS
our $spades_assemblies;
 GetOptions(
	'mac=s'  =>\$MACfile,
	'mic=s'  =>\$MICfile,
	'telo=s'	=>\$telo_file,
	'pi=s'  =>\$Pident_limit,
	'ev=s'  =>\$Evalue_limit,
	'il=s'  =>\$insertion_length_limit,
	'o=s' =>\$outPath,
	'im=s' => \$input_mode,
	'in=s' => \$inPut_File,
	'spades=s' => \$spades_assemblies,  
 );
 
 die "[Usage] perl $0 -mac mac.noTelo.fa -mic mic.fa -telo telomeres.inf.tsv -im [1|2]\n\t-im -- 1: separate files  2: merged file , accompanied with -in if '-im 2' is enabled\n" unless(-s $MACfile and -s $telo_file and -s $MICfile and ((-s $inPut_File and $input_mode==2) or $input_mode==1));
 chop $outPath if($outPath=~/\/$/);

# Check and create working directories
	mkdir $outPath  or die "Fail to creat the directory $outPath\n" unless(-e $outPath);
	mkdir "$outPath/Progess/"  or die "Fail to creat the directory 'Progess'\n" unless( -e "$outPath/Progess/");
	mkdir "$outPath/MAC.annotation/"  or die "Fail to creat the directory 'MAC'\n" unless(-e "$outPath/MAC.annotation/");
	mkdir "$outPath/MIC.annotation/"  or die "Fail to creat the directory 'MIC'\n" unless (-e "$outPath/MIC.annotation/");
	mkdir "$outPath/Updated_genome/"  or die "Fail to creat the directory 'Updated_genome'\n" unless (-e "$outPath/Updated_genome/");
	

# Creat output files
	open Error_log,">","$outPath/Progess/Process.err"  or die $!;
	open Typical_section,">","$outPath/Progess/middle.Typically_Merged_Records.csv" or die $!;
	
	

# Get basename
	use File::Basename;
	my @suffixlist=qw(.fa .fas .fasta);
	my ($MAC_basename, $MAC_inpath, $MAC_suffix)=fileparse($MACfile, @suffixlist);
	my ($MIC_basename, $MIC_inpath, $MIC_suffix)=fileparse($MICfile, @suffixlist);

# time recording
	my $time_start=time();
	local $|=1;

# load MAC genome and MIC genome 
	our %MAChash;
	our %MAClen;
	our %MAC_scaff_codes;
	our %MIChash;
	our %MIClen;
	our %MIC_scaff_codes;	
	&load_genome($MACfile,$MICfile) ;

# stat telomeres of each contig in both Genomes ================================================================================
	our %MAC_tel_5;
	our %MAC_tel_3;
	&stat_MAC_tels($telo_file);

# prepare hash
 our %step01_Blast_result;
 our %step01_MICbased_Blast_result;
 our %MAC_blast_coverage_statistics;
 our %direction_of_MIC;
 our %MICbased_elim1st_Blast_result;
 our %step01_HSP_list;
 our %MACbased_typical_Merge_record;
 our %MICbased_typical_Merge_record;
 our %typical_merged_MAC_MIC_combinations;
 our %step02_Tyically_Merged_HSP_List;
 our %MAChsp_TMR;

# load blast results
	if($input_mode==1){
		opendir MAChsp_dir,"./Blast_result/" or die "Fail to open ./Blast_result/\n" ;  
		my @hsp_files_counts=readdir MAChsp_dir;
		closedir MAChsp_dir;
		print "  Loading Blast Result:\t\t[ ";
		my $process=0;
		open Blast_result,">","$outPath/Progess/Blast.out6" or die $!;
		open Blast_Qualified_result,">","$outPath/Progess/Blast.qualified.out6" or die $!;
		for my $MAC_HSPfile (sort {$a cmp $b } @hsp_files_counts){
			$process++;
			print "=" if( $process%(@hsp_files_counts/20)==0);
			unless ($MAC_HSPfile=~/.*csv/){
				print "\n    Skip $MAC_HSPfile\n" unless ($MAC_HSPfile eq '.' or $MAC_HSPfile eq '..' )		;
				next;
			}
			my ($blastInf,$Qualified_blastInf)=load_Blast_out6($input_mode,"./Blast_result/$MAC_HSPfile",\%step01_Blast_result,\%step01_MICbased_Blast_result,\%direction_of_MIC,\%MAC_blast_coverage_statistics);
			print Blast_result $blastInf;
			print Blast_Qualified_result $Qualified_blastInf;
		}
		close Blast_Qualified_result;
		close Blast_result;
		print " ]  DONE !\n";
	}elsif($input_mode==2){
		print "  Loading Blast Result ...\n";
		#die "\nQualified File has already existed!\n" if(-s "$outPath/Progess/Blast.qualified.out6");		
		open Blast_Qualified_result,">","$outPath/Progess/Blast.qualified.out6" or die $!;	
		my ($tmp,$Qualified_blastInf)=load_Blast_out6($input_mode,$inPut_File,\%step01_Blast_result,\%step01_MICbased_Blast_result,\%direction_of_MIC,\%MAC_blast_coverage_statistics);
		print Blast_Qualified_result $Qualified_blastInf;
	}
	
	
	print "  Building HSP list:\t\t[ ";
	my $process=0;
	for my $MACname (sort keys %step01_Blast_result){
		$process++;
		print "=" if( $process%((keys %step01_Blast_result)/20)==0);
		for my $MICname	(sort keys %{$step01_Blast_result{$MACname}}){			
			# [c1-Generate HSP_list with Primary elimination (e.g. same HSPs and overlapped HSPs)]	
				for my $MACinfo( sort {${$step01_Blast_result{$MACname}{$MICname}}{$a}<=>${$step01_Blast_result{$MACname}{$MICname}}{$b} } keys %{$step01_Blast_result{$MACname}{$MICname}}){
					my @inf=split /,/,$MACinfo;
					my $MIClen=$MIClen{$MICname};
					my ($MACstart,$MACend,$MICstart,$MICend)=@inf[5,6,7,8];					
					anchor_1:
					my $dealted=0;
					# (1) eliminate same hits, may caused by 2-step blast
					if(exists ${$step01_HSP_list{$MACname}{$MICname}}{"$MACstart,$MACend,$MICstart,$MICend"}){
						next;
					}
					# (2) eliminate nested hits, may caused by 2-step blast
					for my $position(sort {$a cmp $b} keys %{$step01_HSP_list{$MACname}{$MICname}}){
						my ($covered_MACstart,$covered_MACend,$covered_MICstart,$covered_MICend)=split /,/,$position;
						if(($covered_MACstart<$MACstart and $MACend<=$covered_MACend) or ($covered_MACstart<=$MACstart and $MACend<$covered_MACend) ){
							&Delete_MAC_MIC_pos($MACname,$MICname,"$MACstart,$MACend,$MICstart,$MICend",\%step01_HSP_list);
							$dealted=1;
							last;
						}
						elsif( ($MACstart<$covered_MACstart and $MACend>=$covered_MACend) or ($MACstart<=$covered_MACstart and $MACend>$covered_MACend)   ){
							&Delete_MAC_MIC_pos($MACname,$MICname,"$covered_MACstart,$covered_MACend,$covered_MICstart,$covered_MICend",\%step01_HSP_list);
							$dealted=2;
							last;
						}
					}
					next if($dealted==1);
					goto anchor_1 if($dealted==2);
					# (3) add this one to the list if it's new
					($dealted==0)? ${$step01_HSP_list{$MACname}{$MICname}}{"$MACstart,$MACend,$MICstart,$MICend"}=$MACinfo : goto Check_again ;  #unnecessary but kept
				}
			# [c1-Complete : Generate step01_HSP_list]

			# [c2-typically merge]					
			for my $coordinate (sort by_MAC_hsp_coordinate keys %{$step01_HSP_list{$MACname}{$MICname}}){
				my $HSP_infos=$step01_HSP_list{$MACname}{$MICname}{$coordinate};
				my $MIClen=$MIClen{$MICname};
				my ($MACstart,$MACend,$MICstart,$MICend,$strand)=(split /,/,$HSP_infos)[5,6,7,8,12];
				$MICbased_elim1st_Blast_result{$MICname}{$step01_HSP_list{$MACname}{$MICname}{$coordinate}}="$MICstart.$MICend";																																						
				anchor_2:
				my $dealted=0;
				for my $covered_position(sort by_MAC_hsp_coordinate keys %{$step02_Tyically_Merged_HSP_List{$MACname}{$MICname}}){
					my $covered_HSP_infos=$step02_Tyically_Merged_HSP_List{$MACname}{$MICname}{$covered_position};
					my @inf=split /,/,$covered_HSP_infos;
					my ($covered_MACstart,$covered_MACend,$covered_MICstart,$covered_MICend,$covered_strand)=@inf[5,6,7,8,12]; #  covered_MACstart 					
					if(  &sections_separate_or_not($MICname,$MACname,$MICstart,$MICend,$MACstart,$MACend,$covered_MICstart,$covered_MICend,$covered_MACstart,$covered_MACend,\%MICbased_typical_Merge_record) eq "Y"  ){
							 	#       		5'___________________________|_______________________/_____|______________________/___________________________________3'
							 	#  MAC:															covered_MACstart			MACstart	covered_MACend				MACend
							 	#
							 	#					5'____________________|_______________________|__________________/_________________________/__________________________3'
							 	#  MIC:												covered_MICstart/end	MICend/start 	covered_MICstart/end			    MICend/start

					 	if($covered_MACstart<$MACstart and -1<=$covered_MACend-$MACstart and $covered_MACend<$MACend){
							my ($new_MICstart,$new_MICend,$new_Strand)=&get_NewStartAndEnd($MICstart,$MICend,$strand,$covered_MICstart,$covered_MICend,$covered_strand);
							unless($new_MICstart>0 and $new_MICend>0){
								my $error_inf="\n    [Error][The posotion is incomplete] adjusted MIC position are $new_MICstart -> $new_MICend\n";
								print $error_inf;
								print Error_log $error_inf;
								next ;
							}
							my $records="$MACname,$MICname,s1,$covered_MACstart,$covered_MACend,$covered_MICstart,$covered_MICend,$covered_strand,s2,$MACstart,$MACend,$MICstart,$MICend,$strand";
							if(exists ${$MACbased_typical_Merge_record{$MACname}{$MICname}}{"$covered_MACstart,$MACend,$new_MICstart,$new_MICend"}){
								print "Conflict mac typical merging: \nnow: ".$records."\n".${$MACbased_typical_Merge_record{$MACname}{$MICname}}{"$covered_MACstart,$MACend,$new_MICstart,$new_MICend"}."\n";
							}elsif(exists ${$MICbased_typical_Merge_record{$MICname}{$MACname}}{"$new_MICstart,$new_MICend,$covered_MACstart,$MACend"}){
								print "Conflict mic typical merging: \nnow: ".$records."\n".${$MICbased_typical_Merge_record{$MICname}{$MACname}}{"$new_MICstart,$new_MICend,$covered_MACstart,$MACend"}."\n";
							}							
							$MACbased_typical_Merge_record{$MACname}{$MICname}{"$covered_MACstart,$MACend,$new_MICstart,$new_MICend"}=$records;
							$MICbased_typical_Merge_record{$MICname}{$MACname}{"$new_MICstart,$new_MICend,$covered_MACstart,$MACend"}=$records;
							$MAChsp_TMR{$MACname}{$MICname}{"$covered_MACstart.$covered_MACend"}=$records;
							$MAChsp_TMR{$MACname}{$MICname}{"$MACstart.$MACend"}=$records;
							$typical_merged_MAC_MIC_combinations{"${MACname}_$MICname"}++;
							print Typical_section ("Merged,$MACname,$MICname,section1,$covered_MACstart,$covered_MACend,$covered_MICstart,$covered_MICend,$covered_strand,section2,$MACstart,$MACend,$MICstart,$MICend,$strand\n");
							&Delete_MAC_MIC_pos($MACname,$MICname,"$covered_MACstart,$covered_MACend,$covered_MICstart,$covered_MICend",\%step02_Tyically_Merged_HSP_List);
							&Delete_MAC_MIC_pos($MACname,$MICname,"$MACstart,$MACend,$MICstart,$MICend",\%step02_Tyically_Merged_HSP_List);

							@inf[3,5,6,7,8,12]=("typical_merge",$covered_MACstart,$MACend,$new_MICstart,$new_MICend,$new_Strand);
							$HSP_infos=join ",",@inf;
							($MACstart,$MACend,$MICstart,$MICend,$strand)=($covered_MACstart,$MACend,$new_MICstart,$new_MICend,$new_Strand);
							$dealted=1;
							last;
					 	}
					}
				}
				goto anchor_2 if($dealted==1);
				${$step02_Tyically_Merged_HSP_List{$MACname}{$MICname}}{"$MACstart,$MACend,$MICstart,$MICend"}=$HSP_infos if ($dealted==0);
			} 
			# [c2-Complete: typical merge]
		}

	}
print " ]  DONE !\n";
	

 # [c3-integrate the merged hits to create a MDSlist]
	$process=0;
	print "  Building MDS list:\t\t[ ";
	my $time_begin=time();
	my $time_here=0;
	for my $MACname( sort { $a cmp $b } keys %step02_Tyically_Merged_HSP_List){
		$process++;
		print "=" if( $process%((keys %step02_Tyically_Merged_HSP_List)/20)==0);		
		#for my $strand( sort { $a cmp $b } keys %{$step02_Tyically_Merged_HSP_List{$MACname}}){
			my %sortlist;
			for my $MICname(sort { $a cmp $b } keys %{$step02_Tyically_Merged_HSP_List{$MACname}}){
				for my $info (sort { $a cmp $b } values %{$step02_Tyically_Merged_HSP_List{$MACname}{$MICname}}){
					my @inf=split /,/,$info; #qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore qcovs (strand)
					unless(@inf==13){
						print "\n    [Error][The contents of %step02_Tyically_Merged_HSP_List is uncomplete or empty!]\n";
						print Error_log ("\n    [Error][The contents of %step02_Tyically_Merged_HSP_List is uncomplete or empty!]\n");
					}
					$sortlist{$info}=$inf[5];
				}
			}
			for my $info( sort { $sortlist{$a} <=> $sortlist{$b} } keys %sortlist  ){  
				my ($MICname,$MACstart,$MACend,$MICstart,$MICend,$bitScore,$strand)=(split /,/,$info)[1,5,6,7,8,10,12]; #qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore qcovs strand
				my $MIClen=$MIClen{$MICname};
				#qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore qcovs				
				anchor_3:
				my $dealted=0;
				goto add_new if(exists ${$MACbased_typical_Merge_record{$MACname}{$MICname}}{"$MACstart,$MACend,$MICstart,$MICend"});
				# (c3.1) merge the same or almostly same hits. only less than 1bp is accepted
				for my $existed_inf(sort {$step03_MDSlist{$MACname}{$a}<=>$step03_MDSlist{$MACname}{$b}} keys %{$step03_MDSlist{$MACname}}){
					my @existed_infos = split /,/,$existed_inf;
					next unless(@existed_infos==7);
					my ($existed_MICname,$existed_MACstart,$existed_MACend,$existed_MICstart,$existed_MICend,$existed_strand)=@existed_infos[1,2,3,4,5,6];					
					my $existed_Score=(split /,/,$step01_HSP_list{$MACname}{$existed_MICname}{"$existed_MACstart,$existed_MACend,$existed_MICstart,$existed_MICend"})[-2];
					if( abs($MACstart-$existed_MACstart)<=1 and abs($MACend-$existed_MACend)<=1 ){
						$dealted=&remove_low_confidence($MACname,$MICname,$MACstart,$MACend,$MICstart,$MICend,$strand,$bitScore,$existed_MICname,$existed_MACstart,$existed_MACend,$existed_MICstart,$existed_MICend,$existed_strand,$existed_Score);
						last;
					}
				}
				next if($dealted==1);
				goto anchor_3 if($dealted==3);
				# (c3.2) merge nested hits after typically merging
				for my $existed_inf(sort {$step03_MDSlist{$MACname}{$a}<=>$step03_MDSlist{$MACname}{$b}} keys %{$step03_MDSlist{$MACname}}){
					my @existed_infos = split /,/,$existed_inf;
					next unless(@existed_infos==7);
					my ($existed_MICname,$existed_MACstart,$existed_MACend,$existed_MICstart,$existed_MICend,$existed_strand)=@existed_infos[1,2,3,4,5,6];
					my $existed_Judgement_MACstart=( abs($MACstart-$existed_MACstart)<=1 ) ? $MACstart : $existed_MACstart;
					my $existed_Judgement_MACend=( abs($MACend-$existed_MACend)<=1 ) ? $MACend : $existed_MACend;
					next unless( ($MACstart<=$existed_Judgement_MACstart and $existed_Judgement_MACend<=$MACend) or ($existed_Judgement_MACstart<=$MACstart and $MACend<=$existed_Judgement_MACend) );
					if($MACend-$MACstart > $existed_MACend-$existed_MACstart and 
						(	!exists ${$MACbased_typical_Merge_record{$MACname}{$existed_MICname}}{"$existed_MACstart,$existed_MACend,$existed_MICstart,$existed_MICend"} or 
							 exists ${$MACbased_typical_Merge_record{$MACname}{$MICname}}{"$MACstart,$MACend,$MICstart,$MICend"}														
						)
					){
						&Delete_MAC_pos($MACname,$existed_inf,\%step03_MDSlist);
						$dealted=3;
					}elsif($MACend-$MACstart < $existed_MACend-$existed_MACstart and 
							(!exists ${$MACbased_typical_Merge_record{$MACname}{$MICname}}{"$MACstart,$MACend,$MICstart,$MICend"} or
							  exists ${$MACbased_typical_Merge_record{$MACname}{$existed_MICname}}{"$existed_MACstart,$existed_MACend,$existed_MICstart,$existed_MICend"}
							)
						  ){
						&Delete_MAC_pos($MACname,"$MACname,$MICname,$MACstart,$MACend,$MICstart,$MICend,$strand",\%step03_MDSlist);
						$dealted=1;
					}
					last;
				}
				next if($dealted==1);
				goto anchor_3 if($dealted==3);
				# (c3.3) add a new one
				add_new:
				$step03_MDSlist{$MACname}{"$MACname,$MICname,$MACstart,$MACend,$MICstart,$MICend,$strand"}="$MACstart.$MACend" if($dealted==0);
			}
	}
	print " ]  DONE !\n";
# [c3-Complete: integrate the merged hits to create a MDSlist]

# [c4- split the typically merged hits]
	for my $MACname(sort {$a cmp $b} keys %step03_MDSlist){
			my $dealted=0;
			for my $MACinfo(sort {$step03_MDSlist{$MACname}{$a} <=> $step03_MDSlist{$MACname}{$b}} keys %{$step03_MDSlist{$MACname}}){
				my %contents;
				my $MACpos=$step03_MDSlist{$MACname}{$MACinfo};
				my ($tmpMACname,$MICname,$MACstart,$MACend,$MICstart,$MICend,$strand)=split /,/,$MACinfo;
				my $key_of_typical_Merge_record="$MACstart,$MACend,$MICstart,$MICend";
				if(exists ${$MACbased_typical_Merge_record{$MACname}{$MICname}}{$key_of_typical_Merge_record}){
					extract_raw_inf_of_Typically_Merged(\%MACbased_typical_Merge_record,$MACname,$MICname,$key_of_typical_Merge_record,\%contents);
					Delete_MAC_pos($MACname,$MACinfo,\%step03_MDSlist);
				}
				for my $pos (sort {$a<=>$b} keys %contents){
					my ($MACstart,$MACend,$MICstart,$MICend)=split /,/,$pos;
					my $info=$contents{$pos};
					$step03_MDSlist{$MACname}{$info}="$MACstart.$MACend";										
				}
			}			
	}
# [c4-Complete: split the typical merged hits]


# transport the data that based on MACname to a new hash  which is based on MICname
	for my $MACname(sort {$a cmp $b} keys %step03_MDSlist){
		for my $MACinfo(sort {$step03_MDSlist{$MACname}{$a} <=> $step03_MDSlist{$MACname}{$b}} keys %{$step03_MDSlist{$MACname}}){
			my ($tmpMACname,$MICname,$MACstart,$MACend,$MICstart,$MICend,$strand)=split /,/,$MACinfo;
			chomp $MICend;			
			if($MICstart<=0 or $MICend<=0){
				print "\nError: line304: $MACname,$MICname,$MACstart,$MACend,$MICstart,$MICend,$strand\n"
			}
			$step05_MIC_MDSlist{$MICname}{"$MACname,$MICname,$MACstart,$MACend,$MICstart,$MICend,$strand"}="$MICstart.$MICend";
		}
	}

# integrate the MACs that mapped to the same reagion of MIC under some specific conditions. These steps are aimed at reducing the disturbance of MAC segments
for my $MICname(sort {$a cmp $b} keys %step05_MIC_MDSlist){
	next if($needless_realMICs{$MICname}>0);
	for my $MICinfo (sort {$step05_MIC_MDSlist{$MICname}{$a} <=> $step05_MIC_MDSlist{$MICname}{$b}} keys %{$step05_MIC_MDSlist{$MICname}}){
		my ($MACname,$MICtmp,$MACstart,$MACend,$MICstart,$MICend,$strand)=split /,/,$MICinfo;
		Check_again_c5_5:
		my $keep_the_behind_one=1;
		for my $MICinfo_front (sort {$step06_MIC_MDSlist{$MICname}{$a} <=> $step06_MIC_MDSlist{$MICname}{$b}} keys %{$step06_MIC_MDSlist{$MICname}}){
			my ($MACname_front,$MICtmp_front,$MACstart_front,$MACend_front,$MICstart_front,$MICend_front,$strand_front)=split /,/,$MICinfo_front;
			next if($MAC_tel_5{$MACname}>0 or $MAC_tel_3{$MACname}>0);
			my $MAC_len=$MAClen{$MACname};
			my $MAC_front_len=$MAClen{$MACname_front};
			if( $MICstart_front<$MICstart and $MICend<=$MICend_front and $MACstart<=2 and abs($MACend-$MAC_len)<=1 ){ # keep the front one
				$keep_the_behind_one=0;
				last
			}elsif($MICstart_front==$MICstart and $MICend_front<=$MICend and abs($MACend_front-$MAC_front_len)<=1 ){ # keep the behind one
				delete ${$step06_MIC_MDSlist{$MICname}}{$MICinfo_front} if(exists ${$step06_MIC_MDSlist{$MICname}}{$MICinfo_front} ); 
				$keep_the_behind_one=2;
				last;
			}
		}
		$step06_MIC_MDSlist{$MICname}{$MICinfo}=$step05_MIC_MDSlist{$MICname}{$MICinfo} if($keep_the_behind_one==1);
		goto Check_again_c5_5 if($keep_the_behind_one==2);
	}	
}

# transport the infos to a MAC-based hash
for my $MICname(sort {$a cmp $b} keys %step06_MIC_MDSlist){
	next if($needless_realMICs{$MICname}>0);
	for my $MICinfo (sort {$step06_MIC_MDSlist{$MICname}{$a} <=> $step06_MIC_MDSlist{$MICname}{$b}} keys %{$step06_MIC_MDSlist{$MICname}}){
		my ($MACname,$MICname,$MACstart,$MACend,$MICstart,$MICend,$strand)=split /,/,$MICinfo;
		$step07_MDSlist{$MACname}{$MICname}{"$MACstart,$MACend,$MICstart,$MICend,$strand"}="$MACstart.$MACend";
	}
}

my $process=0;
print "  Final Checking :\t\t[ ";
for my $MACname (sort {$a cmp $b} keys %step07_MDSlist){
	$process++;
	print "=" if( $process%((keys %step07_MDSlist)/20)==0);
	# (c6) connect the HSPs which may be splited unexpectedly
   	for my $MICname (sort {$a cmp $b} keys %{$step07_MDSlist{$MACname}}){
 		for my $MACinfo (sort {$step07_MDSlist{$MACname}{$MICname}{$a} <=> $step07_MDSlist{$MACname}{$MICname}{$b}} keys %{$step07_MDSlist{$MACname}{$MICname}}){
			my ($MACstart,$MACend,$MICstart,$MICend,$strand)=split /,/,$MACinfo;
			Check_again_c6:
			my $dealted=0;
			for my $covered_info (sort {$step08_MDSlist{$MACname}{$a} <=> $step08_MDSlist{$MACname}{$b}} keys %{$step08_MDSlist{$MACname}}){
				my ($covered_MICname,$covered_MACstart,$covered_MACend,$covered_MICstart,$covered_MICend,$covered_strand)=split /,/,$covered_info;
				next if($covered_strand ne $strand or $covered_MICname ne $MICname ) ;
			    my $limit_MAC=abs($MACend-$covered_MACstart+1)*$insertion_length_limit;
			    my $limit_MIC=abs($MICend-$covered_MICstart+1)*$insertion_length_limit;
			    if(  		$covered_MACstart<=$MACstart and $MACstart<=$covered_MACend and $covered_MACend-$MACstart<=$limit_MAC and $covered_MACend<=$MACend
			    		and $covered_MICstart<=$MICstart and $MICstart<=$covered_MICend and $covered_MICend-$MICstart<=$limit_MIC and $covered_MICend<=$MICend
						and !exists ${$MAChsp_TMR{$MACname}{$MICname}}{"$MACstart.$MACend"} and !exists ${$MAChsp_TMR{$MACname}{$MICname}}{"$covered_MACstart.$covered_MACend"}
			    ){
						&Delete_MAC_pos($MACname,"$MICname,$MACinfo",\%step08_MDSlist);
						&Delete_MAC_pos($MACname,$covered_info,\%step08_MDSlist);
						($MACstart,$MACend,$MICstart,$MICend)=($covered_MACstart,$MACend,$covered_MICstart,$MICend);
						$dealted=2;
						last;
			    }
			}
			goto Check_again_c6 if($dealted==2);
			$step08_MDSlist{$MACname}{"$MICname,$MACstart,$MACend,$MICstart,$MICend,$strand"}="$MACstart.$MACend" if($dealted==0);
		}
	}
}

for my $MACname (sort {$a cmp $b} keys %step08_MDSlist){
	# (c7) Skip the 77bp HSPs locat at the end of both MAC and MIC
	for my $MACinfo (sort {$step08_MDSlist{$MACname}{$a} <=> $step08_MDSlist{$MACname}{$b}} keys %{$step08_MDSlist{$MACname}}){
		my ($MICname,$MACstart,$MACend,$MICstart,$MICend,$strand)=split /,/,$MACinfo;
		my $MIClen=$MIClen{$MICname};
		my $MAClen=$MAClen{$MACname};
		if($spades_assemblies=~/mac/i and $spades_assemblies=~/mic/i ){
			if( abs($MACstart-$MACend)+1==77 and  abs($MICstart-$MICend)+1==77 
					and ( ( ($MACstart==1 and $MAC_tel_5{$MACname}==0) or ($MACend==$MAClen and $MAC_tel_3{$MACname}==0))
					or ( $MICstart==1 or $MICend==$MIClen) )
			){
				next
			}
		}
		$step09_MDSlist{$MACname}{"$MICname,$MACstart,$MACend,$MICstart,$MICend,$strand"}="$MACstart.$MACend";
	}
}

for my $MACname (sort {$a cmp $b} keys %step09_MDSlist){
	# discard the nested MDSs on MAC
	my @infos=sort {$step09_MDSlist{$MACname}{$a} <=> $step09_MDSlist{$MACname}{$b}} keys %{$step09_MDSlist{$MACname}};
	my %wait2delete_MACmds;
		for my $inf (@infos){
			my ($MICname,$MACstart,$MACend,$MICstart,$MICend,$strand)=split /,/,$inf;
			for my $inf2 (@infos){
				next if(exists $wait2delete_MACmds{$inf2});
				my ($existed_MICname,$asRef_MACstart,$asRef_MACend,$asRef_MICstart,$asRef_MICend,$asRef_strand)=split /,/,$inf2;
				if($asRef_MACstart<=$MACstart and $MACend<=$asRef_MACend and $step09_MDSlist{$MACname}{$inf2}!=$step09_MDSlist{$MACname}{$inf} and 
					!(exists ${$MAChsp_TMR{$MACname}{$MICname}}{"$MACstart.$MACend"} and !exists ${$MAChsp_TMR{$MACname}{$existed_MICname}}{"$asRef_MACstart.$asRef_MACend"} )
				){
					$wait2delete_MACmds{$inf}++;
					$inf=~s/,/\t/g;
					$inf2=~s/,/\t/g;
					my $output_inf="#[hypothetical-Nested]\nmac:$MACname\n";
					$output_inf.="mic:$MICname is involved in Typically-Merging\n"  if(exists ${$MAChsp_TMR{$MACname}{$MICname}}{"$MACstart.$MACend"});
					$output_inf.="mic:$existed_MICname is involved in Typically-Merging\n"  if(exists ${$MAChsp_TMR{$MACname}{$existed_MICname}}{"$asRef_MACstart.$asRef_MACend"});
					$output_inf.="inf1:$inf\ninf2:$inf2\n";
					print Error_log $output_inf;
				}
			}
		}
	for my $inf (keys %wait2delete_MACmds){ delete ${$step09_MDSlist{$MACname}}{$inf} }
	
	# discard the insert MDSs between MAC MDSs
	my $dealted=1;
	while($dealted==0){
		$dealted=0;
		my @infos=sort {$step09_MDSlist{$MACname}{$a} <=> $step09_MDSlist{$MACname}{$b}} keys %{$step09_MDSlist{$MACname}};
		for($i=0;$i<@infos;$i++){
			my ($MICname,$MACstart,$MACend,$MICstart,$MICend,$strand)=split /,/,$infos[$i];
			for($k=$i+2;$k<@infos;$k++){
				my ($check_MICname,$check_MACstart,$check_MACend,$check_MICstart,$check_MICend,$check_strand)=split /,/,$infos[$k];
				my $insert_MDSs=$infos[($k-1)];
				my ($insert_MIC,$insert_MACstart,$insert_MACend)=(split /,/,$insert_MDSs)[0,1,2];
				if($check_MACstart<=$MACend and !exists ${$MAChsp_TMR{$MACname}{$insert_MIC}}{"$insert_MACstart.$insert_MACend"}){																					
					my $records="#[hypothetical-Insertion]\nmac:$MACname\nThis MDSs: $insert_MDSs\nis An Insertion Between $infos[$i] \nand $infos[$k]\n";
					print Error_log $records;				
					delete ${$step09_MDSlist{$MACname}}{$infos[($k-1)]} ;
					$dealted=1;					
					last;
				}elsif($check_MACstart>$MACend){
					last;
				}
			}
		}
	}
		
}
print " ]  DONE !\n";

	# 
	for my $MACname (sort {$a cmp $b} keys %step09_MDSlist){
		for my $MACinfo (sort {$step09_MDSlist{$MACname}{$a} <=> $step09_MDSlist{$MACname}{$b}} keys %{$step09_MDSlist{$MACname}}){
			my ($MICname,$MACstart,$MACend,$MICstart,$MICend,$strand)=split /,/,$MACinfo;
			$final_MIC_MDSlist{$MICname}{"$MACname,$MACstart,$MACend,$MICstart,$MICend,$strand"}="$MICstart.$MICend";
			$Preparation_for_IESidentify{"$MICname,$MACname"}{"$MACname,$MACstart,$MACend,$MICstart,$MICend,$strand"}="$MICstart.$MICend";
		}
	}
	
	# Identify IESs from the SAME MAC-MIC pairs
	for my $ctgname (sort keys %Preparation_for_IESidentify){
		my ($MICname,$MACname)=split /,/,$ctgname;		
		my $last_MICstart=0;
		my $last_MICstop=0;
		my $last_MACstart=0;
		my $last_MACend=0;
		for my $MICinfo (sort {$Preparation_for_IESidentify{$ctgname}{$a} <=> $Preparation_for_IESidentify{$ctgname}{$b}} keys %{$Preparation_for_IESidentify{$ctgname}}){			
			my ($MACname,$MACstart,$MACend,$MICstart,$MICend,$strand)=split /,/,$MICinfo;			
			if($last_MICstart>0){
				my $IESstart=$last_MICend+1;
				my $IESstop=$MICstart-1;
				next if($IESstart>$IESstop);
				if( $IESstart<=$IESstop and ($last_MACstart<=$MACstart and $MACstart<=$last_MACend and $last_MACend<=$MACend ) or ($MACstart<=$last_MACstart and $last_MACstart<=$MACend and $MACend<=$last_MACend ) ){
					$IESsource{$MICname}{"$IESstart.$IESstop"}="Homo";
					$MIC_IESs{$MICname}{$MICinfo}="$IESstart.$IESstop";
				}
			}	
			($last_MICstart,$last_MICend,$last_MACstart,$last_MACend)=($MICstart,$MICend,$MACstart,$MACend);
		}
	}
	# Identify IESs according MIC MDSs coord
	for my $MICname (sort keys %final_MIC_MDSlist){
		my $last_MICstart=0;
		my $last_MICend=0;		
		my $last_MACname;
		for my $MICinfo (sort {$final_MIC_MDSlist{$MICname}{$a} <=> $final_MIC_MDSlist{$MICname}{$b}} keys %{$final_MIC_MDSlist{$MICname}}){
			my ($MACname,$MACstart,$MACend,$MICstart,$MICend,$strand)=split /,/,$MICinfo;
			if($last_MICstart>0){
				my $IESstart=$last_MICend+1;
				my $IESstop=$MICstart-1;	
				if($IESstart<=$IESstop){
					$IESsource{$MICname}{"$IESstart.$IESstop"}="Homo-like" if($MACname eq $last_MACname and !exists $IESsource{$MICname}{"$IESstart.$IESstop"});
					$MIC_IESs{$MICname}{$MICinfo}="$IESstart.$IESstop" if(!exists $IESsource{$MICname}{"$IESstart.$IESstop"});
				}
			}
			($last_MACname,$last_MICstart,$last_MICend,$last_MACstart,$last_MACend)=($MACname,$MICstart,$MICend,$MACstart,$MACend);
		}
	}


# =============================================================================================================================================================================================== # 
#																					Output Results																								  #
# =============================================================================================================================================================================================== # 

 # Output Final Results
	print "  Output Final result:\t\t[ ";
	my $total_outFiles=9;
	my $part=int(20/$total_outFiles);
	my $last_part=20-$part*($total_outFiles-1);

	# file1
	outPut_MDSlist_style1("$outPath/Progess/middle.step03_MDSlist","tsv",\%step03_MDSlist,1);
	print "=" for(1..$part);
	
	# file2
	outPut_MDSlist_style2("$outPath/Progess/middle.step07_MDSlist","tsv",\%step07_MDSlist,2);
	print "=" for(1..$part);

	# file3
	outPut_MDSlist_style1("$outPath/Progess/middle.step08_MDSlist","tsv",\%step08_MDSlist,3);
	print "=" for(1..$part);

	# file4
	outPut_MDSlist_style1("$outPath/Progess/middle.step09_MDSlist","tsv",\%step09_MDSlist,3);
	my $mac_out_base="$outPath/MAC.annotation/mac.mds.structures";
	open MACannotationTAB,">","$mac_out_base.tab" or die "Fail to open mac.mds.structures.tab!\n";
	open MACannotationCSV,">","$mac_out_base.csv" or die "Fail to open mac.mds.structures.csv!\n";
	open MACannotationGFF,">","$mac_out_base.gff3" or die "Fail to open mac.mds.structures.gff3!\n";
	my $head_inf="MACname,MICname,type,MACstart,MACstop,MICstart,MICstop,strand\n";
	print MACannotationCSV $head_inf;
	$head_inf=~s/,/\t/g;
	print MACannotationTAB $head_inf;
	$head_inf="#gff-version3\n";
	print MACannotationGFF $head_inf;
	for my $MACname(sort {$a cmp $b} keys %step09_MDSlist){
			#my $MAClen=$MAClen{$MACname};			
			my $circle_time=1;
			my $MACscaff_code=$MAC_scaff_codes{$MACname};
			my $Tel5_len=$MAC_tel_5{$MACname};
			my $Tel3_len=$MAC_tel_3{$MACname};
			my $MAClen=$MAClen{$MACname}+$Tel5_len+$Tel3_len;
			my $tel3_start=$MAClen-$Tel3_len;	
			my $head=&GFFformat($MACname,"Perl","region","1",$MAClen,".","+",".","mac.$MACscaff_code","");
			print MACannotationGFF $head;			 		
			if ($Tel5_len>0){
				my $outinf="$MACname,NA,5tel,1,$Tel5_len,0,0,+\n";
				print MACannotationCSV $outinf;
				$outinf=~s/,/\t/g;
				print MACannotationTAB $outinf;
				$outinf=&GFFformat($MACname,"Perl","tel","1",$Tel5_len,".","+",".","tel.$MACscaff_code.5","");
				print MACannotationGFF $outinf;
			}
			my $tel5_start=$Tel5_len+1 if($Tel5_len>0);
			my $out_CSV;
			my $out_GFF;
			for my $MACinfo(sort {$step09_MDSlist{$MACname}{$a} <=> $step09_MDSlist{$MACname}{$b}} keys %{$step09_MDSlist{$MACname}}){
				my ($MICname,$MACstart,$MACend,$MICstart,$MICend,$strand)=split /,/,$MACinfo;
				   $MACstart+=$Tel5_len;
				   $MACend+=$Tel5_len;
				my $OUT_start;
				my $OUT_stop;
				if($Tel5_len>0 and $MACstart<$Tel5_len and $Tel3_len>0 and $MACend>$tel3_start){							
					$OUT_start=$tel5_start;
					$OUT_stop=$tel3_start;
				}elsif($Tel5_len>0 and $MACstart<$Tel5_len){
					$OUT_start=$tel5_start;
					$OUT_stop=$MACend;
				}elsif($Tel3_len>0 and $MACend>$tel3_start){
					$OUT_start=$MACstart;
					$OUT_stop=$tel3_start;
				}else{
					$OUT_start=$MACstart;
					$OUT_stop=$MACend;
				}
				$out_v4csv.="$MACname,mds_$circle_time,$strand,$OUT_start,$OUT_stop\n";
				$out_CSV.="$MACname,$MICname,mds_$circle_time,$OUT_start,$OUT_stop,$MICstart,$MICend,$strand,mds.$MACscaff_code.$circle_time\n";
				$out_GFF.=&GFFformat($MACname,"Perl","mds",$OUT_start,$OUT_stop,".","$strand",".","mds.$MACscaff_code.$circle_time","Source=$MICname|${MICstart}_$MICend;");								
				$circle_time++; 		
			} 
			print MACannotationCSV $out_CSV;
			$out_CSV=~s/,/\t/g;
			print MACannotationTAB $out_CSV;
			print MACannotationGFF $out_GFF;
			$tel3_start++; 				
			if ($Tel3_len>0){		
				my $outinf="$MACname,NA,3tel,$tel3_start,$MAClen,0,0,+\n";
				print MACannotationCSV $outinf;
				$outinf=~s/,/\t/g;
				print MACannotationTAB $outinf;
				$outinf=&GFFformat($MACname,"Perl","tel",$tel3_start,$MAClen,".","+",".","tel.$MACscaff_code.3","");
				print MACannotationGFF $outinf;
			}
	}
	close MACannotationGFF;
	close MACannotationCSV;
	close MACannotationTAB;
	intergrate_GFF("$mac_out_base.gff3","$mac_out_base.intergrated.gff3");
	print "=" for(1..$part);

	open HSP_lists,">","$outPath/Progess/middle.step01_HSP_list.tsv"  or die $!;
		for my $MACname( sort { $a cmp $b } keys %step01_HSP_list){
			for my $MICname( sort { $a cmp $b } keys %{$step01_HSP_list{$MACname}}){
				for my $coordinate (sort by_MAC_hsp_coordinate keys %{$step01_HSP_list{$MACname}{$MICname}}){
					my $HSP_infos=$step01_HSP_list{$MACname}{$MICname}{$coordinate};
					$HSP_infos=~s/,/\t/g;																						  
					print HSP_lists $HSP_infos."\n";																																											
				}
			}
		}
	print "=" for(1..$part);

	open Tyically_Merged_HSP_List_MACbased,">","$outPath/Progess/middle.step02_Tyically_Merged_HSP_List.MACbased.tsv" or die $!;
		for my $MACname( sort { $a cmp $b } keys %MACbased_typical_Merge_record){
			for my $MICname( sort { $a cmp $b } keys %{$MACbased_typical_Merge_record{$MACname}}){
				for my $coordinate (sort by_MAC_hsp_coordinate keys %{$MACbased_typical_Merge_record{$MACname}{$MICname}}){
					my $HSP_infos=$MACbased_typical_Merge_record{$MACname}{$MICname}{$coordinate};
					$HSP_infos=~s/,/\t/g;
					print Tyically_Merged_HSP_List_MACbased $HSP_infos."\n";																																											
				}
			}	
		}
	open Tyically_merged_HSP_list,">","$outPath/Progess/middle.step02_Tyically_Merged_HSP_List.tsv"  or die $!;		
		for my $MACname( sort { $a cmp $b } keys %step02_Tyically_Merged_HSP_List){
			for my $MICname( sort { $a cmp $b } keys %{$step02_Tyically_Merged_HSP_List{$MACname}}){
				for my $coordinate (sort by_MAC_hsp_coordinate keys %{$step02_Tyically_Merged_HSP_List{$MACname}{$MICname}}){
					my $HSP_infos=$step02_Tyically_Merged_HSP_List{$MACname}{$MICname}{$coordinate};
					$HSP_infos=~s/,/\t/g;
					print Tyically_merged_HSP_list $HSP_infos."\n";																																											
				}
			}	
		}
	print "=" for(1..$part);

	open step05_MIC_MDSlist,">","$outPath/Progess/middle.step05_MIC_MDSlist" or die  $!;
		print step05_MIC_MDSlist "#IESs in this file are just predicted and maybe inaccurate\nMACname,type,strand,start,end\n";
		#my ($MACname,$MICname,$MACstart,$MACend,$MICstart,$MICend,$strand)=split /,/,$step05_MIC_MDSlist{$MICname}{$MICpos};
		for my $MICname(sort { $a cmp $b } keys %step05_MIC_MDSlist){
			my $MIClen=$MIClen{$MICname};
			my $last_MICend=0;
			my $mds_order=1;
			my $ies_order=1;
			for my $MICinfo (sort {$step05_MIC_MDSlist{$MICname}{$a} <=> $step05_MIC_MDSlist{$MICname}{$b}} keys %{$step05_MIC_MDSlist{$MICname}}){
				my ($MACname,$MICname,$MACstart,$MACend,$MICstart,$MICend,$strand)=split /,/,$MICinfo;
				if($last_MICend>0 and $last_MICend<$MICstart ){
			 		my $attribude=( $last_MACname eq $MACname ) ? "Homo" : "Para" ; # Hic means High-Confidence, Prd means Predicited
					my $ies_start=$last_MICend+1;
					my $ies_stop=$MICstart-1;
					print step05_MIC_MDSlist "$MICname,${attribude}ies$ies_order,$strand,$ies_start,$ies_stop\n";
					$ies_order++;
				}
				print step05_MIC_MDSlist "$MICname,mds$mds_order,$strand,$MICstart,$MICend\n";
				$mds_order++;
				$last_MICend=$MICend;
			 	$last_MACname=$MACname;
		 	}
	 	}
	print "=" for(1..$part);


	#  MICname final_MIC_MDSlist
	 	open MIC_MDSlist,">","$outPath/Progess/middle.final_MIC_MDSlist.tsv" or die $!;
		my $out_base="$outPath/MIC.annotation/mic.mds.structures";
		open MICannotationTAB,">","$out_base.tab" or die $!;
		open MICannotationGFF,">","$out_base.gff3" or die $!;
		open MICannotationCSV,">","$out_base.csv" or die $!;
		my $out="MICname,type,strand,start,end\n";
		$out=~s/,/\t/g;	 	 
		print MIC_MDSlist $out;	 	
	 	$out="MACname\tMICname\ttype\tMACstart\tMACstop\tMICstart\tMICstop\tstrand\n";
		print MICannotationTAB $out;	 	
		$out=~s/\t/,/g;
	 	print MICannotationCSV $out;	 	
	 	print MICannotationGFF "#gff-version3\n";
		for my $MICname(sort { $a cmp $b } keys %final_MIC_MDSlist){
			my $MIClen=$MIClen{$MICname};
		 	my $mds_order=1;
		 	my $ies_order=1;
		 	my $MICscaff_code=$MIC_scaff_codes{$MICname};
			my $head=&GFFformat($MICname,"Perl","region","1",$MIClen,".","+",".","mic.$MICscaff_code","");
			print MICannotationGFF $head;
		 	for my $MICinfo(sort {$final_MIC_MDSlist{$MICname}{$a} <=> $final_MIC_MDSlist{$MICname}{$b}} keys %{$final_MIC_MDSlist{$MICname}}){
			 	my ($MACname,$MACstart,$MACend,$MICstart,$MICend,$strand)=split /,/,$MICinfo;	 
				my $MACtelo5_len=$MAC_tel_5{$MACname};
				$MACstart+=$MACtelo5_len;
				$MACend+=$MACtelo5_len;	
			 	my $out="$MICname,mds$mds_order,$strand,$MICstart,$MICend\n";
				$out=~s/,/\t/g;
				print MIC_MDSlist $out;
				$out="$MACname\t$MICname\tmds$mds_order\t$MACstart\t$MACend\t$MICstart\t$MICend\t$strand\tmds.$MICscaff_code.$mds_order\n";				
			 	print MICannotationTAB  $out;
				$out=~s/\t/,/g;
			 	print MICannotationCSV $out;
			 	my $out_GFF=&GFFformat($MICname,"Perl","mds",$MICstart,$MICend,".","$strand",".","mds.$MICscaff_code.$mds_order","Source=$MACname|${MACstart}_$MACend;Type=NA;");
			 	print MICannotationGFF "$out_GFF";
			 	$mds_order++;
		 	}
		 	my $last_MACname="";
		 	for my $MICinfo(sort {$MIC_IESs{$MICname}{$a} <=> $MIC_IESs{$MICname}{$b}} keys %{$MIC_IESs{$MICname}}){			 
				my ($MACname,$MACstart,$MACend,$MICstart,$MICend,$strand)=split /,/,$MICinfo;	
				my $MICpos=$MIC_IESs{$MICname}{$MICinfo};
				my ($IESstart,$IESend)=split /\./,$MICpos;
				my $attribude= (exists ${$IESsource{$MICname}}{$MICpos}) ? $IESsource{$MICname}{$MICpos} : "Para" ;			
				my $out="$MICname,ies$ies_order,$strand,$MICstart,$MICend\n";
				$out=~s/,/\t/g;
				print MIC_MDSlist $out;
				$out="NA\t$MICname\t${attribude}ies$ies_order\tNA\tNA\t$IESstart\t$IESend\tNA\ties.$MICscaff_code.$ies_order";
				print MICannotationTAB "$out\n";
				$out=~s/\t/,/g;
				print MICannotationCSV "$out\n";
				my $out_GFF=&GFFformat($MICname,"Perl","ies",$IESstart,$IESend,".","$strand",".","ies.$MICscaff_code.$ies_order","Source=$MACname|mds;Type=$attribude;");
				print MICannotationGFF "$out_GFF";
				$ies_order++;	
				$last_MACname=$MACname;	
		 	}						
	 	}
		 close MICannotationCSV;
		 close MICannotationTAB;
		 close MICannotationGFF;
	intergrate_GFF("$out_base.gff3","$out_base.intergrated.gff3");
	print "=" for(1..$part);

	


	# Get the infomation of MIC-algnments	
		open MICdetails,">","$outPath/Progess/middle.MICbased_elim1st_Blast_result.tsv" or die $!;
		for my $MICname(sort {$a cmp $b} keys  %MICbased_elim1st_Blast_result){
		 	for my $MICinfo(sort {$MICbased_elim1st_Blast_result{$MICname}{$a}<=>$MICbased_elim1st_Blast_result{$MICname}{$b}} keys %{$MICbased_elim1st_Blast_result{$MICname}}){
				$MICinfo=~s/,/\t/g;
				print MICdetails "$MICinfo\n";
		 	}
		}
		close MICdetails;
	print "=" for(1..$last_part);


	print " ]  DONE !\n";
	$time_end=time();
	$timespending=$time_end-$time_start; #sprintf ("%.1f",($time_end-$time_start)/60);
	print "  All steps cost $timespending s \n";





	sub intergrate_GFF {
		my ($rawGFF,$intergratedGFF)=@_;
		my %all_records;
		open INF,"<",$rawGFF or die $!;
		open OUT,">",$intergratedGFF or die $!;
		while( my $record=<INF>){
			chomp $record;
			my @inf=split /\t/,$record;
			$all_records{$inf[0]}{"$inf[3].$inf[4]"}{"$inf[2]"}{$record}++
		}
		for my $ctgname (sort keys %all_records){
			for my $pos (sort {$a<=>$b} keys %{$all_records{$ctgname}}){
				for my $type (sort keys %{$all_records{$ctgname}{$pos}}){
					my @rec=keys %{$all_records{$ctgname}{$pos}{$type}};
					my @Sources_Types;
					my @Sources;
					#my @IDs;
					if( @rec>1 ){
						for my $rec (@rec){
							#push @IDs,$1 if($rec=~/ID=\w*?\.(\d+\.\d+)/);
							if($rec =~/\t([+-])\t.*?Source=(.*?)\;Type=(.*?)\;/){
								push @Sources_Types,"$2$1\;$3";
							}elsif($rec =~/\t([+-])\t.*?Source=(.*?)\;/){
								push @Sources,"$2$1";
							}
						}
						#my $finalID=max(@IDs);
						my @inf=split /\t/,$rec[0];
						my $ID=$1 if($inf[-1]=~/ID=(.*?)\;/);
						my $i=0;
						$inf[-1]="ID=$ID\;";
						if(@Sources_Types){
							for my $S_T (@Sources_Types){
								$i++;
								$inf[-1].="Source.$i=$1\;Type.$i=$2\;" if($S_T=~/^(.*?)\;(.*?)$/);
							}
						}else{
							for my $S (@Sources){
								$i++;
								$inf[-1].="Source.$i=$S\;";
							}
						}
						my $final=join "\t",@inf;
						print OUT $final."\n";
					}else{
						$rec[0]=~s/Source/Source.1/;
						$rec[0]=~s/Type/Type.1/;
						print OUT $rec[0]."\n";
					}
				}
			}
		}


	}
	sub GFFformat {
		my ($MACname,$source,$type,$start,$end,$score,$strand,$phase,$attribude_ID,$otherINF)=@_;
		my $GFF_format="$MACname\t$source\t$type\t$start\t$end\t$score\t$strand\t$phase\tID=$attribude_ID;$otherINF\n";
		return $GFF_format;
	}

	sub get_reverse_complementry_seq {
		my $Seq=shift;
			 $Seq=tr/ATCG/TAGC/;
			 $Seq=reverse $Seq;
		return $Seq;
	}

	sub min {
			my $Min = shift;
			for (@_) {
					 $Min = $_ if $_ < $Min;
			}
			return $Min;
	}

	sub max {
			my $Max = shift;
			for (@_) {
				 $Max = $_ if $_ > $Max;
			}
			return $Max;
	}

	sub average{
		my $aver=0;
		my $count=0;
		while(my $aver=shift){
			$total+=$aver;
			$count++;
		}
		$aver=$total/$count;
		$aver=sprintf ("%.1f",$aver)
	}



	sub Delete_MAC_pos {
		my $MAC=shift;
		my $position=shift;
		my $hash=shift;
		delete ${$$hash{$MAC}}{$position} if ( defined ${$$hash{$MAC}}{$position} or exists ${$$hash{$MAC}}{$position} );
		delete $$hash{$MAC} if (keys %{$$hash{$MAC}}<1);
	}

	sub Delete_MAC_MIC_pos {
		my $MAC=shift;
		my $MIC=shift;
		my $position=shift;
		my $hash=shift;
		delete ${$$hash{$MAC}{$MIC}}{$position} if ( defined ${$$hash{$MAC}{$MIC}}{$position} or exists ${$$hash{$MAC}{$MIC}}{$position} );
		delete ${$$hash{$MAC}}{$MIC} if (keys %{$$hash{$MAC}{$MIC}}<1);
	}


sub sections_separate_or_not  {
			my ($MICname,$MACname,$section1_MICstart,$section1_MICend,$section1_MACstart,$section1_MACend,$section2_MICstart,$section2_MICend,$section2_MACstart,$section2_MACend,$MACorMIC_based_typical_Merge_record)=@_;

			my %contents1;
			my %contents2;
			my $key_of_typical_Merge_record1="$section1_MICstart,$section1_MICend,$section1_MACstart,$section1_MACend";
			my $key_of_typical_Merge_record2="$section2_MICstart,$section2_MICend,$section2_MACstart,$section2_MACend";
			extract_raw_inf_of_Typically_Merged(\%$MACorMIC_based_typical_Merge_record,$MICname,$MACname,$key_of_typical_Merge_record1,\%contents1) if(exists ${$$MACorMIC_based_typical_Merge_record{$MICname}{$MACname}}{$key_of_typical_Merge_record1});
			extract_raw_inf_of_Typically_Merged(\%$MACorMIC_based_typical_Merge_record,$MICname,$MACname,$key_of_typical_Merge_record2,\%contents2) if(exists ${$$MACorMIC_based_typical_Merge_record{$MICname}{$MACname}}{$key_of_typical_Merge_record2});

			my @querying_Positions;
			if((keys %contents1)>0){
				for my $pos (sort {$a<=>$b} keys %contents1){
					my $qPos="$3.$4" if($pos=~/^(\d+),(\d+),(\d+),(\d+)$/);
					push @querying_Positions,$qPos;
				}
			}else{
				push @querying_Positions,"$section1_MICstart.$section1_MICend"
			}
			if((keys %contents2)>0){
				for my $pos (sort {$a<=>$b} keys %contents2){
					my $qPos="$3.$4" if($pos=~/^(\d+),(\d+),(\d+),(\d+)$/);
					push @querying_Positions,$qPos;
				}
			}else{
				push @querying_Positions,"$section2_MICend.$section2_MICend"
			}

			my $result='Y';
			for(my $i=1;$i<@querying_Positions;$i++){
				my $Previous_start=min(split /\./,$querying_Positions[$i-1]);
				my $Previous_stop=max(split /\./,$querying_Positions[$i-1]);
				my $Behind_start=min(split /\./,$querying_Positions[$i]);
				my $Behind_stop=max(split /\./,$querying_Positions[$i]);
				if( $Previous_start>=$Behind_stop or $Behind_start>=$Previous_stop ){
				}else{
					$result='N';last;
				}				
			}
			return $result;
}

sub extract_raw_inf_of_Typically_Merged {
	my $typically_merged_record=shift;
	my $MACorMIC_name=shift;
	my $corespond_MICorMAC_name=shift;
	my $Pos_4containers=shift;
	my $result=shift;
	my $info=$$typically_merged_record{$MACorMIC_name}{$corespond_MICorMAC_name}{$Pos_4containers};
	unless($info=~/^(.*?),(.*?),s1,(.*?),(.*?),(.*?),(.*?),(.*?),s2,(.*?),(.*?),(.*?),(.*?),(.*)$/){
		my $error_inf="\n[Error] Empty Merge_record! MACname(MICname):$MACorMIC_name\tMICname(MACname):$corespond_MICorMAC_name\tposition:$Pos_4containers\n";
		print $error_inf;
		print Error_log $error_inf;
	}
	$$result{"$3,$4,$5,$6"}="$1,$2,$3,$4,$5,$6,$7";
   	$$result{"$8,$9,$10,$11"}="$1,$2,$8,$9,$10,$11,$12";
	    
	my $dealted=1;
    while($dealted==1){
		$dealted=0;
		for my $Check_pos_4containers (sort {$a cmp $b} keys %{$result}){
			my $info=$$typically_merged_record{$MACorMIC_name}{$corespond_MICorMAC_name}{$Check_pos_4containers};
        	if($info=~/^(.*?),(.*?),s1,(.*?),(.*?),(.*?),(.*?),(.*?),s2,(.*?),(.*?),(.*?),(.*?),(.*)$/){			
				delete $$result{$Check_pos_4containers};
				$$result{"$3,$4,$5,$6"}="$1,$2,$3,$4,$5,$6,$7";
				$$result{"$8,$9,$10,$11"}="$1,$2,$8,$9,$10,$11,$12";
				$dealted++				            
        	}
    	}   
	} 
}

sub get_NewStartAndEnd {
		my $section1_start=shift;
		my $section1_end=shift;
		my $section1_strand=shift;
		my $section2_start=shift;
		my $section2_end=shift;
		my $section2_strand=shift;
		my @NewStartAndEnd;
		$NewStartAndEnd[0]=&min($section1_start,$section1_end,$section2_start,$section2_end);
		$NewStartAndEnd[1]=&max($section1_start,$section1_end,$section2_start,$section2_end);
		my $new_strand="";
		if( $section1_strand eq '-' and $section2_strand eq '-' ){
			$new_strand='-'
		}
		elsif( $section1_strand eq '+' and $section2_strand eq  '+'){
			$new_strand='+'
		}
		elsif( $section1_strand ne  $section2_strand ){
			if($direction_of_MIC{$MACname}{$MICname}{"+"}>$direction_of_MIC{$MACname}{$MICname}{"-"}){
				$new_strand='+'
			}
			else{
				$new_strand='-'
			}
		}
		else{
			my $error_inf= "[Error][sub get_NewStartAndEnd: the position of two sections seems wrong!][Some details for debugging: $section1_start,$section1_end\t|\t$section2_start,$section2_end]\n";
			print $error_inf;
			print Error_log $error_inf;
		}
		$NewStartAndEnd[2]=$new_strand;
		return @NewStartAndEnd;
}

sub by_MAC_hsp_coordinate {
	# format of input: "start1,end1,start2,end2"
	my @a_site=split /,/,$a;
	my @b_site=split /,/,$b;
	my $a_pos1="$a_site[0].$a_site[1]";
	my $b_pos1="$b_site[0].$b_site[1]";
	my $a_pos2="$a_site[2].$a_site[3]";
	my $b_pos2="$b_site[2].$b_site[3]";

	if($a_pos1 < $b_pos1){
		return -1
	}elsif($a_pos1 > $b_pos1){
		return 1
	}elsif( $a_pos2 < $b_pos2){
		return -1
	}elsif( $a_pos2 > $b_pos2){
		return 1		
	}else{
		return -1
	}
}


sub by_HSPs_order {
		$a=~/(.*?),(.*?),(.*?,.*?,.*?),(.*?),(.*?),(.*?),(.*?),(.*?,.*?,.*)$/;
		my ($MACstart_a,$MACend_a,$MICstart_a,$MICend_a)=($4,$5,$6,$7);
		my $MIClen_a=$MIClen{$2};
		$b=~/(.*?),(.*?),(.*?,.*?,.*?),(.*?),(.*?),(.*?),(.*?),(.*?,.*?,.*)$/;
		my ($MACstart_b,$MACend_b,$MICstart_b,$MICend_b)=($4,$5,$6,$7);
		my $MIClen_b=$MIClen{$2};
		#print "a is $a\tb is $b\n";
		if($MACstart_a+1<$MACstart_b){
			-1
		}
		elsif($MACstart_a>1+$MACstart_b){
			1
		}
		elsif($MACend_a<$MACend_b){
			-1
		}
		elsif($MACend_a>$MACend_b){
			1
		}
		elsif($MACstart_a<$MACstart_b){
			-1
		}
		elsif($MACstart_a>$MACstart_b){
			1
		}

}

sub by_len {
	((length $MAChash{$a}) < (length $MAChash{$b})) ? return 1 :  return -1
}

sub load_genome {
	my $MAC_file=shift;
	my $MIC_file=shift;
	local $|=1;
	print "\n  Build a hash for MAC and MIC\t[ ";
	load_fasta($MAC_file,\%MAChash,\%MAClen);
	my $not_standard=0;
	for my $MACname (keys %MAChash){
		unless($MACname=~/NODE_(\d+)_/){
			$not_standard++;
			last;
		}
	}
	if($not_standard==0){
		for my $MACname (keys %MAChash){
			$MACname=~/NODE_(\d+)_/;
			$MAC_scaff_codes{$MACname}=$1;
		}
	}else{
		my $scaff_code=1;
		for my $MACname (sort by_len keys %MAChash){
			$MAC_scaff_codes{$MACname}=$scaff_code;
			$scaff_code++;
		}
	}
	print "==========";

	load_fasta($MIC_file,\%MIChash,\%MIClen);
	my $not_standard=0;
	for my $MICname (keys %MIChash){
		unless($MICname=~/NODE_(\d+)_/){
			$not_standard++;
			last;
		}
	}
	if($not_standard==0){
		for my $MICname (keys %MIChash){
			$MICname=~/NODE_(\d+)_/;
			$MIC_scaff_codes{$MICname}=$1;
		}
	}else{
		my $scaff_code=1;
		for my $MICname (sort by_len keys %MIChash){
			$MIC_scaff_codes{$MICname}=$scaff_code;
			$scaff_code++;
		}
	}
	print "========== ]  DONE !\n";
	my $MACnum=keys %MAChash;
	my $MICnum=keys %MIChash;
	print "\tLoaded $MACnum MACs and $MICnum MICs from genome file\n";
	close MAC;
	close MIC;
}

sub load_fasta {
	my $file=shift;
	my $hash_seq=shift;
	my $hash_len=shift;
	open MAC,"<",$file or die "Fail to load $file!\n";
	{
		local $/=undef;
		my $allSEQs=<MAC>;
		    $allSEQs=~s/\n/!/g;
		my @all=split />/,$allSEQs;
		for my $info(@all){
			next unless($info=~/^(.*?)!(.*)$/);
			my $name=$1;
			my $seq=$2;
			   $seq=~s/!//g;
			$name=$1 if($name=~/^(.*?)(\s|\t)/);
			$$hash_seq{$name}=$seq;
			$$hash_len{$name}=length $seq;
		}
	}
}

sub load_Blast_out6{
	# mode
	my $mode=shift; # 1: separate file 2: merged file
	# files 
	my $infile=shift;
	# hash
	my $Blast_result=shift;
	my $MICbased_Blast_result=shift;
	my $direction_of_sseqID=shift;
	my $blast_coverage_statistics=shift;
	# out inf
	my $blast_out_inf;
	my $qualified_Blast_inf;

	open MAChsp,"<",$infile or die "Fail to load $infile\n";		
	while(my $MACinfo=<MAChsp>){
		#qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore qcovs (strand)
		chomp $MACinfo;
		my @inf=split /,/,$MACinfo;
		next unless(@inf==12);
		my $MACname=$inf[0];
		my $MICname=$inf[1];
		$blast_out_inf.= "$MACinfo\n" if($mode==1);
		next if($inf[2]<$Pident_limit or $inf[9]>$Evalue_limit);
		$$blast_coverage_statistics{$inf[0]}{$inf[11]}++;
		my $strand= ($inf[7]<$inf[8]) ? "+" : "-" ;
		($inf[7],$inf[8])=($inf[8],$inf[7]) if($strand eq '-');
		push @inf,$strand;
		$MACinfo= join ",",@inf;			
		$$direction_of_sseqID{$MACname}{$MICname}{$strand}++;
		$$Blast_result{$MACname}{$MICname}{$MACinfo}="$inf[5].$inf[6]";
		$qualified_Blast_inf.="$MACinfo\n";
		$$MICbased_Blast_result{$MICname}{$MACname}{$MACinfo}="$inf[7].$inf[8]";
	}
	close MAChsp;
	return ($blast_out_inf,$qualified_Blast_inf);
}



sub stat_MAC_tels {
		my $MAC_telo_file=shift;
		local $|=1;
		print "  Mark the telomeres of MACs\t[ ";
		my $i=0;
		my %Tels_stat;
		my $total_lines=$1 if(`wc -l $MAC_telo_file`=~/^(\d+)/);
		open Telo_inf,"<","$MAC_telo_file" or die "Fail to open MAC_Telomeres.inf.tsv\n";		
		while(my $a=<Telo_inf>){			
			$i++;
			print "=" if($total_lines>0 and $i%($total_lines/20)==0);
			chomp $a;
			my ($MACname,$telnum,$tel5len,$tel3len)=split /\t/,$a;
			$MAC_tel_5{$MACname}=$tel5len if($tel5len>0 );
			$MAC_tel_3{$MACname}=$tel3len if($tel3len>0 );
			my $tel_type=($telnum==2) ? 8 :
						($telnum==1 and $tel5len>0) ? 5 :
						($telnum==1 and $tel3len>0) ? 3 : 0 ;
			$Tels_stat{$tel_type}++;
		}
		print " ]  DONE !\t";
		close Telo_inf;
		my ($five_end,$three_end,$both_end,$none)=($Tels_stat{5},$Tels_stat{3},$Tels_stat{8},$Tels_stat{0});
		my $five_part=sprintf("%.1f",($five_end/(keys %MAChash))*100);
		my $three_part=sprintf("%.1f",($three_end/(keys %MAChash))*100);
		my $both_part=sprintf("%.1f",($both_end/(keys %MAChash))*100);
		my $none_part=sprintf("%.1f",($none/(keys %MAChash))*100);
		print "  \n\tSummary:\n\t  both-Tels are ${both_end}[${both_part}%]\t".
						"5-Tels are ${five_end}[${five_part}%]\t".
						"3-Tels are ${three_end}[${three_part}%]\t".
						"0-Tels are ${none}[${none_part}%]\n";
}




sub by_old_MICs_order {
		$sortlist{$a}=~/^(.*?),(\d+),(\d+),(\d+),(\d+)$/;
		my ($MACstart_a,$MACend_a,$MICstart_a,$MICend_a)=($2,$3,$4,$5);
		my $MIClen_a=$MIClen{$1};
		$sortlist{$b}=~/(.*?),(\d+),(\d+),(\d+),(\d+)$/;
		my ($MACstart_b,$MACend_b,$MICstart_b,$MICend_b)=($2,$3,$4,$5);
		my $MIClen_b=$MIClen{$1};	
		if($MACstart_a<$MACstart_b){
			-1
		}
		elsif($MACstart_a>$MACstart_b){
			1
		}
		elsif( ($MICend_a<=2 or $MIClen_a-$MICend_a<=1) and ($MICstart_b<=2 or $MIClen_b-$MICstart_b<=1) ){
			-1
		}
		elsif( ($MICend_b<=2 or $MIClen_b-$MICend_b<=1) and ($MICstart_a<=2 or $MIClen_a-$MICstart_a<=1) ){
			1
		}
		else{
			print "[Error-sort_MIC_order]\n";
			print Error_log ("[Error-sort_MIC_order]\n");
		}
}

sub check_pos_at_end {
		my $end3_pos=shift;
		while(my $pos2check=shift){
			return "Y" if( abs($pos2check-$end3_pos)<=1 or $pos2check<=2);
		}
}

sub fix_minusStrand_position {
  	my $position=shift;
  	my $seq_len=shift;
  	my $new_pos=$seq_len-$position+1;
  	return $new_pos;
}

sub remove_low_confidence {
		my ($MACname,$MICname,$MACstart,$MACend,$MICstart,$MICend,$strand,$bitScore,$existed_MICname,$existed_MACstart,$existed_MACend,$existed_MICstart,$existed_MICend,$existed_strand,$existed_Score)=@_;
		
		my $key_of_existed_in_MACbased_typical_Merge_record="$existed_MACstart,$existed_MACend,$existed_MICstart,$existed_MICend";
		my $key_of_new_in_MACbased_typical_Merge_record="$MACstart,$MACend,$MICstart,$MICend";
		
		my $key_of_new_in_step03_MDSlist="$MACname,$MICname,$MACstart,$MACend,$MICstart,$MICend,$strand";
		my $key_of_existed_in_step03_MDSlist="$MACname,$existed_MICname,$existed_MACstart,$existed_MACend,$existed_MICstart,$existed_MICend,$existed_strand";
		
		my $dealted=0;
		if(exists ${$MACbased_typical_Merge_record{$MACname}{$MICname}}{$key_of_new_in_MACbased_typical_Merge_record} and exists ${$MACbased_typical_Merge_record{$MACname}}{$key_of_existed_in_MACbased_typical_Merge_record}){
			$dealted=0;
		}elsif(exists ${$MACbased_typical_Merge_record{$MACname}{$MICname}}{$key_of_new_in_MACbased_typical_Merge_record} and !exists ${$MACbased_typical_Merge_record{$MACname}}{$key_of_existed_in_MACbased_typical_Merge_record}){
			&Delete_MAC_pos($MACname,$key_of_existed_in_step03_MDSlist,\%step03_MDSlist);
			$dealted=3;
		}
		elsif(!exists ${$MACbased_typical_Merge_record{$MACname}{$MICname}}{$key_of_new_in_MACbased_typical_Merge_record} and exists ${$MACbased_typical_Merge_record{$MACname}}{$key_of_existed_in_MACbased_typical_Merge_record}){
			&Delete_MAC_pos($MACname,$key_of_new_in_step03_MDSlist,\%step03_MDSlist);
			$dealted=1;
		}
		elsif( keys %{$step01_HSP_list{$MACname}{$MICname}} >  keys %{$step01_HSP_list{$MACname}{$existed_MICname}} ){
			&Delete_MAC_pos($MACname,$key_of_existed_in_step03_MDSlist,\%step03_MDSlist);
			$dealted=3;
		}
		elsif(keys %{$step01_HSP_list{$MACname}{$MICname}} < keys %{$step01_HSP_list{$MACname}{$existed_MICname}} ){
			&Delete_MAC_pos($MACname,$key_of_new_in_step03_MDSlist,\%step03_MDSlist);
			$dealted=1;
		}
		elsif( $bitScore>$existed_Score ){
			&Delete_MAC_pos($MACname,$key_of_existed_in_step03_MDSlist,\%step03_MDSlist);
			$dealted=3;
		}
		elsif( $bitScore<$existed_Score ){
			&Delete_MAC_pos($MACname,$key_of_new_in_step03_MDSlist,\%step03_MDSlist);
			$dealted=1;
		}
		elsif( $MIClen{$MICname}>=$MIClen{$existed_MICname} ){
			&Delete_MAC_pos($MACname,$key_of_existed_in_step03_MDSlist,\%step03_MDSlist);
			$dealted=3;
		}
		elsif( $MIClen{$MICname}<$MIClen{$existed_MICname} ){
			&Delete_MAC_pos($MACname,$key_of_new_in_step03_MDSlist,\%step03_MDSlist);
			$dealted=1;
		}
		return $dealted;
}

sub check_site_at_the_end_of_MICcontig {
	my $check_site=shift;
	my $MIC_len=shift;
	my $target_end=shift; # 5 or 3
	my $oligo_length=shift;
	if($target_end==5 and $check_site<=2){
		$$oligo_length++ if($check_site==1);
		return "Y"
	}elsif($target_end==3 and $MIC_len-$check_site<=1){
		$$oligo_length++ if(($MIC_len-$check_site)==1);
		return "Y"
	}else{
		return "N"
	}
}


sub outPut_MDSlist_style1 {
	my $outFile=shift;
	my $outfmt=shift;
	my $outHash=shift;
	my $data_structure=shift;
	$outFile.=($outfmt eq "tsv") ? ".tsv" : ".csv";
	my $head_inf="MACname,type,strand,start,end";
	$head_inf=~s/,/\t/g if($outfmt eq "tsv" );
	open MDSlist,">",$outFile or die "Fail to open $outFile!\n";
	print MDSlist $head_inf."\n";
	my $out;
	for my $MACname(sort {$a cmp $b} keys %{$outHash}){
		my $circle_time=1;
		my $Tel5_len=$MAC_tel_5{$MACname};
		my $Tel3_len=$MAC_tel_3{$MACname};
		my $MAClen=$MAClen{$MACname}+$Tel3_len+$Tel5_len;
		my $tel3_start=$MAClen-$Tel3_len;	
		$out.="$MACname,5tel,+,$Tel5_len"."\n" if ($Tel5_len>0);
		for my $MACinfo (sort {$$outHash{$MACname}{$a} <=> $$outHash{$MACname}{$b}} keys %{$$outHash{$MACname}}){			
			my ($MICname,$MACstart,$MACend,$MICstart,$MICend,$strand);
			if($data_structure ==1 ){ ($MICname,$MACstart,$MACend,$MICstart,$MICend,$strand)=(split /,/,$MACinfo)[1,2,3,4,5,6]	}
			elsif($data_structure ==2){ ($MACstart,$MACend,$MICstart,$MICend,$strand)=(split /,/,$MACinfo) }
			elsif($data_structure ==3){ ($MICname,$MACstart,$MACend,$MICstart,$MICend,$strand)=(split /,/,$MACinfo) }
			$MACstart+=$Tel5_len;
			$MACend+=$Tel5_len;
			$out.=($Tel5_len>0 and $MACstart<$Tel5_len and $Tel3_len>0 and $MACend>$tel3_start) ? "$MACname,mds_$circle_time,$strand,$tel5_start,$tel3_start" :
				  ($Tel5_len>0 and $MACstart<$Tel5_len) ? "$MACname,mds_$circle_time,$strand,$tel5_start,$MACend" :
				  ($Tel3_len>0 and $MACend>$tel3_start) ? "$MACname,mds_$circle_time,$strand,$MACstart,$tel3_start" :
														"$MACname,mds_$circle_time,$strand,$MACstart,$MACend";
			$out.="\n";
			$circle_time++; 
		}
	}
	$tel3_start++;
	$out.="$MACname,3tel,+,$tel3_start,$MAClen\n" if($Tel3_len>0);
	$out=~s/,/\t/g if($outfmt eq "tsv" );
	print MDSlist $out;
	close MDSlist;
}

sub outPut_MDSlist_style2 {
	my $outFile=shift;
	my $outfmt=shift;
	my $outHash=shift;
	my $data_structure=shift;
	$outFile.=($outfmt eq "tsv") ? ".tsv" : ".csv";
	my $head_inf="MACname,type,strand,start,end";
	$head_inf=~s/,/\t/g if($outfmt eq "tsv" );
	open MDSlist,">",$outFile or die "Fail to open $outFile!\n";
	print MDSlist $head_inf."\n";
	my $out;
	for my $MACname(sort {$a cmp $b} keys %{$outHash}){
		my $circle_time=1;
		my $Tel5_len=$MAC_tel_5{$MACname};
		my $Tel3_len=$MAC_tel_3{$MACname};
		my $MAClen=$MAClen{$MACname}+$Tel5_len+$Tel3_len;
		my $tel3_start=$MAClen-$Tel3_len;	
		$out.="$MACname,5tel,+,$Tel5_len"."\n" if ($Tel5_len>0);
		for my $MICname(sort {$a cmp $b} keys %{$$outHash{$MACname}}){
			for my $MACinfo (sort {$$outHash{$MACname}{$MICname}{$a} <=> $$outHash{$MACname}{$MICname}{$b}} keys %{$$outHash{$MACname}{$MICname}}){			
				my ($MICname,$MACstart,$MACend,$MICstart,$MICend,$strand);
				if($data_structure ==1 ){ ($MICname,$MACstart,$MACend,$MICstart,$MICend,$strand)=(split /,/,$MACinfo)[1,2,3,4,5,6]	}
				elsif($data_structure ==2){ ($MACstart,$MACend,$MICstart,$MICend,$strand)=(split /,/,$MACinfo) }
				elsif($data_structure ==3){ ($MICname,$MACstart,$MACend,$MICstart,$MICend,$strand)=(split /,/,$MACinfo) }
				$MACstart+=$Tel5_len;
				$MACend+=$Tel5_len;
				$out.=($Tel5_len>0 and $MACstart<$Tel5_len and $Tel3_len>0 and $MACend>$tel3_start) ? "$MACname,mds_$circle_time,$strand,$tel5_start,$tel3_start" :
					  ($Tel5_len>0 and $MACstart<$Tel5_len) ? "$MACname,mds_$circle_time,$strand,$tel5_start,$MACend" :
					  ($Tel3_len>0 and $MACend>$tel3_start) ? "$MACname,mds_$circle_time,$strand,$MACstart,$tel3_start" :
															"$MACname,mds_$circle_time,$strand,$MACstart,$MACend";
				$out.="\n";
				$circle_time++; 
			}
		}
	}
	$tel3_start++;
	$out.="$MACname,3tel,+,$tel3_start,$MAClen\n" if($Tel3_len>0);
	$out=~s/,/\t/g if($outfmt eq "tsv" );
	print MDSlist $out;
	close MDSlist;
}

