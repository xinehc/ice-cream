#!/usr/bin/perl 
# Meng Liu and Hong-Yu Ou on April-24-2018
# School of Life Sciences & Biotechnology, Shanghai Jiao Tong University
# modified by Yin Xiaole on October 5th, 2020
# University of Hong Kong

use strict;
# use Bio::Perl;
use Bio::SeqIO;
use Getopt::Std;
our ($opt_h, $opt_i, $opt_t, $opt_o, $opt_a) = "";
my $usage = q(
This script is used to detect the signature sequence about integrative and/or conjugative module 
of ICE/IME based on HMM profile and BLASTp.

Usage:	perl script_name.pl <job_id> 

);
die($usage) if ( @ARGV <= 1 ); 
getopts('i:t:o:a:h');

#Declare variables
my $job_id = $opt_a;
my $tmp_path = "$opt_t/$job_id/";
my $seq_fna = $tmp_path."$job_id.fna";
my $seq_faa = $tmp_path."$job_id.gi.faa";
my $seq_ptt = $tmp_path."$job_id.ptt";
my $seq_RNA = $tmp_path."$job_id.rnt";## rnt file containing the tRNA info converted from tmRNA_aragorn.out

#working dir
my $candidate_dir = $tmp_path."candidate/";
system "rm -rf $candidate_dir";
system "mkdir $candidate_dir";
system "chmod 777 $candidate_dir";
my $download_dir = "$opt_o/$job_id/";
system "rm -rf $download_dir";
system "mkdir $download_dir";
system "chmod 777 $download_dir";

my $error_log = $tmp_path."run_region_finder.log";
my $script_dir = "./scripts/";

#database
my $db_dir = "./data/";
my $int_hmm_db = $db_dir."Int.hmm.db";###transposase and integrase hmm
my $ice_hmm_db = $db_dir."ICE.hmm.db";### ICE singature protein hmm
my $threshold_file = $db_dir."threshold";###JUN definded parameters

#cmd
my $hmmscan_cmd = "hmmscan";
my $analyze_cmd = $script_dir."region_analyzer.pl";
#my $ARG_VF_cmd = "./scripts/ARG_VF_scanner.pl";



my $genome_name = "";
my $genome_size = 0;
my $hmm_evalue = 0.0001;
my $rep_tag = 0;
my $dr_max_mismatch = 0;
$| = 1;



#############################################################
#get threshold
#############################################################
my %threshold_evalue = ();##hmmer evalue
my %threshold_cvalue = ();##hmmer cvalue
my %dom_function = ();
open (THRE, $threshold_file) || die "open $threshold_file error: $!\n";
### Example of $threshold_file
#Domain name	E-value	c-Evalue	Function
#Phage_integrase	0.000001	0.0001	Int
#Recombinase	0.000001	0.000001	Int
#rve	0.000001	0.000001	Int
while(<THRE>){
	chomp;
	my @array = split(/\t/, $_);
	next unless $array[0];
	next if $array[0] eq "Domain name";
	$array[0] =~ s/ //g;
	$threshold_evalue{$array[0]} = $array[1];
	$threshold_cvalue{$array[0]} = $array[2];
	$dom_function{$array[0]} = $array[3]; #name is key; function is value
}
close THRE;

##########################################################
##submit integrase hmmscan task and read files.
##########################################################
##check the ptt file
if(! -e $seq_ptt){
	print "ERROR!For $job_id, no protein was found: $seq_ptt file was not found!\n";
	open ERROR, ">>$error_log";
	print ERROR "ERROR!For $job_id, no protein was found: $seq_ptt file was not found!\n";
	close ERROR;
	exit;
}
##check the faa file
if(! -e $seq_faa){
	print "ERROR! For $job_id, no protein was found: $seq_faa file was not found!\n";
	open ERROR, ">>$error_log";
	print ERROR "ERROR!For $job_id, no protein was found: $seq_faa file was not found!\n";
	close ERROR;
	exit;
}

######read the fna file.#####################################
if(! -e $seq_fna){
	print "ERROR!For $job_id, $seq_fna file was not found!\n";
	open ERROR, ">>$error_log";
	print ERROR "ERROR!For $job_id, $seq_fna file was not found!\n";
	close ERROR;
	exit;
}
my $genome_seqio_obj = Bio::SeqIO->new(-file => "$seq_fna", -format => "fasta" );
my $genome_seq_obj = $genome_seqio_obj->next_seq;
my $genome_id = $genome_seq_obj->display_id;
my $genome_desc = $genome_seq_obj->desc;
my $genome_seq = $genome_seq_obj->seq;
my $genome_gc = &calculate_gc($genome_seq); 
my $genome_len = $genome_seq_obj->length;
$genome_size = $genome_len;


#variables
my $left_region_range = 250000;
my $right_region_range = 250000; 
if (-e "$opt_t/$job_id/$job_id.gram"){
	$left_region_range = 42000;
	$right_region_range = 42000;
}

open LOG, ">>$error_log"; 
print LOG "For $job_id:left_region_range:$left_region_range\n";
close LOG;
####submit the integrase hmm task###################
print "Start detecting candidate regions containing integrase flanked by tRNA....\n";

my $int_hmm_out = $tmp_path."int_hmm.out";
if(-e $int_hmm_out){
	system "rm -f $int_hmm_out";
}
my $int_hmm_cmd = "$hmmscan_cmd --domtblout $int_hmm_out --cpu 20 -E $hmm_evalue $int_hmm_db $seq_faa 1> /dev/null 2>>$error_log";
system($int_hmm_cmd);

########read the ptt file#######################
open PTT, "$seq_ptt";
my $protein_num = 0;
my %ptt_gi_line = ();	#gi as key, line as value
my %ptt_gi_coor = ();	#gi as key, left coordinate as value
my %ptt_coor_gi = ();	#left coordinate as key, gi as value
my %ptt_coordinate = (); #left coordinate as key, right as value


my $line1= readline(PTT);
chomp($line1);
$line1 =~ /^(\w+.*\.)\s\-\s(\d+)\.\.(\d+)/; 
$genome_name = $1;


while(<PTT>){
	chomp;
	if(/^(\d+)\.\.(\d+)\t.\t.+\t(\d+)\t.+\t.+\t.+\t.+\t.+/){	#$1 left, $2 right, $3 gi
		$protein_num ++;
		next if $1 > $2;
		$ptt_gi_line{$3} = $_;
		$ptt_gi_coor{$3} = $1;
		$ptt_coor_gi{$1} = $3;
		$ptt_coordinate{$1} = $2;
		}
	}


if($protein_num == 0){
	print "ERROR!For $job_id, no protein was found: No protein was recorded in $seq_ptt file!\n";
	open ERROR, ">>$error_log";
	print ERROR "ERROR!For $job_id, no protein was found: No protein was recorded in $seq_ptt file!\n";
	close ERROR;
	exit;
}
#######read the faa file. generate a hash with gi as key and fasta as value######
my %pro = ();	#gi as key and fasta as value
my $seqio_obj_aa = Bio::SeqIO->new(-file => "$seq_faa", -format => "fasta" );
while ( my $seq = $seqio_obj_aa->next_seq() ) {
	my $dis = $seq->display_id;
	my $des = $seq->desc;
	my $string = $seq->seq;
	if(!$string){
		$string = "";
	}
	my $ref = "";
	if($dis =~ /gi\|(\d+)\|/){	#>gi|\d+|ref|... ...
		$ref = $1;
	}elsif($dis =~ /^(\d+)\_/){	#>\d+_1 ...
		$ref = $1;
	}
	$pro{$ref} = ">$dis $des\n";
	my $k=0;
	while($k<length($string)) {
		$pro{$ref} .= substr($string,$k,70)."\n";
		$k += 70;
	}
}



###########################################################################
##submit the tRNA detection task and read files
##convert the output into rnt file (ptt like format of RNA)
##############################################################################
	if(!-e "$tmp_path/$job_id.rnt"){
	my $tmrna = $tmp_path."/tmRNA_aragorn.out";
#### example of tmRNA_aragorn.out
#>CP003200 CP003200.1 Klebsiella pneumoniae subsp. pneumoniae HS11286, complete genome.
#88 genes found
#1   tRNA-Glu                 [17805,17880]	35  	(ttc)
#2   tRNA-Ile               [122234,122310]	35  	(gat)
#3   tRNA-Ala               [122347,122422]	34  	(tgc)
#4   tRNA-Asp               [125762,125838]	35  	(gtc)
	system("aragorn  -o $tmrna $tmp_path/$job_id.fna");
  open TMRNA, "$tmrna";
  open NEWRNT, "> $seq_RNA";
	my %rna = ();
	my %rna_name = ();
	my $file_rnt="";
	my $trna_num=0; 
	while(<TMRNA>){
	if(/Sequence \[(\d+),(\d+)\]/){
	my $length = $2 - $1 + 1;
	$rna{$1} = $2;
	$rna_name{$1} = "$1..$2\t+\t$length\ttRNA\ttRNA\t-\t-\t-\ttRNA";
	$file_rnt .= $rna_name{$1};
	$file_rnt .= "\n";
	$trna_num++;
	}elsif(/Sequence c\[(\d+),(\d+)\]/){
	$rna{$1} = $2;
	my $length = $2 - $1 + 1;
	$rna_name{$1} = "$1..$2\t-\t$length\ttRNA\ttRNA\t-\t-\t-\ttRNA";
	$file_rnt .= $rna_name{$1};
	$file_rnt .= "\n";
	$trna_num++;
		}
	
	}
	print NEWRNT "Complete RNA information\n";
	print NEWRNT "tRNA: $trna_num\n";
	print NEWRNT join("\t", qw(Location Strand Length PID Gene Synonym Code COG Product)),"\n";
	print NEWRNT $file_rnt;
	close TMRNA;
	close NEWRNT;
}

#######read the rnt file. genetare a hash with start as key and end as value  
my %rna = ();	#rna left coordinate as key, right coordinate as value
my %rna_name = ();	#rna left coordinate as key, line content as value

if(-e $seq_RNA){
	open RNA, "$seq_RNA";
	##### example of a .rnt file
  #Complete RNA information
  #tRNA: 88
  #Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product
  #17805..17880	+	76	tRNA	tRNA	-	-	-	tRNA
  #122234..122310	+	77	tRNA	tRNA	-	-	-	tRNA
  #122347..122422	+	76	tRNA	tRNA	-	-	-	tRNA
	while(<RNA>){
		if(/^\d+\.\./){
			chomp;
			my $line = $_;
			my @line_rna = split (/\t/, $line);
			my $rna_left = $1, my $rna_right = $2 if $line_rna[0] =~ /(\d+)\.\.(\d+)/;
			my $length = $rna_right - $rna_left;
			next if $length < 60;	#select complete rna gene
			if($line_rna[8] eq '-'){
				if($line_rna[4] eq "ssrA"){
					$rna{$rna_left} = $rna_right;
					$rna_name{$rna_left} = "$line\tRNA";
				}
			}else{
				$rna{$rna_left} = $rna_right;
				$rna_name{$rna_left} = "$line\tRNA";
			}
		}
	}
}

my %rna_inverse = reverse %rna;
my $count_rna = 0;
my @rna_left;
my @rna_right;
# save the rna coordanates by order
foreach my $rna_key (sort {$a<=>$b} keys %rna){
	$rna_left[$count_rna] = $rna_key;
	$rna_right[$count_rna] = $rna{$rna_key};
	$count_rna++;
}
##########################################################
##check integrase hmmscan results
##determine candidate regions
##########################################################
##parse hmm result and get unique gi of integrase
if(! -e $int_hmm_out){
	print "ERROR!For $job_id, $int_hmm_out file was not found!\n";
	open ERROR, ">$error_log";
	print ERROR "ERROR!For $job_id, $int_hmm_out file was not found!\n";
	close ERROR;
	exit;
}

my $int_count = 0;
my %integrase_gi;
my $int_hmm_parse_result = &parse_hmm_result($int_hmm_out);
##### Example of $int_hmm_parse_result
## 3873186351	rve
## 3873186602	Phage_integrase	Phage_integrase
## 3873186988	rve
my @int_hmm_parse_line = split(/\n/, $int_hmm_parse_result);
my $ii = 0;
my $i = 0;
my @can_region;
while ($int_hmm_parse_line[$ii]){
	my @c = split(/\t/, $int_hmm_parse_line[$ii]);
	$ii ++;
	my $gi = $c[0];
	$integrase_gi{$gi} = $c[1];
	$int_count ++;
	###looking for the candidate region#######
	my $int_left_end = $ptt_gi_coor{$gi};
	my $int_right_end = $ptt_coordinate{$ptt_gi_coor{$gi}};
	##boundaries by tRNA and tmRNA;
	$count_rna = 0;
	my $can_region_left_rna = 1;
	my $can_region_right_rna = $genome_size;
	while(($rna_right[$count_rna])&&($rna_right[$count_rna] < $int_left_end)){
		$can_region_left_rna = $rna_left[$count_rna];
		$count_rna++;
	}
	if($rna_right[$count_rna]){
		$can_region_right_rna = $rna_right[$count_rna];
	}
	my $can_region_left = $int_left_end - $left_region_range;
	if($can_region_left<0){
		$can_region_left = 1;
	}
	my $can_region_right = $int_right_end + $right_region_range;
	if($can_region_right>$genome_size){
		$can_region_right = $genome_size;
	}
	if($can_region_left_rna > $can_region_left){
		$can_region_left = $can_region_left_rna;
	}
	if($can_region_right_rna < $can_region_right){
		$can_region_right = $can_region_right_rna;
	}
	if($int_count == 1){
		$can_region[$i][0] = $can_region_left;
		$can_region[$i][1] = $can_region_right;
	}elsif(($can_region_left < $can_region[$i][1])&&($can_region_left != $rna_inverse{$can_region[$i][1]})){
		$can_region[$i][1] = $can_region_right;
	}else{
		$i++;
		$can_region[$i][0] = $can_region_left;
		$can_region[$i][1] = $can_region_right;
	}
}
if($int_count == 0){
	print "For $job_id, no integrase was not found! So the job exited!\n\n";
	open ERROR, ">$error_log";
	print ERROR "For $job_id, no integrase was not found! So the job exited!\n\n";
	close ERROR;
	system("rm -rf $download_dir");
	exit;
}

###################################################analysis final candidate region####################################
##############################################################
##extract fna and faa
##############################################################
my $k = 0;
$i = 0;

while($can_region[$k][0]){
	my $candidate_region_fna = $candidate_dir."candidate_region_".$i.".fna";
	my $candidate_region_faa = $candidate_dir."candidate_region_".$i.".faa";
	my $candidate_region_fea = $candidate_dir."candidate_region_".$i.".fea";
	##extract fna
	my $string = $genome_seq_obj->subseq($can_region[$k][0],$can_region[$k][1]);
	my $seq_obj = Bio::Seq->new(-display_id => $genome_id, -seq => $string, -desc =>"candidate_region_$i $can_region[$k][0]..$can_region[$k][1]" );
	my $seqio_obj = Bio::SeqIO->new(-file => ">$candidate_region_fna", -format => 'fasta' );
	$seqio_obj->write_seq($seq_obj);
	##extract faa and annotate information
	open CANFAA, ">$candidate_region_faa";
	open CANFEA, ">$candidate_region_fea";
	print CANFEA "#$genome_name\n".
	             "#Genome GC $genome_gc %\n".
	             "#Candidate region: $can_region[$k][0]..$can_region[$k][1]\n";
######### .fea file example
##test
##Genome GC 0.5 %
##Candidate region: 3433480..3504627
#3433480..3433555	-	76	tRNA	tRNA	-	-	-	Asn tRNA	RNA
#3433717..3434979	+	420	2783483349	-	KPHS_34600	-	-	integrase	INT
#3435173..3436477	-	434	2783483350	-	KPHS_34610	-	-	salicylate synthase Irp9
#3436505..3437785	-	426	2783483351	-	KPHS_34620	-	-	MFS superfamily transporter signal transducer

	print CANFEA "$rna_name{$can_region[$k][0]}\n" if exists $rna{$can_region[$k][0]};
	my $line_num = 0;
	foreach (sort {$a<=>$b} keys %ptt_coor_gi){
		my $pro_left = $_;
		my $pro_right = $ptt_coordinate{$pro_left};
		my $pro_gi = $ptt_coor_gi{$pro_left};
		my $line = $ptt_gi_line{$pro_gi};
		if($pro_left>=$can_region[$k][0] && $pro_right < $can_region[$k][1]){
			if(!exists $pro{$pro_gi}){
			}else{
				$line_num ++;
				print CANFAA "$pro{$pro_gi}\n";
				print CANFEA "$line";
				if(exists $integrase_gi{$pro_gi}){
					print CANFEA "\tINT";
				}
				print CANFEA "\n";
			}
		}
	}
	close CANFAA;
	if(exists $rna_inverse{$can_region[$k][1]}){
		print CANFEA "$rna_name{$rna_inverse{$can_region[$k][1]}}\n";
	}
	close CANFEA;
	if($line_num == 0){
		system("rm -f $candidate_region_fna");
		system("rm -f $candidate_region_faa");
		system("rm -f $candidate_region_fea");
		$k++ ;
		next;
	}

	##########################################
	##submit hmmer jobs
	##########################################
	if(! -e $candidate_region_faa){
	  open LOG, ">>$error_log";
		print LOG "For $job_id, $candidate_region_faa file is not found\n";
		close LOG;
		exit;
	}
	my $ice_hmm_out = $candidate_dir."ICE_hmm_$i.out";
	my $hmm_ice_cmd = "$hmmscan_cmd --domtblout $ice_hmm_out --cpu 20 $ice_hmm_db $candidate_region_faa 1>/dev/null 2>>$error_log";### to dectect ICE specific protein
	`$hmm_ice_cmd`;
	
	my $orit_blastn_cmd = "perl ./scripts/oriT_scanner.pl -a $job_id -n $i -t $opt_t";
	`$orit_blastn_cmd`;
	
	$k++ ;
	$i++;

}


##################################################
##check hmm result
##################################################
print "Start analysing the Releaxase, T4SS, T4CP, Tra, Rep, ARG, VF of the candidate regions...\n";

my $runing_job_count0 = $i; ## for count rep_tag for AICE
my $runing_job_count = $i;
my $runing_job_num0 = $i;  ## for count rep_tag for AICE
my $runing_job_num = $i;


## count rep_tag for AICE
while($runing_job_count0){
	$i = 0;
	while($i < $runing_job_num0){
			$runing_job_count0--;
			my $ice_hmmer_out = $candidate_dir."ICE_hmm_$i.out";
			my $goon_ssICE = $1, my $goon_IME = $2, my $goon_AICE = $3 if &check_candidate_region($ice_hmmer_out) =~ /(\d+)\t(\d+)\t(\d+)/;  # $continue\t$ss #$ss = $hit{Mob} * ($hit{T4SS} + $hit{Tra}); #$continue = ($hit{Mob} * ($hit{T4SS} + $hit{Tra})) + ($hit{Rep} * $hit{Tra});    
      if($goon_AICE >0){$rep_tag += 1;}
      $i++;
	
			if($runing_job_count0 == 0){
			last;
		}
	}
	if($runing_job_count0 == 0){
		last;
	}
}


my %ice_candidate_no;
my $ice_candidate_count = 0; 
while($runing_job_count){
	$i = 0;
	while($i < $runing_job_num){
			$runing_job_count--;
			my $ice_hmmer_out = $candidate_dir."ICE_hmm_$i.out";
			my $goon_ssICE = $1, my $goon_IME = $2, my $goon_AICE = $3 if &check_candidate_region($ice_hmmer_out) =~ /(\d+)\t(\d+)\t(\d+)/;  # $continue\t$ss #$ss = $hit{Mob} * ($hit{T4SS} + $hit{Tra}); #$continue = ($hit{Mob} * ($hit{T4SS} + $hit{Tra})) + ($hit{Rep} * $hit{Tra});    
			## oriT tag for IME 
			my $orit_out = $candidate_dir."oriT_$i.result";  
      my $orit_tag = 0;      
      if(-e $orit_out){      
      	$orit_tag = 1;		     
      }                 
			if(($goon_ssICE + $goon_IME + $goon_AICE + $orit_tag)){
				$ice_candidate_no{$i} = 1;
				$ice_candidate_count++;
				###############################
				##submit vmatch job
				###############################
				my $candidate_region_fna = $candidate_dir."candidate_region_".$i.".fna";
				if(! -e $candidate_region_fna){
					print "ERROR! For $job_id, $candidate_region_fna file is not found!\n";
					open LOG, ">>$error_log";
					print LOG "ERROR! For $job_id, $candidate_region_fna file is not found!\n";
					close LOG;
					exit;
				}
				my $run_analyze_cmd = "perl $analyze_cmd -i $opt_i -a $job_id -n $i -r $rep_tag -d $dr_max_mismatch -g $genome_gc -t $opt_t -o $opt_o "; 
				system($run_analyze_cmd);
			}			
		$i++;
		if($runing_job_count == 0){
			last;
		}
	}
	if($runing_job_count == 0){
		last;
	}

}
################################################
##generate the download file 
################################################
if(!$ice_candidate_count){
	my $complete = $tmp_path."candiate_region_none";
	print  "The ICE/IME detection has been doned. No ICE/IME candidate region was detected!\n\n"; 
	open FIN,">$complete";
	print FIN "The ICE/IME detection has been doned. No ICE/IME candidate region was detected!\n\n";
	close FIN;
	system("rm -rf $download_dir");
	exit;
}
my $ice_count = 0;
$runing_job_count = $ice_candidate_count;
$runing_job_num = $ice_candidate_count;
while($runing_job_count){
	foreach (sort {$a<=>$b} keys %ice_candidate_no){
		$i = $_;
			$runing_job_count--;
			opendir (DIR, "$download_dir");
			for my $file (readdir DIR){
				if($file =~ /^ICEfinder\_$i\_(\d+)$/){
					$ice_count ++;
					open IN, "$download_dir$file";
					my $ice_pro = $download_dir."Protein_$file.fas";
					open ICEFAA, ">$ice_pro";
					while(<IN>){
						if(/^Location: (\d+)\.\.(\d+)/){
							my $ice_fas = $download_dir."DNA_$file.fas";
							my $ice_seq = $genome_seq_obj->subseq($1,$2);
							my $ice_des = $genome_desc.": $1..$2";
							my $ice_seq_obj = Bio::Seq->new(-display_id => $genome_id, -seq => $ice_seq, -desc =>$ice_des );
							my $ice_seqio_obj = Bio::SeqIO->new(-file => ">$ice_fas", -format => 'fasta' );
							$ice_seqio_obj->write_seq($ice_seq_obj);
		
						}elsif(/^\d+\.\.\d+\t.\t\d+\t(\d+)/){
							print ICEFAA "$pro{$1}\n";
						}
					}
					close ICEFAA;
					#my $run_ARG_VF_cmd= "perl $ARG_VF_cmd $job_id $file";
					#system($run_ARG_VF_cmd);
				}
			}
			if($runing_job_count == 0){
				last;
			}
		
		if($runing_job_count == 0){
			last;
		}
	}
}

################################
### result summary
################################
my $ice_all_out = $download_dir.$job_id."_summary.txt\n";
my $left_ice_all; my $right_ice_all;
#my $ARG = "";
#my $VF = "";
#my $ARG_none = "";
#my $VF_none = "";
#my $hit_ARG = 0;
#my $hit_VF = 0;
opendir (DIR, "$download_dir");
#######
#Description: Putative ICE with T4SS
#Location: 3433540..3495705
#oriT: 3479690..3479813
#Length: 62166 bp
#GC content: 52.46 %
#DR: 3433540..3433556	3495689..3495705
	for my $file (readdir DIR){
		 if($file =~ /^ICEfinder\_(\d+)\_(\d+)$/){	
		 	## ARG and VF statistic
		 	#$ARG = "$download_dir$file.ARG";
		 	#$ARG_none = "$download_dir$file.ARG.none";
		 	#$VF = "$download_dir$file.VF";
		 	#$VF_none = "$download_dir$file.VF.none";
		 	#if (-e $ARG){
		 	#	open (RESULT5, $ARG);
		 	#		while(<RESULT5>){
		 	#			next unless /^\d+/;	
			#			chomp;
			#			$hit_ARG ++;
			#	  }
			#	close RESULT5;
			#}elsif(-e $ARG_none){$hit_ARG = 0;}
			#if (-e $VF){
			#	open (RESULT6, $VF) ;
			#	while(<RESULT6>){
			#		next unless /^\d+/;	
			#		chomp;
			#		$hit_VF ++;
			#	}
			#	close RESULT6;	
			#}elsif(-e $VF_none){$hit_VF = 0;}	
			## write summary file	 	
			open IN, "$download_dir$file";		                              
		  while(<IN>){ 
		  	if(/^Description: (Putative.*)/){
		  		chomp;
		  		open ICEALLOUT, ">>$ice_all_out";                           
		    	print ICEALLOUT "$job_id\t$genome_name\t$genome_size\t$file\t$1:\t";
		      close ICEALLOUT;  
		  		}                                                   
		  	elsif(/^Location: (\d+)\.\.(\d+)/){                              
		    	$left_ice_all=$1; $right_ice_all=$2;
		    	my $ice_len = $right_ice_all - $left_ice_all + 1;                        
		    	open ICEALLOUT, ">>$ice_all_out";                           
		    	print ICEALLOUT "$left_ice_all..$right_ice_all\t$ice_len\t";
		      close ICEALLOUT;                                            
				}
				elsif(/^oriT: /){
					open ICEALLOUT, ">>$ice_all_out";  
					chomp;
					my $orit_desc = (split(/\s+/,$_))[1];                         
		    	print ICEALLOUT "$orit_desc\t";
		      close ICEALLOUT;
		    }	
				elsif(/^GC content: (\d+\.\d+) /){
					my $gc_diff = $1- $genome_gc;
					$gc_diff = sprintf("%.2f",abs($gc_diff));
					open ICEALLOUT, ">>$ice_all_out";                           
		    	print ICEALLOUT "$1\t$genome_gc\t$gc_diff\t0\t0\n";
		      close ICEALLOUT;
		    }
											
			}
		}
		close IN;
	}
close DIR;


my $complete = $tmp_path."region_finder_finish";
open FIN,">$complete";
close FIN;


#############################################################
##calculate genomic average GC content
#############################################################
sub calculate_gc{
	my $sequence = shift;
	my $seq_length = length($sequence);
	my $num_a    = $sequence =~ tr/[Aa]//;
	my $num_c    = $sequence =~ tr/[Cc]//;
	my $num_g    = $sequence =~ tr/[Gg]//;
	my $num_t    = $sequence =~ tr/[Tt]//;
	my $num_o    = $seq_length - $num_a - $num_c - $num_g - $num_t;
	my $gc       = sprintf( "%.2f", ( $num_c + $num_g ) / $seq_length * 100 );
	return $gc;
}

#############################################################
##parse HMM --domtblout results
#############################################################
#\#                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
#\# target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
#\#------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
#TIGR02224            TIGR02224    295 2584642392_1         -            301  7.9e-110  355.8   0.0   1   1  4.4e-110  8.9e-110  355.6   0.0     2   295    12   298    11   298 0.96 recomb_XerC: tyrosine recombinase XerC
#TIGR02225            TIGR02225    290 2584642392_1         -            301   2.6e-85  275.2   0.0   1   1   1.5e-85     3e-85  275.0   0.0     2   290    12   296    11   296 0.97 recomb_XerD: tyrosine recombinase XerD
#Phage_integrase      PF00589.17   173 2584642392_1         -            301   1.4e-53  170.4   0.0   1   1   9.5e-54   1.9e-53  170.0   0.0     3   172   116   283   114   284 0.98 Phage integrase family
sub parse_hmm_result{
	my ($hmm) = @_;
	my $content = "";
	my @hit_dom;	#aligned domians of regions within a protein
	open (HMM, "$hmm") || die "open $hmm error: $!\n";
	my $i = 0;
	my $do_num = 0;
	my $gi_tmp = 0;
	while(<HMM>){
		unless(/^#/){
			chomp;
			my @array = split(/\s+/,$_);
			next unless $array[0];

			my $gi = $1 if (($array[3] =~ /^gi\|(\d+)\|/) || ($array[3] =~ /^(\d+)\_\d/));##for select genome or standard parsed gbk#

			next if $array[6] > $threshold_evalue{$array[0]};######evalue#######
			next if $array[11] > $threshold_cvalue{$array[0]};#####cvalue#######
			if($gi != $gi_tmp){
				$do_num = 0;
				$hit_dom[$i][$do_num] = $array[0];
				$i ++;
				if($i-1){
					$content .= "$gi_tmp";
					my $m = 0;
					while($hit_dom[$i-2][$m]){
						$content .= "\t$hit_dom[$i-2][$m]";
						$m ++;
					}
					$content .= "\n";
				}
				$gi_tmp = $gi;
				next;
			}
			$do_num ++;
			$hit_dom[$i-1][$do_num] = $array[0];
		}
	}
	close HMM;
	if(($i-1)>0){ ## note
		$content .= "$gi_tmp";
		my $m = 0;
		while($hit_dom[$i-1][$m]){
			$content .= "\t$hit_dom[$i-1][$m]";
			$m ++;
		}
		$content .= "\n";
	}
	return $content;
}

#############################################################
##Check if a candadate region contains ICE features
#############################################################
sub check_candidate_region{ 
	my $hmm = shift;
	my $continue_ssICE = 0;
	my $continue_IME = 0;
	my $continue_AICE = 0;
	my %hit = (
	           'Rep'  => 0,
	           'Tra'  => 0,
	           'Mob'  => 0,
	           'T4CP' => 0,
	           'T4SS' => 0,
	           );
	open (HMM, "$hmm") || die "open $hmm error: $!\n";
	while(<HMM>){
		unless(/^#/){
			chomp;
			my @array = split(/\s+/,$_);
			next unless $array[0];
			my $gi = $1 if $array[3] =~ /^gi\|(\d+)\|/;
			next if $array[0] eq "Pfam-B_3022";
			next if $array[6] > $threshold_evalue{$array[0]};
			next if $array[11] > $threshold_cvalue{$array[0]};
			$hit{$dom_function{$array[0]}} = 1;
		}
	}
	close HMM;
	if ($hit{Mob}+$hit{T4SS}){
		$continue_ssICE = 1;
		$continue_IME =1;
	}
	elsif($hit{Rep}){
		$continue_AICE =1;
	}
	return "$continue_ssICE\t$continue_IME\t$continue_AICE";
}
