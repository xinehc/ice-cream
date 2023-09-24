#!/usr/bin/perl
# Meng Liu and Hong-Yu Ou on April-24-2018
# School of Life Sciences & Biotechnology, Shanghai Jiao Tong University

use strict;
# use Bio::Perl;
use Bio::SeqIO;
use Getopt::Std;
our ($opt_h, $opt_i, $opt_t, $opt_o, $opt_a, $opt_n, $opt_r, $opt_d, $opt_g) = "";
my $usage = q(
This script is used to analyse and score the releaxase, t4cp, Tra, Rep, and oriT of the candidate regions, 
and further delimit the ICE/IME with direct repeats(DR).

Usage:	perl script_name.pl <job_id> <region_no> <rep_tag> <maxmismatch> <genome_gc>


);
die($usage) if ( @ARGV < 1 ); 
getopts('i:t:o:a:n:r:d:g:h');

###############################################
##declare major variables
###############################################
my $job_id = $opt_a;
my $i = $opt_n; #candidate region ID
my $rep_tag = $opt_r;
my $dr_max_mismatch = $opt_d;
my $genome_gc = $opt_g;
my $tmp_path = "$opt_t/$job_id/";
my $download_dir = "$opt_o/$job_id/";
my $can_region_left;	#candidate region left coordinate
my $can_region_right;	#candidate region right coordinate
my $candidate_dir = $tmp_path."candidate/";
my $seq_fna = $tmp_path."$job_id.fna";
my $db_dir = "./data/";
my $candidate_region_fna = $candidate_dir."candidate_region_".$i.".fna";
my $candidate_region_faa = $candidate_dir."candidate_region_".$i.".faa";
my $candidate_region_fea = $candidate_dir."candidate_region_".$i.".fea";
my $ice_hmmer_out = $candidate_dir."ICE_hmm_$i.out";
my $int_hmmer_out = $tmp_path."int_hmm.out";
my $orit_out = $candidate_dir."oriT_$i.result";
my $direct_repeat_out = $candidate_dir."DR_vmatch_$i.out";
my $Log = $tmp_path."run_region_analyzer.log";
my $RDP_parse_out = $tmp_path.$job_id."_all_blastn_vs_RDP_parsed.txt";
my $ice_tag = 0;
my $aic_tag = 0;
my $at4_tag = 0;
my $ime_tag = 0;
my $gram_postive_tag=0;
if (-e "$opt_t/$job_id/$job_id.gram"){
	$gram_postive_tag=1;
}
open (LOG,">>$Log");

##############################################
##classifiy proteins
##############################################
my %ice_info_designate_no = ( 
                              'Phage_integrase'     => 0,
                              'Pfam-B_1819'					=> 0, 
                              'Recombinase'         => 0,
                              'DDE_Tnp_Tn3 '        => 0,
                              'rve'                 => 0,
                              'RepSAv2'             => 6,
                              'DUF3631'      				=> 6,
                              'Prim-Pol'            => 6,
                              'AAA_10'              => 1,
                              'TraC_F_IV'						=> 1, 
                              'CagE_TrbE_VirB' 			=> 1, 
                              'TrwB_AAD_bind'       => 2,
                              'T4SS-DNA_transf'     => 2,
                              'TraG-D_C'            => 2,
                              'TraD_N'              => 2,
                              'FtsK_SpoIIIE'        => 2,
                              'MobC'                => 3,
                              'Pfam-B_6973'         => 3,
                              'Replic_Relax'        => 3,
                              'TrwC'                => 3,
                              'TraI'                => 3,
                              'TraI_2'              => 3,
                              'Relaxase'            => 3,
                              'MobA_MobL'           => 3,
                              'Mob_Pre'             => 3,
                              'Rep_trans'           => 3,
                              'Pfam-B_4801'         => 3,
                              'MOBL'								=> 3,
                              'TraE'                => 4,
                              'TraH'                => 4,
                              'TraK'                => 4,
                              'TraL'                => 4,
                              'TraN'                => 4,
                              'TraU'                => 4,
                              'TraV'                => 4,
                              'Pfam-B_1474'         => 4,
                              'TrbC_Ftype'          => 4,
                              'Pfam-B_1110'         => 4,
                              'DUF2895'             => 4,
                              'Pfam-B_1690'         => 4,
                              'Pfam-B_724'          => 4,
                              'Pfam-B_2070'         => 4,
                              'traK_tpyeI'          => 4,
                              'traL_tpyeI'          => 4,
                              'IcmL'                => 4,
                              'DUF3625'             => 4,
                              'DUF1525'             => 4,
                              'DUF2976'             => 4,
                              'DUF3438'             => 4,
                              'DUF3487'             => 4,
                              'DUF4400'             => 4,
                              'Plasmid_RAQPRD'      => 4,
                              'T2SSE'             	=> 4,
                              'T4SS'             		=> 4,
                              'Tra_M'             	=> 4,
                              'TraA'             		=> 4,
                              'TraF'             		=> 4,
                              'TraG_N'             	=> 4,
                              'TraQ'            		=> 4,
                              'TraW_N'             	=> 4,
                              'TraX'             		=> 4,
                              'TraY'             		=> 4,
                              'TrbC'             		=> 4,
                              'TrbH'             		=> 4,
                              'TrbI'             		=> 4,
                              'traP_tpyeI'          => 4,
                              'traQ_tpyeI'          => 4,
                              'Pfam-B_13874'        => 4,
                              'Pfam-B_3616'         => 4,
                              'Pfam-B_706'          => 4,
                              'VirB3'               => 4,
                              'TrbL'                => 4,
                              'VirB8'               => 4,
                              'CagX'                => 4
);
my %feature_designate_no_info = ( 0 => 'Integrase',
                                  1 => 'T4SS ATPase',
                                  2 => 'T4CP',
                                  3 => 'Relaxase',
                                  4 => 'T4SS',
                                  5 => '-',
                                  6 => 'Replication'
);

########################################################
##read hmmscan results, generate a hash with gi as key and annotation as value
########################################################
my %annotation = ();
my @mob = ("MobC","Pfam-B_6973","Replic_Relax","TrwC","TraI","TraI_2","Relaxase","MobA_MobL","Mob_Pre","Rep_trans","Pfam-B_4801");
if(-e $ice_hmmer_out){
	open (ICEOUT, "$ice_hmmer_out") || die "open file error: $!\n";
########Example of $ice_hmmer_out
##                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
## target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
##------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
#Pfam-B_3022          PB003022     545 9883795155_1         -            571   3.1e-06   17.5   0.0   1   1   2.7e-07   6.4e-06   16.5   0.0   198   230   342   374   333   385 0.87 -
#AAA_10               PF12846.2    328 9883795155_1         -            571   0.00011   13.4   0.1   1   2   0.00013    0.0032    8.5   0.0     4    23   343   362   340   379 0.84 AAA-like domain	
	my $gi_tmp = 0;
	while(<ICEOUT>){
		my @array = split(/\s+/,$_);
		if(/(\S+)\s+.+\s+gi\|(\d+)\|\s+.+\s+(\d+)/ || /(\S+)\s+.+\s+(\d+)\_\d+\s+\-\s+(\d+)/){
			if($2 == $gi_tmp){
				next;
			}
			$gi_tmp = $2; ## gi is part of query name 
			$annotation{$gi_tmp} = $1; ## target name
			if((grep /^$1$/,@mob) && ($3 >= 199)){ ## the minium length of MOB is set to 200 aa 
				$annotation{$gi_tmp} = $1;
			}elsif((grep /^$1$/,@mob) && ($3 < 199)){
				$annotation{$gi_tmp} = "";
			}	
		}
	}
}
close ICEOUT;

if(-e $int_hmmer_out){
	open (INTOUT, "$int_hmmer_out") || die "open file error: $!\n";
	my $gi_tmp = 0;
	while(<INTOUT>){
		my @array = split(/\s+/,$_);
		if(/(\S+)\s+.+\s+gi\|(\d+)\|\s+.+\s+(\d+)/ || /(\S+)\s+.+\s+(\d+)\_\d+\s+\-\s+(\d+)/){
			if($2 == $gi_tmp){
				next;
			}
			$gi_tmp = $2;
			$annotation{$gi_tmp} = $1;
			if($3 >= 150){ ## the minium length of INT is set to 150 aa ; 
				$annotation{$gi_tmp} = $1;
			}elsif($3 < 150){
				$annotation{$gi_tmp} = "";
			}	
		}
	}
}
close INTOUT;

######################################################
## get the oriT with the max hvalue 
## and generate a hash recoding oriT coordinate 
#################################################
my @cut_off;
my @orit_array;
my $max = 0;
my %orit_coordinate;## key is the left; vaule is the right
my $orit_left = 0;
my $candidate_region_left = 0;
my $candidate_region_right = 0;
open (OUT2,"$candidate_region_fna") ||die "open file $candidate_region_fna error: $!\n";
	my $line1 = readline(OUT2);
	close OUT2;
	chomp($line1);
	if($line1 =~ /.+\s+.+\s+(\d+)\.\.(\d+)$/){
		$candidate_region_left = $1;
		$candidate_region_right = $2;
}
print LOG "***********************************************************\nCandidate region$i:candidate_region_left..candidate_region_right:$candidate_region_left..$candidate_region_right\n";	

my $orit_tag = 0;
if(-e $orit_out){
	open (OUT, "$orit_out") || die "open file $orit_out error: $!\n";
	while(<OUT>){
		next if (/^qseqid/);
		chomp; #qseqid	sseqid	identities	length	coverage	hvalue	qstart	qend	sstart	send
		my @orit_array = split(/\t/, $_); 
		push(@cut_off, $orit_array[5]);
	}
	close OUT;
	foreach (@cut_off) {
		$max = $_ if $_ > $max;
	}
	
	open (OUT2,"$candidate_region_fna") ||die "open file $candidate_region_fna error: $!\n";
	my $line1 = readline(OUT2);
	close OUT2;
	chomp($line1);
	if($line1 =~ /.+\s+.+\s+(\d+)\.\.(\d+)$/){$candidate_region_left = $1;}
	
	open (OUT, "$orit_out") || die "open file $orit_out error: $!\n";
	while(<OUT>){
		next if (/^qseqid/);
		chomp; #qseqid	sseqid	identities	length	coverage	hvalue	qstart	qend	sstart	send
		@orit_array = split(/\t/, $_); 
		next if $orit_array[5] != $max;
		$orit_left = $orit_array[6] + $candidate_region_left;
		$orit_coordinate{$orit_left} = $orit_array[7] + $candidate_region_left; ## orit_right
	}
	close OUT;
		
}




#######################################################
##score orfs and find conjugation region
#######################################################
my $rna_count = 0;	##record number of rna genes
my %rna = ();	##record rna coordinates if there is any. left coordinate as key, right coordinate as value
my %rna_name = ();	##record rna infomation if there is any. left coordinate as key, information as value
my %rna_strand = ();	##record rna strand if there is any. left coordinate as key, strand as value
my $position = 0;	##orf position
my @score = ([0, 0, 0, 0, 0, 5, 0]);	##$socre[$position][0], int; $socre[$position][1], virB4&traU; $socre[$position][2], T4CP; $socre[$position][3], MOB; $socre[$position][4], T4SS; $socre[$position][5], note; $socre[$position][6], Rep
my @orf_score = (0);
my @inter_feature;	##regions between each feature
my $inter_feature_count = 0;	##number of regions between each feature
my @conj_region;	#t4ss and t4cp
my $conj_region_count = 0;
my @inter_conj;
my $inter_conj_count = 0;
my @int_region;
my $int_region_count = 0;
my $distance = 0;	##distance (orf number) between T4SS proteins
my $max_distance = 5;	##allowed max distance (orf number) +2 between conjugation proteins in one ICE
my @pos_coor;	## coordinate of each position
my %ptt; ## protein information
my %ptt2; ## signature protein info
my %orf;	##protein coordinate


open (CANFEA, "$candidate_region_fea") || die "open file error: $!\n";
#\#Streptomyces avermitilis MA-4680 = NBRC 14893 DNA, complete genome.
#\#Genome GC 68.77 %
#\#Candidate region: 4589760..4710527
#4589760..4589835	+	76	tRNA	tRNA	-	-	-	tRNA	RNA
#4590262..4590483	+	73	739792727	-	SAVERM_RS19135	-	-	hypothetical protein
#4590544..4591728	-	394	1373155552	-	SAVERM_RS19140	-	-	site-specific integrase	INT
#4591725..4591925	-	66	499293882	-	SAVERM_RS19145	-	-	DNA-binding protein
#4592396..4593871	-	491	499293884	-	SAVERM_RS19150	-	-	hypothetical protein
while(<CANFEA>){
	if(/^\d+\.\./){
		chomp;
		my @line = split (/\t/, $_);	#0, location; 1, strand; 2, length/aa; 3, PID; 4 ,gene; 5, synonym; 6, code; 7, cog; 8, product; [9] note
		$line[4] = $line[5] if $line[4] eq '-';
		if(($line[9]) && ($line[9] eq "RNA")){	#read rna information
			$rna_count ++;
			my @rna_coordinate = split (/\.\./, $line[0]);
			$rna{$rna_coordinate[0]} = $rna_coordinate[1];
			$rna_name{$rna_coordinate[0]} = "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\tRNA\t$line[8]\n";
			$rna_strand{$rna_coordinate[0]} = $line[1];
			next;
		}else{	#read protein information
			$position ++;
			my $orf_left;
			my $orf_right;
			if($line[0] =~ /(\d+)\.\.(\d+)/){
				$orf_left = $1;
				$orf_right = $2;
				$pos_coor[$position][0] = $1;
				$pos_coor[$position][1] = $2;
			}
			if($distance){
				$distance ++;
				if($distance > $max_distance){
						$distance = 0;
						$conj_region_count ++;
				}
			}
			#initialize values
			$score[$position][0] = 0; #int
			$score[$position][1] = 0; #virB4&traU
			$score[$position][2] = 0; #T4CP
			$score[$position][3] = 0; #MOB
			$score[$position][4] = 0; #T4SS
			$score[$position][5] = 5;	#note
			$score[$position][6] = 0; #Rep
			$orf_score[$position] = $orf_score[$position-1];
			##integrase?
			if(($line[9]) && ($line[9] eq "INT")){
				$score[$position][0] = 1;
				$score[$position][5] = 0;
				$orf_score[$position] += $score[$position][0];
				$inter_feature[$inter_feature_count][1] = $orf_right;
				$inter_feature_count ++;
				$inter_feature[$inter_feature_count][0] = $orf_left;
				$int_region[$int_region_count][0] = $orf_left;
				$int_region[$int_region_count][1] = $orf_right;
				$int_region_count ++;
			}else{
			##conjugation proteins?
				if(exists $annotation{$line[3]}){
					#print "okokok\n";
					$score[$position][$ice_info_designate_no{$annotation{$line[3]}}] = 1;
					$orf_score[$position] += $score[$position][$ice_info_designate_no{$annotation{$line[3]}}];
					$score[$position][5] = $ice_info_designate_no{$annotation{$line[3]}};
					$inter_feature[$inter_feature_count][1] = $orf_right;
					$inter_feature_count ++;
					$inter_feature[$inter_feature_count][0] = $orf_left;
					if($score[$position][5] == 2 || $score[$position][5] == 4 || $score[$position][5] == 3 ||$score[$position][5] == 6){	# T4CP ; T4SS ; add relaxase(for IME) && Rep(for AICE)
						if($distance == 0){
							$distance ++;
							$inter_conj_count ++;
							$inter_conj[$inter_conj_count][0] = $orf_left;
							$conj_region[$conj_region_count][0] = $orf_left;
						}else{
							$distance = 1;
						}
						$inter_conj[$inter_conj_count-1][1] = $orf_right;
						$conj_region[$conj_region_count][1] = $orf_right;
					}
				}
			}
			$orf{$orf_left} = $orf_right; ## orf coordinates are saved in %orf 
			$ptt{$orf_left} = "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$feature_designate_no_info{$score[$position][5]}\t$annotation{$line[3]}\t$line[8]\n";

			if (($feature_designate_no_info{$score[$position][5]} ne "-") && ($feature_designate_no_info{$score[$position][5]} ne "")){
				$ptt2{$orf_left} = "$orf_left..$orf_right\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$feature_designate_no_info{$score[$position][5]}\t$annotation{$line[3]}\n";#\t$line[8]
				print LOG "Candidate region$i:core_protein: $ptt2{$orf_left}";
			}
		}
	}elsif(/^#Candidate region: (\d+)\.\.(\d+)/){ 
		$can_region_left = $1;
		$inter_conj[$inter_conj_count][0] = $can_region_left;
		$inter_feature[$inter_feature_count][0] = $can_region_left;
		$can_region_right = $2;
	}
}
close CANFEA;

### to count the T4SS components in the whole candidate region
		
    my $int_tag = 0;	
    my $mob_tag = 0;
    my $t4ss_com_num = 0; ## except t4ss ATPase
    my %t4ss_pro_left_name;
    my %t4ss_pro_left_right;
   	my @t4ss_region;## $t4ss_region[0][0]..$t4ss_region[0][1] is the fisrt qualied t4ss region in the candidate region
   	my $h = 0; ## the qualified t4ss number in the candidate region
    foreach (sort {$a<=>$b} keys %orf){                            
		    	my $orf_key = chomp($_);     
		    	if(($_ >= $can_region_left) && ($orf{$_} <= $can_region_right)){                                             
		    		my @ptt2_array_tmp = (split (/\t/,$ptt2{$_}));
		    		my $marker_name = $ptt2_array_tmp[5]; 
		    		my $sign_pro_name = $ptt2_array_tmp[6];
		    		if ($marker_name eq "Integrase" ){
		    			$int_tag += 1;
		    		}
		    		if ($marker_name eq "Relaxase" ){
		    			$mob_tag += 1;
		    		}
		    		if ($marker_name eq "T4SS" or $marker_name eq "T4SS ATPase" ){
		    			chomp($sign_pro_name);
		    			$t4ss_pro_left_name{$_} = $sign_pro_name;
		    			$t4ss_pro_left_right{$_} = $orf{$_};
		    			$t4ss_com_num +=1; ## the candidate t4ss components number in the whole candidate region before co-localization
		    		}
		    		  		
		    	}		    				    						 		                                                                                           
		}
	if( $t4ss_com_num >=5){
    ## t4ss co-localization
    my $t4ss_tag = 0; 
		my @region;
		my @region_name;
		my $i_t4ss = 0;
		my $region_num = 0;
		my $mpf_distance = 10000; ## the distance between ecch Mpf gene <= 10 k 
    my $t4ss_core = 5;
	
		foreach (sort {$a<=>$b} keys %t4ss_pro_left_right){
			if($i_t4ss == 0){## $i_t4ss here serve as a tag
						if($t4ss_pro_left_name{$_} ne "AAA_10" && $t4ss_pro_left_name{$_} ne "TraC_F_IV" && $t4ss_pro_left_name{$_} ne "CagE_TrbE_VirB" ){## the first component of T4SS should not be T4SS ATPase
                $region[$region_num][0] = $_;
                $region[$region_num][1] = $t4ss_pro_left_right{$_};
                $region_name[$region_num] = $t4ss_pro_left_name{$_};
                $i_t4ss ++;
            }else{next;}
        }else{
                if($_ <= ($region[$region_num][1]+ $mpf_distance)){#; can not be $distance here
                        $region[$region_num][1] = $t4ss_pro_left_right{$_};
                        $region_name[$region_num] .= "\t$t4ss_pro_left_name{$_}";
                        #push (@region_name,$t4ss_pro_left_name{$_});
                        #$region_name[$region_num] = $t4ss_pro_left_name{$_};
                        $i_t4ss ++;
                }else{
                			if($t4ss_pro_left_name{$_} ne "AAA_10" && $t4ss_pro_left_name{$_} ne "TraC_F_IV" && $t4ss_pro_left_name{$_} ne "CagE_TrbE_VirB" ){#the first component of T4SS should not be T4SS ATPase
                        $region_num ++;## to count the next T4SS region 
                        $region[$region_num][0] = $_;
                        $region[$region_num][1] = $t4ss_pro_left_right{$_};
                        $region_name[$region_num] = $t4ss_pro_left_name{$_};      
                        $i_t4ss = 1;
                       }else{next;}
                }
        }
			
			}
			

my $u = 0;
my $region_content = "";
while($region[$u][0]){
        $u ++;
				my @arr_t4ss_com = split(/\t/, $region_name[$u-1]);
				my $co_t4ss_com_num = $#arr_t4ss_com + 1; ## the core t4ss components number after co-localization
				my $candidate_t4ss_region_desc = "Candidate region$i:#Candidate T4SS Region $u\t$region[$u-1][0]..$region[$u-1][1];co_t4ss_com_num:$co_t4ss_com_num\n $region_name[$u-1]\n";
				#print LOG "$candidate_t4ss_region_desc";
				if($co_t4ss_com_num >=$t4ss_core){
					my $t4ss_region_desc = "Candidate region$i:#T4SS Region $u\t$region[$u-1][0]..$region[$u-1][1];co_t4ss_com_num:$co_t4ss_com_num\n $region_name[$u-1]\n";
					print LOG "$t4ss_region_desc";
					#my $h = 0; ## the qualified t4ss number in the candidate region
					$t4ss_region[$h][0] = $region[$u-1][0];# $t4ss_region[0][0]..$t4ss_region[0][1] is the fisrt qualied t4ss region in the candidate region
					$t4ss_region[$h][1] = $region[$u-1][1];
					print LOG "test:t4ss_region[$h][0]..t4ss_region[$h][1]:$t4ss_region[$h][0]..$t4ss_region[$h][1]";
					$h ++;
				}
				
				
}

}


$inter_conj[$inter_conj_count][1] = $can_region_right;
$inter_feature[$inter_feature_count][1] = $can_region_right;
my %rna_inverse = reverse %rna;
my $x = 0;
while($conj_region[$x][0]){
	print LOG "Candidate region$i:conj_region$x:conj_region[$x][0]..conj_region[$x][1]:$conj_region[$x][0]..$conj_region[$x][1]\n";
	$x ++;
}
$x = 0;
while($inter_conj[$x][0]){
	$x ++;
}
$x = 0;
while($inter_feature[$x][0]){
	$x ++;
}



###############################################
##run vmatch
###############################################
my $mismatch_arg = "-h $dr_max_mismatch";
if($dr_max_mismatch == 0){
	$mismatch_arg = "";
}

my $dr_min_length = 15;
if($gram_postive_tag == 1 && $rep_tag >= 1)
{$dr_min_length = 40;}
my $maktree_cmd = "mkvtree -db $candidate_region_fna -indexname $candidate_region_fna -dna -pl -lcp -suf -tis -ois -bwt -bck -sti1";
my $vmatch_cmd = "vmatch -l $dr_min_length $mismatch_arg $candidate_region_fna >$direct_repeat_out";
system "$maktree_cmd";
system "$vmatch_cmd";

###############################################
##read and process vmatch results
###############################################
##get DR
open DROUT, "$direct_repeat_out";
#\# args=-l ./tmp/NC_003155_Sav_2/candidate/candidate_region_1.fna
#   16    0  32767   D    16    0  41442   0    3.42e+00     32   100.00
#   15    0  31548   D    15    0  32767   0    1.37e+01     30   100.00
#   15    0  31548   D    15    0  41442   0    1.37e+01     30   100.00
#   15    0  31548   D    15    0  31571   0    1.37e+01     30   100.00
my $dr_count = 0;	#number of DR pairs
my %dr_pair = ();	#left coordinate of upstream DR as key while left coordinate of downstream DR as value
my %dr_pair_reverse = ();	#left coordinate of downstream DR as key while left coordinate of upstream DR as value
my %dr;	#dr coordinate, left coordinate as key while right coordinate as value
while(<DROUT>){
	if(/^\s+(\d+)\s+(\d+)\s+(\d+)\s+D\s+(\d+)\s+(\d+)\s+(\d+)\s+\-?(\d+)\s+.+\s+.+\s+.+/){
		my $dr_u_start = $2 + $3 + $can_region_left;
		my $dr_u_end = $1 + $dr_u_start - 1;
		my $dr_d_start = $5 + $6 + $can_region_left;
		my $dr_d_end = $4 + $dr_d_start - 1;
		my $dr_mismatch = $7;
		my $count_inter = 0;
		my $length_u = $1;
		my $length_d = $4;
		## filter DR length > 200
		if($length_u > 200){
			next;
		}
		##filter DR pair that locates in two RNA genes.
		if($rna_count == 2){
			if(($dr_u_start < $rna{$can_region_left})&&($dr_d_end > $rna_inverse{$can_region_right})){
				next;
			}
		}
		##filter DR pair locating in inter-conjugation region
		my $j = 0;
		my $space_trash = 0;
		while($inter_conj[$j][0]){
			if($dr_u_start >= $inter_conj[$j][0] && $dr_d_end <= $inter_conj[$j][1]){
				$space_trash = 1;
				last;
			}
			$j ++;
		}
		##filter DR pair, select most far pairs
		if(exists $dr_pair{$dr_u_start}){
			if(($dr_d_start - $dr_u_start) < ($dr_pair{$dr_u_start} - $dr_u_start)){
				next;
			}
		}elsif(exists $dr_pair_reverse{$dr_d_start}){
			if(($dr_d_start - $dr_u_start) < ($dr_d_start - $dr_pair_reverse{$dr_d_start})){
				next;
			}
		}
		##filter included DR pair
		$dr{$dr_u_start} = $dr_u_end;
		$dr{$dr_d_start} = $dr_d_end;
		$dr_pair{$dr_u_start} = $dr_d_start;
		$dr_pair_reverse{$dr_d_start} = $dr_u_start;
		$dr_count ++;
	}
}

close DROUT;


my @candidate_dr_up;	##candidate upstream DR
my @candidate_dr_down;	##candidate downstream DR
my $candidate_dr_count = 0;	##number of candidate DR pairs
my $dr_num = 0;	##count of DR pairs
my $gap_len = 3;	##gap length allowed to combine two tendem DRs
my @ice_dr;	##record dr of ICEs
my $ice_dr_count = 0;	##record number of dr of ICEs
my @select_dr_score;	##record score of selected dr; #dr_score[$n][0] = candidate dr number, dr_score[$n][1] = feature score,  dr_score[$n][2] = ice score, dr_score[$n][3] = conj score, dr_score[$n][4] = mob score, dr_score[$n][5] = distance to the nearest integrase,
my $select_dr_number = 0;	##record number score of selected dr;
my $length_thred = 30000; ## DR max length; mliu added 
##process DR and score each DR
foreach (sort {$a<=>$b} keys %dr_pair){
	if($dr_num == 0){
		$candidate_dr_up[$candidate_dr_count][0] = $_;
		$candidate_dr_up[$candidate_dr_count][1] = $dr{$_};
		$candidate_dr_down[$candidate_dr_count][0] = $dr_pair{$_};
		$candidate_dr_down[$candidate_dr_count][1] = $dr{$dr_pair{$_}};
	}elsif(($_ <= ($candidate_dr_up[$candidate_dr_count][1] + $gap_len)) && ($dr_pair{$_} <= ($candidate_dr_down[$candidate_dr_count][1] + $gap_len))){
		#combine overlapping or neighbour DRs
		if($dr{$_} > $candidate_dr_up[$candidate_dr_count][1]){
			$candidate_dr_up[$candidate_dr_count][1] = $dr{$_};
			$candidate_dr_down[$candidate_dr_count][1] = $dr{$dr_pair{$_}};
		}
	}else{	#new DR
		##calculate score of the regions in last DR pair
		my $k = 1;	##count of posistion
		my @dr_position;	##temnimal position in the DR pair
		my @score_sum = (0, 0, 0, 0, 0, 0);	##$socre_sum[0], int; $socre_sum[1], virB&4traU; $socre_sum[2], T4CP; $socre_sum[3], MOB; $socre_sum[4], T4SS; $socre_sum[5], total Score; $score_sum[6],Rep
		##find orf positions in the DR pair
		my $dr_region_len = $candidate_dr_down[$candidate_dr_count][1] - $candidate_dr_up[$candidate_dr_count][0] + 1;
		while($pos_coor[$k][0]){
			if(!$dr_position[0]){
				if($candidate_dr_up[$candidate_dr_count][1] < $pos_coor[$k][0]){
					$dr_position[0] = $k;
					$score_sum[0] += $score[$k][0];	#int
					$score_sum[1] += $score[$k][1];	#virB4&traU ,that is T4SS ATPase 
					$score_sum[2] += $score[$k][2];	#T4CP
					$score_sum[3] += $score[$k][3];	#MOB
					$score_sum[4] += $score[$k][4];	#T4SS:4
					$score_sum[6] += $score[$k][6]; #Rep
				}
			}else{
				if($candidate_dr_down[$candidate_dr_count][0] > $pos_coor[$k][1]){
					$dr_position[1] = $k;
					$score_sum[0] += $score[$k][0];	#int
					$score_sum[1] += $score[$k][1];	#virB4&traU ,that is T4SS ATPase 
					$score_sum[2] += $score[$k][2];	#T4CP
					$score_sum[3] += $score[$k][3];	#MOB
					$score_sum[4] += $score[$k][4];	#T4SS
					$score_sum[6] += $score[$k][6]; #Rep
				}else{
					last;
				}
			}
			$k++;
		}
		if($orit_left > $candidate_dr_up[$candidate_dr_count][0] && $orit_coordinate{$orit_left} < $candidate_dr_down[$candidate_dr_count][1])
		{$orit_tag = 1;}else{$orit_tag = 0;}

		##calculate feature sores
		my $ice_score = 0;
		my $ime_score = 0;
		if($gram_postive_tag == 0){
		$ice_score = ($score_sum[1] + $score_sum[4]) * $score_sum[2] * $score_sum[3] * $score_sum[0] if ($score_sum[4] >2 && $score_sum[1]<4 ) ;	#t4ss ice = (VirB4/TraU or T4SS)+ int + mob	+ t4cp (must have t4cp );at least two mpf gene    
    $ime_score = $score_sum[0] * ($orit_tag + $score_sum[3]) if (($score_sum[4] < 2)&&($score_sum[1] == 0)); ## IME = Int + Mob - (virB4 )- full T4SS
    }
    else{ 
		$ice_score = ($score_sum[1] + $score_sum[4]) * $score_sum[2] * $score_sum[3] * $score_sum[0] if ($score_sum[4]>= 1 && $score_sum[1]<4 );	#t4ss ice = (VirB4/TraU or T4SS)+ int + mob; can not identify G+ ICE if $score_sum[4] >= 2
    $ime_score = $score_sum[0] * ($orit_tag + $score_sum[3]) if (($score_sum[4] <= 0)&&($score_sum[1] == 0)); ## IME = Int + Mob - (virB4)- full T4SS
  	}
		my $aic_score = 0;
		$aic_score = $score_sum[0] * $score_sum[2] * $score_sum[6] if ($score_sum[3]<1);	#aice = int + t4cp + Rep - mob
		my $at4_score = 0;
		$at4_score = $score_sum[0] * $score_sum[2] * ($score_sum[1] + $score_sum[4]) if ($dr_region_len < $length_thred && $score_sum[4] <= 2 && $score_sum[4] > 0);	#t4ss-like ice = (VirB4/TraU or T4SS) + int + t4cp - (full t4ss)
		my $conj_score = ($score_sum[1] + $score_sum[4]) * $score_sum[2] * $score_sum[3]; #$conj_score = (virB4/TraU or T4SS) + t4cp + mob 
		my $mob_score = $score_sum[3] * $score_sum[0]; # mob_score = mob + int
		$score_sum[5] = $orf_score[$dr_position[1]] - $orf_score[$dr_position[0]-1];	#total Score
		if ($ice_score >0){$ice_tag = $ice_score;}
		elsif($aic_score >0){$aic_tag = $aic_score;}
		elsif($ime_score >0){$ime_tag = $ime_score;}
		
	
		if(($ice_score + $aic_score + $ime_score) > 0){
		#check dr quality
		print LOG "Candidate region$i:dr_position$k:0int;1VirB4&traU;2T4CP;3MOB;4T4SS;6Rep;orit_tag;total_score: $score_sum[0] $score_sum[1] $score_sum[2] $score_sum[3] $score_sum[4] $score_sum[6] $orit_tag $score_sum[5]\t";
	  print LOG "candidate_dr_count:$candidate_dr_count;candidate_dr_region:$candidate_dr_up[$candidate_dr_count][0]..$candidate_dr_up[$candidate_dr_count][1];candidate_dr_down_region:$candidate_dr_down[$candidate_dr_count][0]..$candidate_dr_down[$candidate_dr_count][1];temnimal position:$dr_position[0]..$dr_position[1]\tice_score,aic_score,ime_score:$ice_score,$aic_score,$ime_score\n";
		
			my $length;
			##calculate distance to nearest integrase
			$k = 0;
			my $kk = 0;
			while($int_region[$k][0]){
				if($int_region[$k][0] > $candidate_dr_up[$candidate_dr_count][1] && $int_region[$k][1] < $candidate_dr_down[$candidate_dr_count][0]){
					my $l1 = $int_region[$k][0] - $candidate_dr_up[$candidate_dr_count][1];
					my $l2 = $candidate_dr_down[$candidate_dr_count][0] - $int_region[$k][1];
					my $l;
					if($l1 > $l2){
						$l = $l2;
					}else{
						$l = $l1;
					}
					if(!$kk){
						$length = $l;
					}elsif($length > $l){
						$length = $l;
					}
					$kk ++;
				}
				$k ++;
			}
			my $dr_record = 1;
			if($rna_count){
				##if dr locates in tRNA/tmRNA
				if((($candidate_dr_up[$candidate_dr_count][0] < $rna{$can_region_left}) && ($rna_strand{$can_region_left} eq "+")) || (($candidate_dr_down[$candidate_dr_count][1] > $rna_inverse{$can_region_right}) && ($rna_strand{$rna_inverse{$can_region_right}} eq "-"))){
					if($ice_dr[$ice_dr_count][1]){
						my $replace = 0;
						##if has overlap with previous determined DR
						if($candidate_dr_down[$ice_dr[$ice_dr_count][0]][1] > $candidate_dr_up[$candidate_dr_count][0]){
							##select by feature score
							if($score_sum[5] > $ice_dr[$ice_dr_count][1]){
								$replace = 1;
							}elsif($score_sum[5] == $ice_dr[$ice_dr_count][1]){
								##select by ice score
								if($ice_score == 0 && $ice_dr[$ice_dr_count][2] > 0){
									$replace = 2;
								}elsif($ice_score > 0 && $ice_dr[$ice_dr_count][2] == 0){
									$replace = 1;
								}else{
									##select by aice_score
									if($aic_score == 0 && $ice_dr[$ice_dr_count][3] > 0){
										$replace = 2;
									}elsif($aic_score > 0 && $ice_dr[$ice_dr_count][3] == 0){
										$replace = 1;
									}else{
										##select by t4ss-related aice
										if($at4_score == 0 && $ice_dr[$ice_dr_count][4] > 0){
											$replace = 2;
										}elsif($at4_score > 0 && $ice_dr[$ice_dr_count][4] == 0){
											$replace = 1;
										}else{
											## select by ime score 
											if($ime_score == 0 && $ice_dr[$ice_dr_count][6] > 0){ # add for ime
												$replace = 2; 
											}elsif($ime_score > 0 && $ice_dr[$ice_dr_count][6] == 0){ # add for ime
												$replace = 1; 
											}else{ 
												if($length >= $ice_dr[$ice_dr_count][5]){
													$replace = 2;
												}else{
													$replace = 1;
												}
											} 
										}
									}#here
								}
							}else{
								$replace = 2;	##trash
							}
						}
						##if new a ice dr
						if($replace == 0){
							$ice_dr[$ice_dr_count][0] = $candidate_dr_count;	##dr number
							$ice_dr[$ice_dr_count][1] = $score_sum[5];	##feature score
							$ice_dr[$ice_dr_count][2] = $ice_score;	##ice score
							$ice_dr[$ice_dr_count][3] = $aic_score;	##aic_score
							$ice_dr[$ice_dr_count][4] = $at4_score;##at4_score 
							$ice_dr[$ice_dr_count][5] = $length;	##distance to the nearest integrase
							$ice_dr[$ice_dr_count][6] = $ime_score;	##ime score 
							$ice_dr_count ++;
						}elsif($replace == 1){
							$ice_dr_count --;
							$ice_dr[$ice_dr_count][0] = $candidate_dr_count;	##dr number
							$ice_dr[$ice_dr_count][1] = $score_sum[5];	##feature score
							$ice_dr[$ice_dr_count][2] = $ice_score;	##ice score
							$ice_dr[$ice_dr_count][3] = $aic_score;	##conj score
							$ice_dr[$ice_dr_count][4] = $at4_score;##mob score
							$ice_dr[$ice_dr_count][5] = $length;	##distance to the nearest integrase
							$ice_dr[$ice_dr_count][6] = $ime_score;	##ime score 
							$ice_dr_count ++;
						}
					}else{
						##record last as ICE dr
						$ice_dr[$ice_dr_count][0] = $candidate_dr_count;	##dr number
						$ice_dr[$ice_dr_count][1] = $score_sum[5];	##feature score
						$ice_dr[$ice_dr_count][2] = $ice_score;	##ice score
						$ice_dr[$ice_dr_count][3] = $aic_score;	##conj score
						$ice_dr[$ice_dr_count][4] = $at4_score;##mob score
						$ice_dr[$ice_dr_count][5] = $length;	##distance to the nearest integrase
						$ice_dr[$ice_dr_count][6] = $ime_score;	##ime score
						$ice_dr_count ++;
					}
					$dr_record = 0;
				}
			}
			if(($dr_record == 1)&&($candidate_dr_up[$candidate_dr_count][1] - $candidate_dr_up[$candidate_dr_count][0]) <= 200){
				##record last dr
				if(!$select_dr_number){
					$select_dr_score[$select_dr_number][0] = $candidate_dr_count;	##dr number
					$select_dr_score[$select_dr_number][1] = $score_sum[5];	##feature score
					$select_dr_score[$select_dr_number][2] = $ice_score;	##ice score
					$select_dr_score[$select_dr_number][3] = $aic_score;	##conj score
					$select_dr_score[$select_dr_number][4] = $at4_score;##mob score
					$select_dr_score[$select_dr_number][5] = $length;	##distance to the nearest integrase
					$select_dr_score[$select_dr_number][6] = $ime_score;	##ime score 
					$select_dr_number ++;
				}else{
					my $count = 0;
					while($select_dr_score[$count][1]){
						##select by feature score
						if($score_sum[5] > $select_dr_score[$count][1]){
							last;
						}elsif($score_sum[5] == $select_dr_score[$count][1]){
							##select by ice score
							if($ice_score > 0 && $select_dr_score[$count][2] == 0){
								last;
							}else{
								##select by conj score ?? aice score
								if($aic_score > 0 && $select_dr_score[$count][3] == 0){
									last;
								}else{
									##select by mob score ?? at4_score
									if($at4_score > 0 && $select_dr_score[$count][4] == 0){
										last;
									}else{# add for ime
										##select by ime score # add for ime
										if ($ime_score > 0 && $select_dr_score[$count][6] == 0){
											last;
										}else{
											##select by length
											if($length < $select_dr_score[$count][5]){
												last;
											}
										}
									}
								}## here
							}
						}
						$count ++;
					}
					while($select_dr_score[$count][1]){
						my $tmp0 = $select_dr_score[$count][0];
						my $tmp1 = $select_dr_score[$count][1];
						my $tmp2 = $select_dr_score[$count][2];
						my $tmp3 = $select_dr_score[$count][3];
						my $tmp4 = $select_dr_score[$count][4];
						my $tmp5 = $select_dr_score[$count][5];
						my $tmp6 = $select_dr_score[$count][6]; 
						$select_dr_score[$count][0] = $candidate_dr_count;	##dr number
						$select_dr_score[$count][1] = $score_sum[5];	##feature score
						$select_dr_score[$count][2] = $ice_score;	##ice score
						$select_dr_score[$count][3] = $aic_score;	##conj score
						$select_dr_score[$count][4] = $at4_score;##mob score
						$select_dr_score[$count][5] = $length;	##distance to the nearest integrase
						$select_dr_score[$count][6] = $ime_score;	##ime score
						$candidate_dr_count = $tmp0;
						$score_sum[5] = $tmp1;
						$ice_score = $tmp2;
						$aic_score = $tmp3;
						$at4_score = $tmp4;
						$length = $tmp5;
						$ime_score = $tmp6;
						$count ++;
					}
					$select_dr_score[$count][0] = $candidate_dr_count;	##dr number
					$select_dr_score[$count][1] = $score_sum[5];	##feature score
					$select_dr_score[$count][2] = $ice_score;	##ice score
					$select_dr_score[$count][3] = $aic_score;	##conj score
					$select_dr_score[$count][4] = $at4_score;##mob score
					$select_dr_score[$count][5] = $length;	##distance to the nearest integrase
					$select_dr_score[$count][6] = $ime_score;	##ime score 
				}
			}
		}

		##new DR
		$candidate_dr_count ++;
		$candidate_dr_up[$candidate_dr_count][0] = $_;
		$candidate_dr_up[$candidate_dr_count][1] = $dr{$_};
		$candidate_dr_down[$candidate_dr_count][0] = $dr_pair{$_};
		$candidate_dr_down[$candidate_dr_count][1] = $dr{$dr_pair{$_}};
	}
	$dr_num ++;

}


if($dr_count >= 0){
		##calculate score of the regions in last DR pair
		my $k = 1;	##count of posistion
		my @dr_position;	##temnimal position in the DR pair
		my @score_sum = (0, 0, 0, 0, 0, 0);	##$socre_sum[0], int; $socre_sum[1], virB&traU; $socre_sum[2], T4CP; $socre_sum[3], MOB; $socre_sum[4], T4SS; $socre_sum[5], total Score
		##find orf positions in the DR pair
		my $dr_region_len = $candidate_dr_down[$candidate_dr_count][1] - $candidate_dr_up[$candidate_dr_count][0] + 1;
		while($pos_coor[$k][0]){
			if(!$dr_position[0]){
				if($candidate_dr_up[$candidate_dr_count][1] < $pos_coor[$k][0]){
					$dr_position[0] = $k;
					$score_sum[0] += $score[$k][0];	#int
					$score_sum[1] += $score[$k][1];	#virB4&traU, that is T4SS ATPase
					$score_sum[2] += $score[$k][2];	#T4CP
					$score_sum[3] += $score[$k][3];	#MOB
					$score_sum[4] += $score[$k][4];	#T4SS
					$score_sum[6] += $score[$k][6]; #Rep
				}
			}else{
				if($candidate_dr_down[$candidate_dr_count][0] > $pos_coor[$k][1]){
					$dr_position[1] = $k;
					$score_sum[0] += $score[$k][0];	#int
					$score_sum[1] += $score[$k][1];	#virB4&traU, that is T4SS ATPase
					$score_sum[2] += $score[$k][2];	#T4CP
					$score_sum[3] += $score[$k][3];	#MOB
					$score_sum[4] += $score[$k][4];	#T4SS
					$score_sum[6] += $score[$k][6]; #Rep
				}else{
					last;
				}
			}
			$k++;
		}
		##calculate feature sores
		if($orit_left > $candidate_dr_up[$candidate_dr_count][0] && $orit_coordinate{$orit_left} < $candidate_dr_down[$candidate_dr_count][1])
		{$orit_tag = 1;}else{$orit_tag = 0;}
		my $ice_score = 0;
		my $ime_score = 0;
		if($gram_postive_tag == 0){
		$ice_score = ($score_sum[1] + $score_sum[4]) * $score_sum[2] * $score_sum[3] * $score_sum[0] if ($score_sum[4] >2 && $score_sum[1]<5);	#t4ss ice = (VirB4/TraU or T4SS)+ int + mob	
    $ime_score = $score_sum[0] * ($orit_tag + $score_sum[3]) if (($score_sum[4] < 2)&&($score_sum[1] == 0)); ## IME = Int + Mob - (virB4/TraU )- full T4SS}
    }
    else{ 
		$ice_score = ($score_sum[1] + $score_sum[4]) * $score_sum[2] * $score_sum[3] * $score_sum[0] if ($score_sum[4] >= 1 && $score_sum[1]<5);	#t4ss ice = (VirB4/TraU or T4SS)+ int + mob; can not identify G+ ICE if $score_sum[4] >= 2
    $ime_score = $score_sum[0] * ($orit_tag + $score_sum[3]) if (($score_sum[4] < 1)&&($score_sum[1] == 0)); ## IME = Int + Mob - (virB4/TraU )- full T4SS}
  	}
    my $aic_score = 0;
		$aic_score = $score_sum[0] * $score_sum[2] * $score_sum[6] if ($score_sum[3]<1);	#aice = int + t4cp + Rep - mob 
		my $at4_score = 0;
		$at4_score = $score_sum[0] * $score_sum[2] * ($score_sum[1] + $score_sum[4]) if ($dr_region_len < $length_thred && $score_sum[4] <= 2 && $score_sum[4] > 0);	#T4ss-related AICE
		$score_sum[5] = $orf_score[$dr_position[1]] - $orf_score[$dr_position[0]-1];	#total Score
		if ($ice_score >0){$ice_tag = $ice_score;}
		elsif($aic_score >0){$aic_tag = $aic_score;}
		elsif($ime_score >0){$ime_tag = $ime_score;}
		
		
		if(($ice_score + $aic_score + $ime_score) > 0){
			#check last dr quality
		print LOG "Candidate region$i:last DR: dr_position$k:0int;1VirB4&traU;2T4CP;3MOB;4T4SS;6Rep;orit_tag: $score_sum[0] $score_sum[1] $score_sum[2] $score_sum[3] $score_sum[4] $score_sum[6] $orit_tag $score_sum[5] \t";
	  print LOG "candidate_dr_count:last DR: $candidate_dr_count;candidate_dr_region:$candidate_dr_up[$candidate_dr_count][0]..$candidate_dr_up[$candidate_dr_count][1];candidate_dr_down_region:$candidate_dr_down[$candidate_dr_count][0]..$candidate_dr_down[$candidate_dr_count][1];temnimal position:$dr_position[0]..$dr_position[1]\tice_score,aic_score,ime_score:$ice_score, $aic_score,$ime_score\n";
	
			my $length;
			##calculate distance to nearest integrase
			$k = 0;
			my $kk = 0;
			while($int_region[$k][0]){
				if($int_region[$k][0] > $candidate_dr_up[$candidate_dr_count][1] && $int_region[$k][1] < $candidate_dr_down[$candidate_dr_count][0]){
					my $l1 = $int_region[$k][0] - $candidate_dr_up[$candidate_dr_count][1];
					my $l2 = $candidate_dr_down[$candidate_dr_count][0] - $int_region[$k][1];
					my $l;
					if($l1 > $l2){
						$l = $l2;
					}else{
						$l = $l1;
					}
					if(!$kk){
						$length = $l;
					}elsif($length > $l){
						$length = $l;
					}
					$kk ++;
				}
				$k ++;
			}
			my $dr_record = 1;
			if($rna_count){
				##if dr locates in tRNA/tmRNA
				if((($candidate_dr_up[$candidate_dr_count][0] < $rna{$can_region_left}) && ($rna_strand{$can_region_left} eq "+")) || (($candidate_dr_down[$candidate_dr_count][1] > $rna_inverse{$can_region_right}) && ($rna_strand{$rna_inverse{$can_region_right}} eq "-"))){
					if($ice_dr[$ice_dr_count][1]){
						my $replace = 0;
						##if has overlap with previous determined DR
						if($candidate_dr_down[$ice_dr[$ice_dr_count][0]][1] > $candidate_dr_up[$candidate_dr_count][0]){
							##select by feature score
							if($score_sum[5] > $ice_dr[$ice_dr_count][1]){
								$replace = 1;
							}elsif($score_sum[5] == $ice_dr[$ice_dr_count][1]){
								##select by ice score
								if($ice_score == 0 && $ice_dr[$ice_dr_count][2] > 0){
									$replace = 2;
								}elsif($ice_score > 0 && $ice_dr[$ice_dr_count][2] == 0){
									$replace = 1;
								}else{
									##select by conj score
									if($aic_score == 0 && $ice_dr[$ice_dr_count][3] > 0){
										$replace = 2;
									}elsif($aic_score > 0 && $ice_dr[$ice_dr_count][3] == 0){
										$replace = 1;
									}else{
										##select by mob score
										if($at4_score == 0 && $ice_dr[$ice_dr_count][4] > 0){
											$replace = 2;
										}elsif($at4_score > 0 && $ice_dr[$ice_dr_count][4] == 0){
											$replace = 1;
										}else{
											## select by ime score 
											if($ime_score == 0 && $ice_dr[$ice_dr_count][6] > 0){ 
												$replace = 2; 
											}elsif($ime_score > 0 && $ice_dr[$ice_dr_count][6] == 0){ 
												$replace = 1; 
											}else{ 
												if($length >= $ice_dr[$ice_dr_count][5]){
													$replace = 2;
												}else{
													$replace = 1;
												}
											}
										}
									}
								}
							}else{
								$replace = 2;	##trash
							}
						}
						##if new a ice dr
						if($replace == 0){
							$ice_dr[$ice_dr_count][0] = $candidate_dr_count;	##dr number
							$ice_dr[$ice_dr_count][1] = $score_sum[5];	##feature score
							$ice_dr[$ice_dr_count][2] = $ice_score;	##ice score
							$ice_dr[$ice_dr_count][3] = $aic_score;	##conj score
							$ice_dr[$ice_dr_count][4] = $at4_score;##mob score
							$ice_dr[$ice_dr_count][5] = $length;	##distance to the nearest integrase
							$ice_dr[$ice_dr_count][6] = $ime_score;	##ime score 
							$ice_dr_count ++;
						}elsif($replace == 1){
							$ice_dr_count --;
							$ice_dr[$ice_dr_count][0] = $candidate_dr_count;	##dr number
							$ice_dr[$ice_dr_count][1] = $score_sum[5];	##feature score
							$ice_dr[$ice_dr_count][2] = $ice_score;	##ice score
							$ice_dr[$ice_dr_count][3] = $aic_score;	##conj score
							$ice_dr[$ice_dr_count][4] = $at4_score;##mob score
							$ice_dr[$ice_dr_count][5] = $length;	##distance to the nearest integrase
							$ice_dr[$ice_dr_count][6] = $ime_score;	##ime score 
							$ice_dr_count ++;
						}
					}else{
						##record last as ICE dr
						$ice_dr[$ice_dr_count][0] = $candidate_dr_count;	##dr number
						$ice_dr[$ice_dr_count][1] = $score_sum[5];	##feature score
						$ice_dr[$ice_dr_count][2] = $ice_score;	##ice score
						$ice_dr[$ice_dr_count][3] = $aic_score;	##conj score
						$ice_dr[$ice_dr_count][4] = $at4_score;##mob score
						$ice_dr[$ice_dr_count][5] = $length;	##distance to the nearest integrase
						$ice_dr[$ice_dr_count][6] = $ime_score;	##ime score 
						$ice_dr_count ++;
					}
					$dr_record = 0;
				}
			}
			if(($dr_record == 1)&&($candidate_dr_up[$candidate_dr_count][1] - $candidate_dr_up[$candidate_dr_count][0]) <= 200){
				##record last dr
				if(!$select_dr_number){
					$select_dr_score[$select_dr_number][0] = $candidate_dr_count;	##dr number
					$select_dr_score[$select_dr_number][1] = $score_sum[5];	##feature score
					$select_dr_score[$select_dr_number][2] = $ice_score;	##ice score
					$select_dr_score[$select_dr_number][3] = $aic_score;	##conj score
					$select_dr_score[$select_dr_number][4] = $at4_score;##mob score
					$select_dr_score[$select_dr_number][5] = $length;	##distance to the nearest integrase
					$select_dr_score[$select_dr_number][6] = $ime_score;	##ime score 
					$select_dr_number ++;
				}else{
					my $count = 0;
					while($select_dr_score[$count][1]){
						##select by feature score
						if($score_sum[5] > $select_dr_score[$count][1]){
							last;
						}elsif($score_sum[5] == $select_dr_score[$count][1]){
							##select by ice score
							if($ice_score > 0 && $select_dr_score[$count][2] == 0){
								last;
							}else{
								##select by conj score
								if($aic_score > 0 && $select_dr_score[$count][3] == 0){
									last;
								}else{
									##select by conj score
									if($at4_score > 0 && $select_dr_score[$count][4] == 0){
										last
									}else{	
										##select by ime score 
										if ($ime_score > 0 && $select_dr_score[$count][6] == 0){ # add for ime
											last;
										}else{
											##select by length
											if($length < $select_dr_score[$count][5]){
												last;
											}
										}
									}
								}
							}
						}
						$count ++;
					}
					while($select_dr_score[$count][1]){
						my $tmp0 = $select_dr_score[$count][0];
						my $tmp1 = $select_dr_score[$count][1];
						my $tmp2 = $select_dr_score[$count][2];
						my $tmp3 = $select_dr_score[$count][3];
						my $tmp4 = $select_dr_score[$count][4];
						my $tmp5 = $select_dr_score[$count][5];
						my $tmp6 = $select_dr_score[$count][6];
						$select_dr_score[$count][0] = $candidate_dr_count;	##dr number
						$select_dr_score[$count][1] = $score_sum[5];	##feature score
						$select_dr_score[$count][2] = $ice_score;	##ice score
						$select_dr_score[$count][3] = $aic_score;	##conj score
						$select_dr_score[$count][4] = $at4_score;##mob score
						$select_dr_score[$count][5] = $length;	##distance to the nearest integrase
						$select_dr_score[$count][6] = $ime_score;	##ime score
						$candidate_dr_count = $tmp0;
						$score_sum[5] = $tmp1;
						$ice_score = $tmp2;
						$aic_score = $tmp3;
						$at4_score = $tmp4;
						$length = $tmp5;
						$ime_score = $tmp6;
						$count ++;
					}
					$select_dr_score[$count][0] = $candidate_dr_count;	##dr number
					$select_dr_score[$count][1] = $score_sum[5];	##feature score
					$select_dr_score[$count][2] = $ice_score;	##ice score
					$select_dr_score[$count][3] = $aic_score;	##conj score
					$select_dr_score[$count][4] = $at4_score;##mob score
					$select_dr_score[$count][5] = $length;	##distance to the nearest integrase
					$select_dr_score[$count][6] = $ime_score;	##ime score
				}
			}
		}
}

my $select_dr_score_count = 0;
while($select_dr_score[$select_dr_score_count][1]){
	my $k = 0;
	my $l1 = $candidate_dr_up[$select_dr_score[$select_dr_score_count][0]][0];
	my $l2 = $candidate_dr_down[$select_dr_score[$select_dr_score_count][0]][1];
	my $trash = 0;
	while($ice_dr[$k][1]){
		my $i1 = $candidate_dr_up[$ice_dr[$k][0]][0];
		my $i2 = $candidate_dr_down[$ice_dr[$k][0]][1];
		if(($l1 > $i1 && $l1 < $i2) || ($l2 > $i1 && $l2 < $i2)){
			$trash = 1;
			last;
		}
		$k ++;
	}
	if(!$trash){
		$ice_dr[$ice_dr_count][0] = $select_dr_score[$select_dr_score_count][0];	##dr number
		$ice_dr[$ice_dr_count][1] = $select_dr_score[$select_dr_score_count][1];	##feature score
		$ice_dr[$ice_dr_count][2] = $select_dr_score[$select_dr_score_count][2];	##ice score
		$ice_dr[$ice_dr_count][3] = $select_dr_score[$select_dr_score_count][3];	##aic score
		$ice_dr[$ice_dr_count][4] = $select_dr_score[$select_dr_score_count][4];	##at4 score
		$ice_dr[$ice_dr_count][5] = $select_dr_score[$select_dr_score_count][5];	##distance to the nearest integrase
		$ice_dr[$ice_dr_count][6] = $select_dr_score[$select_dr_score_count][6];  ##ime score
		$ice_dr_count ++;
	}
		$select_dr_score_count ++;
}

##select inter ice region
my @inter_ice;
my $inter_ice_count = 0;
$inter_ice[$inter_ice_count][0] = $can_region_left;
my $ice_count = 0;
while($ice_count < $ice_dr_count){
	$inter_ice[$inter_ice_count][1] = $candidate_dr_up[$ice_dr[$ice_count][0]][0];
	$inter_ice_count ++; 
	$inter_ice[$inter_ice_count][0] = $candidate_dr_down[$ice_dr[$ice_count][0]][1];
	$ice_count ++;
}
my $e1=0;
my $e2=0;

my $genome_seqio_obj = Bio::SeqIO->new(-file => "$seq_fna", -format => "fasta" );
my $genome_seq_obj = $genome_seqio_obj->next_seq;

## genus is used to exclude some non Actinobacteria
my $genome_desc = $genome_seq_obj->desc;
my $genus_gbk = (split(" ",$genome_desc))[1]; 
my $genus_RDP = "";
if (-e $RDP_parse_out){
	open RDP, "$RDP_parse_out";
	my $RDP_line1= readline(RDP);
	close RDP;
	chomp($RDP_line1);
	$genus_RDP = (split(" ",$RDP_line1))[5];
}

my $string = "";
my $gc_content = 0;

my @ptt2_array;
my $tra_num = 0; ## for AICE: FtsK_SpoIIIE
my $marker_name = "";

##define ICEs
if(($ice_dr_count>0)){ ## with dr;
	my $e1 = 0;
	my $k = 0;
	my @ice_region;
	
	##define ICE with dr
	while($k < $ice_dr_count){
		my $ice_left = $candidate_dr_up[$ice_dr[$k][0]][0];
		my $ice_right = $candidate_dr_down[$ice_dr[$k][0]][1];
		$ice_region[$e1][0] = $ice_left;
		$ice_region[$e1][1] = $ice_right;
		my $ice_desc = "";
		my $dr_desc = "$candidate_dr_up[$ice_dr[$k][0]][0]..$candidate_dr_up[$ice_dr[$k][0]][1]\t$candidate_dr_down[$ice_dr[$k][0]][0]..$candidate_dr_down[$ice_dr[$k][0]][1]";
		my $insert_desc = "-";
		my $orit_desc = "";

		$insert_desc = $rna_name{$can_region_left} if $ice_left <= $rna{$can_region_left};
		$insert_desc = $rna_name{$rna_inverse{$can_region_right}} if $ice_right >= $rna_inverse{$can_region_right};
		if($insert_desc eq "-"){
			foreach (sort {$a<=>$b} keys %orf){
				if(($_ <= $ice_left) && ($orf{$_} >= $ice_left) || ($_ <= $ice_right) && ($orf{$_} >= $ice_right)){
					$insert_desc = $ptt{$_};
					last;
				}
			}
		}
		if($insert_desc eq "-"){
		$insert_desc = "-\n";}
		my $ice_out = $download_dir."ICEfinder_$i"."_"."0\n"; ## final result
		my $ice_out2 = $download_dir."ICEfinder_$i"."_"."0.core\n"; ## core element
		my $ice_len = $ice_right - $ice_left +1;
		my $string = $genome_seq_obj->subseq($ice_left,$ice_right);
    
   
    ## oriT desc
    if( $orit_left >= $ice_left && $orit_coordinate{$orit_left} <= $ice_right ){
    	$orit_desc = "$orit_left..$orit_coordinate{$orit_left}"; 
    	$orit_tag = 1;
    }else{
    	$orit_desc = "-";
    	$orit_tag = 0;}

    
    	
    my $int_tag = 0;	
    my $mob_tag = 0;
    my %t4ss_pro_left_name;
    my %t4ss_pro_left_right;
    foreach (sort {$a<=>$b} keys %orf){                            
		    	my $orf_key = chomp($_);     
		    	if(($_ >= $ice_left) && ($orf{$_} <= $ice_right)){                                             
		    		my @ptt2_array_tmp = (split (/\t/,$ptt2{$_}));
		    		$marker_name = $ptt2_array_tmp[5]; 
		    		my $sign_pro_name = $ptt2_array_tmp[6];
		    		if ($marker_name eq "Integrase" ){
		    			$int_tag += 1;
		    		}
		    		if ($marker_name eq "Relaxase" ){
		    			$mob_tag += 1;
		    		}
		    		if ($marker_name eq "T4SS" ){
		    			$t4ss_pro_left_name{$_} = $sign_pro_name;
		    			$t4ss_pro_left_right{$_} = $orf{$_};
		    		}
		    		  		
		    	}		    				    						 		                                                                                           
		}
    ## t4ss co-localization
    my $t4ss_tag = 0; 
		my @region;
		my @region_name;
		my $i_t4ss = 0;
		my $region_num = 0;
		my $mpf_distance = 10000; ## the distance between ecch Mpf gene <= 10 k 
    my $t4ss_core = 2;
		if($gram_postive_tag == 0){## strict restrication of t4ss core components for G- T4SS 
			$t4ss_core = 2;
			}else{ ## less strict restrication of t4ss core components for G- T4SS
				$t4ss_core = 1;
		}
		if($gram_postive_tag == 1 && $rep_tag >= 1){ ## more strict restrication of t4ss core components for AICE-like bacteria
			$t4ss_core = 5;
		}
		
		foreach (sort {$a<=>$b} keys %t4ss_pro_left_right){
			if($i_t4ss == 0){
                $region[$region_num][0] = $_;
                $region[$region_num][1] = $t4ss_pro_left_right{$_};
                $region_name[$region_num] = $t4ss_pro_left_name{$_};
                $i_t4ss ++;
        }else{
                if($_ <= ($region[$region_num][1]+ $distance)){
                        $region[$region_num][1] = $t4ss_pro_left_right{$_};
                        $region_name[$region_num] .= "\t$t4ss_pro_left_name{$_}";
                        $i_t4ss ++;
                }else{
                        $region_num ++;
                        $region[$region_num][0] = $_;
                        $region[$region_num][1] = $t4ss_pro_left_right{$_};
                        $region_name[$region_num] = $t4ss_pro_left_name{$_};      
                        $i_t4ss = 1;
                }
        }
			
			}
			
		my $core_num = $#region_name + 1;
		if ($core_num >= $t4ss_core){
			$t4ss_tag = 1;
			}
		
		## gc diff
		my $gc_diff_tag = 0;
		my $gc_content = &calculate_gc($string); 
		my $gc_diff = $gc_content - $genome_gc;
		$gc_diff = abs($gc_diff);

		if ($aic_tag >0){ ## actually no gc deviation for AICE
			if($gc_diff >= 1){$gc_diff_tag = 1;}
			else{$gc_diff_tag = 0;}
		}
		if($ice_tag >0){
			if($ice_len <250000){ #strict gc deviation for T4SS-type ICE with the size over 250 kb
				$gc_diff_tag = 1;
			}else{
				if($gc_diff >= 1){
					$gc_diff_tag = 1;
				}else{$gc_diff_tag = 0;}}
		}
		if($ime_tag >0){
			if($ice_len <50000){ #strict gc deviation for IME with the size over 50 kb
				$gc_diff_tag = 1;
			}else{
				if($gc_diff >= 1){
					$gc_diff_tag = 1;
				}else{$gc_diff_tag = 0;}}
		}
	#print LOG "Candidate region$i: withDR:0int;1VirB4;2T4CP;3MOB;4T4SS;6Rep;orit_tag;total_score: $score_sum[0] $score_sum[1] $score_sum[2] $score_sum[3] $score_sum[4] $score_sum[6] $orit_tag $score_sum[5] \n";
	print LOG "\nCandidate region$i: withDR: dr_desc: $dr_desc;ice:$ice_left..$ice_right\n";
	print LOG "e1:$e1	ice_tag;aic_tag;ime_tag:$ice_tag,$aic_tag,$ime_tag\n"; 
	print LOG "Candidate region$i:int_tag;mob_tag;t4ss_tag;gc_diff_tag;ice_len:$int_tag;$mob_tag;$t4ss_tag;$gc_diff_tag;$ice_len\n";				
	if ($int_tag >=1){     
			if($ice_tag > 0 && $mob_tag >0 && $t4ss_tag >0 && $gc_diff_tag >0 && $ice_len < 600000 && $ice_len> 4000){ 
				$ice_desc = "Putative ICE with T4SS";
				open ICEOUT, ">$ice_out";                                      
		    open ICEOUT2,">$ice_out2";                                                          
		    print ICEOUT "Description: $ice_desc\n".                                    
		                 "Location: $ice_left..$ice_right\n".
		                 "oriT: $orit_desc\n".              
		                 "Length: $ice_len bp\n".  
		                 "GC content: $gc_content %\n".                    
		                 "DR: $dr_desc\n".                                 
		                 "Insert: $insert_desc";                           
		    foreach (sort {$a<=>$b} keys %orf){                            
		    	my $orf_key = chomp($_);    
		    		if(($_ >= $ice_left) && ($orf{$_} <= $ice_right)){         
		    			print ICEOUT $ptt{$_};                                     
		    			if($ptt2{$_}){
				 				print ICEOUT2 "$job_id\t$ptt2{$_}";
				 			}                                  
		    		}	                                                           
		    }                                                              
		    close ICEOUT;                                                  
		    close ICEOUT2;                                                                     
			}
			elsif($ime_tag >0 && ($mob_tag >0 ||$orit_tag>0) && $t4ss_tag == 0 && $gc_diff_tag >0 &&  $ice_len < 90000 && $ice_len> 1000){   
			if($h >0){
				my $t4ss_n = 0;
				while($t4ss_region[$t4ss_n][0]){	
				##find neighbour integrase
					my $ice_left_new = $t4ss_region[$t4ss_n][0];
					my $ice_right_new = $t4ss_region[$t4ss_n][1];
					my $int_n = 0;
					while($int_region[$int_n][0]){
						if(($t4ss_region[$t4ss_n][0] - $int_region[$int_n][1]) < 50000){
							$ice_left_new = $int_region[$int_n][0] if $int_region[$int_n][0] < $ice_left_new;
							last;
						}
						$int_n ++;
					}
					while($int_region[$int_n][0]){
						if(($int_region[$int_n][0] - $t4ss_region[$t4ss_n][1]) < 50000){
							$ice_right_new = $int_region[$int_n][1] if $int_region[$int_n][1] > $ice_right_new;
							$int_n ++;
						}else{
							last;
						}
					}
				#}
				if( $ice_left_new < $ice_left || $ice_right_new > $ice_right){## T4SS + int 
					$ice_left = $ice_left_new;
					$ice_right = $ice_right_new;
					$ice_desc = "Putative ICE with T4SS";
					$dr_desc = "-";
					$insert_desc = "-\n";
					open ICEOUT, ">$ice_out";                                      
		    	open ICEOUT2,">$ice_out2";                                                          
		   		print ICEOUT "Description: $ice_desc\n".                                    
		                 "Location: $ice_left..$ice_right\n".
		                 "oriT: $orit_desc\n".              
		                 "Length: $ice_len bp\n".  
		                 "GC content: $gc_content %\n".                    
		                 "DR: $dr_desc\n".                                 
		                 "Insert: $insert_desc";                           
		    	foreach (sort {$a<=>$b} keys %orf){                            
		    	my $orf_key = chomp($_);    
		    		if(($_ >= $ice_left) && ($orf{$_} <= $ice_right)){         
		    			print ICEOUT $ptt{$_};                                     
		    			if($ptt2{$_}){
				 				print ICEOUT2 "$job_id\t$ptt2{$_}";
				 			}                                  
		    		}	                                                           
		    	}                                                              
		    	close ICEOUT;                                                  
		   		close ICEOUT2;  
				}
				else{
					$ice_desc = "Putative IME";
					open ICEOUT, ">$ice_out";
					open ICEOUT2,">$ice_out2";                                  
					print ICEOUT "Description: $ice_desc\n".                       
				              "Location: $ice_left..$ice_right\n".
				              "oriT: $orit_desc\n". 
				              "Length: $ice_len bp\n".  
				              "GC content: $gc_content % \n".           
				              "DR: $dr_desc\n". 				                                 
				              "Insert: $insert_desc";              
				 	foreach (sort {$a<=>$b} keys %orf){                          
				 		my $orf_key = chomp($_);  
				 		if(($_ >= $ice_left) && ($orf{$_} <= $ice_right)){        
				 			print ICEOUT $ptt{$_}; 
				 			if($ptt2{$_}){
				 				print ICEOUT2 "$job_id\t$ptt2{$_}";
				 			}                                                        
				 		}	                                                          
				 	}                                                            
				 	close ICEOUT;
				 	close ICEOUT2;                                                 
					}
				$t4ss_n++;
				}
			}else{## if h>0
				$ice_desc = "Putative IME";
					open ICEOUT, ">$ice_out";
					open ICEOUT2,">$ice_out2";                                  
					print ICEOUT "Description: $ice_desc\n".                       
				              "Location: $ice_left..$ice_right\n".
				              "oriT: $orit_desc\n". 
				              "Length: $ice_len bp\n".  
				              "GC content: $gc_content % \n".           
				              "DR: $dr_desc\n". 				                                 
				              "Insert: $insert_desc";              
				 	foreach (sort {$a<=>$b} keys %orf){                          
				 		my $orf_key = chomp($_);  
				 		if(($_ >= $ice_left) && ($orf{$_} <= $ice_right)){        
				 			print ICEOUT $ptt{$_}; 
				 			if($ptt2{$_}){
				 				print ICEOUT2 "$job_id\t$ptt2{$_}";
				 			}                                                        
				 		}	                                                          
				 	}                                                            
				 	close ICEOUT;
				 	close ICEOUT2;                                                 
					
				}
			}
			elsif($aic_tag > 0 && $ice_len < 60000 && $gram_postive_tag >0){     
				if(($genus_gbk =~ /Anaerolinea|Bradyrhizobium|Burkholderia|Candidatus|Clostridium|Desulfarculus|Desulfomonile|Enterococcus|Erysipelothrix|Erysipelotrichaceae|Lactobacillus|Lactococcus|Mariprofundus|Marivivens|Mesoplasma|Neisseria|Nitrosomonas|Paenibacillus|Paracoccus|Pseudomonas|Ruminiclostridium|Sphingobium|Staphylococcus|Streptococcus|Xanthomonas|Xylella/i)||($genus_RDP =~ /Anaerolinea|Bradyrhizobium|Burkholderia|Candidatus|Clostridium|Desulfarculus|Desulfomonile|Enterococcus|Erysipelothrix|Erysipelotrichaceae|Lactobacillus|Lactococcus|Mariprofundus|Marivivens|Mesoplasma|Neisseria|Nitrosomonas|Paenibacillus|Paracoccus|Pseudomonas|Ruminiclostridium|Sphingobium|Staphylococcus|Streptococcus|Xanthomonas|Xylella/i)){
				print LOG "This is not in Actinobacteria: genus_gbk:$genus_gbk;genus_RDP:$genus_RDP\n";
			}else{
				$ice_desc = "Putative  AICE with Rep and Tra";
				$orit_desc = "-"; 
				$tra_num = 0;
				foreach (sort {$a<=>$b} keys %orf){                                                                 	
				 		my $orf_key = chomp($_);                                       	
				 		if(($_ >= $ice_left) && ($orf{$_} <= $ice_right)){                                           	
	 						@ptt2_array = split(/\t/,$ptt2{$_});
	 						chomp($ptt2_array[6]); 
	 						if($ptt2_array[6] eq "FtsK_SpoIIIE" ){
	 							$tra_num +=1;
	 						}                                                                         	
						}	
				}
				
				if($tra_num >0){                                                                         
				open ICEOUT, ">$ice_out";  
				open ICEOUT2,">$ice_out2";                                                                          	
				print ICEOUT "Description: $ice_desc\n".                                                                          	
				              "Location: $ice_left..$ice_right\n".
				              "oriT: $orit_desc\n".                                                   	
				              "Length: $ice_len bp\n".
				              "GC content: $gc_content %\n".                                                               	
				              "DR: $dr_desc\n".                                                                      	
				              "Insert: $insert_desc";                                                                	
				 foreach (sort {$a<=>$b} keys %orf){                                                                 	
				 		my $orf_key = chomp($_);                                           	
				 		if(($_ >= $ice_left) && ($orf{$_} <= $ice_right)){                                               	
				 		print ICEOUT $ptt{$_}; 
				 		if($ptt2{$_}){
				 		print ICEOUT2 "$job_id\t$ptt2{$_}";
				 		}                                                                          	
				 	}	                                                                                                 	
				 }                                                                                                   	
				 close ICEOUT;
				 close ICEOUT2; 
				}                                                                                       	
			}}                                                                                                    			
	}	
		$e1 ++;
		$k ++;
	
	}
	
}

else{	##no dr
	my $j = 0;
	my $e2 = 0;
	my $dr_desc = "-";
	my $insert_desc = "-\n"; ## no DR && no tRNA
	my $orit_desc = "";
	while($conj_region[$j][0]){	
		##find neighbour integrase
		my $ice_left = $conj_region[$j][0];
		my $ice_right = $conj_region[$j][1];
		my $int_n = 0;
		while($int_region[$int_n][0]){
			if(($conj_region[$j][0] - $int_region[$int_n][1]) < 50000){
				$ice_left = $int_region[$int_n][0] if $int_region[$int_n][0] < $ice_left;
				last;
			}
			$int_n ++;
		}
		while($int_region[$int_n][0]){
			if(($int_region[$int_n][0] - $conj_region[$j][1]) < 50000){
				$ice_right = $int_region[$int_n][1] if $int_region[$int_n][1] > $ice_right;
				$int_n ++;
			}else{
				last;
			}
		}
		##define the ICE
		my $k = 1;
		my $left_position = 0;
		my $right_position = 0;
		my @score_sum = (0, 0, 0, 0, 0, 0);
		while($pos_coor[$k][0]){
			if(!$left_position){
				if($ice_left == $pos_coor[$k][0]){
					$left_position = $k;
					$score_sum[0] += $score[$k][0];	#int
					$score_sum[1] += $score[$k][1];	#virB&traU
					$score_sum[2] += $score[$k][2];	#T4CP
					$score_sum[3] += $score[$k][3];	#MOB
					$score_sum[4] += $score[$k][4];	#T4SS
					$score_sum[6] += $score[$k][6]; #Rep 
				}
			}else{
				if($ice_right >= $pos_coor[$k][1]){
					$right_position = $k;
					$score_sum[0] += $score[$k][0];	#int
					$score_sum[1] += $score[$k][1];	#virB&traU
					$score_sum[2] += $score[$k][2];	#T4CP
					$score_sum[3] += $score[$k][3];	#MOB
					$score_sum[4] += $score[$k][4];	#T4SS
					$score_sum[6] += $score[$k][6]; #Rep 
				}else{
					last;
				}
			}
			$k++;
		}

		if($orit_left >= $ice_left && $orit_coordinate{$orit_left} <= $ice_right)
		{$orit_tag = 1;}else{$orit_tag = 0;} 

		##calculate feature sores
		my $ice_score = 0;
		my $ice_tag = 0;
		$ice_score = ($score_sum[1] + $score_sum[4]) * $score_sum[2] * $score_sum[3] * $score_sum[0] if $score_sum[4] > 2;	#t4ss ice = (VirB4/TraU or T4SS)+ int + t4cp + mob;
		my $aic_score = 0;
		my $aic_tag = 0;
		$aic_score = $score_sum[0] * $score_sum[2] * $score_sum[6] if ($score_sum[3]<1);	#aice = int + t4cp + Rep - mob 
		my $ime_score =0;
		my $ime_tag = 0;
    $ime_score = $score_sum[0] * ($orit_tag + $score_sum[3]) if (($score_sum[4] <= 2)&&($score_sum[1] == 0)); ## IME = Int + Mob - (virB4/TraU )- full T4SS
    if ($ice_score >0){$ice_tag = $ice_score;}
		elsif($aic_score >0){$aic_tag = $aic_score;}
		elsif($ime_score >0){$ime_tag = $ime_score;}		
		my $conj_score = ($score_sum[1] + $score_sum[4]) * $score_sum[2] * $score_sum[3]; #$conj_score = (virB4/TraU or T4SS) + t4cp + mob
		my $mob_score = $score_sum[3] * $score_sum[0];  # mob_score = mob + int
		$score_sum[5] = $score[$right_position][5] - $score[$left_position-1][5];	#total Score
    my $ice_len = $ice_right - $ice_left +1;
    my $string = $genome_seq_obj->subseq($ice_left,$ice_right);
    my $gc_content = &calculate_gc($string); ## liu added 
		## oriT desc
    if($orit_left >= $ice_left && $orit_coordinate{$orit_left} <= $ice_right ){##$orit_tag >0 && $orit_left > $ice_left && $orit_coordinate{$orit_left} < $ice_right 
    	$orit_desc = "$orit_left..$orit_coordinate{$orit_left}"; 
    }else{$orit_desc = "-";}

		my $int_tag = 0;	
    my $mob_tag = 0;
    my %t4ss_pro_left_name;
    my %t4ss_pro_left_right;
    foreach (sort {$a<=>$b} keys %orf){                            
		    	my $orf_key = chomp($_);     
		    	if(($_ >= $ice_left) && ($orf{$_} <= $ice_right)){                                             
		    		my @ptt2_array_tmp = (split (/\t/,$ptt2{$_}));
		    		$marker_name = $ptt2_array_tmp[5]; 
		    		my $sign_pro_name = $ptt2_array_tmp[6]; 
	
		    		if ($marker_name eq "Integrase" ){
		    			$int_tag += 1;
		    		}
		    		if ($marker_name eq "Relaxase" ){
		    			$mob_tag += 1;
		    		}
		    		if ($marker_name eq "T4SS" ){ ## MPF gene cluster ;does not contain T4SS ATPase
		    			$t4ss_pro_left_name{$_} = $sign_pro_name;
		    			$t4ss_pro_left_right{$_} = $orf{$_};
		    		}
		    		
		    	}		    				    						 		                                                                                           
		}
    ## t4ss co-localization
    my $t4ss_tag = 0; 
		my @region;
		my @region_name;
		my $i_t4ss = 0;
		my $region_num = 0;
		my $mpf_distance = 10000; ## the distance between ecch Mpf gene<= 10 k 
		my $t4ss_core = 2; ## 5
		if($gram_postive_tag == 0){## strict restrication of t4ss core components for G- T4SS
			$t4ss_core = 2;
			}else{ ## less strict restrication of t4ss core components for G+ T4SS
				$t4ss_core = 1;
		}
		if($gram_postive_tag == 1 && $rep_tag >= 1){ ## more strict restrication of t4ss core components for AICE-like bacteria
			$t4ss_core = 5;
		}
		foreach (sort {$a<=>$b} keys %t4ss_pro_left_right){
			if($i_t4ss == 0){
                $region[$region_num][0] = $_;
                $region[$region_num][1] = $t4ss_pro_left_right{$_};
                $region_name[$region_num] = $t4ss_pro_left_name{$_};
                $i_t4ss ++;
        }else{
                if($_ <= ($region[$region_num][1]+ $distance)){
                        $region[$region_num][1] = $t4ss_pro_left_right{$_};
                        $region_name[$region_num] .= "\t$t4ss_pro_left_name{$_}";
                        $i_t4ss ++;
                }else{
                        $region_num ++;
                        $region[$region_num][0] = $_;
                        $region[$region_num][1] = $t4ss_pro_left_right{$_};
                        $region_name[$region_num] = $t4ss_pro_left_name{$_};
                        $i_t4ss = 1;
                }
        }
			
			}
			
		my $core_num = $#region_name + 1;
		if ($core_num >= $t4ss_core){
			$t4ss_tag = 1;
			}
		
		
		## gc diff
		my $gc_diff_tag = 0;
		my $gc_content = &calculate_gc($string); 
		my $gc_diff = $gc_content - $genome_gc;
		$gc_diff = abs($gc_diff);
		
		if($ice_tag >0){
			if($ice_len <250000){ #strict gc deviation for T4SS-type ICE with the size over 250 kb
				$gc_diff_tag = 1;
			}else{
				if($gc_diff >= 1){
					$gc_diff_tag = 1;
				}else{$gc_diff_tag = 0;}}
		}
		if($ime_tag >0){
			if($ice_len <50000){ #strict gc deviation for IME with the size over 50 kb
				$gc_diff_tag = 1;
			}else{
				if($gc_diff >= 1){
					$gc_diff_tag = 1;
				}else{$gc_diff_tag = 0;}}
		}
	print LOG "\nCandidate region$i: noDR:conj_region$j:$conj_region[$j][0]..$conj_region[$j][1];int_region:$int_region[$int_n][0]..$int_region[$int_n][1];ice:$ice_left..$ice_right\n";
	print LOG "0int;1VirB4;2T4CP;3MOB;4T4SS;6Rep;orit_tag;total_score: $score_sum[0] $score_sum[1] $score_sum[2] $score_sum[3] $score_sum[4] $score_sum[6] $orit_tag $score_sum[5] \n";
	print LOG "e2:$e2	ice_tag;aic_tag;ime_tag:$ice_tag,$aic_tag,$ime_tag\n"; 
	print LOG "Candidate region$i:int_tag;mob_tag;t4ss_tag;gc_diff_tag;ice_len:$int_tag;$mob_tag;$t4ss_tag;$gc_diff_tag;$ice_len\n";			

	if ($int_tag >=1){
		if($ice_tag > 0 && $mob_tag >0 && $t4ss_tag >0 && $gc_diff_tag >0 && $ice_len < 600000 && $ice_len> 4000){ ## the largest length of ICE is set tot be less than 600 kb 
				my $ice_desc = "Putative ICE without identified DR"; ## ICE-like region
				my $ice_out = $download_dir."ICEfinder_$i"."_"."1\n"; ## final result
	  		my $ice_out2 = $download_dir."ICEfinder_$i"."_"."1.core\n"; ## core element	
				open ICEOUT, ">$ice_out"; 
				open ICEOUT2,">$ice_out2";                                                             
		    print ICEOUT "Description: $ice_desc\n".                                    
		                 "Location: $ice_left..$ice_right\n".
		                 "oriT: $orit_desc\n".              
		                 "Length: $ice_len bp\n".  
		                 "GC content: $gc_content %\n".                        
		                 "DR: $dr_desc\n".                                 
		                 "Insert: $insert_desc";                           
		    foreach (sort {$a<=>$b} keys %orf){                            
		    	my $orf_key = chomp($_);      
		    		if(($_ >= $ice_left) && ($orf{$_} <= $ice_right)){         
		    		print ICEOUT $ptt{$_};                                     
		    		if($ptt2{$_}){
				 		print ICEOUT2 "$job_id\t$ptt2{$_}";
				 		}                                  
		    	}	                                                           
		    }                                                              
		    close ICEOUT;                                                  
		    close ICEOUT2;                      		                                                                                                       
		}
		elsif($ime_tag >0 && ($mob_tag >0 ||$orit_tag>0) && $t4ss_tag == 0 && $gc_diff_tag >0 &&  $ice_len < 50000 && $ice_len> 1000){   ## typical: 18-33k for MGI; typical:5-18k for Strepto; max:< 50 kb (exper) 
				my $ice_desc = "Putative IME without identified DR";
				my $ice_out = $download_dir."ICEfinder_$i"."_"."1\n"; ## final result
	  		my $ice_out2 = $download_dir."ICEfinder_$i"."_"."1.core\n"; ## core element
		    #$e2 ++;
				open ICEOUT, ">$ice_out"; 
				open ICEOUT2,">$ice_out2";                                                             
		    print ICEOUT "Description: $ice_desc\n".                                    
		                 "Location: $ice_left..$ice_right\n".
		                 "oriT: $orit_desc\n".              
		                 "Length: $ice_len bp\n".  
		                 "GC content: $gc_content %\n".                        
		                 "DR: $dr_desc\n".                                 
		                 "Insert: $insert_desc";                           
		    foreach (sort {$a<=>$b} keys %orf){                            
		    	my $orf_key = chomp($_);     
		    		if(($_ >= $ice_left) && ($orf{$_} <= $ice_right)){         
		    		print ICEOUT $ptt{$_};                                     
		    		if($ptt2{$_}){
				 		print ICEOUT2 "$job_id\t$ptt2{$_}";
				 		}                                  
		    	}	                                                           
		    }                                                              
		    close ICEOUT;                                                  
		    close ICEOUT2;                                               
		}
		elsif($aic_tag > 0 && $ice_len < 60000 && $ice_len> 4000 && $gram_postive_tag >0){   ## typical: 9-24k; max: <60 kb; the largest length of AICE  is set to be less than 60 kb 
			if(($genus_gbk =~ /Anaerolinea|Bradyrhizobium|Burkholderia|Candidatus|Clostridium|Desulfarculus|Desulfomonile|Enterococcus|Erysipelothrix|Erysipelotrichaceae|Lactobacillus|Lactococcus|Mariprofundus|Marivivens|Mesoplasma|Neisseria|Nitrosomonas|Paenibacillus|Paracoccus|Pseudomonas|Ruminiclostridium|Sphingobium|Staphylococcus|Streptococcus|Xanthomonas|Xylella/i)||($genus_RDP =~ /Anaerolinea|Bradyrhizobium|Burkholderia|Candidatus|Clostridium|Desulfarculus|Desulfomonile|Enterococcus|Erysipelothrix|Erysipelotrichaceae|Lactobacillus|Lactococcus|Mariprofundus|Marivivens|Mesoplasma|Neisseria|Nitrosomonas|Paenibacillus|Paracoccus|Pseudomonas|Ruminiclostridium|Sphingobium|Staphylococcus|Streptococcus|Xanthomonas|Xylella/i)){
				print LOG "This is not in Actinobacteria: genus_gbk:$genus_gbk;genus_RDP:$genus_RDP\n";
			}else{
				my $ice_desc = "Putative  AICE without identified DR"; 
				my $ice_out = $download_dir."ICEfinder_$i"."_"."1\n"; ## final result
	  		my $ice_out2 = $download_dir."ICEfinder_$i"."_"."1.core\n"; ## core element 
		    $orit_desc = "-"; 
				$tra_num = 0;
				foreach (sort {$a<=>$b} keys %orf){                                                                 	
				 		my $orf_key = chomp($_);                                       	
				 		if(($_ >= $ice_left) && ($orf{$_} <= $ice_right)){                                             	
	 						@ptt2_array = split(/\t/,$ptt2{$_});
	 						chomp($ptt2_array[6]); 
	 						if($ptt2_array[6] eq "FtsK_SpoIIIE" ){
	 							$tra_num +=1;
	 						}                                                                         	
						}	
				}
				
				if($tra_num >0){                                                                         
				open ICEOUT, ">$ice_out";  
				open ICEOUT2,">$ice_out2";                                                                          	
				print ICEOUT "Description: $ice_desc\n".                                                                          	
				              "Location: $ice_left..$ice_right\n".
				              "oriT: $orit_desc\n".                                                   	
				              "Length: $ice_len bp\n".
				              "GC content: $gc_content %\n".                                                               	
				              "DR: $dr_desc\n".                                                                      	
				              "Insert: $insert_desc";                                                                	
				 foreach (sort {$a<=>$b} keys %orf){                                                                 	
				 		my $orf_key = chomp($_);                                       	
				 		if(($_ >= $ice_left) && ($orf{$_} <= $ice_right)){                                               	
				 		print ICEOUT $ptt{$_}; 
				 		if($ptt2{$_}){
				 		print ICEOUT2 "$job_id\t$ptt2{$_}";
				 		}                                                                          	
				 	}	                                                                                                 	
				 }                                                                                                   	
				 close ICEOUT;
				 close ICEOUT2;                                                                                                                                                 	
			}
		}} 
	}
	$e2 ++;                                                                                                    			
	$j++;
}


##############

}
close LOG;
my $complete = $tmp_path."region_analyzer_finish";
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
