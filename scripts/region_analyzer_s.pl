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

This version will ignore the DR region.

Usage:	perl script_name.pl <job_id> <region_no> <rep_tag> <maxmismatch> <genome_gc>

For more details, please see the README.txt file.

);
die($usage) if ( @ARGV <= 1 ); 
getopts('i:t:o:a:n:r:d:g:h');

###############################################
##declare major variables
###############################################
my $job_id = $opt_a;
my $i = $opt_n;	#candidate region ID
my $rep_tag = $opt_r;
my $dr_max_mismatch = $opt_d;
my $genome_gc = $opt_g;
system "echo $job_id $i $rep_tag $genome_gc";

my $can_region_left;	#candidate region left coordinate
my $can_region_right;	#candidate region right coordinate
my $tmp_path = "$opt_t/$job_id";
my $candidate_dir = "${tmp_path}/candidate";
my $seq_fna = "$tmp_path/$job_id.fna";
my $db_dir = "./data/";
my $download_dir = "$opt_o/$job_id/";
my $candidate_region_fna = "${candidate_dir}/candidate_region_${i}.fna";
my $candidate_region_faa = "${candidate_dir}/candidate_region_${i}.faa";
my $candidate_region_fea = "${candidate_dir}/candidate_region_${i}.fea";
my $ice_hmmer_out = "${candidate_dir}/ICE_hmm_${i}.out";
my $int_hmmer_out = "$tmp_path/int_hmm.out";
my $orit_out = "$candidate_dir/oriT_${i}.result";
my $Log = "$tmp_path/run_region_analyzer.log";
my $ice_tag_1 = 0; ## int + mob + t4cp + t4ss
my $ice_tag_2 = 0; ## int + t4ss +(mob)+(t4cp), that is, conjugative region without mob or t4cp
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
	open (ICEOUT, "$ice_hmmer_out") || die "ERROR: can not open $ice_hmmer_out file!\n";
#########Example of $ice_hmmer_out
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
		
		if(/(\S+)\s+.+\s+gi\|(\d+)\|/ || /(\S+)\s+.+\s+(\d+)\_\d\s+/){
			if($2 == $gi_tmp){
				next;
			}
			$gi_tmp = $2;
			$annotation{$gi_tmp} = $1;
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
open (OUT2,"$candidate_region_fna") ||die "ERROR: can not find $candidate_region_fna file!\n";
	my $line1 = readline(OUT2);
	close OUT2;
	chomp($line1);
	if($line1 =~ /.+\s+.+\s+(\d+)\.\.(\d+)$/){
		$candidate_region_left = $1;
		$candidate_region_right = $2;
}
print LOG "Candidate region$i:candidate_region_left..candidate_region_right:$candidate_region_left..$candidate_region_right\n";
	
my $orit_tag = 0;
my @orit_region;
my $orit_count = 0;
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
		$orit_region[$orit_count][0] = $orit_left ;
		$orit_region[$orit_count][1] = $orit_coordinate{$orit_left};
		$orit_count ++;
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
my @conj_region;	#t4ss;t4cp;relaxase
my $conj_region_count = 0;
my @inter_conj;
my $inter_conj_count = 0;
my @int_region;
my $int_region_count = 0;
my $distance = 0;	##distance (orf number) between T4SS proteins
my $max_distance = 17;	## 5 ##allowed max distance (orf number) +2 between conjugation proteins in one ICE
my @pos_coor;	## coordinate of each position
my %ptt; ##protein information
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
					$score[$position][$ice_info_designate_no{$annotation{$line[3]}}] = 1;
					$orf_score[$position] += $score[$position][$ice_info_designate_no{$annotation{$line[3]}}];
					$score[$position][5] = $ice_info_designate_no{$annotation{$line[3]}};
					$inter_feature[$inter_feature_count][1] = $orf_right;
					$inter_feature_count ++;
					$inter_feature[$inter_feature_count][0] = $orf_left;
					if($score[$position][5] == 2 || $score[$position][5] == 4 || $score[$position][5] == 3 ||$score[$position][5] == 6){	# T4CP ; T4SS ; ##mliu add relaxase(for IME) && Rep(for AICE)
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
			$ptt{$orf_left} = "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$feature_designate_no_info{$score[$position][5]}\t$annotation{$line[3]}\t$line[8]\n"; ## add \t$annotation{$line[3]}
			#Location	Strand	Length	PID	Gene	Type Product              
			#1176316..1177098	-	260	25011181	gbs1130	-	hypothetical protein 
			#1177098..1177370	-	90	25011182	gbs1131	-	hypothetical protein 
			#3016293..3018443	-	716	4863220007	VCD_003714	Relaxase	TraI

			if (($feature_designate_no_info{$score[$position][5]} ne "-") && ($feature_designate_no_info{$score[$position][5]} ne "")){
				##Location	Strand	Length	PID	Gene	Type protein_name Product
				$ptt2{$orf_left} = "$orf_left..$orf_right\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$feature_designate_no_info{$score[$position][5]}\t$annotation{$line[3]}\n";#\t$line[8]
				print LOG "Candidate region$i:core_protein: $ptt2{$orf_left}"; ## $line[0] eq  $orf_left..$orf_right
			}
		}
	}elsif(/^#Candidate region: (\d+)\.\.(\d+)/){ #if(/^\d+\.\./){
		$can_region_left = $1;
		$inter_conj[$inter_conj_count][0] = $can_region_left;
		$inter_feature[$inter_feature_count][0] = $can_region_left;
		$can_region_right = $2;
	}
}
close CANFEA;


$inter_conj[$inter_conj_count][1] = $can_region_right;
$inter_feature[$inter_feature_count][1] = $can_region_right;
my %rna_inverse = reverse %rna;
my $x = 0;
my $conj_region_count2 = 0; ## the real nubmer of conj_region
while($conj_region[$x][0]){
	 print LOG "Candidate region$i:conj_region$x:conj_region[$x][0]..conj_region[$x][1]:$conj_region[$x][0]..$conj_region[$x][1]\n";
	$x ++;
	$conj_region_count2 ++;
}
$x = 0;
while($inter_conj[$x][0]){
	#print "$inter_conj[$x][0]..$inter_conj[$x][1]\n";
	$x ++;
}
$x = 0;
#print "inter_feature\n";
while($inter_feature[$x][0]){
	#print "inter_feature: $inter_feature[$x][0]..$inter_feature[$x][1]\n";
	$x ++;
}

my $genome_seqio_obj = Bio::SeqIO->new(-file => "$seq_fna", -format => "fasta" );
my $genome_seq_obj = $genome_seqio_obj->next_seq;
my @ptt2_array;
my $tra_num = 0; ## for AICE:FtsK_SpoIIIE
my $marker_name = "";

my $j = 0;
my $e2 = 0; ## for ICE/IME count
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
					$score_sum[6] += $score[$k][6]; #Rep ## added for AICE
				}
			}else{
				if($ice_right >= $pos_coor[$k][1]){
					$right_position = $k;
					$score_sum[0] += $score[$k][0];	#int
					$score_sum[1] += $score[$k][1];	#virB&traU
					$score_sum[2] += $score[$k][2];	#T4CP
					$score_sum[3] += $score[$k][3];	#MOB
					$score_sum[4] += $score[$k][4];	#T4SS
					$score_sum[6] += $score[$k][6]; #Rep ## added for AICE
				}else{
					last;
				}
			}
			$k++;
		}
		if($orit_left >=$ice_left && $orit_coordinate{$orit_left} <= $ice_right)
		{$orit_tag = 1;}else{$orit_tag = 0;} 
		# print "ice_left..ice_right:$ice_left..$ice_right;orit_left..orit_right:$orit_left..$orit_coordinate{$orit_left};\n";

    my $ice_len = $ice_right - $ice_left +1;
    my $string = $genome_seq_obj->subseq($ice_left,$ice_right);
    my $gc_content = &calculate_gc($string); 
		## oriT desc
    if($orit_left >= $ice_left && $orit_coordinate{$orit_left} <= $ice_right ){
    	$orit_desc = "$orit_left..$orit_coordinate{$orit_left}"; 
    }else{$orit_desc = "-";}
   # print "Candidate region$i:test:dr_position$k:0int;1VirB4&traU;2T4CP;3MOB;4T4SS;6Rep;orit_tag;total_score: $score_sum[0] $score_sum[1] $score_sum[2] $score_sum[3] $score_sum[4] $score_sum[6] $orit_tag $score_sum[5] \t";
	 			
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
		    		#$ptt2{$_} example: 3433717..3434979	+	420	9139168891	KPHS_34600	Integrase	Phage_integrase
						#core file example: HS11286	3469356..3472094	+	912	9139168909	KPHS_34780	T4SS ATPase	AAA_10
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
		my $t4ss_core = 1; ## 5

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
		#print "T4SS: core_num:$core_num; t4ss_score:$t4ss_core; t4ss_tag:$t4ss_tag;\n@region_name\n";
		
	 ##calculate feature sores
	 my $ice_score_1 = 0;
	 my $ice_score_2 = 0;
	 my $ice_score_3 = 0;
	 my $ice_score_4 = 0;
	 if($score_sum[2]>0){
			$ice_score_1 = ($score_sum[1] + $score_sum[4]) * $score_sum[2] * $score_sum[3] * $score_sum[0] if ($score_sum[4] >=1 && $score_sum[1]<4 ) ;	#t4ss ice = (VirB4/TraU or T4SS)+ int + mob	+ t4cp (must have t4cp );at least two mpf gene    
		#print "ice_score_1:$ice_score_1\n";
		}
	 if($score_sum[2]<=0 && $score_sum[3] >0 && $core_num >= 5){
			$ice_score_2 = ($score_sum[1] + $score_sum[4]) * $score_sum[3] * $score_sum[0] if ($score_sum[4] >2 && $score_sum[1]<4 ) ; ## if colocalized T4SS core proteins >=5, keep as ICE even without T4CP -- output as conjugative region
		#print "ice_score_2:$ice_score_2\n";
		}
	 if($score_sum[2]>0 && $score_sum[3]<=0 && $core_num >= 5){
			$ice_score_3 = ($score_sum[1] + $score_sum[4]) * $score_sum[2] * $score_sum[0] if ($score_sum[4] >2 && $score_sum[1]<4 ) ; ## if colocalized T4SS core proteins >=5, keep as ICE even without relaxase-- output as conjugative region 
		#print "ice_score_3:$ice_score_3\n";
		}
		if($score_sum[3]<=0 && $score_sum[2]<=0 && $core_num >= 5){
			$ice_score_3 = ($score_sum[1] + $score_sum[4]) * $score_sum[0] if ($score_sum[4] >2 && $score_sum[1]<4 ) ; ## if colocalized T4SS core proteins >=5, keep as ICE even without relaxase-- output as conjugative region 
		#print "ice_score_4:$ice_score_4\n";
		}
		my $aic_score = 0;
		$aic_score = $score_sum[0] * $score_sum[2] * $score_sum[6] if ($score_sum[3]<1);	#aice = int + t4cp + Rep - mob 
		my $ime_score =0;
    $ime_score = $score_sum[0] * ($orit_tag + $score_sum[3]) if (($score_sum[4] < 1)&&($score_sum[1] == 0)); ## IME = Int + Mob - (virB4/TraU )- full T4SS}

    if ($ice_score_1 >0){$ice_tag_1 = $ice_score_1;}
    if($ice_score_2 >0 ){$ice_tag_2 = $ice_score_2;}
    if($ice_score_3 >0){$ice_tag_2 = $ice_score_3;}
    if($ice_score_4 >0){$ice_tag_2 = $ice_score_4;}
		elsif($aic_score >0){$aic_tag = $aic_score;}
		elsif($ime_score >0){$ime_tag = $ime_score;}		
		my $conj_score = ($score_sum[1] + $score_sum[4]) * $score_sum[2] * $score_sum[3]; #$conj_score = (virB4/TraU or T4SS) + t4cp + mob
		my $mob_score = $score_sum[3] * $score_sum[0];  # mob_score = mob + int
		$score_sum[5] = $score[$right_position][5] - $score[$left_position-1][5];	#total Score
  print LOG "Candidate region$i:noDR:0int;1VirB4;2T4CP;3MOB;4T4SS;6Rep;orit_tag;total_score: $score_sum[0] $score_sum[1] $score_sum[2] $score_sum[3] $score_sum[4] $score_sum[6] $orit_tag $score_sum[5] \n";
	print LOG "e2:$e2	ice_tag_1;ice_tag_2;aic_tag;ime_tag:$ice_tag_1,$ice_tag_2,$aic_tag, $ime_tag\n"; 
	print LOG "Candidate region$i:int_tag;mob_tag;t4ss_tag;ice_len:$int_tag;$mob_tag;$t4ss_tag;$ice_len\n";	
	if ($int_tag >=1){
		if($ice_tag_1 > 0 && $mob_tag >0 && $t4ss_tag >0 && $ice_len < 600000 && $ice_len> 1000){ ## the largest length of ICE is set tot be less than 700000 kb ## 500kb ;800 kb; here, 1 Mb for NCBI_9434 genome scanning
				my $ice_desc = "Putative ICE without identified DR"; ## ICE-like region
				my $ice_out = $download_dir."ICEfinder_$i"."_"."1\n"; ## final result
	  		my $ice_out2 = $download_dir."ICEfinder_$i"."_"."1.core\n"; ## core element
	  		$e2 ++;
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
		if($ice_tag_2 > 0  && $t4ss_tag >0 && $ice_len < 600000 && $ice_len> 1000){ ## the largest length of ICE is set tot be less than 600kb ## 500kb ;800 kb; 
				my $ice_desc = "Putative conjugative region"; ## ICE-like region
				my $ice_out = $download_dir."ICEfinder_$i"."_"."1\n"; ## final result
	  		my $ice_out2 = $download_dir."ICEfinder_$i"."_"."1.core\n"; ## core element	
	  		$e2 ++;	
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
		elsif($ime_tag >0 && ($mob_tag >0 ||$orit_tag>0) && $t4ss_tag == 0 &&  $ice_len < 50000 && $ice_len> 1000){   ## typical: 18-33k for MGI; typical:5-18k for Strepto; max:< 50 kb (exper);## the largest length of IME  is set to be less than 100 kb ; here, 500 kb for NCBI_9434 genome scanning
				my $ice_desc = "Putative IME without identified DR";
				my $ice_out = $download_dir."ICEfinder_$i"."_"."1\n"; ## final result
	  		my $ice_out2 = $download_dir."ICEfinder_$i"."_"."1.core\n"; ## core element
	  		$e2 ++;
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
		elsif($aic_tag > 0 && $ice_len < 60000 && $ice_len> 4000){    ## typical: 9-24k; max: <60 kb; the largest length of AICE  is set to be less than 70 kb  ## 60 kb here, 100 kb for NCBI_9434 genome scanning
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
				
				if($tra_num >0){  ##check Tra protein
				$e2 ++;                                                                       
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
		} 
	}                                                                                                    			
	$j++;
}


if (($conj_region_count2 == 0) && ($int_region_count >0) && ($orit_count >0) ){ ## find MGI: only Int and oriT
	my $ice_left = $int_region[$j][0];
	my $ice_right = $int_region[$j][1];
	my $j=0;
	while($int_region[$j][0]){	
		##find neighbour integrase
		$ice_left = $int_region[$j][0];
		 $ice_right = $int_region[$j][1];
		my $orit_n = 0;
		print LOG "\nTo find MGI\nCandidate region$i:conj_region_count:$conj_region_count;int_region_count:$int_region_count;orit_count:$orit_count;int_region$j: $int_region[$j][0]..$int_region[$j][1];orit_region$orit_n:$orit_region[$orit_n][0]..$orit_region[$orit_n][1]\n";
		while($orit_region[$orit_n][0]){
			if(($int_region[$j][0] - $orit_region[$orit_n][1]) < 40000){
				$ice_left = $orit_region[$orit_n][0] if $orit_region[$orit_n][0] < $ice_left;
				last;
			}
			$orit_n ++;
		}
		while($orit_region[$orit_n][0]){
			if(($orit_region[$orit_n][0] - $int_region[$j][1]) < 40000){
				$ice_right = $orit_region[$orit_n][1] if $orit_region[$orit_n][1] > $ice_right;	
				$orit_n ++;
			
			}else{
				last;
			}
		}
	$j++;
	}
	my $ice_len = $ice_right - $ice_left + 1;
  my $string = $genome_seq_obj->subseq($ice_left,$ice_right);
  my $gc_content = &calculate_gc($string); 
  if($orit_left >= $ice_left && $orit_coordinate{$orit_left} <= $ice_right ){
    	$orit_desc = "$orit_left..$orit_coordinate{$orit_left}"; 
    }else{$orit_desc = "-";}
 print LOG "ice_left..ice_right:$ice_left..$ice_right;orit_left..orit_right:$orit_left..$orit_coordinate{$orit_left};\n";
	 		
	if( $ice_len < 50000 && $ice_len> 1000){   ## typical: 18-33k for MGI; typical:5-18k for Strepto; max:< 50 kb (exper); the largest length of IME  is set to be less than 90 kb ; here, 100 kb for NCBI_9434 genome scanning
				my $ice_desc = "Putative IME without identified DR";
				my $ice_out = $download_dir."ICEfinder_$i"."_"."2\n"; ## final result
	  		my $ice_out2 = $download_dir."ICEfinder_$i"."_"."2.core\n"; ## core element
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
				 	$e2 ++;                                                
			}	
}



close LOG;
my $complete = "$tmp_path/region_analyzer_finish";
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
