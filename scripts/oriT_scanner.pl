#!/usr/bin/perl
# Meng Liu and Hong-Yu Ou on April-24-2018
# School of Life Sciences & Biotechnology, Shanghai Jiao Tong University

use strict;
use Bio::SeqIO;
use Getopt::Std;

our ($opt_a, $opt_n, $opt_t) = "";
my $usage = qq(
This script is used to predict the oriT across the whole candidate region using blastn(H-value >= 0.49)

Usage: script_name.pl <jobid> <region_id>

<oriT_blastn_whole_out>  -outfmt 6 example
(qseqid sseqid  pident  length  mismatch    gapopen qstart  qend    sstart  send    evalue  bitscore )\#<oriT_blastn_whole_out> without this line
NC_004578	100162	91.89	37	3	0	3628526	3628562	202	238	1e-05	52.8
NC_004578	100162	87.18	39	5	0	3627464	3627502	202	240	0.002	45.4

<oriT_blastn_whole_result> example
qseqid	sseqid	identities	length	coverage	hvalue	qstart	qend	sstart	send
NC_004578	100162	87.18	39	94.1	0.82	3627464	3627502	202	240

);
die($usage) if ( @ARGV < 2 );
getopts('a:n:t');
### Declare variables
my $job_id = $opt_a;
my $i = $opt_n;	#candidate region ID
my $tmp_dir = "$opt_t/$job_id/"; #tmp dir
my $candidate_dir = $tmp_dir."/candidate/";
my $db_dir = "./data/";


### input file
my $seq_fna = $candidate_dir."candidate_region_".$i.".fna"; #<seq_fna_file>
### tools or scripts
my $blast = "./tools/blastn";
my $blast_db = "$db_dir/oriT.fna.db";
### parameters
my $evalue = 0.01; 
my $Hvaluecutoff = 0.49; # search across the whole sequence
### intemidate file
my $Log = "$tmp_dir/oriT_blastn.err";
my $blastn_out = $candidate_dir."oriT_$i.blastn"; # <oriT_blastn_whole_out>
### result file
my $final_result = $candidate_dir."oriT_$i.result"; # <oriT_blastn_whole_result>
 
## check <region_fna_file>
if(! -e $seq_fna){
	open ERROR, ">$Log";
	print ERROR "ERROR: For $job_id, during finding the oriT across the whole sequence,	$seq_fna file was not found!\n";
	close ERROR;
	exit;
} 
my $blast_cmd = "$blast -db $blast_db -query $seq_fna -evalue $evalue -word_size 11 -outfmt 6  -num_alignments 1 -num_threads 10 -out $blastn_out";
print "$blast_cmd\n";
system($blast_cmd);

my $content = "";
if(-e $blastn_out){
	open (OUT, $blastn_out); 
	while(<OUT>){
		next unless /.+\s+\d+/;
		chomp;
		my @array = split(/\t/, $_);
		my $acc = $array[1];
		my $identities = sprintf("%.1f",$array [2]);  
		my $length = $array[3];	
		my $oriT_fna = `find data/oriT_seq/ -name $acc.fna`;
		chomp $oriT_fna;
		my $seqio_obj = Bio::SeqIO->new(-file=>"$oriT_fna", -format=>'fasta');
		my $seq_obj = $seqio_obj->next_seq;
		my $seq_length_subject = $seq_obj->length;
		my $coverage = sprintf("%.1f",$length/$seq_length_subject*100);
 	 	my $hvalue= sprintf("%.2f",$identities*($length/$seq_length_subject)/100);
  	#print "$hvalue";
		if($hvalue > $Hvaluecutoff){
			$content .= "$array[0]\t$array[1]\t$identities\t$length\t$coverage\t$hvalue\t$array[6]\t$array[7]\t$array[8]\t$array[9]\n";
		}
	}
	close OUT;
	if($content){
		open RESULT, ">$final_result";
		print RESULT "qseqid\tsseqid\tidentities\tlength\tcoverage\thvalue\tqstart\tqend\tsstart\tsend\n";
		print RESULT $content;	
		close RESULT;
	}
}
