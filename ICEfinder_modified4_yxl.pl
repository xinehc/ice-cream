#!/usr/bin/env perl

## Modified by Yin Xiaole on October 5th, 2020. yinlele99@gmail.com
# University of Hong Kong

## Meng Liu and Hong-Yu Ou on April-24-2018. hyou@sjtu.edu.cn/liumeng94@sjtu.edu.cn
# School of Life Sciences & Biotechnology, Shanghai Jiao Tong University

use warnings;
use strict;
# use Bio::Perl;
use Bio::SeqIO;
use Getopt::Std;
use File::Basename;
my $SCRIPT_DIR = dirname(__FILE__);


our ($opt_h, $opt_i, $opt_t, $opt_o) = "";
my $usage = <<USE;
  perl $0 -i <folder> -t <folder> -o <folder> -h 
       
       -i Input files directory, required
       -t Temporary files directory, required
       -o Output files directory,required
       -h Print this help information

USE

#die($usage) if ( @ARGV < 1 ); 

getopts('i:t:o:h');

if ($opt_h ){
die "$usage\n"
}

my $minlength = 15; ##direct repeat; dr_min_length
printf "\nThe input folder is $opt_i \n ";

if (-e "$opt_i/split.result.summary.noplasmid.txt"){
    printf "\nThe list file is $opt_i/split.result.summary.noplasmid.txt ";
}else{
    opendir(DIR,$opt_i);
    my @files =grep (/\.gbk$/, readdir(DIR));
    closedir(DIR);
    open(FH,">","$opt_i/split.result.summary.noplasmid.txt") or die $!;
    foreach my $file (@files) {
          print FH "$file\n";
     }
    close(FH);
}

mkdir $opt_t unless -d $opt_t;
mkdir $opt_o unless -d $opt_o;

open FILE, "< $opt_i/split.result.summary.noplasmid.txt" or die "ERROR: can't open the list of gbk files!\n";

while(<FILE>){
    chomp; 
    my $seqinput = $_;#query files
		if ($seqinput =~ /(.+)\.gb/){
			my $acc = $1 ;
			print "For $acc: Examining the input file format...\n";
			if (-e "$opt_i/$acc.gb" || -e "$opt_i/$acc.gbk" ){
				my $tmp_path = "$opt_t/$acc";
				if(-e $tmp_path){
					system("rm -rf $tmp_path");
				}
				system("mkdir $tmp_path");
				system("chmod 777 $tmp_path");		     
				my $candidate_gbk = "$opt_i/$acc.gbk";
				my $candidate_fna = "$tmp_path/$acc.fna";
				my $ptt_file = "$tmp_path/$acc.ptt";
				my $tsv_file = "$tmp_path/$acc.ptt.gi.coords";
				my $ffa_file = "$tmp_path/$acc.gi.ffa";
				my $faa_file = "$tmp_path/$acc.gi.faa";
				#!system("$seqret_cmd $candidate_gbk $candidate_fna 1>/dev/null 2>/dev/null") or die("For $acc:\n Error:  The uploaded file of $acc is not a standard GenBank format!\n");# generate fna file
				#!system("$seqret_cmd $candidate_gbk $candidate_fna 1>/dev/null 2>/dev/null") or die("For $acc:\n Error:  The uploaded file of $acc is not a standard GenBank format!\n");# generate fna file
				!system("perl $SCRIPT_DIR/scripts/gbk2ptt.pl -a $acc -i $opt_i -t $opt_t -o $opt_o") or die("For $acc:\n Error:  The uploaded file of $acc is not a standard GenBank format!\n");	# generate $ptt_file
				TransformCDSFile($ptt_file);
				system "python $SCRIPT_DIR/scripts/readCDSseq.py $candidate_gbk $tsv_file $candidate_fna $ffa_file $faa_file";
				#system "./tools/readCDSseq $candidate_fna $tsv_file $ffa_file 1> /dev/null 2>/dev/null"; ## generate $ffa_file ;$ffa_file contains coding fna 
				#system "$transeq_cmd -clean $ffa_file $faa_file	1> /dev/null 2>/dev/null"; ## 	generate $faa_file with $ffa_file ## add -clean to avoid the stop codon.
				## check the ptt file
				if(! -e $ptt_file){
					print "For $acc: \n ERROR: $ptt_file file was not found!\n";
					next;
				}
				## check the fna file
				if(! -e $candidate_fna){
					print "For $acc: \n ERROR: $candidate_fna file was not found!\n";
					next;
				}
				## check the faa file
				if(! -e $faa_file){
					print "For $acc: \n ERROR: $faa_file file was not found!\n";
					next;
				}
			
				my $genome_seqio_obj = Bio::SeqIO->new(-file => "$candidate_fna", -format => "fasta" );
				my $genome_seq_obj = $genome_seqio_obj->next_seq;
				my $genome_len = $genome_seq_obj->length;
				my $genome_desc = $genome_seq_obj->desc;
				
				if ($genome_len <= 1000000){ ## for contig with size smaller than 1 Mb
					my $cmd_ice_s = "perl $SCRIPT_DIR/scripts/region_finder_s.pl -i $opt_i -a $acc -t $opt_t -o $opt_o";
					system($cmd_ice_s);
				}else{	
					#system("perl ./scripts/gram_detector.pl $acc");		
					system "python $SCRIPT_DIR/scripts/gram_detector.py $acc $opt_t";
     			my $cmd_ice ="perl $SCRIPT_DIR/scripts/region_finder.pl -i $opt_i -a $acc -t $opt_t -o $opt_o";
      		system($cmd_ice);
				}
				    
     		if(-e "$tmp_path/region_finder_finish"){
     			if ( -e "$opt_o/$acc/$acc"."_summary.txt"){
     				print "The detection of ICE/IME has been done! \n\n";  	
    			}else{
    				print "For $acc: The ICE/IME detection has been doned. No ICE/IME was found!\n";
    				system ("rm -rf $opt_o/$acc/");
     			}
    		}
 
    	}else{
					print "For $acc: \n ERROR: Can not find the GenBank file of $opt_i/$acc.gbk in $opt_i directory!\n";
			}
	  }else{
			print "For $seqinput: \n ERROR: $seqinput is not a standard Genbank file(*.gbk or *.gb)!\n $usage";
	 	}
}
		
close FILE;




sub TransformCDSFile{
#####################################################################
## Transform the CDS position file format from 'NCBI PTT' into 'TSV' .
#####################################################################

        my ($infile) = @_;
        my $tempfile = "$infile.gi.coords";
        my @start=();
        my @stop=();
        my @strand=();
        my @name=();
        my @description=();
        my $i=0;
        my $k=0;
        my $CDSNum=0;

        open(INPUT,"<".$infile)  or die "Unable open the file ".$infile." .\n";
        my $lineNo=1;
        while (<INPUT>){          # retrieve file, line by line
             $lineNo++;
             if (( /[0-9]+\.\.[0-9]+/ ) and ($lineNo >2)) {
                my @fields = split "\t";
                if ( $fields[0] =~ /([0-9]+)..([0-9]+)/){
                        $start[$CDSNum] = $1;
                        $stop[$CDSNum] = $2;
                        $strand[$CDSNum] =$fields[1];
                        $name[$CDSNum] = $fields[3];
                        if (scalar @fields > 6){$name[$CDSNum]=$name[$CDSNum]."\n";}
                        $CDSNum++;
                }
             }
          }
        close (INPUT);

        open (OUTFILE, "> $tempfile") || die "Unable open  $tempfile!\n" ;
        for($i=0;$i<$CDSNum;$i++){
                print OUTFILE $start[$i], "\t", $stop[$i], "\t", $strand[$i], "\t", $name[$i];
        }
        close (OUTFILE);
      
}
