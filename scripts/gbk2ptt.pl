#!/usr/bin/perl
# Meng Liu and Hong-Yu Ou on April-24-2018
# School of Life Sciences & Biotechnology, Shanghai Jiao Tong University

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Std;

##########################################
#extract ptt file only form GenBank file
###########################################

our ($opt_h, $opt_i, $opt_t, $opt_o, $opt_a) = "";
die "usage: perl gbk2ptt_and_faa.pl -i <gbk_accession>" if @ARGV <= 1;
getopts('i:t:o:a:h');
my $job_id = $opt_a;
my $tmp_path = "$opt_t/$job_id/";
my $gbk_path = "$opt_i/" ;
my $candidate_gbk_1 = $gbk_path."$job_id.gbk";
my $candidate_gbk_2 = $gbk_path."$job_id.gb";
my $ptt = $tmp_path."$job_id.ptt";
my $out = $tmp_path."$job_id.fake_gi";
my $org = "";
my $gbk = "";

if(-e $candidate_gbk_1){
	$gbk = $candidate_gbk_1;
}elsif(-e $candidate_gbk_2){
	$gbk = $candidate_gbk_2;
}else{
	print "ERROR: Can not find the GenBank file of $job_id in $opt_i directory!";
}

system"dos2unix $gbk >/dev/null 2>&1";
open(GBK, $gbk)||die "open $gbk file error: $!\n";
while(<GBK>){
	if(/  ORGANISM  (.+)\n/){
		$org = $1;
		last;
	}
}
close GBK;

my $j = 0;
my $rand = 0;
while($j < 10){
        my $tm  = 10;
        while($tm == 10){
                $tm = int(rand(10)) + 1;
        }
        $rand = $rand * 10 + $tm;
        $j ++;
}
my $fake = 0;
my $seqio_object = Bio::SeqIO->new(-file => "$gbk", -format=>'genbank');       
my $seq_object = $seqio_object->next_seq;
my @cds = grep { $_->primary_tag eq 'CDS' } $seq_object->get_SeqFeatures;


open PTT, ">$ptt";
print  PTT $seq_object->description, " - 1..",$seq_object->length,"\n";
print  PTT scalar(@cds)," proteins\n";
print  PTT join("\t", qw(Location Strand Length PID Gene Synonym Code COG Product)),"\n";
my $i = 0;
for my $f (@cds) {
	$i++;
	my $code = '-'; ##
	my $gi = 0;
	$gi = $1 if tag($f, 'db_xref') =~ m/\bGI:(\d+)\b/;

	if($gi == 0){
		$gi = $i + $rand;
		$fake = 1;
	}

        my $cog = '-';##
        $cog = $1 if tag($f, 'product') =~ m/^(COG\S+)/;##
	my $translation = tag($f, 'translation');
	if($translation ne "-"){
		my $prot_aa = "";
		my $k=0;
		while($k<length($translation)) {
			$prot_aa .= substr($translation,$k,70)."\n";
			$k += 70;
		}

		if($gi == 0){
	                $gi = $i + $rand;
	                $fake = 1;
	        }
		my @col = (
	                $f->start.'..'.$f->end,
	                $f->strand >= 0 ? '+' : '-',
	                ($f->length/3)-1,
	                $gi,
	                tag($f, 'gene'),
	                tag($f, 'locus_tag'),
	                $code,
	                $cog,
	                tag($f, 'product'),
	        );
	        print PTT join("\t", @col), "\n";
	}

}
close PTT;
if($fake){
	open OUT, ">$out";
	close OUT;
}

sub tag {
        my($f, $tag) = @_;
        return '-' unless $f->has_tag($tag);
        return join(' ', $f->get_tag_values($tag));
}

my $complete = $tmp_path."gbk2ptt_finish";
open FIN,">$complete";
close FIN;
exit;
