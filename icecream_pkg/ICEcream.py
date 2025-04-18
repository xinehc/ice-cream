#!/usr/bin/env python3
"""
Author: 
Date: 2024-09-13
Actinomycetota and Bacillota are updated on Oct 6th, 2024
"""

import os
import sys
import argparse
import shutil
from Bio import SeqIO
import subprocess
import re
import glob
import logging
import pandas as pd
import pkg_resources
import csv
import gzip
from datetime import datetime
from .scripts.GbffParser_YXL import parse_gbk_files
from .scripts.gbk2ptt import process_gbk
from .scripts.region_finder_s import region_finder_s
from .scripts.region_finder import region_finder
from .scripts.amendORF import amendORF
from .scripts.amendORF2 import amendORF2
from .organize_orfs import organize_orfs
from .plotting_script import plotting_script

RESOURCE_HUB = os.path.join(os.path.dirname(__file__), 'resource')
os.environ['DATABASE_FOLDER'] = RESOURCE_HUB
bakta_db_path = pkg_resources.resource_filename(__name__, './bakta_db')

def decompress_files(input_dir):
    gz_files = [f for f in os.listdir(input_dir) if f.endswith('.gz')]
    for gz_file in gz_files:
        gz_file_path = os.path.join(input_dir, gz_file)
        with gzip.open(gz_file_path, 'rb') as f_in:
            with open(gz_file_path[:-3], 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(gz_file_path)  # Remove the .gz file after decompression

def split_gbff_to_gbk(input_dir, output_dir):
    gbff_files = [f for f in os.listdir(input_dir) if f.endswith('.gbff')]
        
    for gbff_file in gbff_files:
        gbff_path = os.path.join(input_dir, gbff_file)
        # Parse the .gbff file and write each record to a separate .gbk file
        with open(gbff_path, 'r') as input_handle:
            for i, record in enumerate(SeqIO.parse(input_handle, "genbank"), start=1):
                output_file = os.path.join(output_dir, f"{gbff_file}.{i}.gbk")
                with open(output_file, 'w') as output_handle:
                    SeqIO.write(record, output_handle, "genbank")

def setup_logger_for_folder(folder_name):
    logger = logging.getLogger(folder_name)
    logger.setLevel(logging.INFO)
    log_file = os.path.join(folder_name, "locate_ice.log")
    file_handler = logging.FileHandler(log_file,mode="a")
    formatter = logging.Formatter('%(message)s')
    file_handler.setFormatter(formatter)
    if not logger.hasHandlers():  # Prevent adding multiple handlers
        logger.addHandler(file_handler)
    return logger

def check_gram_file(acc,genus_gbk,tmp,given_name,run_actino=False,run_gram=False):
    db_dir = os.environ.get('DATABASE_FOLDER')
    firmicutes_path=os.path.join(db_dir,"genus_name_in_Bacillota.txt") # last update: 10/06/2024
    firmi_names=[]
    for line in open(firmicutes_path,'r'):
        firmi_names.append(str(line).strip())

    actino=os.path.join(db_dir,"genus_name_in_Actinomycetota.txt") # last update: 10/06/2024
    actino_names=[]
    for line in open(actino,'r'):
        actino_names.append(str(line).strip())

    if given_name != "":
        genus=given_name.strip()
    else:
        genus=genus_gbk.split()[0]
    #print(f"genus : {genus}")
    outputfilegram = os.path.join(tmp,acc,f"{acc}.gram")
    if genus in firmi_names:
        with open (outputfilegram, "a") as fileoutput:
            fileoutput.write('Bacillota\n')
    elif genus in actino_names or run_actino :
        with open (outputfilegram, "a") as fileoutput:
            fileoutput.write('Actinomycetota\n')
    if run_gram and not os.path.exists(os.path.join(tmp,acc,f"{acc}.gram")):
        with open (outputfilegram, "a") as fileoutput:
                fileoutput.write('unknown\n')
            
def main():
    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    default_temp_folder = f"temporary_{current_time}"
    parser = argparse.ArgumentParser(description='To classify ICE/IME')
    parser.add_argument('-i', '--inputF', required=True, help='Input folder containing .gbk files (each .gbk contains one chromosome/scafold/contig)')
    parser.add_argument('-o', '--outputF', required=True, help='Output folder')
    parser.add_argument('--genome_fna', action='store_true', default=False, help='Choose if your input files are .fna in the input folder')
    parser.add_argument('--gbff_gz', action='store_true', default=False, help='Choose if your input files are .gbff.gz in the input folder')
    parser.add_argument('--gbff', action='store_true', default=False, help='Choose if your input files are .gbff in the input folder')
    parser.add_argument('-t', '--tempF', default=default_temp_folder, help='Temporary folder')
    parser.add_argument('-g', '--genus',default="", help='Optional: if not satisfied, try to give genus name identified by yourself to see if something different')
    parser.add_argument('-n', '--cpu', type=int, default=2, help='Threads')
    parser.add_argument('--resistance', action='store_true', default=False, help='To annotate resistance genes')
    parser.add_argument('--bakta', action='store_true', default=False, help='Optional: when you provide gbk and would like to run bakta for complentary ORF annotation')
    parser.add_argument('--plot', action='store_true', default=False, help='To plot')
    parser.add_argument('--grampositive', action='store_true', default=False, help='Choose if gram positive')
    parser.add_argument('--actinomycetota', action='store_true', default=False, help='Choose if Actinomycetota')
    parser.add_argument('--bakta_db', default="", help='Optional: give the path to the bakta database')
    args = parser.parse_args()
    if not os.path.exists(args.inputF):
        print(f"Input folder {args.inputF} does not exist.")
        exit(1)

    terminal_size = shutil.get_terminal_size((80, 20))  # Default width=80 if not available
    terminal_width = terminal_size.columns
    welcome_message = [ "=== Welcome to the script of ICEcream ===",
                    	"Author: Xiaole (Charlotte) YIN",
                     	"E-mail: yinlele99@gmail.com",
                       	"Hong Kong University, Harvard Medical School",
                       	"Version 1.14, last update in September 2024"]
    print("=" * terminal_width) 
    for line in welcome_message:
        print(line.center(terminal_width))
    print("=" * terminal_width)

    input_folder = os.path.abspath(args.inputF)
    temp_folder = os.path.abspath(args.tempF)
    output_folder = os.path.abspath(args.outputF)
    threads = args.cpu
    run_resistance = args.resistance
    run_plot = args.plot
    run_gbff_gz = args.gbff_gz
    run_gbff = args.gbff
    run_fna = args.genome_fna
    run_grampositive = args.grampositive
    run_actinomycetota = args.actinomycetota
    given_name = args.genus
    run_bakta = args.bakta
    self_bakta_db = args.bakta_db
    
    if os.path.isdir(temp_folder) or os.path.isdir(output_folder):
        print('The temporary folder or the output folder already exist, please delete it')
        shutil.rmtree(temp_folder, ignore_errors=True)
        shutil.rmtree(output_folder, ignore_errors=True)
        os.makedirs(temp_folder)
        os.makedirs(output_folder)
        #sys.exit(1)
    else:
        os.makedirs(temp_folder)
        os.makedirs(output_folder)
    
    print(f"\n=== The input folder is {input_folder}")
    if run_gbff_gz:
        print("\n=== Step 0: Unzip compressed files, and split .gbff into .gbk")
        decompress_files(input_folder)
        new_gbk_folder = os.path.join(temp_folder,"gbk_dir")
        os.makedirs(new_gbk_folder)
        split_gbff_to_gbk(input_folder, new_gbk_folder)
        input_folder = new_gbk_folder
    elif run_gbff:
        print("\n=== Step 0: Split .gbff into .gbk")
        new_gbk_folder = os.path.join(temp_folder,"gbk_dir")
        os.makedirs(new_gbk_folder)
        split_gbff_to_gbk(input_folder, new_gbk_folder)
        input_folder = new_gbk_folder
    elif run_fna:
        print("\n=== Step 0: Bakta genome fna to get .gbk")
        new_gbk_folder = os.path.join(temp_folder,"gbk_dir")
        os.makedirs(new_gbk_folder)
        for filename_i in os.listdir(input_folder):
            if filename_i.endswith(".fna"):  
                filepath_i = os.path.join(input_folder,filename_i)
                output_i = os.path.splitext(filename_i)[0]
                if self_bakta_db=="":
                    print("""Bakta should be installed: 
                                mkdir -p bakta_db
                                wget https://zenodo.org/record/10522951/files/db-light.tar.gz bakta_db
                                tar -xzf bakta_db/db-light.tar.gz
                                rm bakta_db/db-light.tar.gz""")
                    #prokka_i = ["bakta", "--db", bakta_db,"--output", new_gbk_folder,"--force", "--prefix", output_i, "--keep-contig-headers","--threads", str(threads),  filepath_i ]
                else: 
                    prokka_i = ["bakta","--skip-plot", "--db", self_bakta_db,"--output", new_gbk_folder,"--force", "--prefix", output_i, "--keep-contig-headers","--threads", str(threads),  filepath_i ]
                subprocess.run(prokka_i, check=True)
        input_folder = new_gbk_folder
                            
    print("\n=== Step 1: Exclude potential plasmid by searching keywords")
    parse_gbk_files(input_folder, '.gbk', temp_folder)
    summary_file = os.path.join(temp_folder,"list.gbk.txt")
    summary_file_gbk = os.path.join(temp_folder,"list.gbk.exclude.plasmid.txt")
    if not os.path.exists(summary_file_gbk):
        with open(summary_file, 'w') as fh:
            for file in os.listdir(input_folder):
                if file.endswith('.gbk'):
                    fh.write(f"{file}\n")
        summary_file_gbk = summary_file
    with open(summary_file_gbk, 'r') as s_f_k:
        num_lines_k = sum(1 for _ in s_f_k)
    print(f"\n=== Step 2: Annotate genomes one by one. The count of genomes: {num_lines_k}")
    min_length = 15  # direct repeat; dr_min_length
    str_df = pd.read_csv(os.path.join(RESOURCE_HUB,'3.SARG_v3.2_20220917_Long_structure.icecream.txt'), sep="\t", header=0)
    id_cut = 70
    evalue_cut = 1e-10
    alignratio_cut = 80
    with open(summary_file_gbk, 'r') as file:
        for line in file:
            seq_input = line.strip()
            print(f"\n\n=== Step 2: Annotate genomes for {seq_input}")
            print(f"=== Step 2.1: Locate ICEs")
            if seq_input.endswith(".gb") or seq_input.endswith(".gbk"):
                acc = os.path.splitext(seq_input)[0]
                if os.path.exists(os.path.join(input_folder, f"{acc}.gb")) or os.path.exists(os.path.join(input_folder, f"{acc}.gbk")):
                    tmp_path = os.path.join(temp_folder, acc)
                    folder_logger = setup_logger_for_folder(tmp_path)
                    candidate_gbk = os.path.join(input_folder, f"{acc}.gbk")
                    candidate_fna = os.path.join(tmp_path, f"{acc}.fna") 
                    ptt_file = os.path.join(tmp_path, f"{acc}.ptt")
                    tsv_file = os.path.join(tmp_path, f"{acc}.ptt.gi.coords")
                    ffa_file = os.path.join(tmp_path, f"{acc}.gi.ffa")
                    faa_file =  os.path.join(tmp_path, f"{acc}.gi.faa")
                    process_gbk(acc, input_folder, temp_folder,ffa_file,faa_file) # Generate ptt file   
                    genome_name = ""
                    with open(ptt_file, 'r') as ptt_tmp:
                        ptt_tmp_line1 = ptt_tmp.readline().strip()
                        genome_name_match = re.match(r'^(.*?)\s\-\s(\d+)\.\.(\d+)', ptt_tmp_line1)
                        if genome_name_match:
                            genome_name = genome_name_match.group(1).strip()
                        else:
                            genome_name = ""
                    #print(f"genome name: {genome_name}")
                    if os.path.getsize(tsv_file) <= 0:
                        if not run_bakta:
                            error_message = f"ERROR! For {acc}, no protein (CDS) was recorded in .gbk file!\n"
                            folder_logger.info(error_message)
                            print(error_message)
                            continue
                        else:
                            error_message = f"No protein was recorded in .gbk file! Thus running bakta for {acc}...\n"
                            folder_logger.info(error_message)                            
                            prokka_output = os.path.join(tmp_path, 'bakta_out')
                            os.makedirs(prokka_output)
                            if self_bakta_db=="":
                                print("""Bakta should be installed: 
                                            mkdir -p bakta_db
                                            wget https://zenodo.org/record/10522951/files/db-light.tar.gz bakta_db
                                            tar -xzf bakta_db/db-light.tar.gz
                                            rm bakta_db/db-light.tar.gz""")
                            else:
                                subprocess.run(["bakta","--skip-plot", "--db", self_bakta_db,"--output", prokka_output, "--force", "--prefix", acc, "--keep-contig-headers","--threads", str(threads),   candidate_fna ],stdout=subprocess.DEVNULL)
                                #subprocess.run(["bakta", "--db", bakta_db,"--output", prokka_output, "--force", "--prefix", acc, "--locus-tag", genome_name.split()[0],"--keep-contig-headers","--threads", str(threads),   candidate_fna ])
                            #prokka_out2 = os.path.join(prokka_output,acc)
                            ptt_file_prokka = os.path.join(prokka_output,acc,f"{acc}.ptt") #rename ptt_file
                            tsv_file_prokka = os.path.join(prokka_output,acc,f"{acc}.ptt.gi.coords") # rename tsv_file
                            process_gbk(acc, prokka_output, prokka_output, ffa_file, faa_file)
                            if os.path.getsize(tsv_file_prokka) > 0:
                                shutil.copy(tsv_file_prokka, tsv_file)
                                shutil.copy(ptt_file_prokka, ptt_file)
                            else:
                                print(f"N protein (CDS) was recorded in .gbff file. even after Bakta! So the job exited!\n")
                                continue

                    for file in [ptt_file, candidate_fna, faa_file]:
                        if not os.path.exists(file):
                            folder_logger.info(f"For {acc}: \n ERROR: {file} file was not found!")
                            continue

                    genome_seq_obj = next(SeqIO.parse(candidate_fna, "fasta")) 
                    genome_len = len(genome_seq_obj) 
                    folder_logger.info(f"genome_len: {genome_len}")
                    if genome_len <= 1000000:  
                        region_finder_s(acc, temp_folder, output_folder,threads)
                    else:
                        check_gram_file(acc, genome_name, temp_folder,given_name,run_actinomycetota,run_grampositive)
                        region_finder(acc, temp_folder, output_folder,threads)                     
                    
                    if os.path.exists(os.path.join(tmp_path, "region_finder_finish")):
                        summary_file = os.path.join(output_folder, acc, f"{acc}_summary.txt")
                        if os.path.exists(summary_file):
                            folder_logger.info("Found ICE/IME! \n")
                            with open(summary_file, 'r') as s_f:
                                num_lines = sum(1 for _ in s_f)
                            print(f"=== Step 2.1: Found ICE/IME! In total: {num_lines} was/were detected. Identification output stored in folder {os.path.join(output_folder,acc)}")
                            print(f"=== Step 2.2: Classify ICEs/IMEs")
                            family_log = open(os.path.join(tmp_path, 'family.log'), 'w')
                            cmd = ['prodigal', '-p', 'meta', '-i', candidate_fna, '-a',  os.path.join(tmp_path,f"{seq_input}.fa_prodigal.faa")]
                            subprocess.run(cmd, stdout=family_log, stderr=family_log)
                            faa_file_simple = os.path.join(tmp_path, f"{acc}.faa")
                            if os.path.exists(faa_file_simple):
                                amendORF(faa_file_simple, os.path.join(tmp_path,f"{seq_input}.fa_prodigal.faa"))
                                with open(os.path.join(tmp_path,f"{seq_input}.fa2.faa"), 'w') as wfd:
                                    with open(faa_file_simple, 'r') as infile1:
                                        wfd.write(infile1.read())
                                    with open(os.path.join(tmp_path,f"{acc}.faa_append.faa"), 'r') as infile2:
                                        wfd.write(infile2.read())
                            else:
                                amendORF2(os.path.join(tmp_path,seq_input))
                            
                            iceclasshmm =  os.path.join(tmp_path,f"{seq_input}.fa2.faa.icefinder.hmmscan")
                            cmd = ['hmmscan', '--tblout', f"{iceclasshmm}.out", '--cpu', str(args.cpu), os.path.join(RESOURCE_HUB,'conju.rough.v3.db'), os.path.join(tmp_path,f"{seq_input}.fa2.faa")]
                            subprocess.run(cmd, stdout=family_log, stderr=family_log)
                            with open(f"{iceclasshmm}.out", 'r') as infile, open(f"{iceclasshmm}.modified.out", 'w') as outfile:
                                for lineline in infile:
                                    if not lineline.startswith('#'):
                                        modified_line = re.sub(r' +', '\t', lineline)
                                        outfile.write(modified_line)

                            iceclasshmm_yxl =  os.path.join(tmp_path,f"{seq_input}.fa2.faa.YXL.hmmscan")
                            cmd = ['hmmscan', '--tblout', f"{iceclasshmm_yxl}.out", '--cpu', str(args.cpu), os.path.join(RESOURCE_HUB,'conjugation','YXL_26family.hmm'), os.path.join(tmp_path,f"{seq_input}.fa2.faa")]
                            subprocess.run(cmd, stdout=family_log, stderr=family_log)
                            with open(f"{iceclasshmm_yxl}.out", 'r') as infile, open(f"{iceclasshmm_yxl}.modified.out", 'w') as outfile:
                                for lineline in infile:
                                    if not lineline.startswith('#'):
                                        modified_line = re.sub(r' +', '\t', lineline)
                                        outfile.write(modified_line)
                            
                            iceclasshmm_merge =  os.path.join(tmp_path,f"{seq_input}.fa2.faa.merge.all.hmmscan.out")
                            with open(iceclasshmm_merge, 'w') as outfile:
                                for fname in [f"{iceclasshmm}.modified.out", f"{iceclasshmm_yxl}.modified.out"]:
                                    with open(fname) as infile:
                                        outfile.write(infile.read())
                            os.makedirs(os.path.join(output_folder,f"{acc}_ICEfamily"))
                            cmd = ['Rscript',os.path.join(RESOURCE_HUB, 'merge.version14.R'),acc,'1e-5',output_folder,tmp_path] 
                            subprocess.run(cmd, stdout=family_log, stderr=family_log)  
                            family_log.close()
                            family_folder = os.path.join(output_folder,f"{acc}_ICEfamily")
                            print(f"=== Step 2.2: Classification output stored in folder {family_folder}")
                            blast_output_filter = os.path.join(tmp_path,f"extracted_classification_{seq_input}.fa2.faa.tab.txt")
                            if run_resistance:
                                print(f"=== Step 2.3: Annotate resistance genes on ICEs/IMEs")
                                resistance_log = open(os.path.join(tmp_path, 'resistance.log'), 'w')
                                blastx_out = os.path.join(tmp_path,f'{seq_input}.fa2.faa.tab')
                                cmd = ['blastp', '-db', os.path.join(RESOURCE_HUB,'3.SARG_v3.2_20220917_Long_subdatabase.fasta'), '-query', os.path.join(tmp_path,f'{seq_input}.fa2.faa'), '-out', blastx_out, '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen', '-max_target_seqs', '5']
                                subprocess.run(cmd, stdout=resistance_log, stderr=resistance_log) 
                                if os.path.getsize(blastx_out) > 0:
                                    blastx = pd.read_csv(blastx_out, sep="\t", header=None)
                                    blastx.columns = ['query', 'subject', 'identity', 'length', 'mismatch', 'gapopen', 'querystart', 'queryend', 'sstart', 'send','evalue', 'bitscore', 'slen', 'qlen']
                                    blastx = blastx.sort_values(by=['query', 'evalue', 'bitscore'], ascending=[True, True, False])
                                    blastx = blastx.drop_duplicates(subset='query')
                                    blastx = blastx[(blastx['identity'] >= id_cut) & (blastx['evalue'] <= evalue_cut)]
                                    blastx['alignlength'] = abs(blastx['queryend'] - blastx['querystart']) + 1
                                    blastx['alignqueryratio'] = 100 * blastx['alignlength'] / blastx['qlen']
                                    blastx['alignsubjectratio'] = 100 * blastx['alignlength'] / blastx['slen']
                                    blastx = blastx[(blastx['alignqueryratio'] >= alignratio_cut) & (blastx['alignsubjectratio'] >= alignratio_cut)]
                                    blastx_classify = pd.merge(blastx, str_df, left_on='subject', right_on='SARG.Seq.ID', how='left')
                                    blastx_classify.to_csv(blast_output_filter, sep="\t", index=False, header=True, quoting=False)
                                
                            if run_plot:
                                print(f"=== Step 2.4: Plot identifed ICEs/IMEs in {seq_input}")
                                result_file_path = os.path.join(output_folder,f"{acc}_ICEfamily",f"{acc}_ICEfamily_result.txt")
                                with open(result_file_path, "r") as result_file:
                                    reader = csv.DictReader(result_file, delimiter="\t")
                                    for row_num, row in enumerate(reader, 1):
                                        detail = row.get("details_file_number", "").strip()
                                        file4 = row.get("icefinder_folder", "")
                                        start = row.get("ICE_start", "")
                                        end = row.get("ICE_end", "")
                                        if not detail:
                                            continue  
                                        parts = detail.split('_')
                                        if len(parts) == 2:
                                            if not run_resistance:
                                                with open(blast_output_filter,"w") as blast_file_w:
                                                    pass
                                            detail_name = parts[1]
                                            organize_orfs(f"{start}..{end}",
                                                        os.path.join(output_folder,f"{acc}_ICEfamily",f"{acc}_classification_summary_details{detail_name}.txt"),
                                                blast_output_filter,
                                                os.path.join(output_folder,acc,file4),
                                                os.path.join(output_folder,acc))
                                            plotting_script(f"{output_folder}/{acc}_combined_orfs.txt",output_folder)
                                        elif len(parts) > 2:
                                            print("Executing another type of command")
                                        else:
                                            continue
                            print(f"=== Step 2: {seq_input} finish!")
                        else:
                            folder_logger.info(f"For {acc}: The ICE/IME detection has been done. No ICE/IME was found!")
                            shutil.rmtree(os.path.join(output_folder, acc))
                            print(f"=== Step 2: No ICE/IME was detected!")
                else:
                    print(f"ERROR: Can not find the GenBank file of {acc}.gbk in {input_folder} directory!")
            else:
                print.info(f"ERROR: {seq_input} is not a standard Genbank file(*.gbk or *.gb)!")
    print("\n=== Great! Task finished!")

if __name__ == "__main__":
    main()
    