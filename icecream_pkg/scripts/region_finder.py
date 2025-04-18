#!/usr/bin/env python3
"""
Author: 
Date: YYYY-MM-DD
"""

import sys
import os
import argparse
from Bio import SeqIO
import subprocess
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import shutil
import logging
from .oriT_scanner import oriT_scanner
from .region_analyzer import judge_ice_ime_aice

def setup_subcode_logger(log_file):
    logger = logging.getLogger('scripts.region_finder')  # Module-specific logger
    logger.setLevel(logging.INFO)
    
    # Create a file handler for subcode-specific logging
    file_handler = logging.FileHandler(log_file, mode='a')
    
    # Create a log format
    formatter = logging.Formatter("%(message)s")
    file_handler.setFormatter(formatter)

    # Add the handler to the logger
    if not logger.hasHandlers():  # Ensure no duplicate handlers are added
        logger.addHandler(file_handler)

    return logger

def calculate_gc(sequence):
    seq_length = len(sequence)
    num_a = sequence.count('A') + sequence.count('a')
    num_c = sequence.count('C') + sequence.count('c')
    num_g = sequence.count('G') + sequence.count('g')
    num_t = sequence.count('T') + sequence.count('t')
    num_o = seq_length - num_a - num_c - num_g - num_t
    gc = round((num_c + num_g) / seq_length * 100, 2)
    return gc

def parse_hmm_result(hmm,threshold_evalue,threshold_cvalue):
    content = ""
    hit_dom = {}  # Aligned domains of regions within a protein
    try:
        with open(hmm, "r") as HMM:
            for line in HMM:
                if not line.startswith("#"):
                    array = line.strip().split()
                    if not array:
                        continue

                    if re.match(r"^gi\|(\d+)\|", array[3]):
                        gi = re.match(r"^gi\|(\d+)\|", array[3]).group(1)
                    elif re.match(r"^(\d+)_\d", array[3]):
                        gi = re.match(r"^(\d+)_\d", array[3]).group(1)
                    else:
                        continue

                    domain_name = array[0].replace(' ', '')
                    evalue = float(array[6])
                    if evalue > threshold_evalue.get(domain_name, float('inf')):
                        continue
                    if float(array[11]) > threshold_cvalue.get(domain_name, float('inf')):
                        continue

                    if gi not in hit_dom:
                        hit_dom[gi]={'domain': domain_name, 'evalue': evalue}
                    else:
                        if evalue < hit_dom[gi]['evalue']:
                            hit_dom[gi] = {'domain': domain_name, 'evalue': evalue}
                                               
    except IOError as e:
        print(f"open {hmm} error: {e}")

    for gi, data in hit_dom.items():
        content += f"{gi}\t{data['domain']}\n"
        
    return content


def check_candidate_region(hmm,threshold_evalue,threshold_cvalue,dom_function):
    continue_ssICE = 0
    continue_IME = 0
    continue_AICE = 0
    
    # Initialize the hit dictionary
    hit = {
        'Rep': 0,
        'Tra': 0,
        'Mob': 0,
        'T4CP': 0,
        'T4SS': 0,
    }
    
    try:
        with open(hmm, "r") as hmm_file:
            for line in hmm_file:
                if not line.startswith("#"):
                    array = line.strip().split()
                    if not array:
                        continue
                    
                    gi = None
                    if re.match(r"^gi\|(\d+)\|", array[3]):
                        gi = re.match(r"^gi\|(\d+)\|", array[3]).group(1)
                    
                    if array[0] == "Pfam-B_3022":
                        continue
                    
                    if float(array[6]) > threshold_evalue.get(array[0], float('inf')):
                        continue
                    if float(array[11]) > threshold_cvalue.get(array[0], float('inf')):
                        continue
                    
                    if array[0] in dom_function:
                        hit[dom_function[array[0]]] = 1
                    
    except IOError as e:
        print(f"open {hmm} error: {e}")
    
    if hit['Mob'] + hit['T4SS']:
        continue_ssICE = 1
        continue_IME = 1
    elif hit['Rep']:
        continue_AICE = 1
    
    return continue_ssICE, continue_IME, continue_AICE


def region_finder(job_id,tmp_folder,output_folder,threads):
    tmp_path = os.path.join(tmp_folder, job_id)
    seq_fna = os.path.join(tmp_path, f"{job_id}.fna")
    seq_faa = os.path.join(tmp_path, f"{job_id}.gi.faa")
    seq_ptt = os.path.join(tmp_path, f"{job_id}.ptt")
    seq_RNA = os.path.join(tmp_path, f"{job_id}.rnt")  # rnt file containing the tRNA info converted from tmRNA_aragorn.out
    
    # Working directories
    candidate_dir = os.path.join(tmp_path, "candidate")
    download_dir = os.path.join(output_folder, job_id)
    
    shutil.rmtree(candidate_dir, ignore_errors=True)
    os.makedirs(candidate_dir, mode=0o777, exist_ok=True)
    os.chmod(candidate_dir, 0o777)
    shutil.rmtree(download_dir, ignore_errors=True)
    os.makedirs(download_dir, mode=0o777, exist_ok=True)
    os.chmod(download_dir, 0o777)
    
    error_log = os.path.join(tmp_path, "run_region_finder.log")
    logger = setup_subcode_logger(error_log)
    
    db_dir = os.environ.get('DATABASE_FOLDER')
    int_hmm_db = os.path.join(db_dir,"integrase.v2.db")  # transposase and integrase hmm
    ice_hmm_db = os.path.join(db_dir,"conju.rough.v3.db")  # ICE signature protein hmm
    threshold_file = os.path.join(db_dir,"threshold")  # JUN defined parameters

    genome_name = ""
    genome_size = 0
    hmm_evalue = 0.0001
    rep_tag = 0
    dr_max_mismatch = 0
    
    # Read threshold file
    threshold_evalue = {}
    threshold_cvalue = {}
    dom_function = {}
    
    ###### Example of $threshold_file
    #Domain name	E-value	c-Evalue	Function
    #Phage_integrase	0.000001	0.0001	Int
    #Recombinase	0.000001	0.000001	Int
    #rve	0.000001	0.000001	Int
    with open(threshold_file, 'r') as THRE:
        next(THRE)  # Skip header
        for line in THRE:
            array = line.strip().split('\t')
            if not array[0] or array[0] == "Domain name":
                continue
            domain_name = array[0].replace(' ', '')
            threshold_evalue[domain_name] = float(array[1])
            threshold_cvalue[domain_name] = float(array[2])
            dom_function[domain_name] = array[3]
    
    # Integrase hmmscan task and read files:
    if not os.path.exists(seq_ptt):
        print(f"ERROR! For {job_id}, no protein was found: {seq_ptt} file was not found!")
        logger.info(f"ERROR! For {job_id}, no protein was found: {seq_ptt} file was not found!")
        return
    
    if not os.path.exists(seq_faa):
        print(f"ERROR! For {job_id}, no protein was found: {seq_faa} file was not found!")
        logger.info(f"ERROR! For {job_id}, no protein was found: {seq_faa} file was not found!")
        return
    
    if not os.path.exists(seq_fna):
        print(f"ERROR! For {job_id}, {seq_fna} file was not found!")
        logger.info(f"ERROR! For {job_id}, {seq_fna} file was not found!")
        return
    
    # Read genome sequence
    genome_seq_obj = next(SeqIO.parse(seq_fna, "fasta"))
    genome_id = genome_seq_obj.id
    genome_desc = genome_seq_obj.description
    genome_seq = str(genome_seq_obj.seq)
    genome_gc = calculate_gc(genome_seq)
    genome_len = len(genome_seq)
    genome_size = genome_len
    
    left_region_range = 250000
    right_region_range = 250000
    
    if os.path.exists(os.path.join(tmp_folder,job_id,f"{job_id}.gram")):
        left_region_range = 43000
        right_region_range = 43000
    
    #logger.info(f"For {job_id}:left_region_range:{left_region_range}")
        
    logger.info("Start detecting candidate regions containing integrase flanked by tRNA....")
    int_hmm_out = os.path.join(tmp_path,"int_hmm.out")
    if os.path.exists(int_hmm_out):
        os.remove(int_hmm_out)
 
    int_hmm_cmd = ['hmmscan', '--domtblout', int_hmm_out, "--cpu", str(threads), "-E", str(hmm_evalue), int_hmm_db, seq_faa]
    #print("Running command:", " ".join(int_hmm_cmd))
    with open(error_log, 'a') as error_log_file:
        subprocess.run(int_hmm_cmd, stdout=subprocess.DEVNULL, stderr=error_log_file)
        #result = subprocess.run(int_hmm_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    ########## Process ptt file
    protein_num = 0
    ptt_gi_line = {}    # GI as key, line as value
    ptt_gi_coor = {}    # GI as key, left coordinate as value
    ptt_coor_gi = {}    # Left coordinate as key, GI as value
    ptt_coordinate = {} # Left coordinate as key, right coordinate as value
    
    try:
        with open(seq_ptt, 'r') as ptt:
            line1 = ptt.readline().strip()
            genome_name_match = re.match(r'^(.*?)\s\-\s(\d+)\.\.(\d+)', line1)
            if genome_name_match:
                genome_name = genome_name_match.group(1).strip()
            else:
                genome_name = ""
            for line in ptt:
                line = line.strip()
                match = re.match(r'^(\d+)\.\.(\d+)\t.\t.+\t(\d+)\t.+\t.+\t.+\t.+\t.+', line)
                if match:
                    left, right = int(match.group(1)), int(match.group(2))
                    gi = match.group(3)
                    
                    protein_num += 1
                    if left > right:
                        continue
                    
                    ptt_gi_line[gi] = line
                    ptt_gi_coor[gi] = left
                    ptt_coor_gi[left] = gi
                    ptt_coordinate[left] = right
    
    except IOError as e:
        print(f"ERROR! For {job_id}, failed to open {seq_ptt}: {e}")
        logger.info(f"ERROR! For {job_id}, failed to open {seq_ptt}: {e}")
        return
    
    if protein_num == 0:
        error_message = f"ERROR! For {job_id}, no protein was found: No protein was recorded in {seq_ptt} file!\n"
        print(error_message)
        logger.info(error_message)
        return      
    ########## Process faa file
    pro = {}  # Dictionary with GI as key and FASTA as value
    for seq_record in SeqIO.parse(seq_faa, "fasta"):
        dis = seq_record.id
        des = seq_record.description
        string = str(seq_record.seq)
        
        if not string:
            string = ""
        
        ref = ""
        if "gi|" in dis:
            ref_match = re.match(r'gi\|(\d+)\|', dis)
            if ref_match:
                ref = ref_match.group(1)
        elif re.match(r'^(\d+)_', dis):
            ref_match = re.match(r'^(\d+)_', dis)
            if ref_match:
                ref = ref_match.group(1)
        
        pro[ref] = f">{dis} {des}\n"
        
        # Add the sequence in lines of 70 characters
        for k in range(0, len(string), 70):
            pro[ref] += string[k:k+70] + "\n"
    #print("test1")
    if not os.path.exists(seq_RNA):
        tmrna = os.path.join(tmp_path, "tmRNA_aragorn.out")
        result = subprocess.run(['aragorn', "-o", tmrna, seq_fna],capture_output=True, text=True)
        with open(tmrna, "r") as TMRNA, open(seq_RNA, "w") as NEWRNT:
            rna = {}
            rna_name = {}
            file_rnt = ""
            trna_num = 0
            
            for line in TMRNA:
                match = re.search(r"Sequence \[(\d+),(\d+)\]", line)
                if match:
                    start, end = int(match.group(1)), int(match.group(2))
                    length = end - start + 1
                    rna[start] = end
                    rna_name[start] = f"{start}..{end}\t+\t{length}\ttRNA\ttRNA\t-\t-\t-\ttRNA"
                    file_rnt += rna_name[start] + "\n"
                    trna_num += 1
                
                match = re.search(r"Sequence c\[(\d+),(\d+)\]", line)
                if match:
                    start, end = int(match.group(1)), int(match.group(2))
                    length = end - start + 1
                    rna[start] = end
                    rna_name[start] = f"{start}..{end}\t-\t{length}\ttRNA\ttRNA\t-\t-\t-\ttRNA"
                    file_rnt += rna_name[start] + "\n"
                    trna_num += 1
            
            # Write output to the new rnt file
            NEWRNT.write("Complete RNA information\n")
            NEWRNT.write(f"tRNA: {trna_num}\n")
            NEWRNT.write("\t".join(["Location", "Strand", "Length", "PID", "Gene", "Synonym", "Code", "COG", "Product"]) + "\n")
            NEWRNT.write(file_rnt)
    
    #print("test2")
    ########## Process rna file
    rna = {}        # RNA left coordinate as key, right coordinate as value
    rna_name = {}   # RNA left coordinate as key, line content as value
    
    if os.path.exists(seq_RNA):
        with open(seq_RNA, 'r') as rna_file:
            for line in rna_file:
                if re.match(r'^\d+\.\.', line):
                    line = line.strip()
                    line_rna = line.split('\t')
                    
                    match = re.match(r'(\d+)\.\.(\d+)', line_rna[0])
                    if match:
                        rna_left, rna_right = int(match.group(1)), int(match.group(2))
                        length = rna_right - rna_left
                        if length < 60:
                            continue  # Skip if the length is less than 60
                        
                        if line_rna[8] == '-':
                            if line_rna[4] == "ssrA":
                                rna[rna_left] = rna_right
                                rna_name[rna_left] = f"{line}\tRNA"
                        else:
                            rna[rna_left] = rna_right
                            rna_name[rna_left] = f"{line}\tRNA"
    
    # Reverse the rna dictionary (right coordinate as key, left coordinate as value)
    rna_inverse = {v: k for k, v in rna.items()}
    
    count_rna = 0
    rna_left = []
    rna_right = []
    for rna_key in sorted(rna.keys()):
        rna_left.append(rna_key)
        rna_right.append(rna[rna_key])
        count_rna += 1

    ############ Process integrase HMM results
    if not os.path.exists(int_hmm_out):
        error_message = f"ERROR! {int_hmm_out} file was not found!\n"
        print(error_message)
        logger.info(error_message)
        return
    
    int_count = 0
    integrase_gi = {}
    #print("test3")

    int_hmm_parse_result = parse_hmm_result(int_hmm_out,threshold_evalue,threshold_cvalue)
    #example: 4779552763_1    Phage_integrase    TIGR02249
    int_hmm_parse_line = int_hmm_parse_result.strip().split('\n')
    can_region = []
    i = 0
    ii = 0
    #logger.info("test5")
    if len(int_hmm_parse_line) > 0 and int_hmm_parse_line != ['']:
        while ii < len(int_hmm_parse_line):
            #logger.info(f"ii is {ii}")
            c = int_hmm_parse_line[ii].split('\t')
            print(f"int_hmm_parse_line: {int_hmm_parse_line}")
            print(f"c: {c}")
            ii += 1
            gi = c[0]
            integrase_gi[gi] = c[1]
            int_count += 1
            
            # Looking for the candidate region
            int_left_end = ptt_gi_coor[gi]
            int_right_end = ptt_coordinate[ptt_gi_coor[gi]]
            #logger.info(f"int_left_end is {int_left_end}")
            #logger.info(f"int_right_end is {int_right_end}")
            # Boundaries by tRNA and tmRNA
            count_rna = 0
            can_region_left_rna = 1
            can_region_right_rna = genome_size
            
            while count_rna < len(rna_right) and rna_right[count_rna] < int_left_end:
                #logger.info(f"count_rna is {count_rna}; rna_right[count_rna] is {rna_right[count_rna]} ")
                can_region_left_rna = rna_left[count_rna]
                count_rna += 1

            #logger.info(f"can_region_left_rna is {can_region_left_rna}")
            if count_rna < len(rna_right):
                can_region_right_rna = rna_right[count_rna]
            
            can_region_left = int_left_end - left_region_range
            if can_region_left < 0:
                can_region_left = 1
            
            can_region_right = int_right_end + right_region_range
            if can_region_right > genome_size:
                can_region_right = genome_size
            #logger.info(f"can_region_left is {can_region_left}")
            #logger.info(f"can_region_right is {can_region_right}")
            if can_region_left_rna > can_region_left:
                can_region_left = can_region_left_rna
                #logger.info("1")
            #logger.info(f"an_region_right_rna is {can_region_right_rna}")
            if can_region_right_rna < can_region_right:
                can_region_right = can_region_right_rna
                #logger.info("2")
            
            if int_count == 1:
                can_region.append([can_region_left, can_region_right])
                #logger.info(f"can_Region : {can_region}")
            elif can_region_left < can_region[i][1] and can_region_left != rna_inverse.get(can_region[i][1], None):
                #logger.info(f"can_region_left is {can_region_left}; can_region[i][1] is {can_region[i][1]}; rna_inverse.get(can_region[i][1], None) is {rna_inverse.get(can_region[i][1], None)}; i is {i}")
                can_region[i][1] = can_region_right
                #logger.info(f"can_Region : {can_region}")
            else:
                i += 1
                can_region.append([can_region_left, can_region_right])
            #logger.info(f"can_Region : {can_region}")
    #print(f"integrase_gi: {integrase_gi}")        
    # If no integrase is found, print an error and exit
    #print(f"int_count is {int_count}")
    if int_count == 0:
        error_message = f"No integrase was found! So the job exited!\n\n"
        #print(error_message)
        logger.info(error_message)
        # Remove the download directory and exit
        shutil.rmtree(download_dir, ignore_errors=True)
        return
    #print("test6")
    ###### Process candidate regions
    k = 0
    i = 0
    while k < len(can_region) and can_region[k][0]:
        candidate_region_fna = os.path.join(candidate_dir, f"candidate_region_{i}.fna")
        candidate_region_faa = os.path.join(candidate_dir, f"candidate_region_{i}.faa")
        candidate_region_fea = os.path.join(candidate_dir, f"candidate_region_{i}.fea")
        #print(f"k is {k}")
        #print(f"can_region[k][0]:{can_region[k][0]}, can_region[k][1]:{can_region[k][1]},")
        # Extract fna
        #print(type(genome_seq_obj.seq[can_region[k][0]-1:can_region[k][1]]))
        string = str(genome_seq_obj.seq[can_region[k][0]-1:can_region[k][1]])
        seq_obj = SeqRecord(Seq(string), id=genome_id, description=f"candidate_region_{i} {can_region[k][0]}..{can_region[k][1]}")
        
        with open(candidate_region_fna, "w") as output_handle:
            SeqIO.write(seq_obj, output_handle, "fasta")
        
        # Extract faa and annotate information
        with open(candidate_region_faa, "w") as canfaa, open(candidate_region_fea, "w") as canfea:
            canfea.write(f"#{genome_name}\n")
            canfea.write(f"#Genome GC {genome_gc} %\n")
            canfea.write(f"#Candidate region: {can_region[k][0]}..{can_region[k][1]}\n")
            
            # Print RNA information if it exists
            if can_region[k][0] in rna:
                canfea.write(f"{rna_name[can_region[k][0]]}\n")
            
            line_num = 0
            for pro_left in sorted(ptt_coor_gi.keys()):
                pro_right = ptt_coordinate[pro_left]
                pro_gi = ptt_coor_gi[pro_left]
                line = ptt_gi_line[pro_gi]
                if pro_left >= can_region[k][0] and pro_right < can_region[k][1]: # this line is different in region_finder_s.py
                    if pro_gi not in pro:
                        continue
                    line_num += 1
                    canfaa.write(f"{pro[pro_gi]}\n")
                    canfea.write(line)
                    if pro_gi in integrase_gi:
                        canfea.write("\tINT")
                    canfea.write("\n")
            
            # Check for RNA information at the end of the region
            if can_region[k][1] in rna_inverse:
                canfea.write(f"{rna_name[rna_inverse[can_region[k][1]]]}\n")
        
        # Remove files if no proteins are found
        if line_num == 0:
            os.remove(candidate_region_fna)
            os.remove(candidate_region_faa)
            os.remove(candidate_region_fea)
            k += 1
            continue
        
        # Submit HMMER jobs
        if not os.path.exists(candidate_region_faa):
            logger.info(f"For {job_id}, {candidate_region_faa} file is not found")
            return
        
        ice_hmm_out = os.path.join(candidate_dir, f"ICE_hmm_{i}.out")
        hmm_ice_cmd = ['hmmscan', '--domtblout', ice_hmm_out,'--cpu', str(threads) ,ice_hmm_db, candidate_region_faa]
        with open(error_log, 'a') as error_log_file:
            subprocess.run(hmm_ice_cmd, stdout=subprocess.DEVNULL, stderr=error_log_file)
        oriT_scanner(job_id,i,tmp_folder,threads)

        conjscan_out = os.path.join(candidate_dir, f"candidate_region_{i}_t4ss")
        macsyfinder_cmd = [
                    "macsyfinder",
                    "--db-type", "ordered_replicon",
                    "--sequence-db", candidate_region_faa,
                    "--models", "CONJScan/Chromosome",
                    "all",
                    "-o", conjscan_out,
                    "--mute"
                ]
        subprocess.run(macsyfinder_cmd, check=True)

        k += 1
        i += 1
    
    # check hmm result
    logger.info("Start analysing the Relaxase, T4SS, T4CP, Tra, Rep, ARG, VF of the candidate regions...\n")
    # Initialize counters
    runing_job_count0 = i  # For counting rep_tag for AICE
    runing_job_count = i
    runing_job_num0 = i  # For counting rep_tag for AICE
    runing_job_num = i
    
    
    # Count rep_tag for AICE
    while runing_job_count0:
        i = 0
        while i < runing_job_num0:
            runing_job_count0 -= 1
            ice_hmmer_out = os.path.join(candidate_dir, f"ICE_hmm_{i}.out")
            goon_ssICE, goon_IME, goon_AICE = check_candidate_region(ice_hmmer_out,threshold_evalue,threshold_cvalue,dom_function)
            
            if goon_AICE > 0:
                rep_tag += 1
            
            i += 1
            
            if runing_job_count0 == 0:
                break
            
        if runing_job_count0 == 0:
            break
    
    # Initialize variables
    ice_candidate_no = {}
    ice_candidate_count = 0
    #print("test6.5")
    # Process each candidate region
    while runing_job_count:
        i = 0
        while i < runing_job_num:
            runing_job_count -= 1
            ice_hmmer_out = os.path.join(candidate_dir, f"ICE_hmm_{i}.out")
            goon_ssICE, goon_IME, goon_AICE = check_candidate_region(ice_hmmer_out, threshold_evalue, threshold_cvalue,dom_function)
            #print(f"goon_ssICE: {goon_ssICE}, goon_IME: {goon_IME}, goon_AICE: {goon_AICE}")
            # Check for oriT tag for IME
            orit_out = os.path.join(candidate_dir, f"oriT_{i}.result")
            orit_tag = 1 if os.path.exists(orit_out) else 0
            
            # Mark this region as a candidate
            if goon_ssICE + goon_IME + goon_AICE + orit_tag:
                ice_candidate_no[i] = 1
                ice_candidate_count += 1
                #print("test7")
                # Submit vmatch job
                candidate_region_fna = os.path.join(candidate_dir, f"candidate_region_{i}.fna")
                #print("test8")
                if not os.path.exists(candidate_region_fna):
                    error_message = f"ERROR! For {job_id}, {candidate_region_fna} file is not found!\n"
                    print(error_message)
                    logger.info(error_message)
                    return
                
                #print("test9")
                judge_ice_ime_aice(job_id,i,rep_tag,dr_max_mismatch,genome_gc,tmp_folder,output_folder) 
            
            i += 1
            
            if runing_job_count == 0:
                break
        
        if runing_job_count == 0:
            break
    
    #print("test12")
    # Generate the download file
    if not ice_candidate_count:
        complete = os.path.join(tmp_path, "candidate_region_none")
        print("The ICE/IME detection has been done. No ICE/IME candidate region was detected!\n\n")
        with open(complete, 'w') as fin:
            fin.write("The ICE/IME detection has been done. No ICE/IME candidate region was detected!\n\n")
        shutil.rmtree(download_dir, ignore_errors=True)
        return
    
    ice_count = 0
    runing_job_count = ice_candidate_count
    runing_job_num = ice_candidate_count
    
    while runing_job_count:
        for i in sorted(ice_candidate_no.keys()):
            runing_job_count -= 1
            with os.scandir(download_dir) as it:
                for entry in it:
                    if entry.is_file():
                        match = re.match(rf"^ICEfinder_{i}_(\d+)$", entry.name)
                        if match:
                            ice_count += 1
                            with open(os.path.join(download_dir, entry.name)) as in_file:
                                ice_pro = os.path.join(download_dir, f"Protein_{entry.name}.fas")
                                with open(ice_pro, "w") as icefaa_file:
                                    for line in in_file:
                                        if line.startswith("Location:"):
                                            loc_match = re.search(r"Location: (\d+)\.\.(\d+)", line)
                                            if loc_match:
                                                start, end = int(loc_match.group(1)), int(loc_match.group(2))
                                                ice_fas = os.path.join(download_dir, f"DNA_{entry.name}.fas")
                                                ice_seq = genome_seq_obj.seq[start-1:end]  # Adjust for 0-based index
                                                ice_des = f"{genome_desc}: {start}..{end}"
                                                ice_seq_obj = SeqRecord(Seq(ice_seq), id=genome_id, description=ice_des)
                                                with open(ice_fas, "w") as output_handle:
                                                    SeqIO.write(ice_seq_obj, output_handle, "fasta")
                                        elif re.match(r"^\d+\.\.\d+\t.\t\d+\t(\d+)", line):
                                            pro_id = re.search(r"^\d+\.\.\d+\t.\t\d+\t(\d+)", line).group(1)
                                            icefaa_file.write(f"{pro[pro_id]}\n")
            
            if runing_job_count == 0:
                break
        if runing_job_count == 0:
            break
    
    # Result summary
    ice_all_out = os.path.join(download_dir, f"{job_id}_summary.txt")
    
    # Initialize variables
    left_ice_all = None
    right_ice_all = None
    #print("test13")
    # Open the directory and iterate over files
    with os.scandir(download_dir) as dir_entries:
        for entry in dir_entries:
            if entry.is_file() and re.match(r"^ICEfinder_(\d+)_(\d+)$", entry.name):
                # Open the current file for reading
                with open(entry.path, "r") as in_file:
                    for line in in_file:
                        if line.startswith("Description:"):
                            match = re.match(r"^Description: (Putative.*)", line)
                            if match:
                                with open(ice_all_out, "a") as iceallout_file:
                                    iceallout_file.write(f"{job_id}\t{genome_name}\t{genome_size}\t{entry.name}\t{match.group(1)}:\t")
                        elif line.startswith("Location:"):
                            match = re.match(r"^Location: (\d+)\.\.(\d+)", line)
                            if match:
                                left_ice_all, right_ice_all = int(match.group(1)), int(match.group(2))
                                ice_len = right_ice_all - left_ice_all + 1
                                with open(ice_all_out, "a") as iceallout_file:
                                    iceallout_file.write(f"{left_ice_all}..{right_ice_all}\t{ice_len}\t")
                        elif line.startswith("oriT:"):
                            orit_desc = line.split()[1]
                            with open(ice_all_out, "a") as iceallout_file:
                                iceallout_file.write(f"{orit_desc}\t")
                        elif line.startswith("GC content:"):
                            match = re.match(r"^GC content: (\d+\.\d+)", line)
                            if match:
                                gc_content = float(match.group(1))
                                gc_diff = abs(gc_content - genome_gc)
                                gc_diff = f"{gc_diff:.2f}"
                                with open(ice_all_out, "a") as iceallout_file:
                                    iceallout_file.write(f"{gc_content}\t{genome_gc}\t{gc_diff}\t0\t0\n")
    
    # Mark completion
    complete = os.path.join(tmp_path, "region_finder_finish")
    with open(complete, 'w') as fin:
        pass
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Region Finder Script")
    parser.add_argument("-a", help="Job ID")
    parser.add_argument("-t",  required=True, help="Temporary path")
    parser.add_argument("-o",  required=True, help="Output path")
    parser.add_argument('-n', type=int, default=2, help='Threads')
    args = parser.parse_args()
    region_finder(args.a, args.t, args.o,args.n)