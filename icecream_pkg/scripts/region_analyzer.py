#!/usr/bin/env python3
"""
Author: 
Date: 2024-09-13
"""

import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import re
import logging
import pandas as pd
from .process_conjscan import process_conjscan

def check_gram_file(tmp_folder, job_id):
    gram_positive_tag = 0
    is_actino = False  # Default value
    gram_file_path = os.path.join(tmp_folder, job_id, f"{job_id}.gram")
    
    if os.path.exists(gram_file_path):
        gram_positive_tag = 1
        with open(gram_file_path, 'r') as gram_file:
            gram_content = gram_file.read().strip()

        if "Bacillota" in gram_content:
            is_actino = False
        elif "Actinomycetota" in gram_content:
            is_actino = True

    return gram_positive_tag, is_actino

def compare_and_insert(candidate_data, select_dr_score):
    candidate_dr_count, score_sum_5, ice_score, aic_score, at4_score, length, ime_score, int_judge = candidate_data
    #candidate_dr_count, score_sum_5, ice_score, aic_score, at4_score, length, ime_score = candidate_data
    score_sum_5_minus_int_judge = score_sum_5 - int_judge

    count = 0
    while count < len(select_dr_score) and select_dr_score[count][1]:
        if score_sum_5_minus_int_judge > (select_dr_score[count][1] - select_dr_score[count][7] ):
        #if score_sum_5 > select_dr_score[count][1]:
            break
        #elif score_sum_5 == select_dr_score[count][1]:
        elif score_sum_5_minus_int_judge == (select_dr_score[count][1] - select_dr_score[count][7] ):
            if int_judge > 0 and select_dr_score[count][7] == 0:
                break
            elif int_judge > 0 and select_dr_score[count][7] > 0:
                if ice_score > 0 and select_dr_score[count][2] == 0:
                    break
                elif aic_score > 0 and select_dr_score[count][3] == 0:
                    break
                elif at4_score > 0 and select_dr_score[count][4] == 0:
                    break
                elif ime_score > 0 and select_dr_score[count][6] == 0:
                    break
                elif length < select_dr_score[count][5]:
                    break
        count += 1

    select_dr_score.insert(count, [candidate_dr_count, score_sum_5, ice_score, aic_score, at4_score, length, ime_score, int_judge])
    #select_dr_score.insert(count, [candidate_dr_count, score_sum_5, ice_score, aic_score, at4_score, length, ime_score])

def write_ice_result_into_file(description, ice_out, ice_out_core, ice_left, ice_right, orit_desc, ice_len, gc_content, dr_desc, insert_desc,orf,ptt,ptt2,job_id):
    with open(ice_out, "w") as ICEOUT, open(ice_out_core, "w") as ICEOUT2:
        ICEOUT.write(f"Description: {description}\n"
                     f"Location: {ice_left}..{ice_right}\n"
                     f"oriT: {orit_desc}\n"
                     f"Length: {ice_len} bp\n"
                     f"GC content: {gc_content} %\n"
                     f"DR: {dr_desc}\n"
                     f"Insert: {insert_desc}\n")
        for key in sorted(orf.keys()):
            key = int(str(key).strip())
            if key >= ice_left and orf[key] <= ice_right:
                ICEOUT.write(ptt[key])
                if ptt2.get(key):
                    ICEOUT2.write(f"{job_id}\t{ptt2[key]}\n")

def judge_is_ice_ime_or_aice(ice_desc_is_ice,ice_desc_is_ime, ice_desc_is_aic,ice_tag, mob_tag, t4ss_tag, gc_diff_tag,ice_len, ice_len_max,ime_tag, orit_tag, ime_len_max, aic_tag, aice_len_min, genus_gbk, orf,download_dir,region_i, ice_left, ice_right, orit_desc,  gc_content, dr_desc, insert_desc,ptt,ptt2,job_id,logger,h,t4ss_region,int_region,dr_TF,gram_positive_tag, is_actino): #dr_TF choose from 0 or 1
    ice_out = os.path.join(download_dir, f"ICEfinder_{region_i}_{dr_TF}")  # exported final result 
    ice_out_core = os.path.join(download_dir, f"ICEfinder_{region_i}_{dr_TF}.core")  # exported core element
    if ice_tag > 0 and mob_tag > 0 and t4ss_tag > 0 and gc_diff_tag > 0 and 4000 < ice_len < ice_len_max:  # is ICE. 元件全都得有才是ICE The largest length of ICE is set to be less than 600 kb
        write_ice_result_into_file(ice_desc_is_ice, ice_out, ice_out_core, ice_left, ice_right, orit_desc, ice_len, gc_content, dr_desc, insert_desc,orf,ptt,ptt2,job_id)
    elif ime_tag > 0 and (mob_tag > 0 or orit_tag > 0) and t4ss_tag == 0 and gc_diff_tag > 0 and 1000 < ice_len < ime_len_max:  # IME characteristics
        if h > 0: # which means len(region) > 0, which means there are more than one conj regions
            t4ss_n = 0
            while t4ss_n < len(t4ss_region) and t4ss_region[t4ss_n][0]:
                # Find neighbor integrase
                ice_left_new = t4ss_region[t4ss_n][0]
                ice_right_new = t4ss_region[t4ss_n][1]
                int_n = 0
                while int_n < len(int_region) and int_region[int_n][0]:
                    if (t4ss_region[t4ss_n][0] - int_region[int_n][1]) < 50000:
                        ice_left_new = int_region[int_n][0] if int_region[int_n][0] < ice_left_new else ice_left_new
                        break
                    int_n += 1
                while int_n < len(int_region) and int_region[int_n][0]:
                    if (int_region[int_n][0] - t4ss_region[t4ss_n][1]) < 50000:
                        ice_right_new = int_region[int_n][1] if int_region[int_n][1] > ice_right_new else ice_right_new
                        int_n += 1
                    else:
                        break
                if ice_left_new < ice_left or ice_right_new > ice_right:  # T4SS + int # is ICE
                    ice_left = ice_left_new
                    ice_right = ice_right_new
                    write_ice_result_into_file(ice_desc_is_ice, ice_out, ice_out_core, ice_left, ice_right, orit_desc, ice_len, gc_content, "-", "-",orf,ptt,ptt2,job_id)
                else: # is IME
                    write_ice_result_into_file(ice_desc_is_ime, ice_out, ice_out_core, ice_left, ice_right, orit_desc, ice_len, gc_content, dr_desc, insert_desc,orf,ptt,ptt2,job_id)
                t4ss_n += 1
        else: # is IME
            write_ice_result_into_file(ice_desc_is_ime, ice_out, ice_out_core, ice_left, ice_right, orit_desc, ice_len, gc_content, dr_desc, insert_desc,orf,ptt,ptt2,job_id)
    elif aic_tag > 0 and aice_len_min < ice_len < 60000 and gram_positive_tag > 0:  # AICE characteristics
        if is_actino == False:
            logger.info(f"This is not in Actinobacteria: genus_gbk:{genus_gbk}")
        else:
            tra_num = 0
            for key in sorted(orf.keys()):
                key = int(str(key).strip())
                if key >= ice_left and orf[key] <= ice_right:
                    ptt2_array = ptt2.get(key,"\t\t\t\t\t\t\t").split('\t')
                    #print(f"key is {key}; ptt2.getkey: {ptt2_array}")
                    if ptt2_array[6].strip() == "FtsK_SpoIIIE":
                        tra_num += 1
            if tra_num > 0:
                write_ice_result_into_file(ice_desc_is_aic, ice_out, ice_out_core, ice_left, ice_right, "-", ice_len, gc_content, dr_desc, insert_desc,orf,ptt,ptt2,job_id)

def setup_subcode_logger(log_file):
    logger = logging.getLogger('scripts.region_analyzer')  # Module-specific logger

    if logger.hasHandlers():
        logger.handlers.clear()  # Clear all handlers
    logger.setLevel(logging.INFO)
    # Create a file handler for subcode-specific logging
    file_handler = logging.FileHandler(log_file, mode='a')
    
    # Create a log format
    formatter = logging.Formatter("%(message)s")
    file_handler.setFormatter(formatter)

    logger.addHandler(file_handler)

    return logger

def safe_access(lst, idx, default):
    try:
        return lst[idx]
    except IndexError:
        return default

def calculate_gc(sequence):
    seq_length = len(sequence)
    num_a = sequence.lower().count('a')
    num_c = sequence.lower().count('c')
    num_g = sequence.lower().count('g')
    num_t = sequence.lower().count('t')
    num_o = seq_length - num_a - num_c - num_g - num_t
    gc = round((num_c + num_g) / seq_length * 100, 2)
    return gc

def run_command(command, shell=False, capture_output=False, output_file=None):
    try:
        if output_file:
            with open(output_file, 'w') as f:
                result = subprocess.run(command, shell=shell, stdout=f, stderr=subprocess.PIPE, text=True)
        else:
            result = subprocess.run(command, shell=shell, capture_output=capture_output, text=True)

        if result.returncode != 0:
            print(f"Command failed with return code: {result.returncode}")
            if not output_file:
                print(f"Standard Output: {result.stdout}")
            print(f"Standard Error: {result.stderr}")
            return False
        else:
            if capture_output and not output_file:
                print(f"Command ran successfully: {result.stdout}")
            return True
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return False

def judge_ice_ime_aice(job_id,i,rep_tag,dr_max_mismatch,genome_gc,tmp_folder,output_folder):
    #print(f"gemoe_gc: {genome_gc}")
    tmp_path = os.path.join(tmp_folder, job_id)
    download_dir = os.path.join(output_folder, job_id)
    can_region_left = None  # candidate region left coordinate
    can_region_right = None  # candidate region right coordinate
    candidate_dir = os.path.join(tmp_path, "candidate")
    seq_fna = os.path.join(tmp_path, f"{job_id}.fna")
    candidate_region_fna = os.path.join(candidate_dir, f"candidate_region_{i}.fna")
    candidate_region_faa = os.path.join(candidate_dir, f"candidate_region_{i}.faa")
    candidate_region_fea = os.path.join(candidate_dir, f"candidate_region_{i}.fea")
    ice_hmmer_out = os.path.join(candidate_dir, f"ICE_hmm_{i}.out")
    int_hmmer_out = os.path.join(tmp_path, "int_hmm.out")
    orit_out = os.path.join(candidate_dir, f"oriT_{i}.result")
    direct_repeat_out = os.path.join(candidate_dir, f"DR_vmatch_{i}.out")
    log_file = os.path.join(tmp_path, "run_region_analyzer.log")
    ice_tag = 0
    aic_tag = 0
    at4_tag = 0
    ime_tag = 0
    logger = setup_subcode_logger(log_file)

    db_dir = os.environ.get('DATABASE_FOLDER')
    locate_db = os.path.join(db_dir,"locate_ice_element_hmm_name.txt")
    locate_db_df = pd.read_csv(locate_db, sep="\t", header=0)
    locate_db_df_dict = locate_db_df.set_index('hmm_name')[['element', 'index']].to_dict(orient='index')
    mob = [key for key, value in locate_db_df_dict.items() if value['element'] == 'Relaxase']
    logger.info(f"mob is {mob}")
    t4ss = [key for key, value in locate_db_df_dict.items() if value['element'] == 'T4SS']

    gram_positive_tag, is_actino = check_gram_file(tmp_folder, job_id)

    conjscan_t4ss,conjscan_mob,annotation_conjscan_andhmm = process_conjscan(tmp_folder,job_id,i)
    annotation = {}

    threshold_file = os.path.join(db_dir,"threshold") 
    threshold_evalue = {}
    threshold_cvalue = {}
    with open(threshold_file, 'r') as THRE:
        next(THRE)  # Skip header
        for line in THRE:
            array = line.strip().split('\t')
            if not array[0] or array[0] == "Domain name":
                continue
            domain_name = array[0].replace(' ', '')
            threshold_evalue[domain_name] = float(array[1])
            threshold_cvalue[domain_name] = float(array[2])
    best_hit = {}
    if os.path.exists(ice_hmmer_out):
        with open(ice_hmmer_out, 'r') as ICEOUT:
            for line in ICEOUT:
                array = re.split(r'\s+', line)  #  matches any sequence of whitespace characters.
                match = re.search(r'(\S+)\s+.+\s+gi\|(\d+)\|\s+.+\s+(\d+)', line) or re.search(r'(\S+)\s+.+\s+(\d+)_\d+\s+\-\s+(\d+)', line)
                if not match:
                    continue
            
                target_name = match.group(1)
                gi_number = match.group(2)
                length = int(match.group(3))
                array = line.strip().split()
                evalue = float(array[6])  
                cvalue = float(array[11])
                if evalue > threshold_evalue.get(target_name, float('inf')):
                    continue
                if cvalue > threshold_cvalue.get(target_name, float('inf')):
                    continue
                
                if gi_number not in best_hit:
                    best_hit[gi_number] = {'domain': target_name, 'evalue': evalue, 'length': length}
                else:
                    if evalue < best_hit[gi_number]['evalue']:
                        best_hit[gi_number] = {'domain': target_name, 'evalue': evalue, 'length': length}

        for gi, data in best_hit.items():
            target_name = data['domain']
            length = data['length']
            matched = next((item for item in annotation_conjscan_andhmm if item['hit_id'] == gi), None)
            if matched:
                annotation[gi] = matched['gene_name']
            else:
                if target_name in mob:
                    if length >= 199:
                        annotation[gi] = target_name
                    else:
                        annotation[gi] = ""
                else:
                    if conjscan_t4ss and target_name in t4ss:
                        annotation[gi] = ""
                    else:
                        annotation[gi] = target_name


    if os.path.exists(int_hmmer_out):
        with open(int_hmmer_out, 'r') as INTOUT:
            gi_tmp = 0
            for line in INTOUT:
                array = re.split(r'\s+', line)
                match = re.search(r'(\S+)\s+.+\s+gi\|(\d+)\|\s+.+\s+(\d+)', line) or re.search(r'(\S+)\s+.+\s+(\d+)_\d+\s+\-\s+(\d+)', line)
                if match:
                    target_name = match.group(1)
                    gi_number = match.group(2)
                    length = int(match.group(3))
                    # if gi_number == gi_tmp:
                    #     continue
                    # gi_tmp = gi_number
                    # annotation[gi_tmp] = target_name
                    # if length >= 150:
                    #     annotation[gi_tmp] = target_name
                    # elif length < 150:
                    #     annotation[gi_tmp] = ""
                    evalue = float(array[6])  
                    cvalue = float(array[11])
                    #if evalue > threshold_evalue.get(target_name, float('inf')):
                    #    continue
                    #if cvalue > threshold_cvalue.get(target_name, float('inf')):
                    #    continue
                    
                    if gi_number not in best_hit:
                        best_hit[gi_number] = {'domain': target_name, 'evalue': evalue, 'length': length}
                    else:
                        if evalue < best_hit[gi_number]['evalue']:
                            best_hit[gi_number] = {'domain': target_name, 'evalue': evalue, 'length': length}

        for gi, data in best_hit.items():
            annotation[gi] = data['domain']

    candidate_region_left = 0
    candidate_region_right = 0
    try:
        with open(candidate_region_fna, "r") as file:
            line1 = file.readline().strip()
            match = re.match(r".+\s+.+\s+(\d+)\.\.(\d+)$", line1)
            if match:
                candidate_region_left, candidate_region_right = int(match.group(1)), int(match.group(2))
    except IOError as e:
        print(f"open file {candidate_region_fna} error: {e}")
        raise
    logger.info(f"***********************************************************\nCandidate region{i}: candidate_region_left..candidate_region_right: {candidate_region_left}..{candidate_region_right}")

    ######################## get the oriT with the max hvalue and generate a hash recoding oriT coordinate 
    cut_off = []
    orit_array = []
    max_value = 0
    orit_coordinate = {}  # key is the left; value is the right
    orit_left = 0
    orit_tag = 0

    if os.path.exists(orit_out):
        with open(orit_out, 'r') as OUT:
            for line in OUT:
                if line.startswith("qseqid"):
                    continue
                orit_array = line.strip().split("\t")
                cut_off.append(float(orit_array[5]))

        max_value = max(cut_off) if cut_off else max_value

        if os.path.exists(candidate_region_fna):
            with open(candidate_region_fna, 'r') as OUT2:
                line1 = OUT2.readline().strip()
                match = re.search(r'.+\s+.+\s+(\d+)\.\.(\d+)$', line1)
                if match:
                    candidate_region_left = int(match.group(1))

        with open(orit_out, 'r') as OUT:
            for line in OUT:
                if line.startswith("qseqid"):
                    continue
                orit_array = line.strip().split("\t")
                if float(orit_array[5]) != max_value:
                    continue
                orit_left = int(orit_array[6]) + candidate_region_left
                orit_coordinate[orit_left] = int(orit_array[7]) + candidate_region_left  # orit_right

    ####################### score ORFs and find conjugation region
    rna_count = 0  # record number of RNA genes
    rna = {}  # record RNA coordinates if there is any. left coordinate as key, right coordinate as value
    rna_name = {}  # record RNA information if there is any. left coordinate as key, information as value
    rna_strand = {}  # record RNA strand if there is any. left coordinate as key, strand as value
    position = 0  # orf position
    ptt = {}  # protein information
    ptt2 = {}  # signature protein info
    orf = {}  # protein coordinate

    score = {}
    #score = [[0, 0, 0, 0, 0, 5, 0, 0]]  # same as feature_designate_no_info # $socre[$position][0], int; $socre[$position][1], virB4&traU; $socre[$position][2], T4CP; $socre[$position][3], MOB; $socre[$position][4], T4SS; $socre[$position][5], note; $socre[$position][6], Rep
    orf_score = [0]
    inter_feature = {}  # regions between each feature
    inter_feature_count = 0  # number of regions between each feature
    conj_region = {}  # t4ss and t4cp
    conj_region_count = 0
    inter_conj = {}
    inter_conj_count = 0
    int_region = {}
    int_region_count = 0
    distance = 0  # distance (orf number) between T4SS proteins
    max_distance = 15  # allowed max distance (orf number) +2 between distant conjugation proteins in one ICE # 
    pos_coor = {}  # coordinate of each position                                                                                                
    with open(candidate_region_fea, "r") as CANFEA:
        for line in CANFEA: #0, location; 1, strand; 2, length/aa; 3, PID; 4 ,gene; 5, synonym; 6, code; 7, cog; 8, product; [9] note
            if re.match(r'^\d+\.\.', line):
                line = line.strip()
                line_data = line.split('\t')  
                if line_data[4] == '-':
                    line_data[4] = line_data[5]
                orf_left, orf_right = 0, 0
                match = re.match(r'(\d+)\.\.(\d+)', line_data[0])
                if match:
                    orf_left, orf_right = int(match.group(1)),int(match.group(2))

                if len(line_data) > 9 and line_data[9] == "RNA":  # read RNA information
                    rna_count += 1
                    rna[orf_left] = orf_right
                    rna_name[orf_left] = f"{line_data[0]}\t{line_data[1]}\t{line_data[2]}\t{line_data[3]}\t{line_data[4]}\tRNA\t{line_data[8]}\n"
                    rna_strand[orf_left] = line_data[1]
                    continue
                else:  
                    position = len(score) + 1
                    if match:
                        pos_coor[position] = [orf_left, orf_right]
                    if distance > 0:
                        distance += 1
                        if distance > max_distance:
                            distance = 0
                            conj_region_count += 1 # every five non ICE orfs cut the region, the conj_region_count add 1
                    score[position] = {
                    'Integrase': 0,
                    'T4SS ATPase': 0,
                    'T4CP': 0,
                    'Relaxase': 0,
                    'T4SS': 0,
                    'note': "-",
                    'Replication': 0,
                    '-':0
                    }
                    if len(orf_score) <= position:
                        orf_score.extend([0] * (position - len(orf_score) + 1))
                    orf_score[position] = orf_score[position - 1]
                    if len(line_data) > 9 and line_data[9] == "INT":                     # Integrase?
                        logger.info(f"this position is INT: {position}")
                        score[position]['Integrase'] = 1
                        score[position]['note'] = "Integrase"
                        orf_score[position] += 1
                        inter_feature[inter_feature_count][1] = orf_right
                        inter_feature_count += 1
                        inter_feature[inter_feature_count] = [orf_left, None]
                        int_region[int_region_count] = [orf_left, orf_right]
                        int_region_count += 1
                    else:
                        # Conjugation proteins?
                        logger.info(f"this position is conj: {position}")
                        logger.info(f"line_data[3]: {line_data[3]}")
                        if line_data[3] in annotation:
                            element_type = locate_db_df_dict.get(annotation[line_data[3]], {}).get("element","-")
                            score[position][element_type] = 1 # should use 7 as default, the origninal perl code is wrongly assign 0 when the line_data[3] is not within annotation, the reason is that it is not meet the threhold (either evalu or cvalue). When I use 7 not 0, the rest resutls are the same with perl result, only difference is the smaller number of 0int (which is score_sum[0])
                            logger.info(f"element_type: {element_type}")
                            if element_type != "-":
                                orf_score[position]  += 1
                            score[position]["note"] = element_type
                            inter_feature[inter_feature_count][1] = orf_right
                            inter_feature_count += 1
                            inter_feature[inter_feature_count] = [orf_left, None]
                            if score[position]["note"] in ["T4CP","Relaxase","T4SS","Replication",'T4SS ATPase']: # t4ss atp 不算 # T4CP ; T4SS ; add relaxase(for IME) && Rep(for AICE)
                                if distance == 0:
                                    distance += 1
                                    inter_conj_count += 1
                                    logger.info(f"inter_conj_count : {inter_conj_count};orf_left : {orf_left};  len(inter_con): {len(inter_conj)}")
                                    inter_conj[inter_conj_count] = [orf_left, None]
                                    conj_region[conj_region_count] = [orf_left, None] # every five non ICE orfs, the conj_region append one element
                                else:
                                    distance = 1
                                inter_conj[inter_conj_count - 1][1] = orf_right
                                conj_region[conj_region_count][1] = orf_right # adjust the right end of this conj_region to the last one's right orientate
                                logger.info(f"inter_conj[inter_conj_count - 1]: {inter_conj[inter_conj_count - 1]}")
                            
                    orf[orf_left] = orf_right  # orf coordinates are saved in orf dictionary
                    ptt[orf_left] = f"{line_data[0]}\t{line_data[1]}\t{line_data[2]}\t{line_data[3]}\t{line_data[4]}\t{score[position]['note']}\t{annotation.get(line_data[3], 'Unknown')}\t{line_data[8]}\n"

                    if score[position]["note"] is not None and score[position]['note'] not in ["-", ""]:
                        ptt2[orf_left] = f"{orf_left}..{orf_right}\t{line_data[1]}\t{line_data[2]}\t{line_data[3]}\t{line_data[4]}\t{score[position]['note']}\t{annotation.get(line_data[3], 'Unknown')}"
                        logger.info(f"Candidate region{i}:core_protein: {ptt2[orf_left]}")
                        #print(f"Candidate region{i}:core_protein: {ptt2[orf_left]}")
            elif re.match(r'^#Candidate region: (\d+)\.\.(\d+)', line):
                match = re.match(r'^#Candidate region: (\d+)\.\.(\d+)', line)
                if match:
                    can_region_left, can_region_right = int(match.group(1)), int(match.group(2))
                    inter_conj[inter_conj_count] = [can_region_left, None]
                    inter_feature[inter_feature_count] = [can_region_left, None]

    ### to count the T4SS components in the whole candidate region
    int_tag = 0
    mob_tag = 0
    t4ss_com_num = 0  # except t4ss ATPase
    t4ss_pro_left_name = {}
    t4ss_pro_left_right = {}
    t4ss_region = []  # t4ss_region[0][0]..t4ss_region[0][1] is the first qualified t4ss region in the candidate region
    h = 0  # the qualified t4ss number in the candidate region

    for orf_key in sorted(orf.keys()):
        orf_key = int(str(orf_key).strip())  # Remove newline character if any
        if orf_key >= can_region_left and orf[orf_key] <= can_region_right:
            ptt2_array_tmp = ptt2.get(orf_key,"unknown\t\t\t\t\t\t\t").split('\t')
            marker_name = ptt2_array_tmp[5]
            sign_pro_name = ptt2_array_tmp[6].strip()

            if marker_name == "Integrase":
                int_tag += 1
            if marker_name == "Relaxase":
                mob_tag += 1
            if marker_name == "T4SS" or marker_name == "T4SS ATPase":
                t4ss_pro_left_name[orf_key] = sign_pro_name
                t4ss_pro_left_right[orf_key] = orf[orf_key]
                t4ss_com_num += 1  # candidate t4ss components number in the whole candidate region before co-localization

    if t4ss_com_num >= 5:  # which means there will be more than one conj regions, which mean conj_region_count over 1
        # t4ss co-localization
        t4ss_tag = 0
        region = []
        region_name = []
        i_t4ss = 0
        region_num = 0
        mpf_distance = 10000  # the distance between each Mpf gene <= 10k
        t4ss_core = 5

        for key in sorted(t4ss_pro_left_right.keys()):
            if i_t4ss == 0:  # i_t4ss serves as a tag
                if t4ss_pro_left_name[key] not in ["AAA_10", "TraC_F_IV", "CagE_TrbE_VirB"]:  # First component should not be T4SS ATPase
                    region.append([key, int(t4ss_pro_left_right[key])])
                    region_name.append(t4ss_pro_left_name[key])
                    i_t4ss += 1
                else:
                    continue
            else:
                if int(key) <= (region[region_num][1] + mpf_distance):
                    region[region_num][1] = int(t4ss_pro_left_right[key])
                    region_name[region_num] += f"\t{t4ss_pro_left_name[key]}"
                    i_t4ss += 1
                else:
                    if t4ss_pro_left_name[key] not in ["AAA_10", "TraC_F_IV", "CagE_TrbE_VirB"]:
                        region_num += 1  # Count the next T4SS region          
                        if len(region) <= region_num:
                            region.extend([None] * (region_num - len(region) + 1))  # Expand the region list
                            region_name.extend([None] * (region_num - len(region_name) + 1))  # Expand the region_name list
                        region[region_num] = [key, int(t4ss_pro_left_right[key])]
                        region_name[region_num] = t4ss_pro_left_name[key]
                        i_t4ss = 1
                    else:
                        continue

        u = 0
        while u < len(region) and region[u][0]:
            u += 1
            arr_t4ss_com = region_name[u-1].split("\t")
            co_t4ss_com_num = len(arr_t4ss_com)
            candidate_t4ss_region_desc = f"Candidate region{i}:#Candidate T4SS Region {u}\t{region[u-1][0]}..{region[u-1][1]};co_t4ss_com_num:{co_t4ss_com_num}\n {region_name[u-1]}\n"
            
            if co_t4ss_com_num >= t4ss_core:
                t4ss_region_desc = f"Candidate region{i}:#T4SS Region {u}\t{region[u-1][0]}..{region[u-1][1]};co_t4ss_com_num:{co_t4ss_com_num}\n {region_name[u-1]}\n"
                #print(t4ss_region_desc)   
                if len(t4ss_region) <= h:
                    t4ss_region.append([None, None])       
                t4ss_region[h][0] = region[u-1][0]
                t4ss_region[h][1] = region[u-1][1]
                #print(f"test:t4ss_region[{h}][0]..t4ss_region[{h}][1]:{t4ss_region[h][0]}..{t4ss_region[h][1]}")
                h += 1


    inter_conj[inter_conj_count][1] = can_region_right
    inter_feature[inter_feature_count][1] = can_region_right

    rna_inverse = {v: k for k, v in rna.items()}
    x = 0
    while x < len(conj_region):
        logger.info(f"Candidate region{i}:conj_region{x}:conj_region[{x}][0]..conj_region[{x}][1]:{conj_region[x][0]}..{conj_region[x][1]}")
        x += 1

    x = 0
    while x < len(inter_conj):
        x += 1

    x = 0
    while x < len(inter_feature):
        x += 1
    ############################################### vmatch
    dr_min_length = 15
    if gram_positive_tag == 1 and rep_tag >= 1:
        dr_min_length = 40

    maktree_cmd = ['mkvtree', "-db", candidate_region_fna,"-indexname", candidate_region_fna, "-dna", "-pl", "-lcp", "-suf", "-tis", "-ois", "-bwt", "-bck", "-sti1"]
    success = run_command(maktree_cmd)
    if not success:
        print("Failed to execute mkvtree.")
        return
    
    if dr_max_mismatch > 0:
        vmatch_cmd = ['vmatch', "-l", str(dr_min_length), '-h', str(dr_max_mismatch), candidate_region_fna]
    else:
        vmatch_cmd = ['vmatch', "-l", str(dr_min_length), candidate_region_fna]
    success = run_command(vmatch_cmd, output_file=direct_repeat_out)
    if not success:
        print("Failed to execute vmatch.")
        return
    
    dr_count = 0  # number of DR pairs
    dr_pair = {}  # left coordinate of upstream DR as key while left coordinate of downstream DR as value
    dr_pair_reverse = {}  # left coordinate of downstream DR as key while left coordinate of upstream DR as value
    dr = {}  # dr coordinate, left coordinate as key while right coordinate as value
    logger.info(f"rna_count: {rna_count}")
    with open(direct_repeat_out, 'r') as dropt:
        for line in dropt:
            match = re.match(r'^\s+(\d+)\s+(\d+)\s+(\d+)\s+D\s+(\d+)\s+(\d+)\s+(\d+)\s+\-?(\d+)\s+.+\s+.+\s+.+', line)
            if match:
                length_u, start_u1, start_u2, length_d,start_d1, start_d2, dr_mismatch = map(int, match.groups())
                dr_u_start = start_u1 + start_u2 + can_region_left
                dr_u_end = length_u + dr_u_start - 1
                dr_d_start = start_d1 + start_d2 + can_region_left
                dr_d_end = length_d + dr_d_start - 1
                count_inter = 0
                logger.info(f"dr_u_start:{dr_u_start}")
                # Filter DR length > 200
                if length_u > 200:
                    continue

                # Filter DR pair that locates in two RNA genes
                if rna_count == 2:
                    if dr_u_start < rna[can_region_left] and dr_d_end > rna_inverse.get(can_region_right,-1):
                        continue

                # Filter DR pair locating in inter-conjugation region
                space_trash = 0
                for inter_n, inter_values in inter_conj.items():
                    inter_start, inter_end  = inter_values
                    if dr_u_start >= inter_start and dr_d_end <= inter_end:
                    #    print("3")
                    #    print(f"inter_start: {inter_start}, inter_end:{inter_end}")
                        space_trash = 1
                        break

                # Filter DR pair, select most far pairs
                if dr_u_start in dr_pair:
                    if (dr_d_start - dr_u_start) < (dr_pair[dr_u_start] - dr_u_start):
                        logger.info("4")
                        continue
                elif dr_d_start in dr_pair_reverse:
                    if (dr_d_start - dr_u_start) < (dr_d_start - dr_pair_reverse[dr_d_start]):
                        logger.info("5")
                        continue

                dr[dr_u_start] = dr_u_end
                dr[dr_d_start] = dr_d_end
                dr_pair[dr_u_start] = dr_d_start
                dr_pair_reverse[dr_d_start] = dr_u_start
                dr_count += 1

    candidate_dr_up = {}  # candidate upstream DR
    candidate_dr_down = {}  # candidate downstream DR
    candidate_dr_count = 0  # number of candidate DR pairs
    dr_num = 0  # count of DR pairs
    gap_len = 3  # gap length allowed to combine two tandem DRs
    ice_dr = []  # record dr of ICEs
    ice_dr_count = 0  # record number of dr of ICEs
    select_dr_score = []  # record score of selected dr
    # select_dr_score[n][0] = candidate dr number, select_dr_score[n][1] = feature score,
    # select_dr_score[n][2] = ice score, select_dr_score[n][3] = conj score,
    # select_dr_score[n][4] = mob score, select_dr_score[n][5] = distance to the nearest integrase,
    select_dr_number = 0  # record number score of selected dr
    length_thred = 30000  # DR max length
    logger.info(f"dr_count: {dr_count}")
    logger.info(f"dr_pair: {dr_pair}")
    logger.info(f"dr: {dr}")
    m = 0
    # Process DR and score each DR
    for dr_key in sorted(dr_pair.keys()):
        logger.info(f"dr_key: {dr_key}")
        logger.info(f"dr_num: {dr_num}")
        logger.info(f"candidate_dr_count: {candidate_dr_count}")
        m += 1
        if dr_num == 0:
            logger.info(f"round1 : {m}")
            candidate_dr_up[candidate_dr_count] = [dr_key, dr[dr_key]]
            candidate_dr_down[candidate_dr_count] = [dr_pair[dr_key], dr[dr_pair[dr_key]]]
            #logger.info(f"candidate_dr_up[candidate_dr_count][0]:{candidate_dr_up[candidate_dr_count][0]}")

        elif (dr_key <= candidate_dr_up[candidate_dr_count][1] + gap_len) and (dr_pair[dr_key] <= candidate_dr_down[candidate_dr_count][1] + gap_len):
            logger.info(f"round2: {m}")
            if dr[dr_key] > candidate_dr_up[candidate_dr_count][0]:
            # combine overlapping or neighboring DRs
                candidate_dr_up[candidate_dr_count][1] =  dr[dr_key]
                candidate_dr_down[candidate_dr_count][1] = dr[dr_pair[dr_key]]

        else:  # new DR
            # Calculate score of the regions in the last DR pair
            dr_position = []  # terminal position in the DR pair
            score_sum = [0] * 7  # int, virB4&traU, T4CP, MOB, T4SS, total Score, Rep # old perl script here is only six values
            dr_region_len = candidate_dr_down[candidate_dr_count][1] - candidate_dr_up[candidate_dr_count][0] + 1
            int_judge = 0
            for k in sorted(pos_coor.keys()): # start from 1
                logger.info(f"k is {k}")
                logger.info(f"pos_coor[k][0]: {pos_coor[k][0]}")
                if not dr_position and candidate_dr_up[candidate_dr_count][1] < pos_coor[k][0]: # up 右 小于 orf 左
                    logger.info(f"round3a : {m}")
                    dr_position.append(k) # find orf from left to right, now get the first orf that is in the region
                    score_sum[0] += score[k]['Integrase']  # int
                    score_sum[1] += score[k]['T4SS ATPase']  # virB4&traU (T4SS ATPase)
                    score_sum[2] += score[k]['T4CP']  # T4CP
                    score_sum[3] += score[k]['Relaxase']  # MOB
                    score_sum[4] += score[k]['T4SS']  # T4SS
                    score_sum[6] += score[k]['Replication']  # Rep
                    if score[k]['note'] != "-":
                        int_judge += 1
                elif dr_position and candidate_dr_down[candidate_dr_count][0] > pos_coor[k][1]: # down 左 大于 orf 右, the loop will get all the orf in this dr region
                    logger.info(f"round3b : {m}")
                    if len(dr_position) == 1:
                        dr_position.append(k)
                    elif len(dr_position) > 1:
                        dr_position[1] = k
                    score_sum[0] += score[k]['Integrase']  # int
                    score_sum[1] += score[k]['T4SS ATPase']  # virB4&traU (T4SS ATPase)
                    score_sum[2] += score[k]['T4CP']  # T4CP
                    score_sum[3] += score[k]['Relaxase']  # MOB
                    score_sum[4] += score[k]['T4SS']  # T4SS
                    score_sum[6] += score[k]['Replication']  # Rep
                    if score[k]['note'] != "-":
                        int_judge += 1
                elif dr_position :
                    break
            # oriT 坐标位于 up 左和 down 右之间
            orit_tag = 1 if orit_left > candidate_dr_up[candidate_dr_count][0] and orit_coordinate[orit_left] < candidate_dr_down[candidate_dr_count][1] else 0
            logger.info(f"score_sum: {score_sum}")
            virB_conjscan_score = t4ss_conjscan_score = ice_score = ime_score = aic_score = at4_score = 0
            score_sum_1_4 = score_sum[1] + score_sum[4]

            if conjscan_t4ss :
                virB_conjscan_score = 1
                t4ss_conjscan_score = 1
                score_sum_1_4 = 1 #score_sum_1_4 = score_sum[1] + score_sum[4]

            if gram_positive_tag == 0:
                if score_sum[4] > 2 and score_sum[1] < 4 or (t4ss_conjscan_score == 1 and virB_conjscan_score == 1):
                    if score_sum[4] > 2 and score_sum[1] < 4:
                        ice_score = (score_sum_1_4) * score_sum[2] * score_sum[3] * score_sum[0]  # T4SS ICE
                    elif t4ss_conjscan_score == 1 and virB_conjscan_score == 1:
                        ice_score = (t4ss_conjscan_score + virB_conjscan_score) * score_sum[2] * score_sum[3] * score_sum[0]
                if score_sum[4] < 2 and score_sum[1] == 0: #没有virB
                    ime_score = score_sum[0] * (orit_tag + score_sum[3])  # IME
            else:
                if (score_sum[4] >= 1 and score_sum[1] < 4) or (t4ss_conjscan_score == 1 and virB_conjscan_score == 1):
                    if (score_sum[4] >= 1 and score_sum[1] < 4):
                        ice_score = (score_sum_1_4) * score_sum[2] * score_sum[3] * score_sum[0]  # T4SS ICE (G+)
                    elif (t4ss_conjscan_score == 1 and virB_conjscan_score == 1):
                        ice_score = (t4ss_conjscan_score + virB_conjscan_score) * score_sum[2] * score_sum[3] * score_sum[0]  # T4SS ICE (G+)
                elif score_sum[4] <= 0 and score_sum[1] == 0:
                    ime_score = score_sum[0] * (orit_tag + score_sum[3])  # IME (G+)

            if score_sum[3] < 1:
                aic_score = score_sum[0] * score_sum[2] * score_sum[6]  # AICE

            if dr_region_len < length_thred and 0 < score_sum[4] <= 2:
                at4_score = score_sum[0] * score_sum[2] * (score_sum_1_4)  # T4SS-like ICE

            conj_score = (score_sum_1_4) * score_sum[2] * score_sum[3]  # Conjugation score
            mob_score = score_sum[3] * score_sum[0]  # Mobility score
            dr_position_0 = dr_position[0] if len(dr_position) > 0 else 0
            dr_position_1 = dr_position[1] if len(dr_position) > 1 else 0
            score_sum[5] = orf_score[dr_position_1] - orf_score[dr_position_0 - 1] if dr_position_0 > 0 else orf_score[dr_position_1]  # Total score (to see if it is over 5)
            if ice_score > 0:
                ice_tag = ice_score
            elif aic_score > 0:
                aic_tag = aic_score
            elif ime_score > 0:
                ime_tag = ime_score
            logger.info(f"ice_score: {ice_score}")
            logger.info(f"aic_score: {aic_score}")
            logger.info(f"ime_score: {ime_score}")
            if (ice_score + aic_score + ime_score) > 0:
                # Check DR quality
                logger.info(f"Candidate region{i}:dr_position{k}:0int;1VirB4&traU;2T4CP;3MOB;4T4SS;6Rep;orit_tag;total_score: {score_sum[0]} {score_sum[1]} {score_sum[2]} {score_sum[3]} {score_sum[4]} {score_sum[6]} {orit_tag} {score_sum[5]}\tcandidate_dr_count:{candidate_dr_count};candidate_dr_region:{candidate_dr_up[candidate_dr_count][0]}..{candidate_dr_up[candidate_dr_count][1]};candidate_dr_down_region:{candidate_dr_down[candidate_dr_count][0]}..{candidate_dr_down[candidate_dr_count][1]};termimal position:{dr_position[0]}..{dr_position[1]}\tice_score,aic_score,ime_score:{ice_score},{aic_score},{ime_score}")

                # Calculate distance to the nearest integrase
                length = 0
                kk = 0
                for k, int_values in int_region.items():
                    int_left, int_right = int_values  
                    if candidate_dr_count < len(candidate_dr_down) and candidate_dr_count < len(candidate_dr_up):
                        if int_left > candidate_dr_up[candidate_dr_count][1] and int_right < candidate_dr_down[candidate_dr_count][0]:
                            l1 = int_left - candidate_dr_up[candidate_dr_count][1]
                            l2 = candidate_dr_down[candidate_dr_count][0] - int_right
                            l = min(l1, l2)

                            if kk == 0:
                                length = l
                            elif length > l: # this line is not clear to me # This line ensures length is updated if a smaller value is found
                                length = l
                            kk += 1

                dr_record = 1
                if rna_count:
                    # If DR locates in tRNA/tmRNA
                    #logger.info(f"candidate_dr_count: {candidate_dr_count}")
                    logger.info(f"can_region_left: {can_region_left}")
                    logger.info(f"can_region_right: {can_region_right}")
                    logger.info(f"ice_dr_count: {ice_dr_count}")
                    logger.info(f"candidate_dr_count: {candidate_dr_count}")
                    logger.info(f"candidate_dr_up[candidate_dr_count][0]:{candidate_dr_up[candidate_dr_count][0]}")
                    logger.info(f"rna.get(can_region_left,-1): {rna.get(can_region_left,-1)}")
                    logger.info(f"rna_strand.get(can_region_left,None): {rna_strand.get(can_region_left,None)}")
                    logger.info(f"candidate_dr_down[candidate_dr_count][1] :{candidate_dr_down[candidate_dr_count][1] }")
                    logger.info(f"rna_inverse.get(can_region_right,-1):{rna_inverse.get(can_region_right,-1)}")
                    logger.info(f"rna_strand.get(rna_inverse.get(can_region_right,-1),None):{rna_strand.get(rna_inverse.get(can_region_right,-1),None)}")
                    ##if dr locates in tRNA/tmRNA
                    if ((candidate_dr_up[candidate_dr_count][0] < rna.get(can_region_left,-1) and rna_strand.get(can_region_left,None) == "+") or (candidate_dr_down[candidate_dr_count][1] > rna_inverse.get(can_region_right,-1) and rna_strand.get(rna_inverse.get(can_region_right,-1),None) == "-")):
                        #print("dr locates in trna true")
                        if ice_dr_count < len(ice_dr) and ice_dr[ice_dr_count][1]:
                            logger.info(f"ice_dr[ice_dr_count]: {ice_dr[ice_dr_count]}")
                            replace = 0
                            # If overlap with previous determined DR
                            if candidate_dr_down[ice_dr[ice_dr_count][0]][1] > candidate_dr_up[candidate_dr_count][0]:
                                # Select by feature score
                                if score_sum[5] > ice_dr[ice_dr_count][1]:
                                    replace = 1
                                elif score_sum[5] == ice_dr[ice_dr_count][1]:
                                    # Select by ice score
                                    if ice_score == 0 and ice_dr[ice_dr_count][2] > 0:
                                        replace = 2
                                    elif ice_score > 0 and ice_dr[ice_dr_count][2] == 0:
                                        replace = 1
                                    else:
                                        # Select by aic_score
                                        if aic_score == 0 and ice_dr[ice_dr_count][3] > 0:
                                            replace = 2
                                        elif aic_score > 0 and ice_dr[ice_dr_count][3] == 0:
                                            replace = 1
                                        else:
                                            # Select by t4ss-related aic_score
                                            if at4_score == 0 and ice_dr[ice_dr_count][4] > 0:
                                                replace = 2
                                            elif at4_score > 0 and ice_dr[ice_dr_count][4] == 0:
                                                replace = 1
                                            else:
                                                # Select by ime score
                                                if ime_score == 0 and ice_dr[ice_dr_count][6] > 0:  # add for ime
                                                    replace = 2
                                                elif ime_score > 0 and ice_dr[ice_dr_count][6] == 0:  # add for ime
                                                    replace = 1
                                                else:
                                                    if length >= ice_dr[ice_dr_count][5]:
                                                        replace = 2
                                                    else:
                                                        replace = 1

                            # If new ICE DR
                            #print(f"replace: {replace}")
                            if replace == 0:
                                if len(ice_dr) <= ice_dr_count:
                                    ice_dr.extend([[None, None, None, None, None, None, None] for _ in range(ice_dr_count - len(ice_dr) + 1)])
                                ice_dr[ice_dr_count] = [candidate_dr_count, score_sum[5], ice_score, aic_score, at4_score, length, ime_score]
                                ice_dr_count += 1
                            elif replace == 1:
                                ice_dr_count -= 1
                                ice_dr[ice_dr_count] = [candidate_dr_count, score_sum[5], ice_score, aic_score, at4_score, length, ime_score]
                                ice_dr_count += 1
                        else:
                            # Record last as ICE DR
                            if len(ice_dr) <= ice_dr_count:
                                    ice_dr.extend([[None, None, None, None, None, None, None] for _ in range(ice_dr_count - len(ice_dr) + 1)])
                            ice_dr[ice_dr_count] = [candidate_dr_count, score_sum[5], ice_score, aic_score, at4_score, length, ime_score]
                            ice_dr_count += 1
                        dr_record = 0
                        #logger.info(f"replace: {replace}")
                        #logger.info(f"ice_dr_count: {ice_dr_count}")
                logger.info(f"candidate_dr_up[candidate_dr_count][1]:{candidate_dr_up[candidate_dr_count][1]}")
                logger.info(f"dr_record: {dr_record}")
                if (dr_record == 1) and ((candidate_dr_up[candidate_dr_count][1] - candidate_dr_up[candidate_dr_count][0]) <= 200):
                    # Record last DR
                    logger.info(f"select_dr_number:{select_dr_number}")
                    candidate_data_element = [candidate_dr_count, score_sum[5], ice_score, aic_score, at4_score, length, ime_score,int_judge]
                    #candidate_data_element = [candidate_dr_count, score_sum[5], ice_score, aic_score, at4_score, length, ime_score]
                    if select_dr_number == 0:
                        select_dr_score.append(candidate_data_element)  # dr number# feature score  # ice score  # aic_s # at4_score   # distance to the nearest integrase # ime score 
                    else:
                        compare_and_insert(candidate_data_element, select_dr_score)
                    select_dr_number += 1
                    logger.info(f"select_dr_score:{select_dr_score}")

            # New DR
            candidate_dr_count += 1
            candidate_dr_up[candidate_dr_count] = [dr_key,  dr[dr_key]]
            candidate_dr_down[candidate_dr_count] = [dr_pair[dr_key],  dr[dr_pair[dr_key]]]

        dr_num += 1
    logger.info("start to test")
    logger.info(f"candidate_dr_count: {candidate_dr_count}")
    m=0
    if dr_count >= 0:        # Calculate score of the regions in the last DR pair
        m += 1
        dr_position = []  # terminal position in the DR pair
        score_sum = [0, 0, 0, 0, 0, 0, 0]  # score_sum[0]: int, score_sum[1]: virB4&traU, etc. # old perl script is wrong to have only six values
        if candidate_dr_count < len(candidate_dr_down) and candidate_dr_count < len(candidate_dr_up):
            if len(candidate_dr_down[candidate_dr_count]) > 1 and len(candidate_dr_up[candidate_dr_count]) > 0:
                dr_region_len = candidate_dr_down[candidate_dr_count][1] - candidate_dr_up[candidate_dr_count][0] + 1
            else:
                dr_region_len = 0  
        else:
            dr_region_len = 0  
        int_judge = 0
        # Find ORF positions in the DR pair
        for k in sorted(pos_coor.keys()): # start from 1
            logger.info("dr_count > 0")
            logger.info(f"k is {k}")
            logger.info(f"pos_coor[k][0]: {pos_coor[k][0]}")
            if not dr_position:
                logger.info(f"round3a+ : {m}")
                if candidate_dr_count < len(candidate_dr_up):
                    if candidate_dr_up[candidate_dr_count][1] < pos_coor[k][0]: # first one that up right is smaller than orf left
                        dr_position.append(k)
                        score_sum[0] += score[k]['Integrase']  # int
                        score_sum[1] += score[k]['T4SS ATPase']  # virB4&traU (T4SS ATPase)
                        score_sum[2] += score[k]['T4CP']  # T4CP
                        score_sum[3] += score[k]['Relaxase']  # MOB
                        score_sum[4] += score[k]['T4SS']  # T4SS
                        score_sum[6] += score[k]['Replication']  # Rep
                        if score[k]['note'] != "-":
                            int_judge += 1
            else:
                logger.info(f"round3b : {m}")
                if candidate_dr_count < len(candidate_dr_down):
                    if candidate_dr_down[candidate_dr_count][0] > pos_coor[k][1]:
                        if len(dr_position) == 1:
                            dr_position.append(k)
                        elif len(dr_position) > 1:
                            dr_position[1] = k
                        score_sum[0] += score[k]['Integrase']  # int
                        score_sum[1] += score[k]['T4SS ATPase']  # virB4&traU (T4SS ATPase)
                        score_sum[2] += score[k]['T4CP']  # T4CP
                        score_sum[3] += score[k]['Relaxase']  # MOB
                        score_sum[4] += score[k]['T4SS']  # T4SS
                        score_sum[6] += score[k]['Replication']  # Rep
                        if score[k]['note'] != "-":
                            int_judge += 1
                    else:
                        break

        logger.info(f"dr_position: {dr_position}")
        logger.info(f"score_sum: {score_sum}")
        # Calculate feature scores
        if candidate_dr_count < len(candidate_dr_down) and candidate_dr_count < len(candidate_dr_up):
            if (orit_left > candidate_dr_up[candidate_dr_count][0] and orit_coordinate[orit_left] < candidate_dr_down[candidate_dr_count][1]):
                orit_tag = 1 
            else:
                orit_tag = 0
        else:
            orit_tag = 0  
            
        ice_score = 0
        ime_score = 0
        t4ss_conjscan_score = 0
        virB_conjscan_score = 0
        score_sum_1_4 = score_sum[1] + score_sum[4]
        if conjscan_t4ss:
            virB_conjscan_score = 1
            t4ss_conjscan_score = 1
            score_sum_1_4 = 1

        if gram_positive_tag == 0:
            if (score_sum[4] > 2 and score_sum[1] < 5) or (t4ss_conjscan_score == 1 and virB_conjscan_score == 1):
                if score_sum[4] > 2 and score_sum[1] < 5 :
                    ice_score = (score_sum_1_4) * score_sum[2] * score_sum[3] * score_sum[0]  # T4SS ICE
                elif t4ss_conjscan_score == 1 and virB_conjscan_score == 1:
                    ice_score = (t4ss_conjscan_score + virB_conjscan_score) * score_sum[2] * score_sum[3] * score_sum[0]
            if score_sum[4] < 2 and score_sum[1] == 0:
                ime_score = score_sum[0] * (orit_tag + score_sum[3])  # IME
        else:
            if (score_sum[4] >= 1 and score_sum[1] < 4) or (t4ss_conjscan_score == 1 and virB_conjscan_score == 1):
                if (score_sum[4] >= 1 and score_sum[1] < 4):
                    ice_score = (score_sum_1_4) * score_sum[2] * score_sum[3] * score_sum[0]  # T4SS ICE (G+)
                elif (t4ss_conjscan_score == 1 and virB_conjscan_score == 1):
                    ice_score = (t4ss_conjscan_score + virB_conjscan_score) * score_sum[2] * score_sum[3] * score_sum[0]  # T4SS ICE (G+)
            elif score_sum[4] < 1 and score_sum[1] == 0:
                ime_score = score_sum[0] * (orit_tag + score_sum[3])  # IME (G+)

        aic_score = 0
        if score_sum[3] < 1:
            aic_score = score_sum[0] * score_sum[2] * score_sum[6]  # AICE

        at4_score = 0
        if dr_region_len < length_thred and 0 < score_sum[4] <= 2:
            at4_score = score_sum[0] * score_sum[2] * (score_sum_1_4)  # T4SS-related AICE

        dr_position_0 = dr_position[0] if len(dr_position) > 0 else 0
        dr_position_1 = dr_position[1] if len(dr_position) > 1 else 0
        score_sum[5] = orf_score[dr_position_1] - orf_score[dr_position_0 - 1] if dr_position_0 > 0 else orf_score[dr_position_1]  # Total score

        if ice_score > 0:
            ice_tag = ice_score
        elif aic_score > 0:
            aic_tag = aic_score
        elif ime_score > 0:
            ime_tag = ime_score
        logger.info(f"ice_score: aic_score; ime_tag: {ice_score}; {aic_score}; {ime_tag}")
        if (ice_score + aic_score + ime_score) > 0:
            # Check last DR quality
            logger.info(f"Candidate region{i}:last DR: dr_position{k}:0int;1VirB4&traU;2T4CP;3MOB;4T4SS;6Rep;orit_tag: {score_sum[0]} {score_sum[1]} {score_sum[2]} {score_sum[3]} {score_sum[4]} {score_sum[6]} {orit_tag} {score_sum[5]}\tcandidate_dr_count:last DR: {candidate_dr_count};candidate_dr_region: {safe_access(candidate_dr_up, candidate_dr_count, [0, 0])[0]}..{safe_access(candidate_dr_up, candidate_dr_count, [0, 0])[1]};candidate_dr_down_region:{safe_access(candidate_dr_down, candidate_dr_count, [0, 0])[0]}..{safe_access(candidate_dr_down, candidate_dr_count, [0, 0])[1]}; temnimal position:{dr_position[0]}..{dr_position[1]}\tice_score,aic_score,ime_score: {ice_score}, {aic_score}, {ime_score}")

            # Calculate distance to the nearest integrase
            length = 0
            kk = 0
            for k, int_values in int_region.items():
                int_left, int_right = int_values 
                
                if candidate_dr_count < len(candidate_dr_down) and candidate_dr_count < len(candidate_dr_up):         # Ensure the integrase region is within the candidate region's boundaries
                    if int_left > candidate_dr_up[candidate_dr_count][1] and int_right < candidate_dr_down[candidate_dr_count][0]:
                        l1 = int_left - candidate_dr_up[candidate_dr_count][1]
                        l2 = candidate_dr_down[candidate_dr_count][0] - int_right
                        l = min(l1, l2)
                        if kk == 0:
                            length = l
                        elif length > l: # this line is not clear to me # This line ensures length is updated if a smaller value is found
                            length = l
                        kk += 1

            dr_record = 1
            if rna_count:
                # If DR locates in tRNA/tmRNA 
                if candidate_dr_count < len(candidate_dr_down) and candidate_dr_count < len(candidate_dr_up):
                    if ((candidate_dr_up[candidate_dr_count][0] < rna.get(can_region_left,-1) and rna_strand.get(can_region_left,None) == "+") or
                            (candidate_dr_down[candidate_dr_count][1] > rna_inverse.get(can_region_right,-1) and rna_strand.get(rna_inverse.get(can_region_right,-1),None) == "-")):
                        if ice_dr and ice_dr[ice_dr_count][1]:
                            replace = 0
                            # If has overlap with previous determined DR
                            if candidate_dr_down[ice_dr[ice_dr_count][0]][1] > candidate_dr_up[candidate_dr_count][0]:
                                # Select by feature score
                                if score_sum[5] > ice_dr[ice_dr_count][1]:
                                    replace = 1
                                elif score_sum[5] == ice_dr[ice_dr_count][1]:
                                    # Select by ICE score
                                    if ice_score == 0 and ice_dr[ice_dr_count][2] > 0:
                                        replace = 2
                                    elif ice_score > 0 and ice_dr[ice_dr_count][2] == 0:
                                        replace = 1
                                    else:
                                        # Select by conj score
                                        if aic_score == 0 and ice_dr[ice_dr_count][3] > 0:
                                            replace = 2
                                        elif aic_score > 0 and ice_dr[ice_dr_count][3] == 0:
                                            replace = 1
                                        else:
                                            # Select by mob score
                                            if at4_score == 0 and ice_dr[ice_dr_count][4] > 0:
                                                replace = 2
                                            elif at4_score > 0 and ice_dr[ice_dr_count][4] == 0:
                                                replace = 1
                                            else:
                                                # Select by IME score
                                                if ime_score == 0 and ice_dr[ice_dr_count][6] > 0:
                                                    replace = 2
                                                elif ime_score > 0 and ice_dr[ice_dr_count][6] == 0:
                                                    replace = 1
                                                else:
                                                    if length >= ice_dr[ice_dr_count][5]:
                                                        replace = 2
                                                    else:
                                                        replace = 1
        
                            # If new ICE DR
                            if replace == 0:
                                if len(ice_dr) <= ice_dr_count:
                                    ice_dr.extend([[None, None, None, None, None, None, None] for _ in range(ice_dr_count - len(ice_dr) + 1)])
                                ice_dr[ice_dr_count] = [candidate_dr_count, score_sum[5], ice_score, aic_score, at4_score, length, ime_score]
                                ice_dr_count += 1
                            elif replace == 1:
                                ice_dr_count -= 1
                                ice_dr[ice_dr_count] = [candidate_dr_count, score_sum[5], ice_score, aic_score, at4_score, length, ime_score]
                                ice_dr_count += 1
                        else:
                            # Record last as ICE DR
                            if len(ice_dr) <= ice_dr_count:
                                    ice_dr.extend([[None, None, None, None, None, None, None] for _ in range(ice_dr_count - len(ice_dr) + 1)])
                            ice_dr[ice_dr_count] = [candidate_dr_count, score_sum[5], ice_score, aic_score, at4_score, length, ime_score]
                            ice_dr_count += 1
                        dr_record = 0

            if dr_record == 1 and candidate_dr_count < len(candidate_dr_down) and candidate_dr_count < len(candidate_dr_up):
                if (candidate_dr_up[candidate_dr_count][1] - candidate_dr_up[candidate_dr_count][0]) <= 200:
                    # Record last DR
                    candidate_data_element = [candidate_dr_count, score_sum[5], ice_score, aic_score, at4_score, length, ime_score,int_judge]
                    #candidate_data_element = [candidate_dr_count, score_sum[5], ice_score, aic_score, at4_score, length, ime_score]
                    if select_dr_number == 0:   
                        select_dr_score.append(candidate_data_element)  # dr number# feature score  # ice score  # aic_s # at4_score   # distance to the nearest integrase # ime score
                    else:
                        compare_and_insert(candidate_data_element, select_dr_score)
                    select_dr_number += 1
        
    select_dr_score_count = 0
    logger.info(f"select_dr_score:{select_dr_score}")
    logger.info(f"ice_dr:{ice_dr}")
    while select_dr_score_count < len(select_dr_score) and select_dr_score[select_dr_score_count][1]:
        k = 0
        logger.info(f"ice_dr_count: {ice_dr_count}")
        logger.info(f"select_dr_score_count: {select_dr_score_count}")
        l1 = candidate_dr_up[select_dr_score[select_dr_score_count][0]][0]
        l2 = candidate_dr_down[select_dr_score[select_dr_score_count][0]][1]
        logger.info(f"l1: {l1}; l2: {l2}")
        trash = False
        logger.info(f"ice_dr:{ice_dr}")
        while k < len(ice_dr) and ice_dr[k][1]:
            i1 = candidate_dr_up[ice_dr[k][0]][0]
            i2 = candidate_dr_down[ice_dr[k][0]][1]
            logger.info(f"i1: {i1}; i2: {i2}")
            logger.info(f"k: {k}")
            if (l1 > i1 and l1 < i2) or (l2 > i1 and l2 < i2): # if overlap with previous one
                logger.info("yes")
                trash = True
                break
            k += 1
        if not trash: # the first one is the the one highest possibility (best alignment region)
            if ice_dr_count == len(ice_dr):
                ice_dr.append([None] * 7)
            logger.info(f"select_dr_score_count: {select_dr_score_count}")
            logger.info(f"ice_dr_count: {ice_dr_count}")
            logger.info(f"Length of ice_dr: {len(ice_dr)}")
            ice_dr[ice_dr_count][0] = select_dr_score[select_dr_score_count][0]  # DR number
            ice_dr[ice_dr_count][1] = select_dr_score[select_dr_score_count][1]  # feature score
            ice_dr[ice_dr_count][2] = select_dr_score[select_dr_score_count][2]  # ICE score
            ice_dr[ice_dr_count][3] = select_dr_score[select_dr_score_count][3]  # AIC score
            ice_dr[ice_dr_count][4] = select_dr_score[select_dr_score_count][4]  # AT4 score
            ice_dr[ice_dr_count][5] = select_dr_score[select_dr_score_count][5]  # Distance to the nearest integrase
            ice_dr[ice_dr_count][6] = select_dr_score[select_dr_score_count][6]  # IME score
            ice_dr_count += 1
        select_dr_score_count += 1

    # Select inter-ICE region
    inter_ice = [[0,0]]
    inter_ice_count = 0
    inter_ice[inter_ice_count][0] = can_region_left
    ice_count = 0
    while ice_count < ice_dr_count:
        if ice_dr and len(ice_dr[ice_count]) > 0:
            inter_ice[inter_ice_count][1] = candidate_dr_up[ice_dr[ice_count][0]][0]
            inter_ice_count += 1
            inter_ice.append([candidate_dr_down[ice_dr[ice_count][0]][1], None])
        ice_count += 1

    e1 = 0
    e2 = 0

    genome_seqio_obj = SeqIO.parse(seq_fna, "fasta")
    genome_seq_obj = next(genome_seqio_obj)

    # Genus is used to exclude some non-Actinobacteria
    genome_desc = genome_seq_obj.description
    genus_gbk = genome_desc.split()[1]

    string = ""
    gc_content = 0

    ptt2_array = []
    tra_num = 0  # for AICE: FtsK_SpoIIIE
    marker_name = ""
    logger.info(f"ice_dr_count: {ice_dr_count}" )
    logger.info(f"ice_dr: {ice_dr}" )
    logger.info(f"candidate_dr_up: {candidate_dr_up}")
    # Define ICEs
    if ice_dr_count > 0:  # with dr
        e1 = 0
        k = 0
        ice_region = []

        # Define ICE with dr
        while k < ice_dr_count:
            ice_left = candidate_dr_up[ice_dr[k][0]][0]
            ice_right = candidate_dr_down[ice_dr[k][0]][1]
            if len(ice_region) <= e1:
                ice_region.extend([[None, None] for _ in range(e1 - len(ice_region) + 1)])
            ice_region[e1][0] = ice_left
            ice_region[e1][1] = ice_right
            ice_desc = ""
            logger.info(f"in ice_dr_count k:{k}")
            dr_desc = f"{candidate_dr_up[ice_dr[k][0]][0]}..{candidate_dr_up[ice_dr[k][0]][1]}\t{candidate_dr_down[ice_dr[k][0]][0]}..{candidate_dr_down[ice_dr[k][0]][1]}"
            insert_desc = "-"
            orit_desc = ""

            if ice_left <= rna.get(can_region_left,-1):
                insert_desc = rna_name.get(can_region_left,None)
            if ice_right >= rna_inverse.get(can_region_right,-1):
                insert_desc = rna_name.get(rna_inverse.get(can_region_right,-1),None)

            if insert_desc == "-":
                for key in sorted(orf.keys()):
                    key = int(str(key).strip())
                    if (key <= ice_left and orf[key] >= ice_left) or (key <= ice_right and orf[key] >= ice_right):
                        insert_desc = ptt[key]
                        break

            if insert_desc == "-":
                insert_desc = "-\n"

            ice_len = ice_right - ice_left + 1
            string = genome_seq_obj.seq[ice_left - 1:ice_right]  # Subtract 1 because Python uses 0-based indexing

            # oriT desc
            if orit_left >= ice_left and orit_coordinate[orit_left] <= ice_right:
                orit_desc = f"{orit_left}..{orit_coordinate[orit_left]}"
                orit_tag = 1
            else:
                orit_desc = "-"
                orit_tag = 0

            #orf_loop(0,0,0,{},{},[],0,orf,ice_left,ice_right)
            int_tag = 0
            mob_tag = 0
            t4ss_pro_left_name = {}
            t4ss_pro_left_right = {}
            for orf_key in sorted(orf.keys()):
                orf_key = int(str(orf_key).strip())
                if orf_key >= ice_left and orf[orf_key] <= ice_right:
                    ptt2_array_tmp = ptt2.get(orf_key,"unknown\t\t\t\t\t\t\t").split("\t")
                    marker_name = ptt2_array_tmp[5]
                    sign_pro_name = ptt2_array_tmp[6].strip()
                    if marker_name == "Integrase":
                        int_tag += 1
                    if marker_name == "Relaxase":
                        mob_tag += 1
                    if marker_name == "T4SS" or marker_name == "T4SS ATPase":
                        t4ss_pro_left_name[orf_key] = sign_pro_name
                        t4ss_pro_left_right[orf_key] = orf[orf_key]

            # T4SS co-localization
            # t4ss_colocalization(0,[],[],0,0,10000,1,2,[''])
            t4ss_tag = 0
            region = []
            region_name = []
            i_t4ss = 0
            region_num = 0
            mpf_distance = 10000  # the distance between each Mpf gene <= 10k
            t4ss_core = 2
            if gram_positive_tag == 0:  # strict restriction of T4SS core components for G- T4SS
                t4ss_core = 2
            else:  # less strict restriction of T4SS core components for G+ T4SS
                t4ss_core = 1
            if gram_positive_tag == 1 and rep_tag >= 1:  # more strict restriction of T4SS core components for AICE-like bacteria
                t4ss_core = 5

            for key in sorted(t4ss_pro_left_right.keys()):
                if i_t4ss == 0:
                    if t4ss_pro_left_name[key] not in ["AAA_10", "TraC_F_IV", "CagE_TrbE_VirB"]:  # First component should not be T4SS ATPase
                        region.append([key, int(t4ss_pro_left_right[key])])
                        region_name.append(t4ss_pro_left_name[key])
                        region[region_num][0] = key
                        region[region_num][1] = t4ss_pro_left_right[key]
                        region_name[region_num] = t4ss_pro_left_name[key]
                        i_t4ss += 1
                    else:
                        continue
                else:
                    if key <= region[region_num][1] + mpf_distance: # this was worng in perl script
                        region[region_num][1] = t4ss_pro_left_right[key]
                        region_name[region_num] += f"\t{t4ss_pro_left_name[key]}"
                        i_t4ss += 1
                    else:
                        if t4ss_pro_left_name[key] not in ["AAA_10", "TraC_F_IV", "CagE_TrbE_VirB"]:
                            region_num += 1
                            if len(region) <= region_num:
                                region.extend([None] * (region_num - len(region) + 1))  # Expand the region list
                                region_name.extend([None] * (region_num - len(region_name) + 1))  # Expand the region_name list
                            region[region_num] = [key, int(t4ss_pro_left_right[key])]
                            region_name[region_num] = t4ss_pro_left_name[key]
                            i_t4ss = 1
                        else:
                            continue

            core_num = len(region_name)
            if core_num >= t4ss_core:
                t4ss_tag = 1
            elif conjscan_t4ss: 
                t4ss_tag = 1
                if conjscan_mob and mob_tag == 0:
                    mob_tag = 1
            else:
                t4ss_tag = 0
            # GC content difference
            gc_diff_tag = 0
            gc_content = calculate_gc(string)
            gc_diff = abs(gc_content - genome_gc)
            logger.info(f"gc_content:{gc_content}")
            logger.info(f"gc_diff:{gc_diff}")
            if aic_tag > 0:  # no GC deviation for AICE
                gc_diff_tag = 1 if gc_diff >= 1 else 0
            if ice_tag > 0:
                if ice_len < 250000:  # strict GC deviation for T4SS-type ICE with size over 250 kb
                    gc_diff_tag = 1
                else:
                    gc_diff_tag = 1 if gc_diff >= 1 else 0
            if ime_tag > 0:
                if ice_len < 50000:  # strict GC deviation for IME with size over 50 kb
                    gc_diff_tag = 1
                else:
                    gc_diff_tag = 1 if gc_diff >= 1 else 0
            
            logger.info(f"\nCandidate region{i}: withDR: dr_desc: {dr_desc};ice:{ice_left}..{ice_right}")
            logger.info(f"e1:{e1}\tice_tag;aic_tag;ime_tag:{ice_tag},{aic_tag},{ime_tag}")
            logger.info(f"Candidate region{i}:int_tag;mob_tag;t4ss_tag;gc_diff_tag;ice_len:{int_tag};{mob_tag};{t4ss_tag};{gc_diff_tag};{ice_len}")
            
            #print(f"\nCandidate region{i}: withDR: dr_desc: {dr_desc};ice:{ice_left}..{ice_right}")
            #print(f"e1:{e1}\tice_tag;aic_tag;ime_tag:{ice_tag},{aic_tag},{ime_tag}")
            #print(f"Candidate region{i}:int_tag;mob_tag;t4ss_tag;gc_diff_tag;ice_len:{int_tag};{mob_tag};{t4ss_tag};{gc_diff_tag};{ice_len}")
            
            if int_tag >= 1:
                judge_is_ice_ime_or_aice("Putative ICE with T4SS","Putative IME","Putative AICE with Rep and Tra",ice_tag, mob_tag, t4ss_tag, gc_diff_tag,ice_len, 600000,ime_tag, orit_tag, 90000, aic_tag, 0, genus_gbk, orf,download_dir,i, ice_left, ice_right, orit_desc, gc_content, dr_desc, insert_desc,ptt,ptt2,job_id,logger,h,t4ss_region,int_region,0,gram_positive_tag, is_actino) 
            e1 += 1
            k += 1
    else:
        e2 = 0
        dr_desc = "-"
        insert_desc = "-"  # no DR && no tRNA
        orit_desc = ""
        for j, conj_values in conj_region.items():
            ice_left, ice_right = conj_values  # Get ice_left and ice_right from conj_region
            match_found = False
            logger.info(f"j: {j}; ice_left: {ice_left}; ice_right: {ice_right}; len(int_region): {len(int_region)}")

            int_n = 0 # First loop: Look for integrase to the left
            for int_n, int_values in int_region.items():
                int_left, int_right = int_values  # Get int_region's start and end coordinates
                logger.info(f"int_n:{int_n}; ice_left: {ice_left}; int_right: {int_right}")
                if (ice_left - int_right) < 50000:
                    ice_left = min(int_left, ice_left)
                    match_found = True
                    break

            if match_found: # Second loop: Look for integrase to the right
                for int_n, int_values in list(int_region.items())[int_n:]:
                    int_left, int_right = int_values  # Get int_region's start and end coordinates
                    logger.info(f"int_n:{int_n}; ice_right: {ice_right}; int_left: {int_left}")
                    if (int_left - ice_right) < 50000:
                        ice_right = max(int_right, ice_right)
                        break
            else:
                int_n = 0                 # Second loop: Look for integrase to the right
                for int_n, int_values in int_region.items():
                    int_left, int_right = int_values  # Get int_region's start and end coordinates
                    logger.info(f"int_n:{int_n}; ice_right: {ice_right}; int_left: {int_left}")
                    if (int_left - ice_right) < 50000:
                        ice_right = max(int_right, ice_right)
                        break
            # Define the ICE
            left_position = 0
            right_position = 0
            score_sum = [0, 0, 0, 0, 0, 0, 0]
            for k in sorted(pos_coor.keys()): # start from 1
                logger.info(f"k is : {k}; left_positino is {left_position}; ice_left: {ice_left}; pos_coor[k][0] is {pos_coor[k][0]}; ice_right:{ice_right}; pos_coor[k][1] is {pos_coor[k][1]}")
                if not left_position:
                    if ice_left == pos_coor[k][0]: 
                        left_position = k  # to find the order number of the first orf within region
                        score_sum[0] += score[k]['Integrase']  # int
                        score_sum[1] += score[k]['T4SS ATPase']  # virB4&traU (T4SS ATPase)
                        score_sum[2] += score[k]['T4CP']  # T4CP
                        score_sum[3] += score[k]['Relaxase']  # MOB
                        score_sum[4] += score[k]['T4SS']  # T4SS
                        score_sum[6] += score[k]['Replication']  # Rep
                else:
                    if ice_right >= pos_coor[k][1]:
                        right_position = k
                        score_sum[0] += score[k]['Integrase']  # int
                        score_sum[1] += score[k]['T4SS ATPase']  # virB4&traU (T4SS ATPase)
                        score_sum[2] += score[k]['T4CP']  # T4CP
                        score_sum[3] += score[k]['Relaxase']  # MOB
                        score_sum[4] += score[k]['T4SS']  # T4SS
                        score_sum[6] += score[k]['Replication']  # Rep
                    else:
                        break
                logger.info(f"k is : {k}; right_position is {right_position}; left_position is {left_position}")
                logger.info(f"score_sum: {score_sum}")
            
            int_tag = 0 # start of from-to that were wrongly in the logic flow in the perl script, there is a end in lines later
            mob_tag = 0
            t4ss_pro_left_name = {}
            t4ss_pro_left_right = {}

            for key in sorted(orf.keys()):
                orf_key = int(str(key).strip())
                if orf_key >= ice_left and orf[orf_key] <= ice_right:
                    #print(f"orf_key: {orf_key}")
                    ptt2_array_tmp = ptt2.get(orf_key,"\t\t\t\t\t\t\t").split("\t")
                    marker_name = ptt2_array_tmp[5]
                    sign_pro_name = ptt2_array_tmp[6]

                    if marker_name == "Integrase":
                        int_tag += 1
                    if marker_name == "Relaxase":
                        mob_tag += 1
                    if marker_name == "T4SS" or marker_name == "T4SS ATPase":  # MPF gene cluster; does not contain T4SS ATPase in the 原版
                        t4ss_pro_left_name[orf_key] = sign_pro_name
                        t4ss_pro_left_right[orf_key] = orf[orf_key]
            #print(f"t4ss_pro_left_name: {t4ss_pro_left_name}")
            # T4SS co-localization
            t4ss_tag = 0
            region = []
            region_name = []
            i_t4ss = 0
            region_num = 0
            mpf_distance = 10000  # The distance between each Mpf gene <= 10k
            t4ss_core = 2  # 5
            if gram_positive_tag == 0:  # Strict restriction of T4SS core components for G- T4SS
                t4ss_core = 2
            else:  # Less strict restriction of T4SS core components for G+ T4SS
                t4ss_core = 1
            if gram_positive_tag == 1 and rep_tag >= 1:  # More strict restriction of T4SS core components for AICE-like bacteria
                t4ss_core = 5

            for key in sorted(t4ss_pro_left_right.keys()):
                if i_t4ss == 0:
                    region.append([key,t4ss_pro_left_right[key]])
                    region[region_num][0] = key
                    region[region_num][1] = t4ss_pro_left_right[key]
                    region_name.append(t4ss_pro_left_name[key])
                    region_name[region_num] = t4ss_pro_left_name[key]
                    i_t4ss += 1
                else:
                    if key <= (region[region_num][1] + mpf_distance): # in perl script, this is wrongly used distance instead of mpf_distance
                        region[region_num][1] = t4ss_pro_left_right[key]
                        region_name[region_num] += f"\t{t4ss_pro_left_name[key]}"
                        i_t4ss += 1
                    else:
                        if t4ss_pro_left_name[key] not in ["AAA_10", "TraC_F_IV", "CagE_TrbE_VirB"]:    
                            region_num += 1
                            if len(region) <= region_num:
                                region.extend([None] * (region_num - len(region) + 1))  # Expand the region list
                                region_name.extend([None] * (region_num - len(region_name) + 1))  # Expand the region_name list
                            region[region_num] = [key, int(t4ss_pro_left_right[key])]
                            region_name[region_num] = t4ss_pro_left_name[key]
                            i_t4ss = 1
                        else:
                            continue

            core_num = len(region_name)
            #print(f"core_num: {core_num}; len(region_name): {len(region_name)}")
            if core_num >= t4ss_core:
                t4ss_tag = 1
            elif conjscan_t4ss: 
                t4ss_tag = 1
            else:
                t4ss_tag = 0 
            if conjscan_mob and mob_tag == 0:
                mob_tag = 1
            #print(f"t4ss_tag: {t4ss_tag}")
            if orit_left >= ice_left and orit_coordinate[orit_left] <= ice_right:
                orit_tag = 1
                orit_desc = f"{orit_left}..{orit_coordinate[orit_left]}"
            else:
                orit_tag = 0
                orit_desc = "-"  # end of from-to that were wrongly in the logic flow in the perl script, there is a start in lines before

            # Calculate feature scores
            ice_score = 0
            ice_tag = 0
            virB_conjscan_score = 0
            t4ss_conjscan_score = 0
            score_sum_1_4 = score_sum[1] + score_sum[4]
            if conjscan_t4ss:
                virB_conjscan_score = 1
                t4ss_conjscan_score = 1
                score_sum_1_4 = 1
            if t4ss_tag > 0 or (t4ss_conjscan_score == 1 and virB_conjscan_score == 1): # use t4ss_tag instead of score_sum[4] which is used in perl
                if t4ss_tag > 0:
                    ice_score = (score_sum_1_4) * score_sum[2] * score_sum[3] * score_sum[0]  # T4SS ICE = (VirB4/TraU or T4SS) + int + T4CP + MOB
                elif t4ss_conjscan_score == 1 and virB_conjscan_score == 1:
                    ice_score = (t4ss_conjscan_score + virB_conjscan_score) * score_sum[2] * score_sum[3] * score_sum[0]
            aic_score = 0
            aic_tag = 0
            if score_sum[3] < 1:
                aic_score = score_sum[0] * score_sum[2] * score_sum[6]  # AICE = int + T4CP + Rep - MOB
            ime_score = 0
            ime_tag = 0
            if score_sum[4] <= t4ss_core  and score_sum[1] == 0: # this should be t4ss_score instead of a fix number 2 which is used in the perl script
                ime_score = score_sum[0] * (orit_tag + score_sum[3])  # IME = Int + Mob - (VirB4/TraU) - full T4SS
            if ice_score > 0:
                ice_tag = ice_score
            elif aic_score > 0:
                aic_tag = aic_score
            elif ime_score > 0:
                ime_tag = ime_score
            conj_score = (score_sum_1_4) * score_sum[2] * score_sum[3]  # Conj_score = (VirB4/TraU or T4SS) + T4CP + MOB
            mob_score = score_sum[3] * score_sum[0]  # MOB score = MOB + int
            score_sum[5] = None  # Total score
            logger.info(f"ice_score: {ice_score}")
            logger.info(f"aic_score: {aic_score}")
            logger.info(f"ime_score: {ime_score}")
            logger.info(f"conj_score: {conj_score}")
            logger.info(f"mob_score: {mob_score}")
            logger.info(f"score_sum[5]: {score_sum[5]}")
            ice_len = ice_right - ice_left + 1
            string = str(genome_seq_obj.seq[ice_left - 1:ice_right])  # Subtract 1 because Python uses 0-based indexing
            gc_content = calculate_gc(string)
            
            # GC diff
            gc_diff_tag = 0
            gc_diff = abs(gc_content - genome_gc)

            if ice_tag > 0:
                if ice_len < 250000:  # Strict GC deviation for T4SS-type ICE with the size over 250 kb
                    gc_diff_tag = 1
                else:
                    gc_diff_tag = 1 if gc_diff >= 1 else 0
            if ime_tag > 0:
                if ice_len < 50000:  # Strict GC deviation for IME with the size over 50 kb
                    gc_diff_tag = 1
                else:
                    gc_diff_tag = 1 if gc_diff >= 1 else 0
            logger.info(f"\nCandidate region{i}: noDR:conj_region{j}:{conj_region[j][0]}..{conj_region[j][1]};int_region:{int_region[int_n][0]}..{int_region[int_n][1]};ice:{ice_left}..{ice_right}")
            logger.info(f"0int;1VirB4;2T4CP;3MOB;4T4SS;6Rep;orit_tag;total_score: {score_sum[0]} {score_sum[1]} {score_sum[2]} {score_sum[3]} {score_sum[4]} {score_sum[6]} {orit_tag} {score_sum[5]}\n")
            logger.info(len(int_region))
            logger.info(int_n)

            #if len(int_region):
            #    #print(f"\nCandidate region{i}: noDR:conj_region{j}:{conj_region[j][0]}..{conj_region[j][1]};int_region:{int_region[int_n][0]}..{int_region[int_n][1]};ice:{ice_left}..{ice_right}")
            #else:
            #    print(f"\nCandidate region{i}: noDR:conj_region{j}:{conj_region[j][0]}..{conj_region[j][1]};int_region:;ice:{ice_left}..{ice_right}")
            #print(f"0int;1VirB4;2T4CP;3MOB;4T4SS;6Rep;orit_tag;total_score: {score_sum[0]} {score_sum[1]} {score_sum[2]} {score_sum[3]} {score_sum[4]} {score_sum[6]} {orit_tag} {score_sum[5]}")
            #print(f"e2:{e2}\tice_tag;aic_tag;ime_tag:{ice_tag},{aic_tag},{ime_tag}")
            #print(f"Candidate region{i}:int_tag;mob_tag;t4ss_tag;gc_diff_tag;ice_len:{int_tag};{mob_tag};{t4ss_tag};{gc_diff_tag};{ice_len}")

            if int_tag >= 1:
                judge_is_ice_ime_or_aice("Putative ICE without identified DR","Putative IME without identified DR","Putative AICE without identified DR",ice_tag, mob_tag, t4ss_tag, gc_diff_tag,ice_len, 600000,ime_tag, orit_tag, 50000, aic_tag, 4000, genus_gbk, orf, download_dir,i, ice_left, ice_right, orit_desc, gc_content, dr_desc, insert_desc,ptt,ptt2,job_id,logger,0,t4ss_region,int_region,1,gram_positive_tag, is_actino) # h=0, ice_len_max = 60000, ime_len_max = 50000, dr_TF = 1, aice_len_min = 4000
            e2 += 1

    complete = os.path.join(tmp_path,"region_analyzer_finish")
    with open(complete, "w") as FIN:
        pass
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze and score the relaxase, t4cp, Tra, Rep, and oriT of candidate regions, and further delimit the ICE/IME with direct repeats (DR).")
    parser.add_argument('-a', help="Job ID")
    parser.add_argument('-n', help="Candidate region ID")
    parser.add_argument('-r', help="Rep tag")
    parser.add_argument('-d', type=int, help="Maximum mismatch for direct repeats")
    parser.add_argument('-g', type=float, help="Genome GC content")
    parser.add_argument('-t', help="Temporary path")
    parser.add_argument('-o', help="Output path")
    args = parser.parse_args()
    judge_ice_ime_aice(args.a,args.n,args.r,args.d,args.g,args.t,args.o)
