#!/usr/bin/env python3
"""
Author: 
Date: YYYY-MM-DD
"""
# this script will ignore the DR region
import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import logging
import subprocess
import pandas as pd
from .process_conjscan import process_conjscan

def setup_subcode_logger(log_file):
    logger = logging.getLogger('scripts.region_analyzer_s')  # Module-specific logger
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

def judge_is_ice_ime_or_aice(ice_tag_1,ice_tag_2, mob_tag, t4ss_tag, ice_len, ime_tag, orit_tag, aic_tag, orf,download_dir,region_i, ice_left, ice_right, orit_desc,  gc_content, dr_desc, insert_desc,ptt,ptt2,job_id,dr_TF): #dr_TF choose from 0 or 1
    ice_out = os.path.join(download_dir, f"ICEfinder_{region_i}_{dr_TF}")  # exported final result 
    ice_out_core = os.path.join(download_dir, f"ICEfinder_{region_i}_{dr_TF}.core")  # exported core element
    if ice_tag_1 > 0 and mob_tag > 0 and t4ss_tag > 0 and  1000 < ice_len < 600000:  # is ICE. 元件全都得有才是ICE The largest length of ICE is set to be less than 600 kb
        write_ice_result_into_file("Putative ICE without identified DR", ice_out, ice_out_core, ice_left, ice_right, orit_desc, ice_len, gc_content, dr_desc, insert_desc,orf,ptt,ptt2,job_id)
    if ice_tag_2 > 0 and t4ss_tag > 0 and  1000 < ice_len < 600000:  # is ICE. 元件全都得有才是ICE The largest length of ICE is set to be less than 600 kb
        write_ice_result_into_file("Putative conjugative region", ice_out, ice_out_core, ice_left, ice_right, orit_desc, ice_len, gc_content, dr_desc, insert_desc,orf,ptt,ptt2,job_id)
    elif ime_tag > 0 and (mob_tag > 0 or orit_tag > 0) and t4ss_tag == 0 and 1000 < ice_len < 50000:  #typical: 18-33k for MGI; typical:5-18k for Strepto; max:< 50 kb (exper);## the largest length of IME  is set to be less than 100 kb ; here, 500 kb for NCBI_9434 genome scanning
        write_ice_result_into_file("Putative IME without identified DR", ice_out, ice_out_core, ice_left, ice_right, orit_desc, ice_len, gc_content, dr_desc, insert_desc,orf,ptt,ptt2,job_id)
    elif aic_tag > 0 and 4000 < ice_len < 60000 :  # #  ## typical: 9-24k; max: <60 kb; the largest length of AICE  is set to be less than 70 kb  ## 60 kb here, 100 kb for NCBI_9434 genome scanning
        tra_num = 0
        for key in sorted(orf.keys()):
            key = int(str(key).strip())
            if key >= ice_left and orf[key] <= ice_right:
                ptt2_array = ptt2.get(key,"\t\t\t\t\t\t\t").split('\t')
                if ptt2_array[6].strip() == "FtsK_SpoIIIE":
                    tra_num += 1
        if tra_num > 0:
            write_ice_result_into_file("Putative AICE without identified DR", ice_out, ice_out_core, ice_left, ice_right, "-", ice_len, gc_content, dr_desc, insert_desc,orf,ptt,ptt2,job_id)

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

def extract_canfea(candidate_region_fea,max_distance_between_cds,annotation,locate_db_df_dict,logger,i): ####################### score ORFs and find conjugation region
    rna_count = 0  # record number of RNA genes
    rna = {}  # record RNA coordinates if there is any. left coordinate as key, right coordinate as value
    rna_name = {}  # record RNA information if there is any. left coordinate as key, information as value
    rna_strand = {}  # record RNA strand if there is any. left coordinate as key, strand as value
    position = 0  # orf position
    score = {}  # score[position][0]: int, score[position][1]: virB4&traU, etc.
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
    max_distance = max_distance_between_cds # 100 or 10  # allowed max distance (orf number) +2 between conjugation proteins in one ICE
    pos_coor = {}  # coordinate of each position
    ptt = {}  # protein information
    ptt2 = {}  # signature protein info
    orf = {}  # protein coordinate
    with open(candidate_region_fea, "r") as CANFEA:
        for line in CANFEA: #0, location; 1, strand; 2, length/aa; 3, PID; 4 ,gene; 5, synonym; 6, code; 7, cog; 8, product; [9] note
            if re.match(r'^\d+\.\.', line):
                line = line.strip()
                line_data = line.split('\t')  # 0: location, 1: strand, 2: length/aa, 3: PID, etc.
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
                            conj_region_count += 1
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
                    # Integrase?
                    if len(line_data) > 9 and line_data[9] == "INT":
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
                        if line_data[3] in annotation:
                            element_type = locate_db_df_dict.get(annotation[line_data[3]], {}).get("element","-")
                            score[position][element_type] = 1
                            if element_type != "-":
                                orf_score[position]  += 1
                            score[position]["note"] = element_type
                            inter_feature[inter_feature_count][1] = orf_right
                            inter_feature_count += 1
                            inter_feature[inter_feature_count] = [orf_left, None]
                            if score[position]["note"] in ["T4CP","Relaxase","T4SS","Replication",'T4SS ATPase']: # T4CP; T4SS; add relaxase(for IME) & Rep(for AICE)
                                if distance == 0:
                                    distance += 1
                                    inter_conj_count += 1
                                    inter_conj[inter_conj_count] = [orf_left, None]
                                    conj_region[conj_region_count] = [orf_left, None]
                                else:
                                    distance = 1
                                inter_conj[inter_conj_count - 1][1] = orf_right
                                conj_region[conj_region_count][1] = orf_right
                    
                    orf[orf_left] = orf_right  # orf coordinates are saved in orf dictionary
                    ptt[orf_left] = f"{line_data[0]}\t{line_data[1]}\t{line_data[2]}\t{line_data[3]}\t{line_data[4]}\t{score[position]['note']}\t{annotation.get(line_data[3], 'Unknown')}\t{line_data[8]}\n"
                    
                    if score[position]["note"] is not None and score[position]['note'] not in ["-", ""]:
                        ptt2[orf_left] = f"{orf_left}..{orf_right}\t{line_data[1]}\t{line_data[2]}\t{line_data[3]}\t{line_data[4]}\t{score[position]['note']}\t{annotation.get(line_data[3], 'Unknown')}"
                        logger.info(f"Candidate region{i}:core_protein: {ptt2[orf_left]}")
            elif re.match(r'^#Candidate region: (\d+)\.\.(\d+)', line):
                match = re.match(r'^#Candidate region: (\d+)\.\.(\d+)', line)
                if match:
                    can_region_left, can_region_right = int(match.group(1)), int(match.group(2))
                    inter_conj[inter_conj_count] = [can_region_left, None]
                    inter_feature[inter_feature_count] = [can_region_left, None]

    inter_conj[inter_conj_count][1] = can_region_right
    inter_feature[inter_feature_count][1] = can_region_right
    x = 0
    conj_region_count2 = 0 ## the real nubmer of conj_region
    while x < len(conj_region):
        logger.info(f"Candidate region{i}:conj_region{x}:conj_region[{x}][0]..conj_region[{x}][1]:{conj_region[x][0]}..{conj_region[x][1]}\n")
        x += 1
        conj_region_count2 += 1

    x = 0
    while x < len(inter_conj):
        x += 1

    x = 0
    while x < len(inter_feature):
        x += 1
    
    rna_inverse = {v: k for k, v in rna.items()}
    return conj_region_count2, rna_inverse, score, orf_score, inter_feature, conj_region, inter_conj, int_region, pos_coor, ptt, ptt2, orf, conj_region_count, int_region_count

def run_for_each_conj_region(score, conj_region, int_region, pos_coor, ptt, ptt2, orf, logger, conjscan_t4ss, conjscan_mob, orit_coordinate, orit_left, genome_seq_obj, dr_desc, insert_desc, job_id, download_dir,i, ime_tag_list):
    ptt2_array = []
    tra_num = 0  # for AICE: FtsK_SpoIIIE
    marker_name = ""
    e2 = 0  # for ICE/IME count
    orit_desc = ""
    for j, conj_values in conj_region.items():    # Find neighbor integrase
        ice_left, ice_right = conj_values
        match_found = False
        int_n = 0
        for int_n, int_values in int_region.items():
            int_left, int_right = int_values
            if (ice_left - int_right) < 50000:
                ice_left = min(int_left, ice_left)
                match_found = True
                break
        
        if match_found: # Second loop: Look for integrase to the right
            for int_n, int_values in list(int_region.items())[int_n:]:
                int_left, int_right = int_values  # Get int_region's start and end coordinates
                #logger.info(f"int_n:{int_n}; ice_right: {ice_right}; int_left: {int_left}")
                if (int_left - ice_right) < 50000:
                    ice_right = max(int_right, ice_right)
                    break
        else:
            int_n = 0                 # Second loop: Look for integrase to the right
            for int_n, int_values in int_region.items():
                int_left, int_right = int_values  # Get int_region's start and end coordinates
                #logger.info(f"int_n:{int_n}; ice_right: {ice_right}; int_left: {int_left}")
                if (int_left - ice_right) < 50000:
                    ice_right = max(int_right, ice_right)
                    break
        
        left_position = 0
        right_position = 0
        score_sum = [0, 0, 0, 0, 0, 0, 0]
        for k in sorted(pos_coor.keys()):
            if not left_position:
                if ice_left == pos_coor[k][0]:
                    left_position = k
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
    
        int_tag = 0
        mob_tag = 0
        t4ss_pro_left_name = {}
        t4ss_pro_left_right = {}
        distance = 0
        max_distance_ime = 10
        ime_distance = False
        for key in sorted(orf.keys()):
            orf_key = int(str(key).strip())
            if orf_key >= ice_left and orf[orf_key] <= ice_right:
                if distance > 0:
                    distance += 1
                    if distance > max_distance_ime:
                        distance = 0
                        ime_distance = True

                ptt2_array_tmp = ptt2.get(orf_key,"\t\t\t\t\t\t\t").split("\t")
                marker_name = ptt2_array_tmp[5]
                sign_pro_name = ptt2_array_tmp[6]

                if marker_name == "Integrase":
                    int_tag += 1
                if marker_name == "Relaxase":
                    mob_tag += 1
                if marker_name == "T4SS":  # MPF gene cluster; does not contain T4SS ATPase
                    t4ss_pro_left_name[orf_key] = sign_pro_name
                    t4ss_pro_left_right[orf_key] = orf[orf_key]

                ptt_array_tmp = ptt.get(orf_key,"\t\t\t\t\t\t\t").split("\t")
                if ptt_array_tmp[5] in ["T4CP","Relaxase","T4SS","Replication",'T4SS ATPase', "Integrase"]: # T4CP; T4SS; add relaxase(for IME) & Rep(for AICE)
                    distance = 1
    
        # T4SS co-localization
        t4ss_tag = 0
        region = []
        region_name = []
        i_t4ss = 0
        region_num = 0
        mpf_distance = 10000  # The distance between each Mpf gene <= 10k
        t4ss_core = 1  # Minimum required core number # this value is different from region_analyzer.py
        
        for key in sorted(t4ss_pro_left_right.keys()):
            if i_t4ss == 0:
                region.append([key,t4ss_pro_left_right[key]])
                region[region_num][0] = key
                region[region_num][1] = t4ss_pro_left_right[key]
                region_name.append(t4ss_pro_left_name[int(key)])
                region_name[region_num] = t4ss_pro_left_name[key]
                i_t4ss += 1
            else:
                if key <= (region[region_num][1] + mpf_distance):
                    region[region_num][1] = t4ss_pro_left_right[int(key)]
                    region_name[region_num] += f"\t{t4ss_pro_left_name[int(key)]}"
                    i_t4ss += 1
                else:
                    region_num += 1
                    if len(region) <= region_num:
                        region.extend([None] * (region_num - len(region) + 1))  # Expand the region list
                        region_name.extend([None] * (region_num - len(region_name) + 1))  # Expand the region_name list
                    region[region_num] = [key, int(t4ss_pro_left_right[key])]
                    region_name[region_num] = t4ss_pro_left_name[key]
                    i_t4ss = 1
        ice_tag_1 =  ice_tag_2 = aic_tag =  at4_tag =  ime_tag = 0
        core_num = len(region_name)
        if core_num >= t4ss_core:
            t4ss_tag = 1
        elif score_sum[4] >= t4ss_core:
            t4ss_tag = 1
        elif conjscan_t4ss: 
            t4ss_tag = 1
        else:
            t4ss_tag = 0 
        if conjscan_mob and mob_tag == 0:
            mob_tag = 1
        if orit_left >= ice_left and orit_coordinate[orit_left] <= ice_right:
            orit_tag = 1
            orit_desc = f"{orit_left}..{orit_coordinate[orit_left]}"
        else:
            orit_tag = 0
            orit_desc = "-"

        virB_conjscan_score = 0
        t4ss_conjscan_score = 0
        score_sum_1_4 = score_sum[1] + score_sum[4]
        if conjscan_t4ss:
            virB_conjscan_score = 1
            t4ss_conjscan_score = 1
            score_sum_1_4 = 1

        ice_len = ice_right - ice_left + 1
        string = ""
        gc_content = 0
        #logger.info(f"ice_left:{ice_left}; ice_right: {ice_right}")
        string = str(genome_seq_obj.seq[ice_left - 1:ice_right])
        gc_content = calculate_gc(string)
        # Calculate feature scores
        ice_score_1 = 0
        ice_score_2 = 0
        ice_score_3 = 0
        ice_score_4 = 0

        if score_sum[2] > 0:
            if t4ss_tag > 0 and score_sum[1] < 4:
                if score_sum_1_4 > 0:
                    ice_score_1 = score_sum_1_4 * score_sum[2] * score_sum[3] * score_sum[0]  ##t4ss ice = (VirB4/TraU or T4SS)+ int + mob	+ t4cp (must have t4cp );at least two mpf gene    
                else:
                    ice_score_1 = (t4ss_conjscan_score + virB_conjscan_score) * score_sum[2] * score_sum[3] * score_sum[0]
        
        if score_sum[2] <= 0 and score_sum[3] > 0 and core_num >= 5:
            if t4ss_tag > 0 and score_sum[1] < 4:
                if score_sum_1_4 > 0:
                    ice_score_2 = score_sum_1_4 * score_sum[3] * score_sum[0]  ### if colocalized T4SS core proteins >=5, keep as ICE even without T4CP -- output as conjugative region
                else:
                    ice_score_2 = (t4ss_conjscan_score + virB_conjscan_score) * score_sum[3] * score_sum[0]
        
        if score_sum[2] > 0 and score_sum[3] <= 0 and core_num >= 5:
            if t4ss_tag > 0 and score_sum[1] < 4:
                if score_sum_1_4 > 0:
                    ice_score_3 = score_sum_1_4 * score_sum[2] * score_sum[0]  ### if colocalized T4SS core proteins >=5, keep as ICE even without relaxase-- output as conjugative region
                else:
                    ice_score_3 = (t4ss_conjscan_score + virB_conjscan_score) * score_sum[2] * score_sum[0]
        
        if score_sum[3] <= 0 and score_sum[2] <= 0 and core_num >= 5:
            if t4ss_tag > 0 and score_sum[1] < 4:
                if score_sum_1_4 > 0:
                    ice_score_4 = score_sum_1_4 * score_sum[0]  ### if colocalized T4SS core proteins >=5, keep as ICE even without relaxase-- output as conjugative region
                else:
                    ice_score_4 = (t4ss_conjscan_score + virB_conjscan_score)  * score_sum[0] 
        # here is a typo in the old perl script that ice_score_4 shoulde be here, but ice_score_3 was defined again. 
        aic_score = ime_score = 0
        if score_sum[3] < 1:
            aic_score = score_sum[0] * score_sum[2] * score_sum[6]  ##aice = int + t4cp + Rep - mob 

        if score_sum[4] < t4ss_core and score_sum[1] == 0:
            ime_score = score_sum[0] * (orit_tag + score_sum[3])  ### IME = Int + Mob - (virB4/TraU )- full T4SS}
        
        if ice_score_1 > 0:
            ice_tag_1 = ice_score_1
        elif ice_score_2 > 0 or ice_score_3 > 0 or ice_score_4 > 0:
            ice_tag_2 = max(ice_score_2, ice_score_3, ice_score_4)
        elif aic_score > 0:
            aic_tag = aic_score
        elif ime_score > 0 and ime_distance == False:
            ime_tag = ime_score
        elif ime_score > 0 and ime_distance == True:
            ime_tag = 0
            ime_tag_list.append(-1)

        conj_score = score_sum_1_4 * score_sum[2] * score_sum[3]  ##$conj_score = (virB4/TraU or T4SS) + t4cp + mob
        mob_score = score_sum[3] * score_sum[0]  # mob_score = mob + int
        #if score[right_position]["note"] is not None and score[left_position - 1]["note"] is not None:
        #    score_sum[5] = score[right_position][5] - score[left_position - 1][5]
        #else:
        score_sum[5] = None  # Total score
        
        logger.info(f"Candidate region{i}: noDR: 0int; 1VirB4; 2T4CP; 3MOB; 4T4SS; 6Rep; orit_tag; total_score: {score_sum[0]} {score_sum[1]} {score_sum[2]} {score_sum[3]} {score_sum[4]} {score_sum[6]} {orit_tag} {score_sum[5]}")
        logger.info(f"e2:{e2}\tice_tag_1; ice_tag_2; aic_tag; ime_tag: {ice_tag_1}, {ice_tag_2}, {aic_tag}, {ime_tag}")
        logger.info(f"Candidate region{i}: int_tag; mob_tag; t4ss_tag; ice_len: {int_tag}; {mob_tag}; {t4ss_tag}; {ice_len}")
        
        if int_tag >= 1:
            judge_is_ice_ime_or_aice(ice_tag_1, ice_tag_2, mob_tag, t4ss_tag,ice_len, ime_tag, orit_tag, aic_tag,  orf, download_dir,i, ice_left, ice_right, orit_desc, gc_content, dr_desc, insert_desc,ptt,ptt2,job_id,1)  
        e2 += 1     
        return ime_tag_list

def judge_ice_ime_aice_s(job_id,i,rep_tag,dr_max_mismatch,genome_gc,tmp_folder,output_folder):
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
    log_file = os.path.join(tmp_path, "run_region_analyzer.log")
    logger = setup_subcode_logger(log_file)

    db_dir = os.environ.get('DATABASE_FOLDER')
    locate_db = os.path.join(db_dir,"locate_ice_element_hmm_name.txt")
    locate_db_df = pd.read_csv(locate_db, sep="\t", header=0)
    locate_db_df_dict = locate_db_df.set_index('hmm_name')[['element', 'index']].to_dict(orient='index')
    mob = [key for key, value in locate_db_df_dict.items() if value['element'] == 'Relaxase']
    #logger.info(f"mob is {mob}")
    t4ss = [key for key, value in locate_db_df_dict.items() if value['element'] == 'T4SS']
    
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
            gi_tmp = 0
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
            target_name = data['domain']
            length = data['length']
            if target_name in mob and length >= 199:
                annotation[gi] = target_name
            elif target_name in mob and length < 199:
                annotation[gi_tmp] = ""   
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
                if not match:
                    continue
            
                target_name = match.group(1)
                gi_number = match.group(2)
                length = int(match.group(3))
                array = line.strip().split()
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

    logger.info(f"Candidate region{i}:candidate_region_left..candidate_region_right:{candidate_region_left}..{candidate_region_right}")

    ######################## get the oriT with the max hvalue and generate a hash recoding oriT coordinate 
    cut_off = []
    orit_array = []
    max_value = 0
    orit_coordinate = {}  # key is the left; value is the right
    orit_left = 0

    orit_tag = 0
    orit_region = {}
    orit_count = 0

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
                orit_region[orit_count]=[orit_left, orit_coordinate[orit_left]]
                orit_count += 1
    genome_seqio_obj = SeqIO.parse(seq_fna, "fasta")
    genome_seq_obj = next(genome_seqio_obj)                
    dr_desc = "-"
    insert_desc = "-"  # no DR && no tRNA

    conj_region_count2, rna_inverse, score, orf_score, inter_feature, conj_region, inter_conj, int_region, pos_coor, ptt, ptt2, orf, conj_region_count, int_region_count = extract_canfea(candidate_region_fea,100,annotation,locate_db_df_dict,logger,i)
    ime_tag_list = [0]

    if  len(conj_region) > 0 :
        ime_tag_list = run_for_each_conj_region(score, conj_region, int_region, pos_coor, ptt, ptt2, orf, logger, conjscan_t4ss, conjscan_mob, orit_coordinate, orit_left, genome_seq_obj, dr_desc, insert_desc, job_id, download_dir,i, ime_tag_list)

    elif conj_region_count2 == 0 and int_region_count > 0 and orit_count > 0:  # Find MGI: only Int and oriT
        ice_left = 0
        ice_right = 0

        for j, int_values in int_region.items():
            ice_left, ice_right = int_values
            orit_n = 0
            logger.info(f"\nTo find MGI\nCandidate region {i}: conj_region_count: {conj_region_count}; int_region_count: {int_region_count}; orit_count: {orit_count}; int_region {j}: {int_region[j][0]}..{int_region[j][1]}; orit_region {orit_n}: {orit_region[orit_n][0]}..{orit_region[orit_n][1]}")
            
            for orit_n, orit_values in orit_region.items():
                orit_left, orit_right = orit_values
                if (ice_left - orit_right) < 40000:
                    ice_left = min(orit_left, ice_left)
                    break
                orit_n += 1
            
            for orit_n, orit_values in orit_region.items():
                orit_left, orit_right = orit_values
                if (orit_left - ice_right) < 40000:
                    ice_right = max(orit_right, ice_right)
                    orit_n += 1
                else:
                    break
        
        ice_len = ice_right - ice_left + 1
        string = ""
        gc_content = 0
        string = str(genome_seq_obj.seq[ice_left-1:ice_right])
        gc_content = calculate_gc(string)
        
        if orit_left >= ice_left and orit_coordinate[orit_left] <= ice_right:
            orit_desc = f"{orit_left}..{orit_coordinate[orit_left]}"
        else:
            orit_desc = "-"
        
        #logger.info(f"ice_left..ice_right: {ice_left}..{ice_right}; orit_left..orit_right: {orit_left}..{orit_coordinate[orit_left]};")
        
        if 1000 < ice_len < 50000:  # typical: 18-33k for MGI; typical: 5-18k for Strepto; max: <50 kb
            ice_out = os.path.join(download_dir, f"ICEfinder_{i}_2")  # final result
            ice_out_core = os.path.join(download_dir, f"ICEfinder_{i}_2.core")  # core element
            write_ice_result_into_file("Putative IME without identified DR", ice_out, ice_out_core, ice_left, ice_right, orit_desc, ice_len, gc_content, dr_desc, insert_desc,orf,ptt,ptt2,job_id)

    if -1 in ime_tag_list:
        logger.info("IME was there but the element is too long, thus redo the split into multi-conj-regions")
        conj_region_count2, rna_inverse, score, orf_score, inter_feature, conj_region, inter_conj, int_region, pos_coor, ptt, ptt2, orf, conj_region_count, int_region_count = extract_canfea(candidate_region_fea,10,annotation,locate_db_df_dict,logger,i)
        ime_tag_list = run_for_each_conj_region(score, conj_region, int_region, pos_coor, ptt, ptt2, orf, logger, conjscan_t4ss, conjscan_mob, orit_coordinate, orit_left, genome_seq_obj, dr_desc, insert_desc, job_id, download_dir,i, ime_tag_list)


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
    judge_ice_ime_aice_s(args.a,args.n,args.r,args.d,args.g,args.t,args.o)


