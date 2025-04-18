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

import os
import re

def process_file1(tmp_folder,job_id,region_id):
    mob_over_zero = False
    file1 = os.path.join(tmp_folder,job_id,"candidate",f"candidate_region_{region_id}_t4ss","best_solution_summary.tsv")
    #print(f"{file1}")
    with open(file1, 'r') as f:
        # Skip lines starting with '#'
        lines = [line for line in f if not line.startswith('#')]

        if len(lines) != 2:
            print("Unexpected number of lines after removing '#' lines.")
            return False

        # Extract the header and data lines
        header = lines[0].strip().split('\t')
        data_row = lines[1].strip().split('\t')

        # Ensure column names are processed in a non-ordered way
        header_to_index = {name: idx for idx, name in enumerate(header)}

        # Check values starting from the second column, excluding "CONJScan/Chromosome/MOB"
        for col_name, col_index in header_to_index.items():
            if col_index == 0:
                continue
            if 'dCONJ' in col_name:
                continue
            if col_name == "CONJScan/Chromosome/MOB":
                if int(data_row[col_index]) > 0:
                    mob_over_zero = True  # Mark as True if MOB value is over 1
                continue  # Skip further processing of this column
            if int(data_row[col_index]) > 0:  # If any value is greater than 0
                return True, mob_over_zero  # Indicate we should process File 2
    return False, mob_over_zero


def process_file2(tmp_folder,job_id,region_id):
    best_hits = []
    file2 = os.path.join(tmp_folder,job_id,"candidate",f"candidate_region_{region_id}_t4ss","best_solution.tsv")
    
    with open(file2, 'r') as f:
        lines = [line for line in f if not line.startswith('#')]

        # Skip the first line, which is the header
        data_lines = lines[1:]

        for line in f:
            # Split the line into fields
            fields = re.split(r'\s+', line.strip())
            hit_id = fields[1]
            gene_name = fields[2]
            model_fqn = fields[4]
            evalue = float(fields[14])

            # Store the best hit for each hit_id based on evalue along with gene_name and model_fqn
            best_hits.append({
                'hit_id': hit_id,
                'gene_name': gene_name,
                'model_fqn': model_fqn,
                'evalue': evalue,
                'line': line.strip()
            })
    
    return best_hits

def process_conjscan(tmp_folder,job_id,region_id):
    file1_result, mob_over_zero = process_file1(tmp_folder,job_id,region_id)
    best_hits = []
    if file1_result:
        best_hits = process_file2(tmp_folder,job_id,region_id)
    
    return file1_result, mob_over_zero, best_hits

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze and score the relaxase, t4cp, Tra, Rep, and oriT of candidate regions, and further delimit the ICE/IME with direct repeats (DR).")
    parser.add_argument('-t', help="Temporary path")
    parser.add_argument('-a', help="Job ID")
    parser.add_argument('-n', help="Candidate region ID")
    args = parser.parse_args()
    process_conjscan(args.t, args.a, args.n)
    
    
