#!/usr/bin/env python3

import sys
import os
import argparse
import subprocess
from Bio import SeqIO
import glob

def oriT_scanner(job_id,i,tmp_folder,threads):
    tmp_path = os.path.join(tmp_folder, job_id)
    candidate_dir = os.path.join(tmp_path, "candidate")
    db_dir = os.environ.get('DATABASE_FOLDER')
    seq_fna = os.path.join(candidate_dir, f"candidate_region_{i}.fna")  # <seq_fna_file>
    
    # Tools or scripts
    blast_db = os.path.join(db_dir,"oriT_seq.db")
    # Parameters
    Hvaluecutoff = 0.49  # search across the whole sequence
    
    # Intermediate file
    log_file = os.path.join(tmp_path, "oriT_blastn.err")
    blastn_out = os.path.join(candidate_dir, f"oriT_{i}.blastn")
    # Result file
    final_result = os.path.join(candidate_dir, f"oriT_{i}.result")  # <oriT_blastn_whole_result>

    # Check <region_fna_file>
    if not os.path.exists(seq_fna):
        with open(log_file, "w") as error_file:
            error_file.write(f"ERROR: For {job_id}, during finding the oriT across the whole sequence, {seq_fna} file was not found!\n")
        return

    blast_cmd = ['blastn', "-db", blast_db, "-query", seq_fna, "-evalue", "0.01", "-word_size", "11", "-outfmt", "6", "-num_alignments", "5", "-num_threads", str(threads), "-out", blastn_out]
    result = subprocess.run(blast_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        # Print the error message
        print(f"Error running command: {blast_cmd}")
        print(f"Return code: {result.returncode}")
        print(f"Standard Output: {result.stdout}")
        print(f"Standard Error: {result.stderr}")

    content = []
    if os.path.exists(blastn_out):
        with open(blastn_out, "r") as out_file:
            current_best_hit = None  
            current_query = None
            for line in out_file:
                if not line.strip() or not line.split()[1].isdigit():
                    continue
                array = line.strip().split("\t")
                qseqid = array[0]
                acc = array[1]
                identities = float(array[2])
                length = int(array[3])
                bit_score = float(array[11])
                qstart = array[6]
                qend = array[7]
                sstart = array[8]
                send = array[9]

                if current_query is None or qseqid != current_query:
                    if current_best_hit:
                        content.append(current_best_hit)
                    current_best_hit = None
                    current_query = qseqid
                
                if current_best_hit is None or bit_score > current_best_hit['bit_score']:
                    oriT_fna_files = glob.glob(os.path.join(db_dir,"oriT_seq",f"{acc}.fna"))
                    oriT_fna = oriT_fna_files[0]  # Get the first match
                    seq_obj = next(SeqIO.parse(oriT_fna, "fasta"))
                    seq_length_subject = len(seq_obj)
                    coverage = round(length / seq_length_subject * 100, 1)
                    hvalue = round(identities * (length / seq_length_subject) / 100, 2)
                    
                    if hvalue > Hvaluecutoff:
                        current_best_hit = {
                            'qseqid': qseqid,
                            'sseqid': acc,
                            'identities': identities,
                            'length': length,
                            'coverage': coverage,
                            'hvalue': hvalue,
                            'qstart': qstart,
                            'qend': qend,
                            'sstart': sstart,
                            'send': send,
                            'bit_score': bit_score
                        }

            if current_best_hit:
                content.append(current_best_hit)

        if content:
            with open(final_result, "w") as result_file:
                result_file.write("qseqid\tsseqid\tidentities\tlength\tcoverage\thvalue\tqstart\tqend\tsstart\tsend\n")
                for hit in content:
                    result_file.write(f"{hit['qseqid']}\t{hit['sseqid']}\t{hit['identities']:.1f}\t"
                                      f"{hit['length']}\t{hit['coverage']}\t{hit['hvalue']:.2f}\t"
                                      f"{hit['qstart']}\t{hit['qend']}\t{hit['sstart']}\t{hit['send']}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script is used to predict the oriT across the whole candidate region using blastn(H-value >= 0.49)")
    parser.add_argument("-a", dest="opt_a", default="", help="Job id")
    parser.add_argument("-n", dest="opt_n", default="", help="Candidate region id")
    parser.add_argument("-t", dest="opt_t", default="", help="Temporary path")
    parser.add_argument('-c', type=int, default=2, help='Threads')
    args = parser.parse_args()
    oriT_scanner(args.a,args.n,args.t,args.c)
