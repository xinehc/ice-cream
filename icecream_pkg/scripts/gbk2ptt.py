import os
import random
import re
import argparse
import logging
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation

def process_gbk(job_id, input_folder, tmp_folder,ffa_file,faa_file):
    tmp_path = os.path.join(tmp_folder, job_id)
    gbk_path = input_folder
    candidate_gbk_1 = os.path.join(gbk_path, f"{job_id}.gbk")
    candidate_gbk_2 = os.path.join(gbk_path, f"{job_id}.gb")
    candidate_gbk_3 = os.path.join(gbk_path, f"{job_id}.gbff")
    ptt_file = os.path.join(tmp_path, f"{job_id}.ptt")
    tsv_file = os.path.join(tmp_path, f"{job_id}.ptt.gi.coords")
    out = os.path.join(tmp_path, f"{job_id}.fake_gi")

    if not os.path.exists(tmp_path):
        os.makedirs(tmp_path)

    def tag(feature, tag_name):
        if tag_name in feature.qualifiers:
            return ' '.join(feature.qualifiers[tag_name])
        return '-'

    gbk = ""

    if os.path.exists(candidate_gbk_1):
        gbk = candidate_gbk_1
    elif os.path.exists(candidate_gbk_2):
        gbk = candidate_gbk_2
    elif os.path.exists(candidate_gbk_3):
        gbk = candidate_gbk_3
    else:
        #print(f"ERROR: Cannot find the GenBank file of {job_id} in {input_folder} directory!")
        return

    rand_id = int(''.join([str(random.randint(1, 9)) for _ in range(10)]))
    fake = False

    seqio_object = SeqIO.parse(gbk, "genbank")
    seq_object = next(seqio_object)
    cds_features = [f for f in seq_object.features if f.type == "CDS"]

    with open(ptt_file, "w") as ptt_output, open(tsv_file, "w") as tsv_output, open(ffa_file, 'w') as ffa_output:
        ptt_output.write(f"{seq_object.description} - 1..{len(seq_object.seq)}\n")
        ptt_output.write(f"{len(cds_features)} proteins\n")
        ptt_output.write("\t".join(["Location", "Strand", "Length", "PID", "Gene", "Synonym", "Code", "COG", "Product"]) + "\n")

        for i, f in enumerate(cds_features, start=1):
            code = '-'
            gi = 0
            db_xref = tag(f, 'db_xref')
            if 'GI:' in db_xref:
                gi = int(db_xref.split('GI:')[1].split()[0])
            
            if gi == 0:
                gi = i + rand_id
                fake = True
            
            cog = '-'
            product = tag(f, 'product')
            if "COG" in product:
                cog = product.split()[0]

            translation = tag(f, 'translation')            
            if isinstance(f.location, CompoundLocation):
                start = min(part.start for part in f.location.parts)
                end = max(part.end for part in f.location.parts)
                location = f"{start + 1}..{end}"
            else:
                start = f.location.start + 1
                end = f.location.end
                # Remove any "<" or ">" from the start and end positions
                location = f"{str(start).replace('<', '').replace('>', '')}..{str(end).replace('<', '').replace('>', '')}"

            strand = '+' if f.location.strand >= 0 else '-'
            length = len(f.location) // 3 - 1
            gene = tag(f, 'gene')
            locus_tag = tag(f, 'locus_tag')
            
            col = [location, strand, str(length), str(gi), gene, locus_tag, code, cog, product]
            ptt_output.write("\t".join(col) + "\n")
            tsv_output.write(f"{location.split('..')[0]}\t{location.split('..')[1]}\t{strand}\t{str(gi)}\n")
            if translation != "-":
                ffa_output.write(f">{str(gi)}_1\n{f.qualifiers['translation'][0]}\n")

    if fake:
        with open(out, "w") as out_file:
            out_file.write("")

    # Convert ffa to faa
    SeqIO.convert(ffa_file, "fasta", faa_file, "fasta")

    # Create completion file
    complete = os.path.join(tmp_path, "gbk2ptt_finish")
    with open(complete, "w") as fin_file:
        fin_file.write("")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process GenBank file to generate PTT, TSV, FFA, and FAA files.')
    parser.add_argument('-i', required=True, help='Input files directory, required')
    parser.add_argument('-t', required=True, help='Temporary files directory, required')
    parser.add_argument('-a', required=True, help='Job ID, required')
    parser.add_argument('-f', required=True, help='Job ID, required')
    parser.add_argument('-b', required=True, help='Job ID, required')
    args = parser.parse_args()
    process_gbk(args.a, args.i, args.t,args.f,args.b)