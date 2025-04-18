"""
Author: 
Date: 2024-09-13
"""
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import glob
import shutil
import logging
import gzip
from Bio.SeqFeature import FeatureLocation

def setup_subcode_logger(log_file):
    logger = logging.getLogger('parse_gbff')  # Module-specific logger
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

def clean_coordinate(coord):
    """ Clean non-numeric characters from a coordinate and handle FeatureLocation objects. """
    if isinstance(coord, FeatureLocation):
        coord_str = str(coord)
    else:
        coord_str = str(coord)
    
    # Remove <, >, and ^ characters
    coord_str = coord_str.replace('<', '').replace('>', '').replace('^', '')
    try:
        return int(coord_str)  # Try to convert to int
    except ValueError:
        return None

def parse_gbk_files(in_dir, file_extension, out_dir):
    log_file = os.path.join(out_dir, "parser.log")
    logger = setup_subcode_logger(log_file)

    if file_extension == '.gbk.gz':
        gz_files = glob.glob(os.path.join(in_dir, '*.gbk.gz'))
        for gz_file in gz_files: 
            out_file_gz = gz_file[:-3]
            with gzip.open(gz_file, 'rb') as f_in:
                with open(out_file_gz, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        input_extension = '.gbk'
    else:
        input_extension = file_extension
    list_gbff = glob.glob(os.path.join(in_dir,'*'+input_extension))
    summary = open(os.path.join(out_dir,"list.gbk.txt"),"a")
    summary_gbk = open(os.path.join(out_dir,"list.gbk.exclude.plasmid.txt"),"a")
    if list_gbff == []:
        logger.info(f"Find no genbank files (.gbff or .gbff.gz) in {in_dir}")
    else:
        Filenumber = 0
        for file_name1 in list_gbff:
            base_name = os.path.basename(file_name1) 
            acc = os.path.splitext(base_name)[0]
            tmp_path = os.path.join(out_dir, acc)
            if os.path.exists(tmp_path):
                shutil.rmtree(tmp_path)
            os.makedirs(tmp_path, mode=0o755, exist_ok=True)
            os.chmod(tmp_path, 0o755)
            gb = SeqIO.parse(os.path.join(in_dir, base_name), "genbank")
            #tmp1 = file_name1.split("_")[0]
            #tmp2 = file_name1.split("_")[1]
            #assem_name = "_".join([tmp1, tmp2])
            #logger.info(f"Processing GenBank chr record {assem_name}")
            Filenumber += 1
            if Filenumber % 10000 == 0:
                print(f"Processing the {Filenumber} record")

            with open(os.path.join(tmp_path, f"{acc}.fna"),'w') as output_handle1, open(os.path.join(tmp_path, f"{acc}.faa"),"w") as output_handle3:
                for record in gb:
                    if "plasmid" in record.description :  ### newly added
                        summary.write(base_name+"\t"+ record.id+"\tplasmid"+"\t"+record.description+"\n")   ### newly added
                        print (f"Skip! {os.path.join(in_dir, base_name)} is a plasmid.")
                    else:   ### newly added
                        summary.write(base_name+"\t"+ record.id+"\tchromosome_or_unsure"+"\t"+record.description+"\n")   ### newly added
                        summary_gbk.write(base_name+"\n")  ### newly added
                        SeqIO.write(record, output_handle1, "fasta") 
                        AA_true = 0                     # fetch CDS sequence
                        for feature in record.features:
                            if feature.type == "CDS":
                                AA_true = 1
                                ID = feature.qualifiers.get('db_xref', 'None')[0]
                                desc = feature.qualifiers.get('protein_id', 'None')[0]
                                product = feature.qualifiers.get('product', 'None')[0]
                                locus = feature.qualifiers.get('locus_tag', 'None')[0]
                                type = feature.type
                                start = clean_coordinate(feature.location.start)
                                end = clean_coordinate(feature.location.end)
                                if start is None or end is None:
                                    #print(f"Skipping feature due to invalid coordinates: start={start}, end={end}")
                                    continue
                                try:
                                    assert len(feature.qualifiers['translation']) == 1
                                    aa_seq = feature.qualifiers['translation'][0]
                                except KeyError:
                                    aa_seq = ''
                                if AA_true == 1:
                                    output_handle3.write(
                                            #">%s_%s_%s_%s_%s %s %s %s %s %s\n%s\n" %
                                            #(base_name, assem_name, record.id, start, end, desc, locus, ID, product, type, aa_seq))
                                            ">%s_%s_%s %s %s %s %s %s\n%s\n" %
                                            (record.id, start, end, desc, locus, ID, product, type, aa_seq))
                                    AA_true = 1

            logger.info('Retrieving whole genome sequences!')
            logger.info('Done!')
    summary.close()
    summary_gbk.close()
    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i",
                        help="input directory or folder of your sequences.",
                        type=str, default='input',metavar='input')
    parser.add_argument("-f",
                        help="file type or filename extension of your gbk.\n \
                            To input genbank file, set \"-f .gbk\" or \"-f .gbk.gz\"",
                        type=str, default='.gbk',metavar='.gbk or .gbk.gz')
    parser.add_argument("-o",
                        help="output directory or folder of your faa.",
                        type=str, default='output',metavar='output')    
    args = parser.parse_args()
    parse_gbk_files(args.i, args.f, args.o)

