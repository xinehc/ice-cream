import shutil
from Bio import SeqIO
import os
import argparse

def main():
    parser = argparse.ArgumentParser(description="gbff to gbk")
    parser.add_argument("-i", help="input")
    parser.add_argument("-o", help="output")
    args = parser.parse_args()
    input_dir = args.i
    output_dir = args.o
    os.makedirs(output_dir)
    gbff_files = [f for f in os.listdir(input_dir) if f.endswith('.gbff')]
    
    for gbff_file in gbff_files:
        gbff_path = os.path.join(input_dir, gbff_file)
        with open(gbff_path, 'r') as input_handle:
            for i, record in enumerate(SeqIO.parse(input_handle, "genbank"), start=1):
                output_file = os.path.join(output_dir, f"{gbff_file}.{i}.gbk")
                with open(output_file, 'w') as output_handle:
                    SeqIO.write(record, output_handle, "genbank")


if __name__ == "__main__":
    main()