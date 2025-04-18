import sys
from Bio import SeqIO 
import argparse

def amendORF2(faa):
    input=str(faa)+".fa_prodigal.faa"
    with open(str(faa)+'.fa2.faa','w') as output:
        for record in SeqIO.parse(input,'fasta'):
            genbank=str(record.id).split("_")[:-1]
            genbank.extend([str(int(record.description.split(" ")[2])-1),record.description.split(" ")[4]])
            record.id='_'.join(genbank)
            output.write('>'+str(record.id)+'\n')
            output.write(str(record.seq.strip())+'\n')

    print("finish")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-a", help="faa", type=str, default='input',metavar='input')
    args = parser.parse_args()
    amendORF2(args.a)