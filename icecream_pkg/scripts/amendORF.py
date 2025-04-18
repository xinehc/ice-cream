import sys
import argparse
from Bio import SeqIO 

def amendORF(faa,faa2) :
    file_name1=faa
    file_name2=faa2

    prodigal={}
    for record in SeqIO.parse(file_name2,'fasta'):
        genbank=str(record.id).split("_")[:-1]
        genbank.extend([str(int(record.description.split(" ")[2])-1),record.description.split(" ")[4]])
        record.id='_'.join(genbank)
        prodigal[str(record.id)]=str(record.seq.strip())

    start=int(0)
    end=int(0)
    with open(f"{file_name1}_append.faa",'w') as output:
        for record in SeqIO.parse(file_name1,'fasta'):
            if (len(record.seq) > 0):
                start=int(str(record.id).split("_")[-2]) + 1
                if (start - end) >= 100:
                    for key in prodigal:
                        if int(key.split("_")[-2]) >= end and int(key.split("_")[-1]) <= start:
                            output.write('>'+str(key)+'\n')
                            output.write(prodigal[key]+'\n')
                end=int(str(record.id).split("_")[-1])
        for key in prodigal:
            if int(key.split("_")[-2]) >= end:
                output.write('>'+str(key)+'\n')
                output.write(prodigal[key]+'\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-a", help="faa", type=str, default='input',metavar='input')
    parser.add_argument("-p", help="prodigal faa", type=str, default='output',metavar='output')
    args = parser.parse_args()
    amendORF(args.a, args.p)
