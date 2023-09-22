import sys

from Bio import SeqIO 
file_name2=str(sys.argv[1])+"_prodigal.faa"
output=open(str(sys.argv[1])+'2.faa','w')


for record in SeqIO.parse(file_name2,'fasta'):
    genbank=str(record.id).split("_")[:-1]
    genbank.extend([str(int(record.description.split(" ")[2])-1),record.description.split(" ")[4]])
    record.id='_'.join(genbank)
    output.write('>'+str(record.id)+'\n')
    output.write(str(record.seq.strip())+'\n')

output.close()

