import sys

from Bio import SeqIO 
file_name1=sys.argv[1]
file_name2=sys.argv[2]
output=open(file_name1+'_append.faa','w')


prodigal={}
for record in SeqIO.parse(file_name2,'fasta'):
    genbank=str(record.id).split("_")[:-1]
    genbank.extend([str(int(record.description.split(" ")[2])-1),record.description.split(" ")[4]])
    record.id='_'.join(genbank)
    prodigal[str(record.id)]=str(record.seq.strip())
    

start=int(0)
end=int(0)
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

output.close()

