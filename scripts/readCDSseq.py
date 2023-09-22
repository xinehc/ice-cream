import sys
from Bio import SeqIO 

file_name1=sys.argv[1]
file_name2=sys.argv[2]
file_name3=sys.argv[3]
file_name4=sys.argv[4]
file_name5=sys.argv[5]

a={}
for line in open(file_name2,'r'):
    a[str(line).strip().split()[0]+".."+str(line).strip().split()[1]]=str(line).strip().split()[3]
	
fileoutput2=open(file_name4,'w')

gb = SeqIO.parse(file_name1, "genbank")
for record in gb:
    with open(file_name3,'w') as fileoutput:
         SeqIO.write(record, fileoutput, "fasta")
    for feature in record.features:
        if feature.type == "CDS":
            start = int(feature.location.start)
            end = int(feature.location.end)
            start_end= str(int(start)+1)+".."+str(end)
            try:
                assert len(feature.qualifiers['translation']) == 1
                aa_seq = feature.qualifiers['translation'][0]
                for k in a:
                    if start_end==k:
                        # output CDS sequence
                        fileoutput2.write(">"+a[k]+"_1\n"+aa_seq+"\n") 
            except KeyError:
                pass

fileoutput.close()
fileoutput2.close()

fileoutput3=open(file_name5,'w')

input_seq_iterator = SeqIO.parse(file_name4, "fasta")
short_seq_iterator = (record for record in input_seq_iterator )

SeqIO.write(short_seq_iterator, fileoutput3, "fasta")
fileoutput3.close()

