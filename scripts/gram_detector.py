import sys
from Bio import SeqIO
import os
PATH = os.path.dirname(os.path.realpath(__file__))

ids = os.path.join(PATH, '..', 'data', 'all.gram.positive.gca.txt')
b=[]
for line in open(ids,'r'):
    b.append(str(line).strip())

actino = os.path.join(PATH, '..', 'data', 'all.gram.positive.gca.actinobacteria.txt')
c=[]
for line in open(actino,'r'):
    c.append(str(line).strip())


acc=sys.argv[1]
tmp=sys.argv[2]
a="_".join(acc.split("_")[:2])
if a in b:
    with open (tmp+"/"+acc+'/'+acc+'.gram', "a") as fileoutput:
        fileoutput.write('a\n')

fileoutput2 = open(tmp+"/"+acc+'/'+acc+'_all_blastn_vs_RDP_parsed.txt', "a")
if a in c:
    fileoutput2.write(acc+"     Anaerolinea\n")
else:
    fileoutput2.write(acc+"     other\n")

fileoutput2.close()
