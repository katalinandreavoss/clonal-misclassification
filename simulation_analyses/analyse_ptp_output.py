import sys
import yaml
import statistics
path1 = sys.argv[1]  
path2 = sys.argv[2] 

file1 = open(path1, 'r')
file2 = open(path2, 'w')
Lines = file1.readlines()

species=False
size=0
sizes=[]

for line in Lines:
    line=line.strip()
    if line.startswith("Species"):
        species=True
        if size!=0 and size!=1:
            sizes.append(size)
            species+=1
        size=0
    elif species and len(line) != 0:
        size+=1
        
if size!=0 and size!=1:
    sizes.append(size)
    species+=1

    
file2.write("fam_num"+"\t"+str(len(sizes))+"\n")
file2.write("median"+"\t"+str(statistics.median(sizes)))
file2.flush()
file2.close()
