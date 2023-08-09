import sys
import yaml
import statistics
path = sys.argv[1]  # .yaml file from partis partition as input

file1 = open(path+"mega_mptp.txt", 'r')
file2 = open(path+"mptp_data_singletons.txt", 'w')
Lines = file1.readlines()

species=False
size=0
sizes=[]

for line in Lines:
    line=line.strip()
    if line.startswith("Species"):
        species=True
        if size!=0:
            sizes.append(size)
            species+=1
        size=0
    elif species and len(line) != 0:
        size+=1
        
if size!=0:
    sizes.append(size)
    species+=1

    
file2.write("fam_num"+"\t"+str(len(sizes))+"\n")
file2.write("median"+"\t"+str(statistics.median(sizes)))
file2.flush()
file2.close()
