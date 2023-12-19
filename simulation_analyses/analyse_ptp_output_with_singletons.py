import sys
import statistics
path = sys.argv[1]  # .yaml file from partis partition as input

file1 = open(path+"mega_mptp.txt", 'r')
file2 = open(path+"mptp_data_singletons.txt", 'w')
Lines = file1.readlines()

species=False
size=0
sizes=[]
clone_id=0

file2.write("clone_id"+"\t"+"sequence_id"+"\n")

for line in Lines:
    line=line.strip()
    if line.startswith("Species"):
        clone_id+=1
        species=True
        if size!=0:
            sizes.append(size)
            species+=1
        size=0
    elif species and len(line) != 0:
        size+=1
        file2.write(str(clone_id) + "\t" + line + "\n")
        
if size!=0:
    sizes.append(size)
    species+=1


file2.flush()
file2.close()