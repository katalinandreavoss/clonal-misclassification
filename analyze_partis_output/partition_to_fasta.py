import sys
a=open(sys.argv[1],"r") #input is .txt file containing sequences of all clonal families (output of yaml_to_families.py)
b= a.readlines()
output=sys.argv[2]
l= len(b)
for i in range(0,l,1):
	if "partition" in b[i]:
		x= output+str(b[i]) + ".fasta"
		c=open(x,"w") #creating .fasta file for each partition
	else:
		if "seq" in b[i]:
			c.write(">"+b[i])
		else:
			c.write(b[i])


