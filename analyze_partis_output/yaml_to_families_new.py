import sys
import yaml

yaml_file = sys.argv[1]  # .yaml file from partis partition as input
output_dir = sys.argv[2]

conf = yaml.load(open(yaml_file), Loader=yaml.Loader)
naive= open(output_dir+"naive.fasta","w") #file for naive sequences
total= open(output_dir+"all.fasta","w")

for key,val in conf.items():
    if key == "partitions":
        n_partitions = len(val)
        if n_partitions != 1: #there should only be one partition for simulations i think
            print("more than one partition")
            break
    if key == "events":
        n = 1
        for event in val:
            family = open(output_dir + "family_"+str(n)+".fasta", "w")
            for element,value in event.items():
                if element == "naive_seq":
                    naive.write(">family_"+str(n) + "\n")
                    naive.write(value+ "\n")
                    n += 1
                if element == "input_seqs":
                    clone = 1
                    for seq in value:
                        family.write(">family_" + str(n) +"_clone_"+ str(clone) + "\n")
                        family.write(seq + "\n")
                        total.write(">family_" + str(n) +"_clone_"+ str(clone) + "\n")
                        total.write(seq + "\n")
                        clone += 1
            family.flush()
            family.close()

naive.flush()
total.flush()
naive.close()
total.close()
#script from non scratch simulations
# if key == "partitions":
#     n_partitions = len(val)
#     for i in range(n_partitions):
            #result = subprocess.Popen(["python", partis_dir+"/bin/parse-output.py", yaml_file, output_dir+"_partition_"+str(i)+".fasta", "--partition-index", str(i)])
            #stdout, stderr = result.communicate()
            #print(stdout)