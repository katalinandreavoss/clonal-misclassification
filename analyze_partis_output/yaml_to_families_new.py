import sys
import subprocess
import yaml

yaml_file = sys.argv[1]  # .yaml file from partis partition as input
partis_dir = sys.argv[2]
output_dir = sys.argv[3]

#conf = yaml.safe_load(Path(yaml_file).read_text())
conf = yaml.load(open(yaml_file), Loader=yaml.Loader)

n=0
for key,val in conf.items():
    if key == "partitions":
        n_partitions = len(key)
        for i in range(n_partitions):
            result = subprocess.run(["python", partis_dir+"/bin/parse-output.py", yaml_file, output_dir+"/partition_"+str(i)+".fasta"  "--partition-index", i], capture_output=True, text=True, check=False)
            print(result.stdout)
            #print("partition"+str(i))

        #for e in val:
        #    n = n+1
        #    c1.write("partition" + "_" + str(n) + "\n")
        #    c2.write(">partition" + "_" + str(n) + "\n")
        #    for seq, v in e.items():
        #        if seq == "input_seqs":
        #            l1 = len(v)
        #            for k in range(0, l1, 1):
        #                c1.write("fam" + str(n) + "seq" + str(k + 1) + "\n" + v[k] + "\n")
        #        if seq == "naive_seq":
        #            c2.write(v+"\n")

#c1.flush()
#c2.flush()
#c1.close()
#c2.close()