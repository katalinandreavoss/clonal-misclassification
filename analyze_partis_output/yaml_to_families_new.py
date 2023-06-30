import sys
import subprocess
import yaml

yaml_file = sys.argv[1]  # .yaml file from partis partition as input
partis_dir = sys.argv[2]
output_dir = sys.argv[3]

conf = yaml.load(open(yaml_file), Loader=yaml.Loader)

for key,val in conf.items():
    if key == "partitions":
        n_partitions = len(val)
        for i in range(n_partitions):
            result = subprocess.Popen(["python", partis_dir+"/bin/parse-output.py", yaml_file, output_dir+"/partition_"+str(i)+".fasta", "--partition-index", str(i)])
            stdout, stderr = result.communicate()
            print(stdout)