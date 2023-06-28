import glob
import yaml
import platform
import os

configfile: "configs/config.yaml"

#specify DATADIR where data is saved and OUTPUT where results are saved
DATADIR=config['datadir']
OUTPUT=config['output']
PARTIS=config['partis']

data = glob.glob(DATADIR+"*/*.tsv")

data = [d.split('.')[0].split('/')[-2] for d in data]

wildcard_constraints:
    d = "[a-zA-Z]+\d"

rule result:
    input:
        expand(OUTPUT + "{d}/sim_1.yaml", d=data)
        
#turn tsv into fasta     
rule tsv_to_fasta:
    resources:
        mem="10G",
    threads: 2
    log: os.path.join(DATADIR, "logs", "tsv_to_fasta_{d}.log")
    input:
        tsv = DATADIR+"{d}/airr-covid-19.tsv", #TODO: change to flexible tsv name
        script = 'Download_data/tsv_to_fasta.py',
    output:
        fasta= OUTPUT + "{d}.fasta"
    shell:
        "echo " + platform.node() + " >> {log} && \
        python {input.script} -i {input.tsv} -o {output.fasta} &&>> {log}"


#cache_parameters partis
rule cache_parameters:
    resources:
        mem="20G",
    threads: 10
    log: os.path.join(DATADIR, "logs", "cache_parameters_{d}.log")
    input:
        fasta = OUTPUT + "{d}.fasta",
        script = 'partis/partis_simulation/cache_parameters.sh',
        partis = PARTIS+"bin/partis"
    output:
        out = directory(OUTPUT + "{d}/")
    shell:
        "module purge && \
        module load gcc/8.3.0 && \
        module load gsl/2.5 && \
        module load git && \
        export PATH=/home1/kavoss/anaconda2/bin:$PATH && \
        echo " + platform.node() + " >> {log} && \
        mkdir {output.out} && \
        sh {input.script} -f {input.fasta} -p {input.partis} -o {output.out} &&>> {log}"

#partition partis
rule partition:
    resources:
        mem="100G",
    threads: 10
    log: os.path.join(DATADIR, "logs", "partition_{d}.log")
    input:
        fasta = OUTPUT + "{d}.fasta",
        script = 'partis/partis_simulation/partition.sh',
        partis = PARTIS+"bin/partis",
        out = OUTPUT + "{d}/"
    output:
        out= OUTPUT + "{d}/pd.yaml"
    shell:
        "module purge && \
        module load gcc/8.3.0 && \
        module load gsl/2.5 && \
        module load git && \
        export PATH=/home1/kavoss/anaconda2/bin:$PATH && \
        echo " + platform.node() + " >> {log} && \
        sh {input.script} -f {input.fasta} -p {input.partis} -o {input.out} &&>> {log}"

#simulation partis
rule simulate:
    resources:
        mem="100G",
    threads: 10
    log: os.path.join(DATADIR, "logs", "simulate_{d}.log")
    input:
        fasta = OUTPUT + "{d}.fasta",
        script = 'partis/partis_simulation/simulate.sh',
        partis = PARTIS+"bin/partis",
        out = OUTPUT + "{d}/",
        yaml = OUTPUT + "{d}/pd.yaml"
    output:
        out= OUTPUT + "{d}/sim_1.yaml"
    shell:
        "module purge && \
        module load gcc/8.3.0 && \
        module load gsl/2.5 && \
        module load git && \
        export PATH=/home1/kavoss/anaconda2/bin:$PATH && \
        export LD_LIBRARY_PATH={input.partis}/packages/bpp-newlik/_build/lib64:$LD_LIBRARY_PATH && \
        echo " + platform.node() + " >> {log} && \
        sh {input.script} -p {input.partis} -o {input.out} &&>> {log}"

#analyze_partis_output
#rule analyze_partis_output:
#     resources:
#            mem="50G",
#     threads: 10
#     log: os.path.join(DATADIR, "logs", "analyze_partis_output_{d}.log")
#     input:
#        yaml = OUTPUT + "{d}.yaml",
#        script1 = 'analyze_partis_simulation_output/yaml_to_families.py',
#        script2 = 'analyze_partis_simulation_output/partition_to_fasta.py',
#        out = OUTPUT + "{d}/*.yaml"
#    output:
#        partitions = OUTPUT + "{d}/partitions.txt",
#        naive = OUTPUT + "{d}/naive.txt"
#    shell:
#        "echo " + platform.node() + " >> {log} && \
#        python {input.script1} {input.yaml} {output.partitions} {output.naive} &&>> {log}" && \
#        python {input.script2} {input.yaml} {output.partitions} &&>> {log}"


