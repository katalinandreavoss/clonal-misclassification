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
        expand(OUTPUT + "{d}/*.yaml", d=data)
        
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


#simulation partis    
rule simulate:
    resources:
        mem="100G",
    threads: 10
    log: os.path.join(DATADIR, "logs", "simulate_{d}.log")
    input:
        fasta = OUTPUT + "{d}.fasta",
        script = 'partis/partis_simulation/simulate.sh',
        partis = PARTIS
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

#simulation partis    
rule partion_simulate:
    resources:
        mem="50G",
    threads: 10
    log: os.path.join(DATADIR, "logs", "partition_simulate_{d}.log")
    input:
        fasta = OUTPUT + "{d}.fasta",
        script = 'partis/partis_simulation/partion_simulate.sh',
        partis = PARTIS,
        out = directory(OUTPUT + "{d}/")
    output:
        out= OUTPUT + "{d}/*.yaml"
    shell:
        "module purge && \
        module load gcc/8.3.0 && \
        module load gsl/2.5 && \
        module load git && \
        export PATH=/home1/kavoss/anaconda2/bin:$PATH && \
        echo " + platform.node() + " >> {log} && \
        mkdir {output.out} && \
        sh {input.script} -f {input.fasta} -p {input.partis} -o {input.out} &&>> {log}"