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


rule result:
    input:
        expand(OUTPUT + "{d}/", d=data)
        
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
#TODO: make flexible for multiple loci    
rule simulate:
    resources:
        mem="10G",
    threads: 2
    log: os.path.join(DATADIR, "logs", "simulate_{d}.log")
    input:
        fasta = OUTPUT + "{d}.fasta",
        script = 'partis/partis_simulation/simulate.sh',
    output:
        out= OUTPUT + "{d}/"
    shell:
        "echo " + platform.node() + " >> {log} && \
        mkdir {output.out} && \
        sh {input.script} -f {input.fasta} -o {output.out} &&>> {log}"
