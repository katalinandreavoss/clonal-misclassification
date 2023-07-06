import glob
import yaml
import platform
import os

configfile: "configs/config.yaml"

#specify DATADIR where data is saved and OUTPUT where results are saved
DATADIR=config['datadir']
OUTPUT=config['output']
PARTIS=config['partis']
HAMPARTIS=config['ham_partis']
RAXML=config['raxml']

data = glob.glob(DATADIR+"*/*.tsv")

data = [d.split('.')[0].split('/')[-2] for d in data]

wildcard_constraints:
    d = "[a-zA-Z]+\d"

rule result:
    input:
        expand(OUTPUT + "{d}/tree_files/", d=data)
        
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
        mem_mb="100G",
	mem="180G",
    threads: 100
    log: os.path.join(DATADIR, "logs", "partition_{d}.log")
    input:
        fasta = OUTPUT + "{d}.fasta",
        script = 'partis/partis_simulation/partition.sh',
        partis = HAMPARTIS+"bin/partis",
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
        script = 'partis/partis_simulation/simulate.sh',
        dir= OUTPUT + "{d}/",
        partis = PARTIS+"bin/partis"
    output:
        out_dir = directory(OUTPUT + "{d}/simulations/"),
        out= OUTPUT + "{d}/simulations/sim_5.yaml"
    shell:
        "module purge && \
        module load gcc/8.3.0 && \
        module load gsl/2.5 && \
        module load git && \
        export PATH=/home1/kavoss/anaconda2/bin:$PATH && \
        echo " + platform.node() + " >> {log} && \
        mkdir {output.out_dir} && \
        export LD_LIBRARY_PATH=/home1/kavoss/partis_with_simulation/partis/packages/bpp-newlik/_build/lib64:$LD_LIBRARY_PATH && \
        sh {input.script} -p {input.partis} -o {input.dir} &&>> {log}"

#analyze_partis_output
rule analyze_partis_output:
     resources:
        mem="50G",
     threads: 10
     log: os.path.join(DATADIR, "logs", "analyze_partis_output_{d}.log")
     input:
        dir= OUTPUT + "{d}/simulations/",
        script = 'analyze_partis_output/simulation_analysis.sh',
        partis = PARTIS,
        sim_check = OUTPUT + "{d}/simulations/sim_5.yaml"
     output:
        out = directory(OUTPUT + "{d}/partitions/"),
        partition = OUTPUT + "{d}/partitions/sim_5_partition_0.fasta"
     shell:
        "echo " + platform.node() + " >> {log} && \
        mkdir  {output.out} && \
        sh {input.script} -d {input.dir} -p {input.partis} -o {output.out} >> {log}"

rule align_partitions:
     resources:
        mem="50G",
     threads: 10
     log: os.path.join(DATADIR, "logs", "align_partitions_{d}.log")
     input:
        partitions = OUTPUT + "{d}/partitions/",
        script = 'tree_building/align_partitions.sh',
        partition_check = OUTPUT + "{d}/partitions/sim_5_partition_0.fasta"
     output:
        out = directory(OUTPUT + "{d}/partitions_aligned/"),
        align = OUTPUT + "{d}/partitions_aligned/sim_5_partition_0_aligned.fasta"
     shell:
        "echo " + platform.node() + " >> {log} && \
        mkdir  {output.out} && \
        sh {input.script} -d {input.partitions} -o {output.out} &&>> {log}"


rule build_tree:
     resources:
        mem="50G",
     threads: 10
     log: os.path.join(DATADIR, "logs", "build_tree_{d}.log")
     input:
        partitions_aligned = OUTPUT + "{d}/partitions_aligned/",
        script = 'tree_building/build_tree.sh',
        raxml  = RAXML+"raxml-ng",
        align_check = OUTPUT + "{d}/partitions_aligned/sim_5_partition_0_aligned.fasta"
     output:
        out = directory(OUTPUT + "{d}/tree_files/"),
        tree = OUTPUT + "{d}/tree_files/sim_5_partition_0_tree_.raxml.bestTree"
     shell:
        "echo " + platform.node() + " >> {log} && \
        mkdir  {output.out} && \
        sh {input.script} -r {input.raxml} -d {input.partitions_aligned} -o {output.out} &&>> {log}"









