import glob
import yaml
import platform
import os

configfile: "configs/config.yaml"

#specify DATADIR where data is saved and OUTPUT where results are saved
#DATADIR=config['datadir']
OUTPUT=config['output']
PARTIS=config['partis']
HAMPARTIS=config['ham_partis']
RAXML=config['raxml']
VQUEST=config['vquest']
RTG=config['RevertToGermline']
data = glob.glob(DATADIR+"*/*.tsv")

#data = [d.split('.')[0].split('/')[-2] for d in data]
sims = range(2,21,2)
data = ["sim_"+str(s) for s in sims]
shm = ["0_01","0_05","0_1","0_2","0_3"] #this has to match the numbers in simulate.sh

wildcard_constraints:
    d = "sim_\d+"
    s = "0_\d+"
   # d = "[a-zA-Z]+\d"

rule result:
    input:
        expand(OUTPUT + "{d}/partitions/sim_{s}_partition_0.fasta", d=data, s=shm)
        #expand(OUTPUT + "{d}/tree_files/", d=data),
        #expand(OUTPUT + "{d}/germline_search/partition_0/germline.fasta", d=data)
        
        
############################ 
#simulate from existing data
#rule tsv_to_fasta:
#    resources:
#        mem="10G",
#    threads: 2
#    log: os.path.join(DATADIR, "logs", "tsv_to_fasta_{d}.log")
#    input:
#        tsv = DATADIR+"{d}/airr-covid-19.tsv", 
#        script = 'Download_data/tsv_to_fasta.py',
#    output:
#        fasta= OUTPUT + "{d}.fasta"
#    shell:
#        "echo " + platform.node() + " &>> {log} && \
#        python {input.script} -i {input.tsv} -o {output.fasta} &&>> {log}"



#cache_parameters partis
#rule cache_parameters:
#    resources:
#        mem="20G",
#    threads: 10
#    log: os.path.join(DATADIR, "logs", "cache_parameters_{d}.log")
#    input:
#        fasta = OUTPUT + "{d}.fasta",
#        script = 'partis/partis_simulation/cache_parameters.sh',
 #       partis = PARTIS+"bin/partis"
 #   output:
 #       out = directory(OUTPUT + "{d}/")
 #   shell:
 #       "module purge &>> {log} && \
 #       module load gcc/8.3.0 &>> {log} && \
 #       module load gsl/2.5 &>> {log} && \
 #       module load git &>> {log} && \
 #       export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log} && \
 #       echo " + platform.node() + " &>> {log} && \
 #       mkdir {output.out} &>> {log} && \
 #       sh {input.script} -f {input.fasta} -p {input.partis} -o {output.out} &>> {log}"

#partition partis
#rule partition:
#    resources:
#        mem_mb="100G",
#	mem="180G",
#    threads: 100
#    log: os.path.join(DATADIR, "logs", "partition_{d}.log")
#    input:
#        fasta = OUTPUT + "{d}.fasta",
#        script = 'partis/partis_simulation/partition.sh',
#        partis = HAMPARTIS+"bin/partis",
#        out = OUTPUT + "{d}/"
#    output:
#        out= OUTPUT + "{d}/pd.yaml"
#    shell:
#        "module purge &>> {log} && \
#        module load gcc/8.3.0 &>> {log} && \
#        module load gsl/2.5 &>> {log} && \
#        module load git &>> {log} && \
#        export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log} && \
#        echo " + platform.node() + " &>> {log} && \
#        sh {input.script} -f {input.fasta} -p {input.partis} -o {input.out} &>> {log}"

############################
#simulate from scratch
#simulation partis

rule simulate:
    resources:
        mem="50G",
    threads: 10
    log: os.path.join(DATADIR, "logs", "simulate_{d}.log")
    input:
        script = 'partis/partis_simulation/simulate.sh',
        partis = PARTIS+"bin/partis"
    output:
        out_dir = directory(OUTPUT + "{d}/simulations/"),
        out= OUTPUT + "{d}/simulations/sim_{s}.yaml"
    shell:
        "module purge &>> {log} && \
        module load gcc/8.3.0 &>> {log} && \
        module load gsl/2.5 &>> {log} && \
        module load git &>> {log} && \
        export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log} && \
        echo " + platform.node() + " &>> {log} && \
        export LD_LIBRARY_PATH=/home1/kavoss/partis_with_simulation/partis/packages/bpp-newlik/_build/lib64:$LD_LIBRARY_PATH &>> {log} && \
        sh {input.script} -p {input.partis} -o {output.dir} &>> {log}"

#analyze_partis_output
rule analyze_partis_output:
     resources:
        mem="50G",
     threads: 10
     log: os.path.join(DATADIR, "logs", "analyze_partis_output_{d}.log")
     input:
        dir = OUTPUT + "{d}/simulations/",
        script = 'analyze_partis_output/simulation_analysis.sh',
        partis = PARTIS,
        sim_check = OUTPUT + "{d}/simulations/sim_{s}.yaml"
     output:
        out = directory(OUTPUT + "{d}/{s}/"),
        naive = OUTPUT + "{d}/{s}/naive.fasta",
        all = OUTPUT + "{d}/{s}/all.fasta",
        clonal_families = OUTPUT + "{d}/{s}/family_1.fasta"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log} && \
        sh {input.script} -d {input.dir} -p {input.partis} -o {output.out} &>> {log}"


rule findVDJ:
     resources:
        mem="10G",
     threads: 10
     log: os.path.join(DATADIR, "logs", "findVDJ_{d}.log")
     input:
        dir = OUTPUT + "{d}/partitions/",
        script = 'germline_search/IMGT_vrequest.sh',
        vquest = VQUEST,
        partition_check = OUTPUT + "{d}/partitions/sim_5_partition_0.fasta"
     output:
        out = directory(OUTPUT + "{d}/germline_search/"),
        dir_check = directory(OUTPUT + "{d}/germline_search/sim_5_partition_0"),
        seq = OUTPUT + "{d}/germline_search/sim_5_partition_0/3_Nt-sequences.txt"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        sh {input.script} -d {input.dir} -v {input.vquest} -o {output.out} &>> {log}"

rule find_germline:
     resources:
        mem="10G",
     threads: 10
     log: os.path.join(DATADIR, "logs", "find_germline{d}.log")
     input:
        dir = OUTPUT + "{d}/germline_search/",
        script = 'germline_search/find_germline_RevertToGermline.sh',
        RevertToGermline = RTG,
        seq_check = OUTPUT + "{d}/germline_search/sim_5_partition_0/3_Nt-sequences.txt"
     output:
        seq = OUTPUT + "{d}/germline_search/sim_5_partition_0/germline.fasta"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log} && \
        sh {input.script} -d {input.dir} -r {input.RevertToGermline} -o {input.dir} &>> {log}"

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
        "echo " + platform.node() + " &>> {log} && \
        sh {input.script} -d {input.partitions} -o {output.out} &>> {log}"


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
        "echo " + platform.node() + " &>> {log} && \
        sh {input.script} -r {input.raxml} -d {input.partitions_aligned} -o {output.out} &>> {log}"









