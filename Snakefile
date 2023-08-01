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
PTP=config['PTP']


clones = range(2,21,2)
shm = ["0_01","0_05","0_1","0_2","0_3"] #this has to match the numbers in simulate.sh
leaves = ["10","20","50","100"]
sims = range(1,11)


wildcard_constraints:
    d = "\d+",
    s = "0_\d+",
    l = "\d+",
    i = "\d+"

rule all:
    input:
        expand(OUTPUT + "{d}/{s}/{l}/{i}/all.fasta", d=clones, s=shm,l=leaves, i=sims)
   

#simulation partis
rule simulate:
    resources:
        mem="20G",
    threads: 10
    log: os.path.join(OUTPUT, "logs", "simulate_{d}_{s}_{l}_{i}.log")
    input:
        script = 'partis/partis_simulation/simulate.sh',
        partis = PARTIS+"bin/partis"
    params:
        out_dir = OUTPUT + "{d}/{s}/{l}/{i}",
        clones = "{d}",
        shm = "{s}",
        leaves = "{l}",
        sim = "{i}"
    output:
        out = OUTPUT + "{d}/{s}/{l}/{i}/clones_{d}_shm_{s}_leaves_{l}_sim_{i}.yaml"
    shell:
        "module purge &>> {log} && \
        module load gcc/8.3.0 &>> {log} && \
        module load gsl/2.5 &>> {log} && \
        module load git &>> {log} && \
        export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log} && \
        echo " + platform.node() + " &>> {log} && \
        export LD_LIBRARY_PATH=/home1/kavoss/partis_with_simulation/partis/packages/bpp-newlik/_build/lib64:$LD_LIBRARY_PATH &>> {log} && \
        sh {input.script} -p {input.partis} -o {params.out_dir} -c {params.clones} -s {params.shm} -l {params.leaves} -i {params.sim} &>> {log}"



#analyze_partis_output
rule analyze_partis_output:
     resources:
        mem="10G",
     threads: 10
     log: os.path.join(OUTPUT, "logs", "analyze_partis_output_{d}_{s}_{l}_{i}.log")
     input:
        script = 'analyze_partis_output/yaml_to_families_new.py',
        partis = PARTIS,
        simulation = OUTPUT + "{d}/{s}/{l}/{i}/clones_{d}_shm_{s}_leaves_{l}_sim_{i}.yaml"
     params:
        out = OUTPUT + "{d}/{s}/{l}/{i}/"
     output:
        naive = OUTPUT + "{d}/{s}/{l}/{i}/naive.fasta",
        all = OUTPUT + "{d}/{s}/{l}/{i}/all.fasta",
        clonal_families = OUTPUT + "{d}/{s}/{l}/{i}/family_1.fasta"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log} && \
        python {input.script} {input.simulation} {params.out} &>> {log}"


rule align_partitions:
     resources:
        mem="50G",
     threads: 10
     log: os.path.join(OUTPUT, "logs", "align_partitions_{d}_{s}.log")
     input:
        all = OUTPUT + "{d}/{s}/all.fasta",
        clonal_families = OUTPUT + "{d}/{s}/family_1.fasta",
        script = 'tree_building/align_partitions.sh'
     params:
         out = OUTPUT + "{d}/{s}/"
     output:
        align = OUTPUT + "{d}/{s}/all_aligned.fasta",
        align_fam = OUTPUT + "{d}/{s}/family_1_aligned.fasta"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        sh {input.script} -d {params.out} &>> {log}"


rule build_tree:
     resources:
        mem="500G",
     threads: 100
     log: os.path.join(OUTPUT, "logs", "build_tree_{d}_{s}.log")
     input:
        script = 'tree_building/build_tree.sh',
        raxml  = RAXML+"raxml-ng",
        align_check = OUTPUT + "{d}/{s}/all_aligned.fasta"
     params:
         out = OUTPUT + "{d}/{s}/"
     output:
        out = directory(OUTPUT + "{d}/{s}/tree_files/"),
        tree = OUTPUT + "{d}/{s}/tree_files/all_tree_.raxml.bestTree"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        sh {input.script} -r {input.raxml} -d {params.out} -o {output.out} &>> {log}"

#this rule does not run on cluster because it needs X11 forwarding: do ssh with -X flag and then run snakemake without running it on the cluster
rule cut_tree:
     resources:
        mem="500G",
     threads: 10
     log: os.path.join(OUTPUT, "logs", "cut_tree_{d}_{s}.log")
     input:
        ptp  = PTP,
        tree = OUTPUT + "{d}/{s}/tree_files/all_tree_.raxml.bestTree"
     params:
         out = OUTPUT + "{d}/{s}/all_PTP"
     output:
        tree = OUTPUT + "{d}/{s}/all_PTP.PTPhSupportPartition.txt.sh.tre",
        summary = OUTPUT + "{d}/{s}/all_PTP.PTPPartitonSummary.txt",
        partitions = OUTPUT+ "{d}/{s}/all_PTP.PTPPartitions.txt",
        png = OUTPUT+ "{d}/{s}/all_PTP.PTPhSupportPartition.txt.png"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log} && \
        python {input.ptp} -t {input.tree} -o {params.out} &>> {log}"



rule get_family_sizes:
     resources:
        mem="150G",
     threads: 10
     log: os.path.join(OUTPUT, "logs", "plot_family_sizes_{d}_{s}.log")
     input:
        script = 'simulation_analyses/get_real_family_sizes.sh',
        sequences = OUTPUT + "{d}/{s}/family_1.fasta"
     params:
         out = OUTPUT + "{d}/{s}"
     output:
        family_sizes = OUTPUT + "{d}/{s}/family_sizes.txt"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log}&& \
        sh {input.script} -d {params.out} &>> {log}"

rule compare_ancestral_seq:
     resources:
        mem="50G",
     threads: 10
     log: os.path.join(OUTPUT, "logs", "compare_ancestral_seq_{d}_{s}.log")
     input:
        script = 'germline_search/combine_naive_discerned.R',
        script_align = 'germline_search/align_combined.sh',
        naive = OUTPUT + "{d}/{s}/naive.fasta",
        discerned = OUTPUT + "{d}/{s}/ancestral_sequences/"
     params:
         out = OUTPUT + "{d}/{s}"
     output:
         out = directory(OUTPUT + "{d}/{s}/combined/")
     shell:
        "echo " + platform.node() + " &>> {log} && \
        export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log}&& \
        mkdir -p {output.out} &>> {log}&& \
        Rscript {input.script} -p {params.out} -n {input.naive} &>> {log}&& \
        sh {input.script_align} -d {output.out} &>> {log}"



