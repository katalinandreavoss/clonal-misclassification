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
MPTP=config['MPTP']


clones = range(4,21,2)
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
        expand(OUTPUT + "{d}/{s}/{l}/{i}/family_sizes.txt", d=clones, s=shm,l=leaves, i=sims),
        expand(OUTPUT+ "{d}/{s}/{l}/{i}/mptp_data.txt", d=clones, s=shm,l=leaves, i=sims),
        expand(OUTPUT+ "{d}/{s}/{l}/{i}/mptp_data_singletons.txt", d=clones, s=shm,l=leaves, i=sims),
        expand(OUTPUT + "{d}/{s}/{l}/{i}/clean.fasta.vdjca.clns_IGH.tsv", d=clones, s=shm,l=leaves, i=sims)
       
   

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


rule remove_singletons:
     resources:
        mem="5G",
     threads: 10
     log: os.path.join(OUTPUT, "logs", "remove_singletons_{d}_{s}_{l}_{i}.log")
     input:
        script = 'tree_building/remove_singletons.R',
        all = OUTPUT + "{d}/{s}/{l}/{i}/all.fasta"
     params:
        out = OUTPUT + "{d}/{s}/{l}/{i}/"
     output:
        clean = OUTPUT + "{d}/{s}/{l}/{i}/clean.fasta",
     shell:
        "echo " + platform.node() + " &>> {log} && \
        Rscript {input.script} -f {input.all} -o {params.out} &>> {log}"

#align
rule align:
     resources:
        mem="5G",
     threads: 10
     log: os.path.join(OUTPUT, "logs", "align_{d}_{s}_{l}_{i}.log")
     input:
        all = OUTPUT + "{d}/{s}/{l}/{i}/clean.fasta",
        script = 'tree_building/align_partitions.sh'
     params:
         out = OUTPUT + "{d}/{s}/{l}/{i}/"
     output:
        align = OUTPUT + "{d}/{s}/{l}/{i}/clean_aligned.fasta",
     shell:
        "echo " + platform.node() + " &>> {log} && \
        sh {input.script} -d {params.out} &>> {log}"

#build megatree
rule build_tree:
     resources:
        mem="10G",
     threads: 16
     log: os.path.join(OUTPUT, "logs", "build_tree_{d}_{s}_{l}_{i}.log")
     input:
        script = 'tree_building/build_tree.sh',
        raxml  = RAXML+"raxml-ng",
        align_check = OUTPUT + "{d}/{s}/{l}/{i}/clean_aligned.fasta"
     params:
         out = OUTPUT + "{d}/{s}/{l}/{i}/"
     output:
        out = directory(OUTPUT + "{d}/{s}/{l}/{i}/tree_files/"),
        tree = OUTPUT + "{d}/{s}/{l}/{i}/tree_files/mega_tree_.raxml.bestTree"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        sh {input.script} -r {input.raxml} -d {params.out} -o {output.out} &>> {log}"

#this rule does not run on cluster because it needs X11 forwarding: do ssh with -X flag and then run snakemake without running it on the cluster
rule cut_tree:
     resources:
        mem="200G",
     threads: 10
     log: os.path.join(OUTPUT, "logs", "cut_tree_{d}_{s}_{l}_{i}.log")
     input:
        ptp  = PTP,
        tree = OUTPUT + "{d}/{s}/{l}/{i}/tree_files/mega_tree_.raxml.bestTree"
     params:
         out = OUTPUT + "{d}/{s}/{l}/{i}/mega"
     output:
        tree = OUTPUT + "{d}/{s}/{l}/{i}/mega.PTPhSupportPartition.txt.sh.tre",
        summary = OUTPUT + "{d}/{s}/{l}/{i}/mega.PTPPartitonSummary.txt",
        partitions = OUTPUT+ "{d}/{s}/{l}/{i}/mega.PTPPartitions.txt",
        png = OUTPUT+ "{d}/{s}/{l}/{i}/mega.PTPhSupportPartition.txt.png"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log} && \
        python {input.ptp} -t {input.tree} -o {params.out} &>> {log}"

#this rule does not run on cluster because it needs X11 forwarding: do ssh with -X flag and then run snakemake without running it on the cluster
rule cut_tree_mptp:
     resources:
        mem="200G",
     threads: 10
     log: os.path.join(OUTPUT, "logs", "cut_tree_mptp_{d}_{s}_{l}_{i}.log")
     input:
        mptp  = MPTP+"mptp",
        tree = OUTPUT + "{d}/{s}/{l}/{i}/tree_files/mega_tree_.raxml.bestTree"
     params:
         out = OUTPUT + "{d}/{s}/{l}/{i}/mega_mptp"
     output:
        partitions = OUTPUT+ "{d}/{s}/{l}/{i}/mega_mptp.txt",
        svg = OUTPUT+ "{d}/{s}/{l}/{i}/mega_mptp.svg"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log} && \
        {input.mptp} --ml --single --tree_file {input.tree} --output_file {params.out} &>> {log}"


rule get_mptp_values:
     resources:
        mem="2G",
     threads: 10
     log: os.path.join(OUTPUT, "logs", "get_mptp_values_{d}_{s}_{l}_{i}.log")
     input:
        script = 'simulation_analyses/analyse_ptp_output.py',
        partitions = OUTPUT+ "{d}/{s}/{l}/{i}/mega_mptp.txt"
     params:
        out = OUTPUT + "{d}/{s}/{l}/{i}/"
     output:
        out = OUTPUT+ "{d}/{s}/{l}/{i}/mptp_data.txt"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log} && \
        python {input.script} {params.out} &>> {log}"


rule get_mptp_values_singletons:
     resources:
        mem="2G",
     threads: 10
     log: os.path.join(OUTPUT, "logs", "get_mptp_values_singletons_{d}_{s}_{l}_{i}.log")
     input:
        script = 'simulation_analyses/analyse_ptp_output_with_singletons.py',
        partitions = OUTPUT+ "{d}/{s}/{l}/{i}/mega_mptp.txt"
     params:
        out = OUTPUT + "{d}/{s}/{l}/{i}/"
     output:
        out = OUTPUT+ "{d}/{s}/{l}/{i}/mptp_data_singletons.txt"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log} && \
        python {input.script} {params.out} &>> {log}"

rule get_family_sizes:
     resources:
        mem="1G",
     threads: 10
     log: os.path.join(OUTPUT, "logs", "plot_family_sizes_{d}_{s}_{l}_{i}.log")
     input:
        script = 'simulation_analyses/get_real_family_sizes.sh',
        sequences = OUTPUT + "{d}/{s}/{l}/{i}/family_1.fasta"
     params:
         out = OUTPUT + "{d}/{s}/{l}/{i}/",
         clones = "{d}",
         shm = "{s}",
         leaves = "{l}",
         sim = "{i}"
     output:
        family_sizes = OUTPUT + "{d}/{s}/{l}/{i}/family_sizes.txt"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log}&& \
        sh {input.script} -d {params.out} -c {params.clones} -s {params.shm} -l {params.leaves} -i {params.sim} &>> {log}"


rule  mixcr:
     resources:
        mem="20G",
     threads: 10
     log: os.path.join(OUTPUT, "logs", "mixcr_{d}_{s}_{l}_{i}.log")
     input:
        fasta = OUTPUT + "{d}/{s}/{l}/{i}/clean.fasta"
     params:
         out = OUTPUT + "{d}/{s}/{l}/{i}/clean.fasta.vdjca.clns.tsv"
     output:
        aligned = OUTPUT + "{d}/{s}/{l}/{i}/clean.fasta.vdjca",
        clones = OUTPUT + "{d}/{s}/{l}/{i}/clean.fasta.vdjca.clns",
        IGH = OUTPUT + "{d}/{s}/{l}/{i}/clean.fasta.vdjca.clns_IGH.tsv"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log} && \
        mixcr align --preset rnaseq-bcr-full-length --species HomoSapiens {input.fasta} {output.aligned} &>> {log} && \
        mixcr assemble {output.aligned} {output.clones} &>> {log} && \
        mixcr exportClones {output.clones} {params.out} &>> {log}"



