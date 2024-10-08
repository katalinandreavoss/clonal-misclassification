import glob
import yaml
import platform
import os

configfile: "configs/config.yaml"


OUTPUT=config['output']
PARTIS=config['partis']
RAXML=config['raxml']
VQUEST=config['vquest']
MPTP=config['MPTP']
VDJ=config['VDJ']

#clones = range(10,11,5)
clones = range(16,21,10)
#clones = [10,20,50]
#clones = ["cattle"]
shm = ["0_001","0_005", "0_01","0_05","0_1","0_2"] 
#shm = ["0_005", "0_01","0_1","0_2"] 
#shm = ["0_05","0_1","0_2","0_3"] 
#shm = ["0_0005", "0_00075", "0_001","0_0015","0_002","0_0025","0_003","0_004","0_005","0_01","0_05","0_1","0_2"] 
#shm = ["0_01"]
#shm = [14363]
leaves = ["10","20","50","100"]
#leaves = ["20","50"]
#leaves = [1]
balance = ["0_0"]
#balance = [1]
#junction_length = ["6","12","18","24","30"]
#junction_length = ["20","30","40"]
sims = range(1,51)
#sims = range(1,11)
#sims = ["day0_z_subset","day28_z_subset"]



wildcard_constraints:
    d = "\d+",
    s = "0_\d+",
    l = "\d+",
    #b = "\d+",
    b = "0_\d+",
    i = "\d+"

# wildcard_constraints:
#     d = "\w+",
#     s = "\d+",
#     l = "\d+",
#     #b = "\d+",
#     b = "\d+",
   # i = "\w+\d"


rule all:
    input:
        expand(OUTPUT + "{d}/{s}/{l}/{b}/{i}/sensitivity_precision.tsv", d=clones, s=shm,l=leaves, b=balance, i=sims)

#simulation partis
# rule simulate:
#     resources:
#         mem="20G",
#     threads: 10
#     log: os.path.join(OUTPUT, "logs", "simulate_{d}_{s}_{l}_{b}_{i}.log")
#     input:
#         script = 'partis/partis_simulation/simulate.sh',
#         partis = PARTIS+"bin/partis"
#     params:
#         out_dir = OUTPUT + "{d}/{s}/{l}/{b}/{i}",
#         clones = "{d}",
#         shm = "{s}",
#         leaves = "{l}",
#         balance = "{b}",
#         sim = "{i}"
#     output:
#         out =  OUTPUT + "{d}/{s}/{l}/{b}/{i}/simu.yaml"
#     shell:
#         "module purge &>> {log} && \
#         module load gcc/8.3.0 &>> {log} && \
#         module load gsl/2.5 &>> {log} && \
#         module load git &>> {log} && \
#         export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log} && \
#         echo " + platform.node() + " &>> {log} && \
#         export LD_LIBRARY_PATH=/home1/kavoss/partis_with_simulation/partis/packages/bpp-newlik/_build/lib64:$LD_LIBRARY_PATH &>> {log} && \
#         sh {input.script} -p {input.partis} -o {params.out_dir} -c {params.clones} -s {params.shm} -l {params.leaves} -b {params.balance} -i {params.sim} &>> {log}"

# rule simulate_restrict_genes:
#     resources:
#         mem="20G",
#     threads: 10
#     log: os.path.join(OUTPUT, "logs", "simulate_restrict_genes_{d}_{s}_{l}_{b}_{i}.log")
#     input:
#         script = 'partis/partis_simulation/simulate_restrict_genes.sh',
#         partis = PARTIS+"bin/partis"
#     params:
#         out_dir = OUTPUT + "{d}/{s}/{l}/{b}/{i}",
#         clones = "{d}",
#         shm = "{s}",
#         leaves = "{l}",
#         balance = "{b}",
#         sim = "{i}"
#     output:
#         out =  OUTPUT + "{d}/{s}/{l}/{b}/{i}/simu.yaml"
#     shell:
#         "module purge &>> {log} && \
#         module load gcc/8.3.0 &>> {log} && \
#         module load gsl/2.5 &>> {log} && \
#         module load git &>> {log} && \
#         export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log} && \
#         echo " + platform.node() + " &>> {log} && \
#         export LD_LIBRARY_PATH=/home1/kavoss/partis_with_simulation/partis/packages/bpp-newlik/_build/lib64:$LD_LIBRARY_PATH &>> {log} && \
#         sh {input.script} -p {input.partis} -o {params.out_dir} -c {params.clones} -s {params.shm} -l {params.leaves} -b {params.balance} -i {params.sim} &>> {log}"


#simulation new
# rule simulate_own:
#     resources:
#         mem="2G",
#     threads: 10
#     log: os.path.join(OUTPUT, "logs", "simulate_{d}_{s}_{l}_{b}_{i}.log")
#     input:
#         script = 'simulation/simulate.R',
#         vdj_dir = VDJ
#     params:
#         out_dir = OUTPUT + "{d}/{s}/{l}/{b}/{i}",
#         clones = "{d}",
#         shm = "{s}",
#         leaves = "{l}",
#         junction = "{b}",
#     output:
#         out =  OUTPUT + "{d}/{s}/{l}/{b}/{i}/clean.fasta",
#         fam = OUTPUT + "{d}/{s}/{l}/{b}/{i}/family_1.fasta"
#     shell:
#         "echo " + platform.node() + " &>> {log} && \
#         Rscript {input.script} -d {params.out_dir} -r {params.shm} -p {input.vdj_dir} -l {params.leaves} -j {params.junction} &>> {log}"


# rule simulate_fake_V:
#     resources:
#         mem="2G",
#     threads: 10
#     log: os.path.join(OUTPUT, "logs", "simulate_{d}_{s}_{l}_{b}_{i}.log")
#     input:
#         script = 'simulation/simulate_new_V_genes.R',
#         vdj_dir = VDJ
#     params:
#         out_dir = OUTPUT + "{d}/{s}/{l}/{b}/{i}",
#         clones = "{d}",
#         shm = "{s}",
#         leaves = "{l}",
#         junction = "{b}",
#     output:
#         out =  OUTPUT + "{d}/{s}/{l}/{b}/{i}/clean.fasta",
#         fam = OUTPUT + "{d}/{s}/{l}/{b}/{i}/family_1.fasta"
#     shell:
#         "echo " + platform.node() + " &>> {log} && \
#         Rscript {input.script} -d {params.out_dir} -r {params.shm} -p {input.vdj_dir} -l {params.leaves} -j {params.junction} &>> {log}"


#analyze_partis_output
rule analyze_partis_output:
     resources:
        mem="2G",
     threads: 10
     log: os.path.join(OUTPUT, "logs", "analyze_partis_output_{d}_{s}_{l}_{b}_{i}.log")
     input:
        script = 'analyze_partis_output/yaml_to_families_new.py',
        partis = PARTIS,
        simulation = OUTPUT + "{d}/{s}/{l}/{b}/{i}/simu.yaml" 
     params:
        out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/"
     output:
        naive = OUTPUT + "{d}/{s}/{l}/{b}/{i}/naive.fasta",
        all = OUTPUT + "{d}/{s}/{l}/{b}/{i}/all.fasta",
        clonal_families = OUTPUT + "{d}/{s}/{l}/{b}/{i}/family_1.fasta"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log} && \
        python {input.script} {input.simulation} {params.out} &>> {log}"


rule remove_singletons:
     resources:
        mem="2G",
     threads: 10
     log: os.path.join(OUTPUT, "logs", "remove_singletons_{d}_{s}_{l}_{b}_{i}.log")
     input:
        script = 'tree_building/remove_singletons.R',
        all = OUTPUT + "{d}/{s}/{l}/{b}/{i}/all.fasta"
     params:
        out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/"
     output:
        clean = OUTPUT + "{d}/{s}/{l}/{b}/{i}/clean.fasta"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        Rscript {input.script} -f {input.all} -o {params.out} &>> {log}"


rule align:
     resources:
        mem="100G",
     threads: 100
     log: os.path.join(OUTPUT, "logs", "align_{d}_{s}_{l}_{b}_{i}.log")
     input:
        all = OUTPUT + "{d}/{s}/{l}/{b}/{i}/clean.fasta",
        script = 'tree_building/align_partitions.sh'
     params:
         out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/",
         filename = OUTPUT + "{d}/{s}/{l}/{b}/{i}/clean"
     output:
        align = OUTPUT + "{d}/{s}/{l}/{b}/{i}/clean_aligned.fasta"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        sh {input.script} -d {params.out} -f {params.filename} &>> {log}"

# rule align_families:
#      resources:
#         mem="2G",
#      threads: 10
#      log: os.path.join(OUTPUT, "logs", "align_families_{d}_{s}_{l}_{b}_{i}.log")
#      input:
#         all = OUTPUT + "{d}/{s}/{l}/{b}/{i}/family_1.fasta",
#         script = 'tree_building/align_families.sh'
#      params:
#          out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/"
#      output:
#         align = OUTPUT + "{d}/{s}/{l}/{b}/{i}/family_1_aligned.fasta"
#      shell:
#         "echo " + platform.node() + " &>> {log} && \
#         sh {input.script} -d {params.out} &>> {log}"


# #important for this rule: you must specify at least --cpus-per-task=16 in the snakemake command, otherwise it takes forever
rule build_megatree:
     resources:
        mem="20G"
     log: os.path.join(OUTPUT, "logs", "build_megatree_{d}_{s}_{l}_{b}_{i}.log")
     input:
        script = 'tree_building/build_tree.sh',
        raxml  = RAXML+"raxml-ng",
        align_check = OUTPUT + "{d}/{s}/{l}/{b}/{i}/clean_aligned.fasta"
     params:
         dir = OUTPUT + "{d}/{s}/{l}/{b}/{i}/",
         out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/tree_files/"
     output:
        tree = OUTPUT + "{d}/{s}/{l}/{b}/{i}/tree_files/mega_tree_.raxml.bestTree"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        sh {input.script} -r {input.raxml} -f {input.align_check} -o {params.out} -f &>> {log}"

rule build_upgma_tree:
     resources:
        mem="20G"
     log: os.path.join(OUTPUT, "logs", "build_upgma_tree_{d}_{s}_{l}_{b}_{i}.log")
     input:
        script = 'tree_building/distance_matrix_UPGMA.R',
        fasta = OUTPUT + "{d}/{s}/{l}/{b}/{i}/clean.fasta"
     params:
         out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/tree_files/"
     output:
        tree = OUTPUT + "{d}/{s}/{l}/{b}/{i}/tree_files/upgma.tre"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        Rscript {input.script} -f {input.fasta} -o {params.out} &>> {log}"


# # build subtrees
# rule build_tree:
#      resources:
#         mem="5G",
#      threads: 32
#      log: os.path.join(OUTPUT, "logs", "build_sub_trees_{d}_{s}_{l}_{b}_{i}.log")
#      input:
#         script = 'tree_building/build_sub_trees.sh',
#         raxml  = RAXML+"raxml-ng",
#         align_check = OUTPUT + "{d}/{s}/{l}/{b}/{i}/family_1_aligned.fasta"
#      params:
#          out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/"
#      output:
#         out = directory(OUTPUT + "{d}/{s}/{l}/{b}/{i}/tree_files/"),
#         tree = OUTPUT + "{d}/{s}/{l}/{b}/{i}/tree_files/build_tree.txt"
#      shell:
#         "echo " + platform.node() + " &>> {log} && \
#         sh {input.script} -r {input.raxml} -d {params.out} -o {output.out} &>> {log}"





# # #this rule does not run on cluster because it needs X11 forwarding: do ssh with -X flag and then run snakemake without running it on the cluster
rule cut_tree_mptp:
     resources:
        mem="200G",
     threads: 10
     log: os.path.join(OUTPUT, "logs", "cut_tree_mptp_{d}/{s}/{l}/{b}/{i}/.log")
     input:
        mptp  = MPTP+"mptp",
        tree = OUTPUT + "{d}/{s}/{l}/{b}/{i}/tree_files/mega_tree_.raxml.bestTree"
     params:
         out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mega_mptp"
     output:
        partitions = OUTPUT+ "{d}/{s}/{l}/{b}/{i}/mega_mptp.txt",
        svg = OUTPUT+ "{d}/{s}/{l}/{b}/{i}/mega_mptp.svg"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        {input.mptp} --ml --single --tree_file {input.tree} --output_file {params.out} &>> {log}"

rule cut_upgma_tree_mptp:
     resources:
        mem="200G",
     threads: 10
     log: os.path.join(OUTPUT, "logs", "cut_upgma_tree_mptp_{d}/{s}/{l}/{b}/{i}/.log")
     input:
        mptp  = MPTP+"mptp",
        tree = OUTPUT + "{d}/{s}/{l}/{b}/{i}/tree_files/upgma.tre"
     params:
         out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mega_upgma_mptp"
     output:
        partitions = OUTPUT+ "{d}/{s}/{l}/{b}/{i}/mega_upgma_mptp.txt",
        svg = OUTPUT+ "{d}/{s}/{l}/{b}/{i}/mega_upgma_mptp.svg"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        {input.mptp} --ml --single --tree_file {input.tree} --output_file {params.out} &>> {log}"


rule get_mptp_values:
     resources:
        mem="2G",
     threads: 10
     log: os.path.join(OUTPUT, "logs", "get_mptp_values_{d}_{s}_{l}_{b}_{i}.log")
     input:
        script = 'simulation_analyses/analyse_ptp_output.py',
        partitions = OUTPUT+ "{d}/{s}/{l}/{b}/{i}/mega_mptp.txt"
     params:
        out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/"
     output:
        out = OUTPUT+ "{d}/{s}/{l}/{b}/{i}/mptp_data.txt"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log} && \
        python {input.script} {input.partitions} {output.out}&>> {log}"

rule get_upgma_mptp_values:
     resources:
        mem="2G",
     threads: 10
     log: os.path.join(OUTPUT, "logs", "get_upgma_mptp_values_{d}_{s}_{l}_{b}_{i}.log")
     input:
        script = 'simulation_analyses/analyse_ptp_output.py',
        partitions = OUTPUT+ "{d}/{s}/{l}/{b}/{i}/mega_upgma_mptp.txt"
     params:
        out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/"
     output:
        out = OUTPUT+ "{d}/{s}/{l}/{b}/{i}/mega_upgma_mptp_data.txt"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log} && \
        python {input.script} {input.partitions} {output.out}&>> {log}"


rule get_mptp_values_singletons:
     resources:
        mem="2G",
     threads: 10
     log: os.path.join(OUTPUT, "logs", "get_mptp_values_{d}_{s}_{l}_{b}_{i}.log")
     input:
        script = 'simulation_analyses/analyse_ptp_output_with_singletons.py',
        partitions = OUTPUT+ "{d}/{s}/{l}/{b}/{i}/mega_mptp.txt"
     params:
        out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/"
     output:
        out = OUTPUT+ "{d}/{s}/{l}/{b}/{i}/mptp_data_singletons.txt"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log} && \
        python {input.script} {input.partitions} {output.out}&>> {log}"


rule get_upgma_mptp_values_singletons:
     resources:
        mem="2G",
     threads: 10
     log: os.path.join(OUTPUT, "logs", "get_upgma_mptp_values_singletons_{d}_{s}_{l}_{b}_{i}.log")
     input:
        script = 'simulation_analyses/analyse_ptp_output_with_singletons.py',
        partitions = OUTPUT+ "{d}/{s}/{l}/{b}/{i}/mega_upgma_mptp.txt"
     params:
        out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/"
     output:
        out = OUTPUT+ "{d}/{s}/{l}/{b}/{i}/mega_upgma_mptp_data_singletons.txt"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log} && \
        python {input.script} {input.partitions} {output.out}&>> {log}"

# rule gmyc:
#      resources:
#         mem="5G",
#      threads: 10
#      log: os.path.join(OUTPUT, "logs", "gmyc_{d}_{s}_{l}_{b}_{i}.log")
#      input:
#         script = 'GMYC/gmyc_mega.R',
#         tree = OUTPUT + "{d}/{s}/{l}/{b}/{i}/tree_files/mega_tree_.raxml.bestTree"
#      params:
#         out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/"
#      output:
#         out = OUTPUT+ "{d}/{s}/{l}/{b}/{i}/gmyc.tsv"
#      shell:
#         "echo " + platform.node() + " &>> {log} && \
#         Rscript {input.script} -d {params.out} &>> {log}"

rule get_family_sizes:
     resources:
        mem="1G",
     threads: 10
     log: os.path.join(OUTPUT, "logs", "plot_family_sizes_{d}_{s}_{l}_{b}_{i}.log")
     input:
        script = 'simulation_analyses/get_real_family_sizes.sh',
        sequences = OUTPUT + "{d}/{s}/{l}/{b}/{i}/family_1.fasta"
     params:
         out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/",
         clones = "{d}",
         shm = "{s}",
         leaves = "{l}",
         sim = "{i}"
     output:
        family_sizes = OUTPUT + "{d}/{s}/{l}/{b}/{i}/family_sizes.txt"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log}&& \
        sh {input.script} -d {params.out} -c {params.clones} -s {params.shm} -l {params.leaves} -i {params.sim} &>> {log}"


# rule mixcr:
#      resources:
#         mem="20G",
#      threads: 10
#      log: os.path.join(OUTPUT, "logs", "mixcr_{d}_{s}_{l}_{b}_{i}.log")
#      input:
#         fasta = OUTPUT + "{d}/{s}/{l}/{b}/{i}/clean.fasta"
#      params:
#          out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/clean.fasta.vdjca.clns.tsv"
#      output:
#         aligned = OUTPUT + "{d}/{s}/{l}/{b}/{i}/clean.fasta.vdjca",
#         clones = OUTPUT + "{d}/{s}/{l}/{b}/{i}/clean.fasta.vdjca.clns",
#         clna = OUTPUT + "{d}/{s}/{l}/{b}/{i}/clean.fasta.vdjca.clna",
#         IGH = OUTPUT + "{d}/{s}/{l}/{b}/{i}/clean.fasta.vdjca.clns_IGH.tsv",
#         mixcr = directory(OUTPUT + "{d}/{s}/{l}/{b}/{i}/mixcr_fastas/")
#      shell:
#         "echo " + platform.node() + " &>> {log} && \
#         mixcr align --preset rnaseq-bcr-full-length --species HomoSapiens -OsaveOriginalReads=true {input.fasta} {output.aligned} >> {log} 2>&1 && \
#         mixcr assemble {output.aligned} {output.clones} &>> {log} && \
#         mixcr assemble {output.aligned} {output.clna} &>> {log} && \
#         mixcr exportClones {output.clones} {params.out} &>> {log} && \
#         mkdir {output.mixcr} && \
#         var=$(cut -f1 {output.IGH} | tail -n +2 | paste -s -d' ') && \
#         mixcr exportReadsForClones --id $var -s {output.clna} {output.mixcr}/clean.fastq.gz &>> {log}"


rule findVDJ:
     resources:
        mem="10G",
     threads: 10
     log: os.path.join(OUTPUT, "logs", "findVDJ_{d}_{s}_{l}_{b}_{i}.log")
     input:
        script = 'germline_search/IMGT_vrequest.sh',
        vquest = VQUEST,
        fasta = OUTPUT + "{d}/{s}/{l}/{b}/{i}/clean.fasta"
     params:
         out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/"
     output:
        out = directory(OUTPUT + "{d}/{s}/{l}/{b}/{i}/vquest_files/"),
        dir_check = directory(OUTPUT + "{d}/{s}/{l}/{b}/{i}/vquest_files/001"),
        seq = OUTPUT + "{d}/{s}/{l}/{b}/{i}/vquest_files/001/3_Nt-sequences.txt"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        sh {input.script} -f {input.fasta} -v {input.vquest} -o {output.out} &>> {log}"

# rule combine_vquest:
#      resources:
#         mem="10G",
#      threads: 10
#      log: os.path.join(OUTPUT, "logs", "combine_vquest_{d}_{s}_{l}_{b}_{i}.log")
#      input:
#         script = 'germline_search/combine_vquest.py',
#         seq = OUTPUT + "{d}/{s}/{l}/{b}/{i}/vquest_files/001/3_Nt-sequences.txt",
#         fasta = OUTPUT + "{d}/{s}/{l}/{b}/{i}/clean.fasta",
#         dir = OUTPUT + "{d}/{s}/{l}/{b}/{i}/vquest_files/"
#      params:
#          out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/"
#      output:
#         out = directory(OUTPUT + "{d}/{s}/{l}/{b}/{i}/vquest_files/combined/"),
#         Nt = OUTPUT + "{d}/{s}/{l}/{b}/{i}/vquest_files/combined/3_Nt-sequences.txt",
#         summary = OUTPUT + "{d}/{s}/{l}/{b}/{i}/vquest_files/combined/1_Summary.txt"
#      shell:
#         "echo " + platform.node() + " &>> {log} && \
#         mkdir -p {output.out} &>> {log} && \
#         python {input.script} {input.dir} &>> {log}"

# rule changeo:
#      resources:
#         mem="2G",
#      threads: 10
#      log: os.path.join(OUTPUT, "logs", "changeo_{d}_{s}_{l}_{b}_{i}.log")
#      input:
#         script = 'simulation_analyses/change-o.sh',
#         summary = OUTPUT + "{d}/{s}/{l}/{b}/{i}/vquest_files/combined/1_Summary.txt",
#         dir = OUTPUT + "{d}/{s}/{l}/{b}/{i}/vquest_files/combined/",
#         fasta = OUTPUT + "{d}/{s}/{l}/{b}/{i}/clean.fasta",
#         vdj = VDJ
#      output:
#         db = OUTPUT + "{d}/{s}/{l}/{b}/{i}/vquest_files/combined_db-pass.tsv",
#         clones = OUTPUT + "{d}/{s}/{l}/{b}/{i}/vquest_files/combined_db-pass_clone-pass.tsv",
#         germline = OUTPUT + "{d}/{s}/{l}/{b}/{i}/vquest_files/combined_db-pass_germ-pass.tsv"
#      shell:
#         "echo " + platform.node() + " &>> {log} && \
#         export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log} && \
#         sh {input.script} -d {input.dir} -f {input.fasta} -v {input.vdj} &>> {log}"

# rule scoper:
#      resources:
#         mem="2G",
#      threads: 10
#      log: os.path.join(OUTPUT, "logs", "scoper_{d}_{s}_{l}_{b}_{i}.log")
#      input:
#         script = 'simulation_analyses/scoper.R',
#         db = OUTPUT + "{d}/{s}/{l}/{b}/{i}/vquest_files/combined_db-pass_clone-pass.tsv"
#      params:
#         dir = OUTPUT + "{d}/{s}/{l}/{b}/{i}/"
#      output:
#         hierarchical = OUTPUT + "{d}/{s}/{l}/{b}/{i}/results_hierClones.tsv",
#         spectral = OUTPUT + "{d}/{s}/{l}/{b}/{i}/results_specClones.tsv"
#      shell:
#         "echo " + platform.node() + " &>> {log} && \
#         Rscript {input.script} -d {input.db} -o {params.dir}&>> {log}"

# # rule scoper_spec:
# #      resources:
# #         mem="10G",
# #      threads: 10
# #      log: os.path.join(OUTPUT, "logs", "scoper_spec_{d}_{s}_{l}_{b}_{i}.log")
# #      input:
# #         script = 'simulation_analyses/scoper_spec.R',
# #         db = OUTPUT + "{d}/{s}/{l}/{b}/{i}/vquest_files/combined_db-pass_germ-pass.tsv"
# #      params:
# #         dir = OUTPUT + "{d}/{s}/{l}/{b}/{i}/"
# #      output:
# #         spectral = OUTPUT + "{d}/{s}/{l}/{b}/{i}/results_specClones_vj.tsv"
# #      shell:
# #         "echo " + platform.node() + " &>> {log} && \
# #         Rscript {input.script} -d {input.db} -o {params.dir}&>> {log}"



rule combine_analysis:
    resources:
       mem="2G",
    threads: 10
    log: os.path.join(OUTPUT, "logs", "combine_analysis_{d}_{s}_{l}_{b}_{i}.log")
    input:
       script = 'simulation_analyses/calculate_measures.py',
       spectral = OUTPUT + "{d}/{s}/{l}/{b}/{i}/results_specClones.tsv",
       hierarchical = OUTPUT + "{d}/{s}/{l}/{b}/{i}/results_hierClones.tsv",
       real = OUTPUT + "{d}/{s}/{l}/{b}/{i}/family_sizes.txt",
       mptp = OUTPUT+ "{d}/{s}/{l}/{b}/{i}/mptp_data_singletons.txt",
       changeo = OUTPUT + "{d}/{s}/{l}/{b}/{i}/vquest_files/combined_db-pass_clone-pass.tsv",
       IGH = OUTPUT + "{d}/{s}/{l}/{b}/{i}/clean.fasta.vdjca.clns_IGH.tsv"
    params:
       dir = OUTPUT + "{d}/{s}/{l}/{b}/{i}/"
    output:
       complete = OUTPUT + "{d}/{s}/{l}/{b}/{i}/analysis_complete.tsv",
       without_singletons = OUTPUT + "{d}/{s}/{l}/{b}/{i}/analysis_no_singletons.tsv"
    shell:
       "echo " + platform.node() + " >> {log} 2>&1 && \
       export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log}&& \
       echo {params.dir} >> {log} 2>&1 && \
       python {input.script} {params.dir} >> {log} 2>&1"

# use this if not all tools return a result (eg for fake V gene sims)
# rule combine_analysis:
#     resources:
#        mem="2G",
#     threads: 10
#     log: os.path.join(OUTPUT, "logs", "combine_analysis_{d}_{s}_{l}_{b}_{i}.log")
#     input:
#        script = 'simulation_analyses/calculate_measures.py',
#        real = OUTPUT + "{d}/{s}/{l}/{b}/{i}/family_sizes.txt",
#        mptp = OUTPUT+ "{d}/{s}/{l}/{b}/{i}/mptp_data_singletons.txt",
#     params:
#        dir = OUTPUT + "{d}/{s}/{l}/{b}/{i}/"
#     output:
#        complete = OUTPUT + "{d}/{s}/{l}/{b}/{i}/analysis_complete.tsv",
#        without_singletons = OUTPUT + "{d}/{s}/{l}/{b}/{i}/analysis_no_singletons.tsv"
#     shell:
#        "echo " + platform.node() + " >> {log} 2>&1 && \
#        export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log}&& \
#        echo {params.dir} >> {log} 2>&1 && \
#        python {input.script} {params.dir} >> {log} 2>&1"


#out = directory(OUTPUT + "{d}/{s}/{l}/{b}/{i}/ancestral_sequences/"),
# rule ancestral_sequence_true:
#     resources:
#        mem="2G",
#     threads: 10
#     log: os.path.join(OUTPUT, "logs", "ancestral_sequence_{d}_{s}_{l}_{b}_{i}.log")
#     input:
#        script = 'germline_search/ancestral_seq_raxml.sh',
#        raxml  = RAXML+"raxml-ng",
#        tree = OUTPUT + "{d}/{s}/{l}/{b}/{i}/tree_files/reroot_tree.txt"
#     params:
#         out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/",
#         dir = OUTPUT + "{d}/{s}/{l}/{b}/{i}/ancestral_sequences/"
#     output:
#        all_naive = OUTPUT + "{d}/{s}/{l}/{b}/{i}/ancestral_sequences/root_naive_rerooted.txt"
#     shell:
#        "echo " + platform.node() + " &>> {log} && \
#        export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log}&& \
#        sh {input.script} -r {input.raxml} -d {params.out} -o {params.dir} &>> {log}"



# rule scoper_hier:
#       resources:
#          mem="2G",
#       threads: 10
#       log: os.path.join(OUTPUT, "logs", "scoper_hier_{d}_{s}_{l}_{b}_{i}.log")
#       input:
#          script = 'simulation_analyses/scoper_hier.R',
#          db = OUTPUT + "{d}/{s}/{l}/{b}/{i}/vquest_files/combined_db-pass_clone-pass.tsv"
#       params:
#          dir = OUTPUT + "{d}/{s}/{l}/{b}/{i}/"
#       output:
#          hierarchical = OUTPUT + "{d}/{s}/{l}/{b}/{i}/results_hierClones_0_05.tsv",
#          hierarchical2 = OUTPUT + "{d}/{s}/{l}/{b}/{i}/results_hierClones_0_5.tsv"
#       shell:
#          "echo " + platform.node() + " &>> {log} && \
#          Rscript {input.script} -d {input.db} -o {params.dir}&>> {log}"

# rule extract_fasta:
#       resources:
#          mem="2G",
#       threads: 10
#       log: os.path.join(OUTPUT, "logs", "extract_fasta_{d}_{s}_{l}_{b}_{i}.log")
#       input:
#          script = 'simulation_analyses/extract_fasta.py',
#          spectral = OUTPUT + "{d}/{s}/{l}/{b}/{i}/results_specClones.tsv",
#          hierarchical = OUTPUT + "{d}/{s}/{l}/{b}/{i}/results_hierClones.tsv",
#          mptp = OUTPUT+ "{d}/{s}/{l}/{b}/{i}/mptp_data_singletons.txt",
#          changeo = OUTPUT + "{d}/{s}/{l}/{b}/{i}/vquest_files/combined_db-pass_clone-pass.tsv",
#       params:
#          dir = OUTPUT + "{d}/{s}/{l}/{b}/{i}/",
#          scoper_hier = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_hier/",
#          scoper_sp = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_sp/",
#          changeo = OUTPUT + "{d}/{s}/{l}/{b}/{i}/changeo/",
#          mptp = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mptp/"
#       output:
#          extract_check = OUTPUT + "{d}/{s}/{l}/{b}/{i}/extract.txt"
#       shell:
#          "echo " + platform.node() + " &>> {log} && \
#          export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log} && \
#          mkdir {params.scoper_hier} {params.scoper_sp} {params.mptp} {params.changeo}
#          echo extracted > {output.extract_check} && \
#          python {input.script} {params.dir} &>> {log}"


rule mixcr_extract_fastas:
     resources:
        mem="2G",
     threads: 32
     log: os.path.join(OUTPUT, "logs", "mixcr_extract_fastas_{d}_{s}_{l}_{b}_{i}.log")
     input:
        script = 'simulation_analyses/create_mixcr_fasta.sh',
        IGH = OUTPUT + "{d}/{s}/{l}/{b}/{i}/clean.fasta.vdjca.clns_IGH.tsv"
     params:
        dir = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mixcr_fastas/"
     output:
        check = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mixcr_fastas/mixcr_fastas.txt"
     shell:
        "echo " + platform.node() + " &>> {log} && \
        sh {input.script} -d {params.dir} &>> {log}"


# rule align_discerned_families:
#      resources:
#         mem="2G",
#      threads: 10
#      log: os.path.join(OUTPUT, "logs", "align_discerned_families_{d}_{s}_{l}_{b}_{i}.log")
#      input:
#         script = 'tree_building/align_families.sh',
#         extract_check = OUTPUT + "{d}/{s}/{l}/{b}/{i}/extract.txt",
#         check = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mixcr_fastas/mixcr_fastas.txt"
#      output:
#         mixcr = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mixcr_fastas/align.txt",
#         mptp = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mptp/align.txt",
#         changeo = OUTPUT + "{d}/{s}/{l}/{b}/{i}/changeo/align.txt",
#         scoper_hier = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_hier/align.txt",
#         scoper_sp = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_sp/align.txt"
#      params:
#          mixcr = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mixcr_fastas/",
#          mptp = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mptp/",
#          changeo = OUTPUT + "{d}/{s}/{l}/{b}/{i}/changeo/",
#          scoper_hier = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_hier/",
#          scoper_sp = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_sp/"
#      shell:
#         "echo " + platform.node() + " &>> {log} && \
#         sh {input.script} -d {params.mixcr} && \
#         sh {input.script} -d {params.mptp} && \
#         sh {input.script} -d {params.changeo} && \
#         sh {input.script} -d {params.scoper_hier} && \
#         sh {input.script} -d {params.scoper_sp} &>> {log}"


# rule build_tree_discerned:
#      resources:
#         mem="5G",
#      threads: 32
#      log: os.path.join(OUTPUT, "logs", "build_tree_discerned_{d}_{s}_{l}_{b}_{i}.log")
#      input:
#         script = 'tree_building/build_sub_trees.sh',
#         raxml  = RAXML+"raxml-ng",
#         mixcr = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mixcr_fastas/align.txt",
#         mptp = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mptp/align.txt",
#         changeo = OUTPUT + "{d}/{s}/{l}/{b}/{i}/changeo/align.txt",
#         scoper_hier = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_hier/align.txt",
#         scoper_sp = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_sp/align.txt"
#      params:
#          mixcr = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mixcr_fastas/",
#          mptp = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mptp/",
#          changeo = OUTPUT + "{d}/{s}/{l}/{b}/{i}/changeo/",
#          scoper_hier = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_hier/",
#          scoper_sp = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_sp/"
#      output:
#          mixcr = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mixcr_fastas/build_tree.txt",
#          mptp = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mptp/build_tree.txt",
#          changeo = OUTPUT + "{d}/{s}/{l}/{b}/{i}/changeo/build_tree.txt",
#          scoper_hier = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_hier/build_tree.txt",
#          scoper_sp = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_sp/build_tree.txt"
#      shell:
#         "echo " + platform.node() + " &>> {log} && \
#         sh {input.script} -r {input.raxml} -d {params.mixcr} -o {params.mixcr} && \
#         sh {input.script} -r {input.raxml} -d {params.mptp} -o {params.mptp} && \
#         sh {input.script} -r {input.raxml} -d {params.changeo} -o {params.changeo} && \
#         sh {input.script} -r {input.raxml} -d {params.scoper_hier} -o {params.scoper_hier} && \
#         sh {input.script} -r {input.raxml} -d {params.scoper_sp} -o {params.scoper_sp} &>> {log}"

# rule reroot_tree:
#      resources:
#         mem="5G",
#      threads: 32
#      log: os.path.join(OUTPUT, "logs", "reroot_tree_{d}_{s}_{l}_{b}_{i}.log")
#      input:
#         script = 'tree_building/reroot_midpoint.R',
#         tree = OUTPUT + "{d}/{s}/{l}/{b}/{i}/tree_files/build_tree.txt",
#         mixcr = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mixcr_fastas/build_tree.txt",
#         mixcr = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mixcr_fastas/build_tree.txt",
#         mptp = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mptp/build_tree.txt",
#         changeo = OUTPUT + "{d}/{s}/{l}/{b}/{i}/changeo/build_tree.txt",
#         scoper_hier = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_hier/build_tree.txt",
#         scoper_sp = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_sp/build_tree.txt"
#      params:
#          out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/tree_files/",
#          mixcr = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mixcr_fastas/",
#          mptp = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mptp/",
#          changeo = OUTPUT + "{d}/{s}/{l}/{b}/{i}/changeo/",
#          scoper_hier = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_hier/",
#          scoper_sp = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_sp/"
#      output:
#         check = OUTPUT + "{d}/{s}/{l}/{b}/{i}/tree_files/reroot_tree.txt",
#         mixcr = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mixcr_fastas/reroot_tree.txt",
#         mptp = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mptp/reroot_tree.txt",
#         changeo = OUTPUT + "{d}/{s}/{l}/{b}/{i}/changeo/reroot_tree.txt",
#         scoper_hier = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_hier/reroot_tree.txt",
#         scoper_sp = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_sp/reroot_tree.txt"
#      shell:
#         "echo " + platform.node() + " &>> {log} && \
#         Rscript {input.script} -d {params.out} &>> {log} && \
#         Rscript {input.script} -d {params.mixcr} &>> {log} && \
#         Rscript {input.script} -d {params.mptp} &>> {log} && \
#         Rscript {input.script} -d {params.changeo} &>> {log} && \
#         Rscript {input.script} -d {params.scoper_hier} &>> {log} && \
#         Rscript {input.script} -d {params.scoper_sp} &>> {log} && \
#         echo rerooted done > {output.check} &>> {log} && \
#         echo rerooted done > {output.mixcr} &>> {log} && \
#         echo rerooted done > {output.mptp} &>> {log} && \
#         echo rerooted done > {output.changeo} &>> {log} && \
#         echo rerooted done > {output.scoper_hier} &>> {log} && \
#         echo rerooted done > {output.scoper_sp} &>> {log} && \
#         &>> {log}"


# rule ancestral_sequence:
#     resources:
#        mem="10G",
#     threads: 10
#     log: os.path.join(OUTPUT, "logs", "ancestral_sequence_{d}_{s}_{l}_{b}_{i}.log")
#     input:
#        script = 'germline_search/ancestral_seq_raxml.sh',
#        raxml  = RAXML+"raxml-ng",
#        mixcr = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mixcr_fastas/build_tree.txt",
#        mptp = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mptp/build_tree.txt",
#        changeo = OUTPUT + "{d}/{s}/{l}/{b}/{i}/changeo/build_tree.txt",
#        scoper_hier = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_hier/build_tree.txt",
#        scoper_sp = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_sp/build_tree.txt"
#     params:
#        mixcr = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mixcr_fastas/",
#        mptp = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mptp/",
#        changeo = OUTPUT + "{d}/{s}/{l}/{b}/{i}/changeo/",
#        scoper_hier = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_hier/",
#        scoper_sp = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_sp/"
#     output:
#        mixcr = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mixcr_fastas/root_naive.txt",
#        mptp = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mptp/root_naive.txt",
#        changeo = OUTPUT + "{d}/{s}/{l}/{b}/{i}/changeo/root_naive.txt",
#        scoper_hier = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_hier/root_naive.txt",
#        scoper_sp = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_sp/root_naive.txt",
#        mixcr_r = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mixcr_fastas/root_naive_rerooted.txt",
#        mptp_r = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mptp/root_naive_rerooted.txt",
#        changeo_r = OUTPUT + "{d}/{s}/{l}/{b}/{i}/changeo/root_naive_rerooted.txt",
#        scoper_hier_r = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_hier/root_naive_rerooted.txt",
#        scoper_sp_r = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_sp/root_naive_rerooted.txt"
#     shell:
#        "echo " + platform.node() + " &>> {log} && \
#        export PATH=/home1/kavoss/anaconda2/bin:$PATH &>> {log} && \
#        sh {input.script} -r {input.raxml} -d {params.mixcr} -o {params.mixcr} &>> {log} && \
#        sh {input.script} -r {input.raxml} -d {params.mptp} -o {params.mptp} &>> {log} && \
#        sh {input.script} -r {input.raxml} -d {params.changeo} -o {params.changeo} &>> {log} && \
#        sh {input.script} -r {input.raxml} -d {params.scoper_hier} -o {params.scoper_hier} &>> {log} && \
#        sh {input.script} -r {input.raxml} -d {params.scoper_sp} -o {params.scoper_sp} &>> {log}"

# rule seq_similarity:
#     resources:
#        mem="2G",
#     threads: 10
#     log: os.path.join(OUTPUT, "logs", "seq_similarity_{d}_{s}_{l}_{b}_{i}.log")
#     input:
#        script = 'simulation_analyses/seq_similarity.py',
#        mixcr = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mixcr_fastas/root_naive.txt",
#        mptp = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mptp/root_naive.txt",
#        mptp_r = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mptp/root_naive_rerooted.txt",
#        changeo = OUTPUT + "{d}/{s}/{l}/{b}/{i}/changeo/root_naive.txt",
#        scoper_hier = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_hier/root_naive.txt",
#        scoper_sp = OUTPUT + "{d}/{s}/{l}/{b}/{i}/scoper_sp/root_naive.txt",
#        all_naive = OUTPUT + "{d}/{s}/{l}/{b}/{i}/ancestral_sequences/root_naive_rerooted.txt"
#     params:
#         out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/"
#     output:
#        similarity = OUTPUT + "{d}/{s}/{l}/{b}/{i}/seq_similarity.tsv"
#     shell:
#        "echo " + platform.node() + " &>> {log} && \
#        python {input.script} {params.out} &>> {log}"


rule mixcr_tsv:
    resources:
       mem="2G",
    threads: 10
    log: os.path.join(OUTPUT, "logs", "mixcr_tsv_{d}_{s}_{l}_{b}_{i}.log")
    input:
       script = 'simulation_analyses/create_mixcr_tsv.py',
       check = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mixcr_fastas/mixcr_fastas.txt"
    params:
        out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/"
    output:
       mixcr = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mixcr.tsv"
    shell:
       "echo " + platform.node() + " &>> {log} && \
       python {input.script} {params.out} &>> {log}"

rule sensitivity_precision:
    resources:
       mem="2G",
    threads: 10
    log: os.path.join(OUTPUT, "logs", "sensitivity_precision_{d}_{s}_{l}_{b}_{i}.log")
    input:
       script = 'simulation_analyses/sensitivity_precision.py',
       check = OUTPUT + "{d}/{s}/{l}/{b}/{i}/results_specClones.tsv"
    params:
        out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/"
    output:
       similarity = OUTPUT + "{d}/{s}/{l}/{b}/{i}/sensitivity_precision.tsv"
    shell:
       "echo " + platform.node() + " &>> {log} && \
       python {input.script} {params.out} &>> {log}"

# rule f1_all:
#     resources:
#        mem="2G",
#     threads: 10
#     log: os.path.join(OUTPUT, "logs", "f1_all_{d}_{s}_{l}_{b}_{i}.log")
#     input:
#        script = 'simulation_analyses/f1_all_sequences.py'
#     params:
#         out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/"
#     output:
#        similarity = OUTPUT + "{d}/{s}/{l}/{b}/{i}/f1_all.tsv"
#     shell:
#        "echo " + platform.node() + " &>> {log} && \
#        python {input.script} {params.out} &>> {log}"

# rule count_mutations:
#     resources:
#        mem="2G",
#     threads: 10
#     log: os.path.join(OUTPUT, "logs", "count_mutations_{d}_{s}_{l}_{b}_{i}.log")
#     input:
#        script = 'simulation_analyses/count_mutations.py',
#        check = OUTPUT + "{d}/{s}/{l}/{b}/{i}/naive.fasta"
#     params:
#         out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/"
#     output:
#        out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/mutations.tsv"
#     shell:
#        "echo " + platform.node() + " &>> {log} && \
#        python {input.script} {params.out} &>> {log}"


# rule ancestral_phangorn:
#     resources:
#        mem="2G",
#     threads: 10
#     log: os.path.join(OUTPUT, "logs", "count_mutations_{d}_{s}_{l}_{b}_{i}.log")
#     input:
#        script = 'germline_search/ancestral_seq_phangorn.R',
#        check = OUTPUT + "{d}/{s}/{l}/{b}/{i}/tree_files"
#     params:
#         out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/"
#     output:
#        out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/ancestral_sequences/ancestral_seqs_parsim.tsv"
#     shell:
#        "echo " + platform.node() + " &>> {log} && \
#        Rscript {input.script} -d {params.out} &>> {log} && \
#        rm -rf  {params.out}*_ancestral_parsimony.fasta &>> {log} && \
#        rm -rf  {params.out}*_ancestral_ml.fasta &>> {log}"


# rule ham_distance:
#     resources:
#        mem="2G",
#     threads: 100
#     log: os.path.join(OUTPUT, "logs", "ham_distance_{d}_{s}_{l}_{b}_{i}.log")
#     input:
#        script = 'simulation_analyses/ham_dist_figure.py',
#        check = OUTPUT + "{d}/{s}/{l}/{b}/{i}/clean.fasta"
#     params:
#         out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/"
#     output:
#        out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/ham_distance.tsv"
#     shell:
#        "echo " + platform.node() + " &>> {log} && \
#        python {input.script} {params.out} &>> {log}"


# rule distance_naive:
#     resources:
#        mem="2G",
#     threads: 100
#     log: os.path.join(OUTPUT, "logs", "ham_distance_{d}_{s}_{l}_{b}_{i}.log")
#     input:
#        script = 'simulation_analyses/ham_dist_naive.py',
#        check = OUTPUT + "{d}/{s}/{l}/{b}/{i}/clean.fasta",
#        naive = OUTPUT + "{d}/{s}/{l}/{b}/{i}/naive.fasta"
#     params:
#         out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/"
#     output:
#        out = OUTPUT + "{d}/{s}/{l}/{b}/{i}/distance.tsv"
#     shell:
#        "echo " + platform.node() + " &>> {log} && \
#        python {input.script} {params.out} &>> {log}"



