# clonal-misclassification

![Pipeline](https://github.com/katalinandreavoss/clonal-misclassification/blob/main/dag.svg?raw=true)

This is a guide to the pipeline with everything you need to know. It is in chronological order.


## Installation
clone this github repo
```
git clone https://github.com/katalinandreavoss/clonal-misclassification.git
```

Make sure you have installed all the dependencies:

partis: https://github.com/psathyrella/partis/tree/88580ef77689d73faa4ae6197a03abd41926e53c
I installed it from scratch. Make sure you also have installed all the dependencies of partis (specifically the ones you need to simulate)

raxml: https://github.com/amkozlov/raxml-ng

vquest: this one is more complicated because they do not have an API. I use https://github.com/ShawHahnLab/vquest with a few tweaks to make sure I get the correct output: https://github.com/katalinandreavoss/v-quest-excel


mptp: https://github.com/Pas-Kapli/mptp

VDJ: download IGHV.fasta, IGHD.fasta, IGHDJ.fasta from https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/ and put in a folder


Once you have completed all those steps you need so specify the locations of all these tools and files in the config file in ``configs/config.yaml``

### Additional installations:

snakemake,clustalo,mixcr,changeo,scoper,seqtk,...

a complete list of dependencies can be found in ``requirements.txt``


## Simulations
### partis simulations

``partis/partis_simulation/simulate.sh`` simulates B-cell data from scratch with partis  
parameters:  
        p) partis directory  
        o) output directory  
        c) number of simulated clones  
        s) SHM rate  
        l) avg number of leaves per clone  
        b) balance of tree (do 0_0, anything else is not supported by partis)  
        i) number of simulations per parameter set  

``analyze_partis_output/yaml_to_families.py`` creates fasta files from yaml returned by partis
parameters:
        1st parameter: .yaml file from partis
        2nd parameter: output directory

### our simulations
``simulation/simulate.R`` currently always simulates 10 clones  
parameters:  
        d) output directory  
        r) SHM rate  
        p) directory of VDJ folder  
        l) avg number of leaves per clone  
        j) junction length

``simulation/simulate_new_V_genes.R`` currently always simulates 10 clones and junction length=12
This script creates "fake" V genes by introducing mutations and indels into the sequences from IMGT. The rest of the simulation is the same as ``simulation/simulate.R``
parameters:  
        d) output directory  
        r) SHM rate  
        p) directory of VDJ folder  
        l) avg number of leaves per clone  
        j) amount of mutations introduced

## Processing
### remove singletons
``tree_building/remove_singletons.R`` remove singletons from fasta
parameters:
        f) fasta file
        o) output directory

### align
``tree_building/align_partitions.sh`` align clean.fasta
parameters:
        d) directory where clean.fasta is located

### build tree
``tree_building/build_tree.sh`` build tree
parameters:
        r) raxml location
        d) directory where aligned file is located
        o) output directory

### get correct family sizes
``simulation_analyses/get_real_family_sizes.sh`` get correct family sizes from simulation
parameters:
        d) input directory
        o) output directory  
        c) number of simulated clones  
        s) SHM rate  
        l) avg number of leaves per clone 
        i) number of simulations per parameter set  


## Apply Tools
### mPTP
rule cut_tree_mptp applies mPTP to the megatree

``simulation_analyses/analyse_ptp_output.py`` extract mPTP results (not counting singletons)
parameters: 1st parameter directory where mPTP results are located

``simulation_analyses/analyse_ptp_output_with_singletons.py`` extract mPTP results (counting singletons)
parameters: 1st parameter directory where mPTP results are located

### MiXCR
rule mixcr applies mixcr to clean.fasta

### Change-O and SCOPer

``germline_search/IMGT_vrequest.sh`` use IMGT vrequest API on clean.fasta
parameters: 
        f) fasta file
        v) path to VQUEST API
        o) output path

``germline_search/combine_vquest.py`` combine VQUEST results (they get split up because it can only handle files with 50 sequences at a time)
parameters: 1st parameter: directory where vquest files are saved

``simulation_analyses/change-o.sh`` script for applying Change-O to data
parameters: 
        d) input directory
        f) clean.fasta
        v) path to directory with VDJ data from IMGT

``simulation_analyses/scoper.R`` script for applying SCOPer (both models) to data
parameters:
        d) db from Change-O
        o) output directory


## Evaluation

``simulation_analyses/calculate_measures.py`` calculate MSE for all tools
parameters: 1st param: input directory

``simulation_analyses/extract_fasta.py`` extract fastas from resuls of all tools except MiXCR (1 fasta file per delimited family)
parameters: 1st param: input directory

``simulation_analyses/create_mixcr_fasta.sh`` extract fastas from resuls of MiXCR (1 fasta file per delimited family)
parameters: d) input directory

``simulation_analyses/create_mixcr_tsv.py`` create tsv for MiXCR results so it matches the rest of the tools' output
parameters: 1st param: input directory

``simulation_analyses/sensitivity_precision.py`` calculate measures dependent on TP,TN, FP,FN for all tools (remove singletons before)
parameters: 1st param: input directory

``simulation_analyses/f1_all_sequences.py`` calculate measures dependent on TP,TN, FP,FN for all tools (include singletons)
parameters: 1st param: input directory

## Ancestral Sequence

``tree_building/align_families.sh`` align all sequences from the families discerned by a tool (need to extract fastas first)
parameters: d) directory of the extracted fastas

``tree_building/build_sub_trees.sh`` build trees for the discerned families
parameters:
        d) directory of aligned sequences
        r) location of raxml
        o) output directory

``tree_building/reroot_midpoint.R`` reroot trees at the midpoint
parameters: d) directory of the trees

``germline_search/ancestral_seq_raxml.sh`` reconstruct ancestral sequence based on tree
parameters: 
        d) directory of the trees
        r) raxml
        o) output directory

``simulation_analyses/seq_similarity.py`` compare naive to reconstructed ancestral sequences
parameters: 1st param: input directory


## Extra Evaluations
``simulation_analyses/count_mutations.py`` count number of mutations introduced by SHM in simulations
parameters: 1st param: directory where naive.fasta and clean.fasta are located

``simulation_analyses/ham_dist_figure.py`` calculate ham distance between all simulated sequences, special focus on max. dist. within a family and min. dist. between families
parameters: 1st param: directory where naive.fasta and clean.fasta are located

``simulation_analyses/ham_dist_naive.py`` calculate ham dist. from sequences to naive germline sequence
parameters: 1st param: directory where naive.fasta and clean.fasta are located
