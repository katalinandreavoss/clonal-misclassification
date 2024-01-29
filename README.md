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

vquest: this one is more complicated because they do not have an API. I use https://github.com/ShawHahnLab/vquest with a few tweaks to make sure I get the correct output:


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

### our simulations
``simulation/simulate.R`` currently always simulates 10 clones
parameters:  
        d) output directory  
        r) SHM rate 
        p) directory of VDJ folder 
        l) avg number of leaves per clone  
        j) junction length

## GMYC

``gmyc_clusters.R`` goes through each step in running the GMYC

``gmyc_mega.R`` runs the GMYC function on the mega trees as a loop

``gmyc_output_formatting.R`` finds the clusters in each tree

## Pseudocode

**part1**
  
``part_1_code.R`` project information

**part2**

``part_2_code.R``project information

## analyze_partis_output

This folder is to partition non-simulated sequences.

``partition_to_fasta.py`` moves each _partis_ partition into a fasta file

``python_partitioning.txt`` instructions on using ``partition_to_fasta.py``  and ``yaml_to_families.py``

``yaml_to_families.py`` converts each yaml file to family files

## analyze_partis_simulation_output

This folder is to partition simulated sequences.

``fasta_per_partition.py`` moves all of the same numbered partitions into a file

``partition_to_fasta.py``  moves each _partis_ partition into a fasta file

``yaml_to_families.py`` converts each yaml file to family files

## germline_search

``find_germline_phangorn.R`` use phangorn for germline sequence search

``germline_search_RevertToGermline`` instructions on using RevertToGermline

## partis

**partis_simulation**

``simulate_partition_steps.txt`` instructions on using partis

## resampling

``resample_families.R`` resample family steps

``resample_families_looped.R`` resample loops to run on all families

## simulation_analyses

This folder corresponds to additional family assignment methods and figures for manuscript:

``ancestral_similarity_figure.R`` figure for different reshuffling levels ancestral sequence similarity

``change-o.sh`` run CHANGE-O

``change-o_analysis.R`` converting CHANGE-O output

``diptest_phylogenetic_distances.R`` calculate phylogenetic distances and distributions

``family_distribution_figure.R`` figure for different methods

``family_distribution_figure_monophyletic_only.R`` figure for different methods on monophyletic families only

``format_families.sh`` steps for data processing

``medusa_analysis.R`` running MEDUSA

``megatree_summary.R`` investigate the trees that are monophyletic

``mixcr.sh`` run MiXCR

``mixcr_analysis.R`` convert MiXCR output

``outputs_trees_5_families.R`` analyzing only monophyletic families

``partis_distribution_output.R`` convert _partis_ output

## tree_building

``ancestral_search_RAxML.txt`` instructions for ancestral sequence search

``build_tree.txt`` instructions for tree building 

``simulations_asr.sh`` ancestral sequence search on simulated data

``simulations_build_family_trees.sh`` family level tree building for simulated data

``simulations_build_mega_trees.sh`` mega tree building for simulated data
