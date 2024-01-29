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

