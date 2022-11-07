##### Code to simulate sequences #####
## simulate sequences in partis


########## Partis Original ##########

##### Step 1: Parameter directory ####

## creating parameter directory from a fasta file 
while getopts p:f:o: flag
do
    case "${flag}" in
        p) partis=${OPTARG};;
        f) fasta=${OPTARG};;
        o) output=${OPTARG};;
    esac
done

$partis cache-parameters --infname $fasta  --parameter-dir $output --outfname $output/simu.yaml

# partition the fasta file used in cache-parameters to create multi-hmm folder in the parameter directory
$partis partition --infname $fasta --outfname $output/pd.yaml --parameter-dir $output --count-parameters

##### Step 2: Simulate Sequences #####

## simulate and rearrange from scratch:

#set x for number of desired simulations $(seq 1 x)
for i in $(seq 1 5); do echo “$i”; $partis simulate --parameter-dir $output --n-sim-events integer --outfname $i.yaml; done
 

#########################################################################################################################

########## Partis Second Version ##########
# run the simulation on Docker
# code files are in a folder called python_test
# change the input for n-sim-events, n-leaves, and shm-freq
 
# autosimulator function
#python python_test/autosimulator2.py --species human --n-procs 1 --yaml-folder run_1_yaml --n-sim-events 5 --n-leaves 5 --fasta-folder run_1_fasta --shm-freq 0.05 --generate-germline 

# move the yaml files from the sub folders in run_1_yaml into a separate folder called run_yaml (there should be no subfolders)

# partition function
# change the n-max-final-cluster and n-partitions-to-write (default is 10 but it may write less)
#python python_test/auto_partition.py --n-procs 1 --n-max-final-cluster <input value> --n-partitions-to-write  <input value> --input-folder run_yaml --output-folder run_parition  --is-simu
 
# partition summarizer function
#python python_test/partition_summarizer.py --input-folder run_parition --outfname run_summarizer --is-simu


