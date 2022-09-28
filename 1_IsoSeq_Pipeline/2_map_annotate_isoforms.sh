#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=20:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=2_map_annotate_isoforms.o
#SBATCH --error=2_map_annotate_isoforms.e


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
export SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/UCL_Collab
source $SC_ROOT/1_IsoSeq_Pipeline/fecd_isoseq.config
source $SC_ROOT/1_IsoSeq_Pipeline/01_source_functions.sh


##-------------------------------------------------------------------------
echo "#************************************* Isoseq3 and Post_Isoseq3 [Function 5, 7, 8]"
## 5) merging_at_refine <output_name> <samples.....>
merging_at_refine $NAME ${SAMPLE_NAMES[@]}

## 7) run_map_cupcakecollapse <output_name> 
run_map_cupcakecollapse $NAME 

## 8) demux <output_name>
demux $NAME


##-------------------------------------------------------------------------
echo "#************************************* SQANTI3 [Function 12]"
## 12) run_sqanti3 <sample> <mode=basic/full/nokallisto/lncrna/basic_exp>
run_sqanti3 $NAME basic_exp


##-------------------------------------------------------------------------
echo "#************************************* TAMA filter [Function 13,14,16]"
## 13) TAMA_remove_fragments <sample> <mode=basic/full/nokallisto/lncrna/basic_exp>
TAMA_remove_fragments $NAME basic_exp

## 14) TAMA_sqanti_filter <sample> <mode=basic/full/nokallisto/lncrna>
TAMA_sqanti_filter $NAME basic_exp

