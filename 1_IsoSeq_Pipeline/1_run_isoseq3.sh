#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the SAMPLE job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=144:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --array=0-8 #9 samples
#SBATCH --output=1_run_isoseq3-%A_%a.o
#SBATCH --error=1_run_isoseq3-%A_%a.e

## print start date and time
echo Job started on:
date -u
echo

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
export SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/UCL_FECD
source $SC_ROOT/1_IsoSeq_Pipeline/fecd_isoseq.config
source $SC_ROOT/1_IsoSeq_Pipeline/01_source_functions.sh


##-------------------------------------------------------------------------

# run as array (defined in config file)
SAMPLE=${SAMPLE_NAMES[${SLURM_ARRAY_TASK_ID}]}
BAM_FILE=${BAM_FILES[${SLURM_ARRAY_TASK_ID}]}

##-------------------------------------------------------------------------
echo "#*************************************  Isoseq3 [Function 1, 2]"
# Isoseq3.4.0
# run_CCS <sample>
# run_LIMA <sample> 
# run_REFINE <sample> 
# run_CLUSTER <sample>
run_CCS ${BAM_FILE} ${SAMPLE} 
run_LIMA ${SAMPLE} 
run_REFINE ${SAMPLE} 
run_CLUSTER ${SAMPLE} 
