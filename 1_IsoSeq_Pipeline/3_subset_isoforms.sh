#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=1:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=3_subset_isoforms.o
#SBATCH --error=3_subset_isoforms.e

## print start date and time
echo Job started on:
date -u


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
export SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/UCL_FECD
source $SC_ROOT/1_IsoSeq_Pipeline/fecd_isoseq.config
source $SC_ROOT/1_IsoSeq_Pipeline/01_source_functions.sh


##-------------------------------------------------------------------------

## parameters 
# loop through each folder, sqanti files and condition
# ** Important ** parameters are in order!!

# output directory 
dir=(
  $WKD_ROOT/9_subset/1_sqanti3 
  $WKD_ROOT/9_subset/2_sqanti3_filtered 
  $WKD_ROOT/9_subset/3_sqanti3_tama_filtered)

# input sqanti classification file
sqfiles=(
  $WKD_ROOT/8_sqanti3/basic_exp/Fecd.collapsed_classification.txt
  $WKD_ROOT/8_sqanti3/basic_exp/Fecd.collapsed_classification.filtered_lite_classification.txt
  $WKD_ROOT/8b_filter_cont/tama/basic_exp/Fecd_sqantitamafiltered_classification.txt
)

# Conditions to split the files 
cond=(Case Control Case_Control)

# parameter for filtering 
# essential for internal script to denote gtf file name
# 0 == not sqanti filtered
# 1 == sqanti filtered
filter_stage=(0 1 1)


##-------------------------------------------------------------------------

# loop through sqanti datasets
for i in {0..2}; do 
  
  source activate nanopore
  
  # make output directory if not already present
  mkdir -p ${dir[i]}; cd ${dir[i]}
  
  ### Subset by condition and by sample
  # subset cases and control by counts
  Rscript $SQCOUNT -f ${sqfiles[i]} -m $WKD_ROOT/0_metadata/sample_metadata.txt -o ${dir[i]}
  
  # subset for each sample by counts 
  Rscript $SQCOUNT_SAMPLE -f ${sqfiles[i]} -m $WKD_ROOT/0_metadata/sample_metadata.txt -o ${dir[i]}
  
  ### generate classification.txt, gtf file, and corresponding SQANTI report
  sq_dir=$(dirname ${sqfiles[i]})
  sqname=$(basename ${sqfiles[i]})
  sqname_prefix=${sqname//"_classification.txt"/}
  
  # for each condition (case, control, case_control)
  for c in ${cond[@]}; do 
    
    echo "Processing $c"
    
    Rscript $SQSUBSET -i $c"_ID.txt" -d $sq_dir -s $sqname_prefix -n $c -o ${dir[i]} -f ${filter_stage[i]}
    Rscript $SQ_Report ${dir[i]}/$c"_classification.txt" ${dir[i]}/$c"_junctions.txt"
    
    if [ $c == "Case_Control" ] ; then
      # for Case and control, convert gtf to bed file 
      convert_gtf_bed12 ${dir[i]}/$c.gtf 
      source activate nanopore
      
      # colour the bed file by abundance 
      python $ISOCOL --bed ${dir[i]}/$c"_sorted.bed12" --a ${dir[i]}/$c"_CaseAbundance.csv" --o ${dir[i]}/$c"sorted_Case_Coloured.bed12"
      python $ISOCOL --bed ${dir[i]}/$c"_sorted.bed12" --a ${dir[i]}/$c"_ControlAbundance.csv" --o ${dir[i]}/$c"sorted_Control_Coloured.bed12"
    fi
  
  done
  
  # for each sample 
  for s in *Sample_ID.txt*; do 
    echo $s
    output_prefix=${s//"_Sample_ID.txt"/}
    Rscript $SQSUBSET -i $s -d $sq_dir -s $sqname_prefix -n $output_prefix -o ${dir[i]} -f ${filter_stage[i]}
  done
  
done


