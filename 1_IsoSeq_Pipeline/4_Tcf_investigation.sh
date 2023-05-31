#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the SAMPLE job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=5:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=4_Tcf.o
#SBATCH --error=4_Tcf.e

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
export SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/UCL_FECD
source $SC_ROOT/1_IsoSeq_Pipeline/fecd_isoseq.config
source $SC_ROOT/1_IsoSeq_Pipeline/01_source_functions.sh

source activate nanopore
# Convert sam file into paf format, which contains information about matches, mismatches,insertions and deletions for every single read, also provides information on directionality
#htsbox samview -pS ${WKD_ROOT}/6_minimap/${NAME}.sorted.sam > ${WKD_ROOT}/6_minimap/${NAME}.paf

# Remove reads with no alignments from the paf file
# Colum 6 = chromosome number, * are no chromosomal aligments
#awk -F'\t' '{if ($6!="*") {print $0}}' ${WKD_ROOT}/6_minimap/${NAME}.paf > ${WKD_ROOT}/6_minimap/${NAME}_filtered.paf

# convert bam to fasta 
#refine2fasta ${WKD_ROOT}/3_refine

python ${EXTRACTPAF} --c_fa=${WKD_ROOT}/5_merged_cluster/${NAME}.clustered.hq.fasta \
  --r_fa=${WKD_ROOT}/3_refine/All.flnc.fasta \
  --paf=${WKD_ROOT}/6_minimap/${NAME}.paf \
  --o_name=${NAME}_Tcf4 \
  --o_dir=${WKD_ROOT}/10_characterisation/1_TCF
  
#run_minimap2 <input_fasta> <output_dir>
#run_minimap2 ${WKD_ROOT}/10_characterisation/1_TCF/${NAME}_Tcf4_transcripts.fasta ${WKD_ROOT}/10_characterisation/1_TCF
run_minimap2 ${WKD_ROOT}/10_characterisation/1_TCF/${NAME}_Tcf4_ccs.fasta ${WKD_ROOT}/10_characterisation/1_TCF


# subsetting TCF novel transcripts
source activate sqanti2_py3
Rscript $SC_ROOT/1_IsoSeq_Pipeline/4_Tcf_investigation.R
for i in Control Case_control Case; do 
  echo $i
  ID=${WKD_ROOT}/10_characterisation/1_TCF/${i}_classification_TCF4NovelTranscripts.txt
  python ${LOGEN}/miscellaneous/subset_fasta_gtf.py \
    ${WKD_ROOT}/8_sqanti3/basic_exp/Fecd.collapsed_classification.filtered_lite.gtf \
    -i ${ID} --gtf -o ${i} -d ${WKD_ROOT}/10_characterisation/1_TCF
done
