################################################################################################
## Define functions for All_Demultiplex_Final_Functions.sh
# 1) run_CCS <input_ccs_bam> <prefix_output_name> 
# 2) run_LIMA <sample> <"no_multiplex"/"multiplex">
# 3) run_REFINE <sample> 
# 5) merging_at_refine <output_name> <samples.....>
# 6) run_CLUSTER $Sample 
# 7) run_map_cupcakecollapse <sample_prefix_input/output_name> 
# 8) demux_targeted <sample>
# 9) run_star <list_of_samples> <J20/Tg4510_input_directory> <output_dir>
# 10) rnaseq_merge_fastq <sample>
# 11) run_kallisto <sample_prefix_output_name> <input_tofu_fasta> <merged_fastq_input_dir> <output_dir>
# 12) run_sqanti3 <sample> <mode=basic/full/nokallisto/lncrna>
# 13) TAMA_remove_fragments <sample> <mode=basic/full/nokallisto/lncrna>
# 14) TAMA_sqanti_filter <sample> <mode=basic/full/nokallisto/lncrna>
# 15) parse_stats_per_sample <sample>
# 16) on_target_rate <input_probe_bed.file> <input_polished.hq.fasta> <input_mappped.fastq.sam> <output_file_name>
# 17) run_target_rate 
# 18) counts_subset_4tappas <input_class> <output_class> <type_genes>
# 19) TAMA_tappas_input <sample> <mode=basic/full/nokallisto/lncrna>
################################################################################################

#************************************* DEFINE VARIABLES
module load Miniconda2/4.3.21

# Listing versions
source activate isoseq3
ccs --version #ccs 5.0.0 (commit v5.0.0)
lima --version #lima 2.0.0 (commit v2.0.0)
isoseq3 --version #isoseq3 3.4.0 (commit v3.4.0)

source activate sqanti2_py3
echo "Minimap2 version:" $(minimap2 --version) # echo version
echo "ToFU Cupcake Version"
head $CUPCAKE/README.md

echo "FASTA SEQUENCE (CLONTECH PRIMERS) FOR NON-MULTIPLEXING"
cat $FASTA
source deactivate

################################################################################################
#*************************************  Isoseq3 [Function 1,2,3,4,5,6]
# run_CCS <input_ccs_bam> <prefix_output_name> 
run_CCS(){
  source activate isoseq3_paper
  
  cd $WKD_ROOT/1_ccs
  echo "Processing Sample $2 from Functions script"
  echo "Processing $1 with old chemistry"
  ccs --version
  
  # ccs <input.subreads.bam> <output.ccs.bam>
  # --minPredictedAccuracy 0.9 = default
  time ccs $1 $2.ccs.bam --minPasses=1 --reportFile $2_ccs_report.txt
  echo "CCS for Sample $2 successful"
  ls *$2*    
  
  source deactivate
}

# run_LIMA <sample> 
# removed --no pbi as this is needed for downstream polishing
run_LIMA(){
  source activate isoseq3
  
  echo "Processing $1 file for demultiplexing"
  
  if [ -f $WKD_ROOT/2_lima/$1.fl.json ]; then
    echo "$1.fl.bam file already exists; LIMA no need to be processed on Sample $1"
  else
    cd $WKD_ROOT/2_lima
    time lima $WKD_ROOT/1_ccs/$1.ccs.bam $FASTA $1.fl.bam --isoseq --dump-clips --dump-removed
    echo "lima $1 successful"
    ls $1.fl*
  fi
  source deactivate
}

# run_REFINE <sample>
run_REFINE(){
  source activate isoseq3
  
  cd $WKD_ROOT/3_refine
  
  echo "Processing $1 file for refine"
  if [ -f $1.flnc.bam ]; then
  echo "$1.flnc bam file already exists; Refine no need to be processed on Sample $1"
  else
    #refine --require-polya <input.lima.consensusreadset.xml> <input.primer.fasta> <output.flnc.bam>
    time isoseq3 refine $WKD_ROOT/2_lima/$1.fl.primer_5p--primer_3p.bam $FASTA $1.flnc.bam --require-polya
  echo "refine $1 successful"
  ls $1.flnc*
    fi
  
  source deactivate
}


# merging_at_refine <output_name> <samples.....>
# aim: merging bam files from refine onwards (similar to run_isoseq3_2_1_merge, but no need to rerun from ccs)
# <output_name> = output name for merged files
# samples.... = list of the sample names
merging_at_refine(){
  
  source activate isoseq3
  
  ###********************* Merging at REFINE
  # Define variable "Merge_Samples" as a list of all samples, in order to find the specified flnc.bam (for dataset create ConsensusReadSet)
  # Define variable "all_flnc_bams" for merged path-directory of all flnc samples (for dataset create ConsensusReadSet)
  Merge_Samples=$(echo "${@:2}")
  
  echo "Merging flnc of samples $Merge_Samples"
  all_flnc_bams=$(
    for i in ${Merge_Samples[@]}; do
    flnc_bam_name=$(find $WKD_ROOT/3_refine -name "*.flnc.bam" -exec basename \{} \; | grep ^$i.flnc )
    flnc_bam=$(find $WKD_ROOT/3_refine -name "*.flnc.bam" | grep "$flnc_bam_name" )
    echo $flnc_bam
    done
  )
  
  cd $WKD_ROOT/5_merged_cluster
  printf '%s\n' "${all_flnc_bams[@]}" > $1.flnc.fofn
  cat $1.flnc.fofn
  
  ###*********************
  
  isoseq3 cluster $1.flnc.fofn $1.clustered.bam --verbose --use-qvs
  gunzip *.gz*
    
    source deactivate
}


# run_CLUSTER <sample>
run_CLUSTER(){
  source activate isoseq3
  
  cd $WKD_ROOT/4_cluster
  echo "Processing $1 file for cluster"
  if [ -f $1.unpolished.bam ]; then
    echo "$1.unpolished.bam file already exists; Cluster no need to be processed on Sample $1"
  else
    # cluster <input.flnc.bam> <output.unpolished.bam>
    time isoseq3 cluster $WKD_ROOT/3_refine/$1.flnc.bam $1.clustered.bam --verbose --use-qvs 2> $1.cluster.log
    echo "cluster $1 successful"
    ls $1.clustered*
  fi
  
  source deactivate
}


################################################################################################
#************************************* Post_Isoseq3 (Minimap2, Cupcake, Demultiplex) [Function 7,8]
# run_map_cupcakecollapse <output_name> 
# Prerequisite: mm10 cage peak
run_map_cupcakecollapse(){
  
  
  source activate sqanti2_py3
  
  # convert between fasta and fastq for downstream process
  echo "fasta2fastq conversion"
  python $SEQUENCE/fa2fq.py $WKD_ROOT/5_merged_cluster/$1.clustered.hq.fasta
  
  echo "Minimap2 version:" $(minimap2 --version) # echo version
  
  echo "Processing Sample $1 for Minimap2 and sort"
  cd $WKD_ROOT/6_minimap 
  minimap2 -t 30 -ax splice -uf --secondary=no -C5 -O6,24 -B4 $GENOME_FASTA $WKD_ROOT/5_merged_cluster/$1.clustered.hq.fastq > $1.sam 2> $1.map.log
  source activate nanopore
  samtools sort -O SAM $1.sam > $1.sorted.sam
  
  echo "Processing Sample $1 for TOFU, with coverage 85% and identity 95%"
  source activate cupcake
  cd $WKD_ROOT/7_tofu
  collapse_isoforms_by_sam.py -c 0.85 -i 0.95 --input $WKD_ROOT/5_merged_cluster/$1.clustered.hq.fastq --fq -s $WKD_ROOT/6_minimap/$1.sorted.sam --dun-merge-5-shorter -o $1 &>> $1.collapse.log
  get_abundance_post_collapse.py $1.collapsed $WKD_ROOT/5_merged_cluster/$1.clustered.cluster_report.csv &>> $1.abundance.log
  
  source activate sqanti2_py3
  # convert rep.fq to rep.fa for SQANTI2 input
  seqtk seq -a $1.collapsed.rep.fq > $1.collapsed.rep.fa
  echo "Processing Sample $1 for TOFU successful"
  source deactivate
}



# demux <output_name>
# read in read.stat file from cupcake collapse output and count abundance of each sample based on CCS ID obtained from each refine report of the sample
demux(){
  source activate nanopore
  
  # script.R <refine_dir> <input_cluster_report> <input_tofu_readstat> <output_path_file>
  Rscript $DEMUXFUNCTIONS $WKD_ROOT/3_refine \
  $WKD_ROOT/5_merged_cluster/$1.clustered.cluster_report.csv \
  $WKD_ROOT/7_tofu/$1.collapsed.read_stat.txt \
  $WKD_ROOT/7_tofu/$1.Demultiplexed_Abundance.txt
}



################################################################################################
#************************************* SQANTI3 [Function 12]

# run_sqanti3 <sample> <mode=basic/full/nokallisto/lncrna/basic_exp>
run_sqanti3(){
  
  source activate sqanti2_py3
  
  # variable 
  sample=$1.collapsed
  gtf=$WKD_ROOT/7_tofu/$1.collapsed.gff
  abundance=$WKD_ROOT/7_tofu/$1.Demultiplexed_Abundance.txt
  
  
  # create directory
  mkdir -p $WKD_ROOT/8_sqanti3/$2; cd $WKD_ROOT/8_sqanti3/$2
  
  # copy STAR output SJ.out.bed files
  SJ_OUT_BED=($(
    for rnaseq in ${RNASEQ_SAMPLES_NAMES[@]}; do
    name=$(find $RNASEQ_MAPPED_DIR -name "*SJ.out.bed" -exec basename \{} \; | grep ^$rnaseq)
    File=$(find $RNASEQ_MAPPED_DIR -name "$name")
    echo $File
    done
  ))
  for file in ${SJ_OUT_BED[@]}; do cp $file $WKD_ROOT/8_sqanti3/with_junc ;done
  
  # prepare kallisto file 
  if [ $2 == "full" ] || [ $2 == "lncrna" ]; then 
  # create tab separated file from kallisto
  # Rscript script.R <input.file_fromkallisto> <output.file_intosqanti>
  Rscript $KALLSTOINPUT $WKD_ROOT/8_kallisto $WKD_ROOT/8_sqanti3/$2/$sample".mod2.abundance.tsv"
  kallisto_expfile=$WKD_ROOT/8_sqanti3/$2/$sample".mod2.abundance.tsv"
  echo "Using $kallisto_expfile"
  fi
  
  # sqanti qc
  echo "Processing Sample $sample for SQANTI3 QC"
  python $SQANTI3_DIR/sqanti3_qc.py -v
  echo $GENOME_GTF
  echo $GENOME_FASTA
  
  if [ $2 == "basic" ]; then
  echo "Processing basic commands"
  python $SQANTI3_DIR/sqanti3_qc.py -t 30 --gtf $gtf $GENOME_GTF $GENOME_FASTA --cage_peak $CAGE_PEAK --polyA_motif_list $POLYA \
  --genename --isoAnnotLite --gff3 $GFF3 --skipORF --report pdf &> $sample.sqanti.qc.log
  
   elif [ $2 == "basic_exp" ]; then
  echo "Processing full commands"
  python $SQANTI3_DIR/sqanti3_qc.py -t 30 --gtf $gtf $GENOME_GTF $GENOME_FASTA --cage_peak $CAGE_PEAK --polyA_motif_list $POLYA \
  --fl_count $abundance --skipORF \
  --genename --isoAnnotLite --gff3 $GFF3 --report pdf &> $sample.sqanti.qc.log
  
  
  elif [ $2 == "full" ]; then
  echo "Processing full commands"
  python $SQANTI3_DIR/sqanti3_qc.py -t 30 --gtf $gtf $GENOME_GTF $GENOME_FASTA --cage_peak $CAGE_PEAK --polyA_motif_list $POLYA \
  -c ./"*SJ.out.bed" --expression $kallisto_expfile --fl_count $abundance \
  --genename --isoAnnotLite --gff3 $GFF3 --report pdf &> $sample.sqanti.qc.log
  
  elif [ $2 == "nokallisto" ]; then 
  echo "Processing basic commands with RNA-Seq bed files but not kallisto"
  python $SQANTI3_DIR/sqanti3_qc.py -t 30 --gtf $gtf $GENOME_GTF $GENOME_FASTA --cage_peak $CAGE_PEAK --polyA_motif_list $POLYA \
  -c "./*SJ.out.bed" \
  --genename --isoAnnotLite --gff3 $GFF3 --report pdf &> $sample.sqanti.qc.log
  
  elif [ $2 == "lncrna" ]; then
  echo "Processing with $LNCRNA_GTF for genome annotation "
  python $SQANTI3_DIR/sqanti3_qc.py -t 30 --gtf $gtf $LNCRNA_GTF $GENOME_FASTA --cage_peak $CAGE_PEAK --polyA_motif_list $POLYA \
  -c "./*SJ.out.bed" --skipORF --fl_count $abundance \
  --genename --isoAnnotLite --gff3 $GFF3 --report pdf &> $sample.sqanti.qc.log
  
  else
    echo "2nd argument required"
  fi
  
  echo "Processing Sample $sample for SQANTI2 filter"
  python $SQANTI3_DIR/sqanti3_RulesFilter.py $sample"_classification.txt" $sample"_corrected.fasta" $sample"_corrected.gtf" -a 0.6 -c 3 &> $1.sqanti.filter.log
  
  if [ $2 != "basic" || $2 != "basic_exp" ]; then 
    # remove temp SJ.out bed files
    rm *SJ.out.bed
  fi
  
  source deactivate
}



################################################################################################
#************************************* TAMA [Function 13,14]
# TAMA_remove_fragments <sample> <mode=basic/full/nokallisto/lncrna>
# remove short fragments from post tofu
# Prerequisite: Require TAMA_prepare.R to change column 4 of bed12 to gene_name: transcript_name for correct file TAMA format
#TAMA_remove_fragments <sample> <mode=basic/full/nokallisto/lncrna>
TAMA_remove_fragments(){
  
  mkdir -p $WKD_ROOT/8b_filter_cont $WKD_ROOT/8b_filter_cont/tama
  mkdir -p $WKD_ROOT/8b_filter_cont/tama/$2; cd $WKD_ROOT/8b_filter_cont/tama/$2
  
  # variables 
  sample=$1.collapsed
  io_dir=$WKD_ROOT/8b_filter_cont/tama/$2
  
  source activate sqanti2
  
  # convert gtf to bed12
  gtfToGenePred $WKD_ROOT/8_sqanti3/$2/$sample"_classification.filtered_lite.gtf" $1.genepred
  genePredToBed $1.genepred $1.bed12
  awk -F'\t' '{print $1,$2,$3,$4,"40",$6,$7,$8,"255,0,0",$10,$11,$12}' $1.bed12| sed s/" "/"\t"/g|sed s/",\t"/"\t"/g|sed s/",$"/""/g > Tama_$1.bed12
  
  # Rscript script.R <name>.bed12 <input_dir>
  Rscript $TAMAMERGE Tama_$1 $io_dir
  python $TAMA_DIR/tama_remove_fragment_models.py -f Tama_$1"_mod.bed12" -o $1
    
  echo "Number of isoforms filtered by TAMA:"
  wc -l $sample"_discarded.txt"
  
  source deactivate
}

# TAMA_sqanti_filter <sample> <mode=basic/full/nokallisto/lncrna>
TAMA_sqanti_filter(){
  
  # variables
  sqname=$1.collapsed_classification.filtered_lite
  sq_dir=$WKD_ROOT/8_sqanti3/$2
  io_dir=$WKD_ROOT/8b_filter_cont/tama/$2
  
  cd $WKD_ROOT/8b_filter_cont/tama/$2
  
  source activate sqanti2_py3
  # Rscript .R <path/tama_filtered_output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_junc_txt> <output_prefix_name> <output_dir>
  Rscript $TAMASUBSET $io_dir/$1.bed $sq_dir $sqname"_classification.txt" $sqname".gtf" $sqname"_junctions.txt" $1 $io_dir
  
  # extract fasta sequence based on the pbid (https://www.biostars.org/p/319099/)
  # script.py <path/sqanti_filtered.fasta> <path/retained_pbid_tama.txt> <path/output.fasta>
  awk '{ print $4 }' $io_dir/$1.bed| cut -d ";" -f 2  > $io_dir/tama_retained_pbid.txt
  python $TAMASUBSETFASTA $sq_dir/$sqname".fasta" $io_dir/tama_retained_pbid.txt $io_dir/$1"_sqantifiltered_tamafiltered_classification.fasta"
  
  Rscript $SQ_Report $WKD_ROOT/8b_filter_cont/tama/$2/$1"_sqantitamafiltered.classification.txt" $WKD_ROOT/8b_filter_cont/tama/$2/$1"_sqantitamafiltered.junction.txt"
}

# remove_3ISM <sample> <mode=basic/full/nokallisto/lncrna>
remove_3ISM(){
  
  mkdir -p $WKD_ROOT/8b_filter_cont $WKD_ROOT/8b_filter_cont/no3ISM
  mkdir -p $WKD_ROOT/8b_filter_cont/no3ISM/$2; cd $WKD_ROOT/8b_filter_cont/no3ISM/$2
  
  # variables 
  sq_dir=$WKD_ROOT/8_sqanti3/$2
  
  # Rscript script.R <input.classfile> <input.gtf> <input.junc> <output.dir> <prefix>
  Rscript $ISMREMOVE $WKD_ROOT/8_sqanti3/$2/$1.collapsed_classification.filtered_lite_classification.txt \
  $sq_dir/$1.collapsed_classification.filtered_lite.gtf \
  $sq_dir/$1.collapsed_classification.filtered_lite_junctions.txt $WKD_ROOT/8b_filter_cont/no3ISM/$2 $1
  
}
