# subset TCF4 novel transcripts into .txt files by case, control, case and control

# input files
dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Collabs/UCL_FECD/10_characterisation/1_TCF/"
TCF4 = read.csv(paste0(dir,"/TCF4_noveltranscripts.csv"))

# loop to extract the isoform column, and write to file
for(i in unique(TCF4$From)){
  
  cat("Extracting isoforms from", i, "\n")
  write.table(TCF4[TCF4$From %in% i,"Isoform"], file = paste0(dir,"/",i,"_TCF4NovelTranscripts.txt"), 
              quote = F, row.names = F, col.names = F, sep = "\t")
  
}