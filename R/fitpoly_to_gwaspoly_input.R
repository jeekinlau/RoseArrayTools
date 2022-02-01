# fitpoly output to GWASpoly input


fitPoly_to_GWASpoly_input <- function(input_file){
  input_file_name = gsub(".csv","",input_file)
  file = read.csv(input_file)
  genomic_pos = read.table("https://raw.githubusercontent.com/jeekinlau/RoseArrayTools/master/docs/for_add_genomic_positions_saintoyant.txt", header = T, sep="\t")
  
  temp_names = names(file)
  temp_names[1] = "Marker"
  colnames(file) = temp_names
  
  merged = merge(file,genomic_pos, by="Marker", all=T)
  
  merged_genotypes = merged[,-which(colnames(merged)%in%c("Marker","LG","Position"))]
  
  final = cbind(merged$Marker,  merged$LG, merged$Position, merged_genotypes)
  final_colnames <- colnames(final)
  final_colnames[1:3] = c("Marker", "Chrom", "Position")
  colnames(final) = final_colnames
  final = final[order(final$LG, final$Position),]
  write.csv(final,paste0(input_file_name,"_GWAS_ready.csv"),row.names = F)
  print("done")
}