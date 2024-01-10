#' Prepare data for MAPpoly
#' Function takes the output from compare_probes() and prepares it for input to use in MAPpoly
#' Function reads in annotation files and matches the genome position of Rosa chinensis genome (Hibrand-Saint Oyant et al., 2018)
#' NOTE: must be connected to internet to read in annotation files
#'
#' @param input_file is a .csv file that has the dosages and Axiom marker names. Most likely an output file of compare_probes() function
#' @param genome specifies which genome order to use "saintoyant" or "raymond" (default is "saintoyant")
#' @param p1 define parent 1 name
#' @param p2 define parent 2 name
#' @return a .csv file that has the added columns of genome position and LG of markers. Will contain markers with no known location at bottom of file.
#'
#' @export fitPoly_to_MAPpoly2_input
#'
#'
#'
fitPoly_to_MAPpoly2_input<-function(input_file, genome=NULL, p1, p2){
  
  
  input_file_name = gsub(".csv","",input_file)
  file = read.csv(input_file)
  
  if(is.null(genome)){
    genomic_pos = read.csv("https://raw.githubusercontent.com/jeekinlau/RoseArrayTools/master/docs/perfect_fits_saintoyant.csv", header = T)
  }else if(genome=="saintoyant"){
    genomic_pos = read.csv("https://raw.githubusercontent.com/jeekinlau/RoseArrayTools/master/docs/perfect_fits_saintoyant.csv", header = T)
  }else if(genome=="raymond"){
    genomic_pos = read.csv("https://raw.githubusercontent.com/jeekinlau/RoseArrayTools/master/docs/perfect_fits_raymond.csv", header = T)
  }else print("please select genome")
  
  
  temp_names = names(file)
  temp_names[1] = "Marker"
  
  colnames(file) = temp_names
  
  
  merged = merge(file,genomic_pos[,c(1,2,14,15,16)], by="Marker", all=T)
  
  p1_geno = merged[,which(colnames(merged)==p1)]
  p2_geno = merged[,which(colnames(merged)==p2)]
  
  merged_progeny = merged[,-which(colnames(merged)%in%c(p1,p2,"Marker","Chrom","snp_pos","Allele_A","Allele_B"))]
  
  final = cbind(merged$Marker, p1_geno, p2_geno, merged$LG, merged$Position, merged_progeny)
  final_colnames <- colnames(final)
  final_colnames[1:7] = c("Marker", "P1", "P2", "chrom", "genome_position","ref","alt")
  colnames(final) = final_colnames
  
  final = final[order(final$LG, final$Position),]
  final = final[-which(is.na(final$P1)),]
  final = final[-which(is.na(final$P2)),]
  
  write.csv(final,paste0(input_file_name,"_mappoly2_ready.csv"),row.names = F)
  print("done")
}
