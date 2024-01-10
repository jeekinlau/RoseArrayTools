#' Function takes the output from compare_probes() and prepares it for input to use in GWASpoly.
#' Function reads in annotation files and matches the genome position of Rosa chinensis genome (Hibrand-Saint Oyant et al., 2018)
#' NOTE: must be connected to internet to read in annotation files
#'
#' @param input_file is a .csv file that has the dosages and Axiom marker names. Most likely an output file of compare_probes() function
#'
#' @return a .csv file that has the added columns of genome position and LG of markers. Will contain markers with no known location at bottom of file.
#'
#' @export fitPoly_to_GWASpoly_input
#'
#'
#'
#'
#'
#'


fitPoly_to_GWASpoly_input <- function(input_file){
  input_file_name = gsub(".csv","",input_file)
  file = read.csv(input_file)
  genomic_pos = read.table("https://raw.githubusercontent.com/jeekinlau/RoseArrayTools/master/data/for_add_genomic_positions_saintoyant2.txt", header = T, sep="\t")

  temp_names = names(file)
  temp_names[1] = "Marker"
  colnames(file) = temp_names

  merged = merge(file,genomic_pos, by="Marker", all=T)

  merged_genotypes = merged[,-which(colnames(merged)%in%c("Marker","LG","Position"))]

  final = cbind(merged$Marker,  merged$LG, merged$Position, merged_genotypes)
  final_colnames <- colnames(final)
  final_colnames[1:3] = c("Marker", "Chrom", "Position")
  colnames(final) = final_colnames
  final = final[order(final$Chrom, final$Position),]
  write.csv(final,paste0(input_file_name,"_GWAS_ready.csv"),row.names = F)
  print("done")
}
