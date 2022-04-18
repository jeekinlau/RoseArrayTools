
#' Function takes the output from compare_probes() and prepares it for input to use in polyOrigin
#' Function reads in annotation files and matches the genome position of Rosa chinensis genome (Hibrand-Saint Oyant et al., 2018)
#' NOTE: must be connected to internet to read in annotation files
#'
#' @param input_file is a .csv file that has the dosages and Axiom marker names. Most likely an output file of compare_probes() function
#'
#' @return a .csv file that has the added columns of genome position and LG of markers. Will contain markers with no known location at bottom of file.
#'
#' @export fitPoly_to_polyOrigin_input
#'
#'
#'


fitPoly_to_polyOrigin_input<-function(input_file){


  input_file_name = gsub(".csv","",input_file)
  file = read.csv(input_file, header = F)
  genomic_pos = read.table("https://raw.githubusercontent.com/jeekinlau/RoseArrayTools/master/docs/for_add_genomic_positions_saintoyant.txt", header = T, sep="\t")

  temp_names = file[1,]
  temp_names[1] = "Marker"

  file = file[-1,]
  colnames(file) = temp_names


  merged = merge(file,genomic_pos, by="Marker", all=T)

  merged_genos = merged[,-which(colnames(merged)%in%c("Marker","LG","Position"))]

  final = cbind(merged$Marker, merged$LG, merged$Position, merged_genos)
  final_colnames <- colnames(final)
  final_colnames[1:3] = c("marker",  "chromosome", "position")
  colnames(final) = final_colnames

  final = final[order(final$chromosome, final$position),]

  write.csv(final,paste0(input_file_name,"_polyOrigin_ready.csv"),row.names = F)
  print("done")
}
