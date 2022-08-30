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
#' @export fitPoly_to_MAPpoly_input
#'
#'
#'
fitPoly_to_MAPpoly_input<-function(input_file, genome=NULL, p1, p2){


input_file_name = gsub(".csv","",input_file)
file = read.csv(input_file)

if(genome==NULL){
genomic_pos = read.table("docs/for_add_genomic_positions_saintoyant.txt", header = T, sep="\t")
}else if(genome=="saintoyant"){
    genomic_pos = read.table("docs/for_add_genomic_positions_saintoyant.txt", header = T, sep="\t")
}else if(genome=="raymond"){
    genomic_pos = read.table("docs/for_add_genomic_positions_raymond.txt", header = T, sep="\t")
}else print("please select genome")


temp_names = names(file)
temp_names[1] = "Marker"

colnames(file) = temp_names


merged = merge(file,genomic_pos, by="Marker", all=T)

p1_geno = merged[,which(colnames(merged)==p1)]
p2_geno = merged[,which(colnames(merged)==p2)]

merged_progeny = merged[,-which(colnames(merged)%in%c(p1,p2,"Marker","LG","Position"))]

final = cbind(merged$Marker, p1_geno, p2_geno, merged$LG, merged$Position, merged_progeny)
final_colnames <- colnames(final)
final_colnames[1:5] = c("Marker", "P1", "P2", "LG", "Position")
colnames(final) = final_colnames

final = final[order(final$LG, final$Position),]
final = final[-which(is.na(final$P1)),]
final = final[-which(is.na(final$P2)),]

write.csv(final,paste0(input_file_name,"_mappoly_ready.csv"),row.names = F)
print("done")
}
