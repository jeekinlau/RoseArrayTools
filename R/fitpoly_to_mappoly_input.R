
fitPoly_to_MAPpoly_input<-function(input_file, p1, p2){


input_file_name = gsub(".csv","",input_file)
file = read.csv(input_file)
genomic_pos = read.table("https://raw.githubusercontent.com/jeekinlau/RoseArrayTools/master/docs/for_add_genomic_positions_saintoyant.txt", header = T, sep="\t")

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
