input_file = "96_SWxBE_scores_compared_calls.csv"



file = read.csv(input_file)
p1= "Brite_Eyes_8"
p2= "Stormy_Weather_3"

p1_geno = file[,which(colnames(file)==p1)]
p2_geno = file[,which(colnames(file)==p2)]
