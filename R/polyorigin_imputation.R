#'
#' this function will take the polyorigin postdoseprob.csv output which contains conditional probabilities and outputs
#'
#' @param input.file name of the file you want to convert from polyorigin
#' @param number.parents number or parental genotypes in the file
#' @param ploidy number as an integer
#' @param probability probability you want to use as a cutoff. recommend to use 0.9
#'
#'
#'
#'
#' @importFrom splitstackshape cSplit
#'
#' @export polyorigin_imputation
#'
#
#
#
#
polyorigin_imputation<-function(input.file, number.parents, ploidy, probability){


  number.parents<-number.parents
  data_original<-read.csv(input.file,header=TRUE,row.names =1)
  data<-data_original[,(number.parents+3):ncol(data_original)]
  file.name<-gsub(".csv","_",input.file)
  #softcode to make it easier to incorporate into a package function
  num.markers <- nrow(data)
  num.ind <- ncol(data)
  ploidy <- ploidy
  dosage.classes <- ploidy+1
  probability <- probability
  error=integer(0)
  #create 3d array to store data
  data.array <- array( NA, dim=c(num.markers, num.ind, dosage.classes))

  #install.packages("splitstackshape")
  #library(splitstackshape)
  data.split<-as.matrix(cSplit(data, splitCols = 1:ncol(data), "|"))


  #fill in array with data
  for(i in 1:num.markers){
    for(j in 1:num.ind){
      for(k in 1:dosage.classes){
        data.array[i,j,1:dosage.classes]<-as.vector(data.split[i,(dosage.classes*j-ploidy):(dosage.classes*j)])
      }}}

  imputed.matrix<-matrix(NA,num.markers,num.ind)


  for(l in 1:num.markers){
    for(m in 1:num.ind){
      imputed.matrix[l,m]<-ifelse(identical(error,which(data.array[l,m,1:dosage.classes]>=probability)),NA, which(data.array[l,m,1:dosage.classes]>=probability)-1)
    }}
  colnames(imputed.matrix)<-colnames(data)
  rownames(imputed.matrix)<-rownames(data)

  final<-cbind(data_original[,1:(number.parents+2)],imputed.matrix)
  write.csv(final,paste0(file.name,"imputed.csv"))
}
