#' Compares both Probes on Axiom Arrays
#'
#' Takes *score.dat output from fitpoly and compares the probes and outputs a CSV with consensus csv of both probes.
#' This function needs internet access as it look for a reference file that is used to differentiate which two probes
#' belong to which SNPs. The reference file is located in /doc folder (if interested which file look in the function's code).
#' The input file is a file exported from fitpoly which has a extension of scores.dat and output will be two .csv files.
#' this is the compare_probes2 function (renamed and made other clearly defunct)
#'
#'
#' @param data_dat a text file with the end *score.dat that has the genotypic dosage call output from fitPoly
#'
#' @param progress if \code{TRUE} progress shown in console; if \code{FALSE}, no output produced default is TRUE
#'
#' @return Two CSV files. One has the consensus dosage call, other has description of the call. S = same, D = different, O = one, NA = not called
#'
#' @author Jeekin Lau, \email{jeekinlau@gmail.com}
#'
#' @importFrom data.table fread
#'
#'
#' @export compare_probes
#'
#'
#'
#'
#'
#'

compare_probes<-function(data_dat,progress=NULL){
  if(is.null(progress)){progress<-T}
  data_dat_name<-gsub(".dat","_",data_dat)

  calls<- as.matrix(fread(data_dat, select = c("marker","MarkerName","SampleName","geno")))

  if(progress==T){print("Done Importing .dat file")}

  ind_unique<-unique(calls[,3])
  num_ind<-as.numeric(length(ind_unique))
  ind<-calls[1:num_ind,3]
  num_markers<-nrow(calls)/num_ind

  markers<-matrix(,num_markers,1)
  colnames(markers)<-"Probes_ID"

  markers = as.matrix(unique(calls[,2]),num_markers,1)

  genocalls<-matrix(, num_markers, num_ind+1)
  colnames(genocalls)<-c("Probes_ID",ind)
  genocalls[,1]<-markers

  stringcall<-as.numeric(calls[,4])

  genocalls<-t(genocalls)
  genocalls[2:nrow(genocalls),1:ncol(genocalls)]<-stringcall
  genocalls<-t(genocalls)


  genocall_order<-as.matrix(read.csv("https://raw.githubusercontent.com/jeekinlau/RoseArrayTools/master/docs/array_snps_flanking_order.csv"))

  marker_col<-matrix(,nrow(genocalls),1)



  colnames(genocall_order)=c("Probes_ID","Affy.SNP.ID")
  temp=merge(genocalls,genocall_order,by="Probes_ID",all=T)
  genocalls2=cbind(temp[,ncol(temp)],temp[,1:(ncol(temp)-1)])


  genocalls3<-genocalls2[order(genocalls2[,1]),]

  genocalls<-genocalls3

  compare_probes<-array(,dim=c(nrow(genocalls)/2,ncol(genocalls),2))


  ####################################################################################################################

  compare_probes[,,1] = as.matrix(genocalls[seq(1,nrow(genocalls),2),])
  compare_probes[,,2] = as.matrix(genocalls[seq(2,nrow(genocalls),2),])


  compared_calls<-matrix(,nrow(compare_probes),ncol(compare_probes))

  compare_probes[is.na(compare_probes)]<-9


  same = which(compare_probes[,,1]==compare_probes[,,2]) #isolates the same calls and places them on new matrix
  compared_calls[same] = compare_probes[,,1][same]

  slice_1_single_probe = which(compare_probes[,,2]==9 & compare_probes[,,1]!=9)
  slice_2_single_probe = which(compare_probes[,,1]==9 & compare_probes[,,2]!=9)

  compared_calls[slice_1_single_probe] = compare_probes[,,1][slice_1_single_probe]
  compared_calls[slice_2_single_probe] = compare_probes[,,2][slice_2_single_probe]
  colnames(compared_calls)<-c("markers","probes",ind)
  NA_ = which(compared_calls==9)
  compared_calls[compared_calls==9]=NA
  compared_calls[,1]=compare_probes[,1,1]



  marker_stats<-matrix(,nrow(compare_probes),ncol(compare_probes))

  same = which(compare_probes[,,1]==compare_probes[,,2])
  different = which(compare_probes[,,1]!=compare_probes[,,2])
  one_s1 = which(compare_probes[,,2]==9 & compare_probes[,,1]!=9)
  one_s2 = which(compare_probes[,,1]==9 & compare_probes[,,2]!=9)

  marker_stats[same]="S"
  marker_stats[different]="D"
  marker_stats[one_s1]="O"
  marker_stats[one_s2]="O"
  marker_stats[NA_]=NA

  colnames(marker_stats)=c("markers","probes",ind)
  marker_stats[,1]=compared_calls[,1]
  if(progress==T){print("Done comparing probes")}

  write.csv(compared_calls[,-2], paste0(data_dat_name,"compared_calls_v2.csv"), row.names = F)
  write.csv(marker_stats[,-2], paste0(data_dat_name,"compared_calls_kind_counts_v2.csv"),row.names = F)


  print("FINISHED")
}


#
#
# call_specs<-function(kind_counts_file, select_genotypes=NULL){
#   if(is.null(select_genotypes)){
#     kinds_of_calls<-as.matrix(read.csv(kind_counts_file,header = T,row.names = 1))
#     percent_same<-sum(kinds_of_calls=="S",na.rm = T)/(nrow(kinds_of_calls)*(ncol(kinds_of_calls)))
#     percent_one<-sum(kinds_of_calls=="O",na.rm = T)/(nrow(kinds_of_calls)*(ncol(kinds_of_calls)))
#     percent_different<-sum(kinds_of_calls=="D",na.rm = T)/(nrow(kinds_of_calls)*(ncol(kinds_of_calls)))
#     percent_NA<-sum(is.na(kinds_of_calls))/(nrow(kinds_of_calls)*(ncol(kinds_of_calls)))
#
#     print(paste("percent______same",percent_same, round((percent_same*nrow(kinds_of_calls)),digits = 0)))
#     print(paste("percent_______one",percent_one, round((percent_one*nrow(kinds_of_calls)),digits = 0)))
#     print(paste("percent_different",percent_different, round((percent_different*nrow(kinds_of_calls)),digits=0)))
#     print(paste("percent________NA",percent_NA, round((percent_NA*nrow(kinds_of_calls)),digits = 0)))}
#
#
#  else{
#   genotypes<-read.csv(select_genotypes,header = F)
#   genotypes<-as.character(genotypes[,1])
#   genotypes<-gsub("-",".",genotypes)
#   kinds_of_calls<-as.matrix(read.csv(kind_counts_file,header = T,row.names = 1))
#   kinds_of_calls<-kinds_of_calls[,colnames(kinds_of_calls)%in%genotypes]
#   percent_same<-sum(kinds_of_calls=="S",na.rm = T)/(nrow(kinds_of_calls)*(ncol(kinds_of_calls)))
#   percent_one<-sum(kinds_of_calls=="O",na.rm = T)/(nrow(kinds_of_calls)*(ncol(kinds_of_calls)))
#   percent_different<-sum(kinds_of_calls=="D",na.rm = T)/(nrow(kinds_of_calls)*(ncol(kinds_of_calls)))
#   percent_NA<-sum(is.na(kinds_of_calls))/(nrow(kinds_of_calls)*(ncol(kinds_of_calls)))
#
#   print(paste("percent______same",percent_same, round((percent_same*nrow(kinds_of_calls)),digits = 0)))
#   print(paste("percent_______one",percent_one, round((percent_one*nrow(kinds_of_calls)),digits = 0)))
#   print(paste("percent_different",percent_different, round((percent_different*nrow(kinds_of_calls)),digits=0)))
#   print(paste("percent________NA",percent_NA, round((percent_NA*nrow(kinds_of_calls)),digits = 0)))}
# }
