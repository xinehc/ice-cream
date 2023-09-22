###################################################################################################
##################################### definition ##################################################
###################################################################################################
args <- commandArgs(trailingOnly=FALSE)
scriptPath<-dirname(sub("--file=","",args[grep("--file=",args)]))
library("reshape2")
options(stringsAsFactors = FALSE)
classify_ice_family_yxl<-function(familyname,orfname=c("orf1.txt","orf2.txt","orf3.txt"),idealORFcount,neighor,startORFname,secondstartORFname,endORFnamme,secondendORFname,filename_order,j,folder){
  conju<-orfname
  m=idealORFcount
  hmmout_cut_filter<-hmmout_cut[which(hmmout_cut$model_name %in% conju),]
  l=nrow(hmmout_cut_filter)
  k=0
  for (i in conju){ if( i %in% hmmout_cut$model_name ){k=k+1}}
  if (icefinder_result[j,"end"]!=1e+07){
    if ( m-k<=1 ){
      w<<-w+1
      write.table(hmmout_cut_filter,file=paste(folder,"/",filename_order,"_classification_summary_details",w,".txt",sep=""),sep="\t",quote=F,row.names=F,col.names=T)
      icefinder_result[j,"ice_family"]<<-paste(icefinder_result[j,"ice_family"],familyname,sep="_")
      icefinder_result[j,"ice_family_conju_ORF_count"]<<-paste(icefinder_result[j,"ice_family_conju_ORF_count"],m,sep="_")
      icefinder_result[j,"detected_conju_ORF_count"]<<-paste(icefinder_result[j,"detected_conju_ORF_count"],l,sep="_")
      icefinder_result[j,"count"]<<-paste(icefinder_result[j,"count"],floor(l/m),sep="_")
      icefinder_result[j,"details_file_number"]<<-paste(icefinder_result[j,"details_file_number"],w,sep="_")
      } 
  } else {
      if (k==m){
        w<<-w+1
        rownames(hmmout_cut_filter)<-1:nrow(hmmout_cut_filter)
        write.table(hmmout_cut_filter,file=paste(folder,"/",filename_order,"_classification_summary_details",w,".txt",sep=""),sep="\t",quote=F,row.names=F,col.names=T)
        tmp_start<-hmmout_cut_filter[which(hmmout_cut_filter$model_name==startORFname),]
        tmp_end<-hmmout_cut_filter[which(hmmout_cut_filter$model_name==endORFnamme),]
        tmp_start2<-data.frame()
        tmp_end2<-data.frame()
        for (z in 1:nrow(tmp_start)){
          if (TRUE %in% (c((as.numeric(rownames(tmp_start[z,]))+m-2):(as.numeric(rownames(tmp_start[z,]))+m+1)) %in% rownames(tmp_end))|
              TRUE %in% (c((as.numeric(rownames(tmp_start[z,]))-m-2):(as.numeric(rownames(tmp_start[z,]))-m+1)) %in% rownames(tmp_end))){
            tmp_start2<-rbind(tmp_start2,tmp_start[z,])}}
        for (z in 1:nrow(tmp_end)){
          if (TRUE %in% (c((as.numeric(rownames(tmp_end[z,]))+m-2):(as.numeric(rownames(tmp_end[z,]))+m+1)) %in% rownames(tmp_start))|
              TRUE %in% (c((as.numeric(rownames(tmp_end[z,]))-m-2):(as.numeric(rownames(tmp_end[z,]))-m+1)) %in% rownames(tmp_start))){
            tmp_end2<-rbind(tmp_end2,tmp_end[z,])}}
        q=nrow(icefinder_result)
     if (nrow(tmp_end2) !=0 & nrow(tmp_start2) !=0 ){
      for ( i in 1:min(nrow(tmp_end2),nrow(tmp_start2),floor(l/m))){
          icefinder_result[q+i,"V1"]<<-filename_order
          icefinder_result[q+i,"ice_family"]<<-familyname
          icefinder_result[q+i,"ice_family_conju_ORF_count"]<<-m
          icefinder_result[q+i,"detected_conju_ORF_count"]<<-l
          icefinder_result[q+i,"count"]<<-floor(l/m)
          icefinder_result[q+i,"details_file_number"]<<-w
        if (neighor=="reverse"){
          tmp_end<- tmp_end2[i,]
          if (c(hmmout_cut_filter[which(rownames(hmmout_cut_filter)==(as.numeric(rownames(tmp_end))-1)),"model_name"]==secondendORFname,"FALSE")[1]){
          tmp_start<-tmp_start2
          tmp_start$dis<-tmp_end$start-tmp_start$end
          if (nrow(tmp_start[which(tmp_start$dis > 0),])!=0) {
              tmp_start<-tmp_start[which(tmp_start$dis > 0),]
              tmp_start<-tmp_start[which(tmp_start$dis==min(tmp_start$dis)),]
              icefinder_result[q+i,"start"]<<-tmp_start[,"start"]
              icefinder_result[q+i,"end"]<<-tmp_end[,"end"] 
          }else if (c(hmmout_cut_filter[which(rownames(hmmout_cut_filter)==(as.numeric(rownames(tmp_end))+1)),"model_name"]==secondendORFname,"FALSE")[1]){
          tmp_start$dis<-tmp_start$start-tmp_end$end
          if (nrow(tmp_start[which(tmp_start$dis > 0),])!=0) {
              tmp_start<-tmp_start[which(tmp_start$dis > 0),]
              tmp_start<-tmp_start[which(tmp_start$dis==min(tmp_start$dis)),]
              icefinder_result[q+i,"start"]<<-tmp_end[,"start"]
              icefinder_result[q+i,"end"]<<-tmp_start[,"end"] 
          }}
         } else if (c(hmmout_cut_filter[which(rownames(hmmout_cut_filter)==(as.numeric(rownames(tmp_end))+1)),"model_name"]==secondendORFname,"FALSE")[1]){
          tmp_start$dis<-tmp_start$start-tmp_end$end
          if (nrow(tmp_start[which(tmp_start$dis > 0),])!=0) {
              tmp_start<-tmp_start[which(tmp_start$dis > 0),]
              tmp_start<-tmp_start[which(tmp_start$dis==min(tmp_start$dis)),]
              icefinder_result[q+i,"start"]<<-tmp_end[,"start"]
              icefinder_result[q+i,"end"]<<-tmp_start[,"end"] 
          }}
        } else if (neighor=="forward"){
          tmp_start<- tmp_start2[i,]
          tmp_end<-tmp_end2
          if (c(hmmout_cut_filter[which(rownames(hmmout_cut_filter)==(as.numeric(rownames(tmp_start))-1)),"model_name"]==secondstartORFname,"FALSE")[1]){
          tmp_end$dis<-tmp_start$start-tmp_end$end
          if (nrow(tmp_end[which(tmp_end$dis > 0),])!=0) {
              tmp_end<-tmp_end[which(tmp_end$dis > 0),]
              tmp_end<-tmp_end[which(tmp_end$dis==min(tmp_end$dis) ),]
              icefinder_result[q+i,"start"]<<-tmp_end[,"start"]
              icefinder_result[q+i,"end"]<<-tmp_start[,"end"] 
          }else if (c(hmmout_cut_filter[which(rownames(hmmout_cut_filter)==(as.numeric(rownames(tmp_start))+1)),"model_name"]==secondstartORFname,"FALSE")[1]){
              tmp_end$dis<-tmp_end$start-tmp_start$end
              if (nrow(tmp_end[which(tmp_end$dis > 0),])!=0) {
                  tmp_end<-tmp_end[which(tmp_end$dis > 0),]
                  tmp_end<-tmp_end[which(tmp_end$dis==min(tmp_end$dis)),]
                  icefinder_result[q+i,"start"]<<-tmp_start[,"start"]
                  icefinder_result[q+i,"end"]<<-tmp_end[,"end"] 
          }}
          } else if (c(hmmout_cut_filter[which(rownames(hmmout_cut_filter)==(as.numeric(rownames(tmp_start))+1)),"model_name"]==secondstartORFname,"FALSE")[1]){
          tmp_end$dis<-tmp_end$start-tmp_start$end
          if (nrow(tmp_end[which(tmp_end$dis > 0),])!=0) {
              tmp_end<-tmp_end[which(tmp_end$dis > 0),]
              tmp_end<-tmp_end[which(tmp_end$dis==min(tmp_end$dis)),]
              icefinder_result[q+i,"start"]<<-tmp_start[,"start"]
              icefinder_result[q+i,"end"]<<-tmp_end[,"end"] 
          }}
        }
          }
     }
      }
  }
}
#library(rlist)
list_ice_family<-readRDS(file.path(scriptPath,"list_ice_family.rds"))
filename<-args[6]
evalue<-as.numeric(args[7])
icefinder_result<-read.delim(file=args[8],sep="\t",header=F)
icefinder_result<-cbind(icefinder_result,colsplit(icefinder_result$V6, "\\.", c("start","N","end")))
superF<-args[9]
tempF<-args[10]


###################################################################################################
##################################### programe ##################################################
###################################################################################################
if(file.size(filename) >0){
no_col <- max(count.fields(filename, sep = "\t"))
hmmout<-read.delim(file=filename,sep="\t",heade=F,fill=TRUE,col.names=1:no_col)
names(hmmout)[1:6]<-c("model_name","model_accession","query_name","query_accession","E_value","score")
names(hmmout)[19]<-"description of target"
hmmout<-hmmout[which(hmmout$E_value<=evalue),]
hmmout<-hmmout[order(hmmout$query_name,hmmout$E_value,-hmmout$score),]
hmmout<-hmmout[!duplicated(hmmout$query_name),]
spl <- strsplit(as.character(hmmout$query_name), "_")
hmmout$query_start_end<-sapply(lapply(spl, tail, 2), paste, collapse="_")
hmmout<-cbind(hmmout,colsplit(hmmout$query_start_end, "\\_", c("start","end")))
hmmout$start<-as.numeric(hmmout$start)
hmmout$end<-as.numeric(hmmout$end)
hmmout$genbank_accession<-sapply(lapply(spl, head, 3), paste, collapse="_")
hmmout<-hmmout[order(hmmout$start),]
hmmout<-hmmout[,c("start","end","model_name","model_accession","E_value","score","description of target","genbank_accession")]
filename_order<-gsub(".gbk.fa2.faa.merge.all.hmmscan.out","",basename(filename))
write.table(hmmout,file=paste(tempF,"/",filename_order,"_dereplicated.txt",sep=""),quote=F,col.names=T,row.names=F,sep="\t")
icefinder_result<-icefinder_result[which(icefinder_result$V1==filename_order),]
icefinder_result[nrow(icefinder_result)+1,c("V1","start","end")]<-c(filename_order,as.numeric(0),as.numeric(10000000))
icefinder_result$start<-as.numeric(icefinder_result$start)
icefinder_result$end<-as.numeric(icefinder_result$end)

w=0
for ( j in 1:nrow(icefinder_result)){
icefinder_result[j,"ice_family"]<-""
icefinder_result[j,"ice_family_conju_ORF_count"]<-""
icefinder_result[j,"detected_conju_ORF_count"]<-""
icefinder_result[j,"count"]<-""
icefinder_result[j,"miss_ORF"]<-""
icefinder_result[j,"details_file_number"]<-""
tmp_start<-icefinder_result[j,"start"]-50000
tmp_end<-icefinder_result[j,"end"]+50000
hmmout_cut<-hmmout[which(hmmout$start>=tmp_start & hmmout$end <=tmp_end),]
hmmout_cut$model_name<-gsub("traC.protein.18","TraC",hmmout_cut$model_name)
hmmout_cut$model_name<-gsub("traB.protein.12","TraB",hmmout_cut$model_name)
hmmout_cut$model_name<-gsub("traA.protein.5","TraA",hmmout_cut$model_name)
hmmout_cut$model_name<-gsub("traW.protein.14","TraW",hmmout_cut$model_name)
hmmout_cut$model_name<-gsub("traD.protein.15","TraD",hmmout_cut$model_name)
hmmout_cut$model_name<-gsub("TraI_2","TraI",hmmout_cut$model_name)
hmmout_cut$model_name<-gsub("TraG_N","TraG",hmmout_cut$model_name)
hmmout_cut$model_name<-gsub("DUF4400","TraJ",hmmout_cut$model_name)
hmmout_cut$model_name<-gsub("TrwB_AAD_bind","TraD",hmmout_cut$model_name)
newfolder<-paste(superF,"/",filename_order,"_ICEfamily",sep="")

for (x in 1:length(list_ice_family)){
classify_ice_family_yxl(list_ice_family[[x]][["familyname"]],list_ice_family[[x]][["orfname"]],list_ice_family[[x]][["idealORFcount"]],list_ice_family[[x]][["neighor"]],list_ice_family[[x]][["startORFname"]],list_ice_family[[x]][["secondstartORFname"]],list_ice_family[[x]][["endORFnamme"]],list_ice_family[[x]][["secondendORFname"]],filename_order,j,newfolder)
}
}
icefinder_result$ice_family<-gsub("^_","",icefinder_result$ice_family)
icefinder_result_modifi<-data.frame()
tmp_yxl<-icefinder_result[is.na(icefinder_result$V5),]
  tmp <- icefinder_result[!is.na(icefinder_result$V5),]
icefinder_result_modifi<-rbind(icefinder_result_modifi,tmp)
for (j in 1:nrow(tmp_yxl)){
  tmp_ice_family<-tmp_yxl[j,"ice_family"]
  if (tmp_ice_family !="" & !is.na(tmp_yxl[j,"start"])){
      if (nrow(tmp[grepl(tmp_ice_family,tmp$ice_family),]) ==0 ){
          icefinder_result_modifi<-rbind(icefinder_result_modifi,tmp_yxl[j,])
      } else if (tmp_yxl[j,"start"]<tmp[grepl(tmp_ice_family,tmp$ice_family),"start"][1] | 
                 tmp_yxl[j,"end"] > tmp[grepl(tmp_ice_family,tmp$ice_family),"end"][1]){
                 if (nrow(tmp[grepl(tmp_ice_family,tmp$ice_family),]) > 1){
                     if( tmp_yxl[j,"start"]<tmp[grepl(tmp_ice_family,tmp$ice_family),"start"][2] | 
                         tmp_yxl[j,"end"] > tmp[grepl(tmp_ice_family,tmp$ice_family),"end"][2]){
                         if (nrow(tmp[grepl(tmp_ice_family,tmp$ice_family),]) > 2){
                             if( tmp_yxl[j,"start"]<tmp[grepl(tmp_ice_family,tmp$ice_family),"start"][3] | 
                                 tmp_yxl[j,"end"] > tmp[grepl(tmp_ice_family,tmp$ice_family),"end"][3]){
                                 if (nrow(tmp[grepl(tmp_ice_family,tmp$ice_family),]) > 3){
                                    if( tmp_yxl[j,"start"]<tmp[grepl(tmp_ice_family,tmp$ice_family),"start"][4] | 
                                        tmp_yxl[j,"end"] > tmp[grepl(tmp_ice_family,tmp$ice_family),"end"][4]){
                                        if (nrow(tmp[grepl(tmp_ice_family,tmp$ice_family),]) > 4){
                                           if( tmp_yxl[j,"start"]<tmp[grepl(tmp_ice_family,tmp$ice_family),"start"][5] | 
                                               tmp_yxl[j,"end"] > tmp[grepl(tmp_ice_family,tmp$ice_family),"end"][5]){
                                               icefinder_result_modifi<-rbind(icefinder_result_modifi,tmp_yxl[j,])}
                                        } else {icefinder_result_modifi<-rbind(icefinder_result_modifi,tmp_yxl[j,])
                                        }}
                                 }else{icefinder_result_modifi<-rbind(icefinder_result_modifi,tmp_yxl[j,])
                                 }}
                         }else{icefinder_result_modifi<-rbind(icefinder_result_modifi,tmp_yxl[j,])
                         }}
                 }else{icefinder_result_modifi<-rbind(icefinder_result_modifi,tmp_yxl[j,])
                 }
      }
  }
}
names(icefinder_result_modifi)[1:16]<-c("id","genome_info","genome_length","icefinder_folder","ICE_description","ICE_start..end","ICE_length","oriT","GC","Genome_GC","|Delta_GC|","ARG","VF","ICE_start","V15","ICE_end")
icefinder_result_modifi2<-icefinder_result_modifi[,c("id","ice_family","ice_family_conju_ORF_count","detected_conju_ORF_count","count",	"miss_ORF","details_file_number","genome_info","genome_length","icefinder_folder","ICE_description","ICE_start","ICE_end","ICE_length")]
write.table(icefinder_result,file=paste(tempF,"/",filename_order,"_classification_summary.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=T)
write.table(icefinder_result_modifi2,file=paste(superF,"/",filename_order,"_ICEfamily","/",filename_order,"_ICEfamily_result.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=T)
}
