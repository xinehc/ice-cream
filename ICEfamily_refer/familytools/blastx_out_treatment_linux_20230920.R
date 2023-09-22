args <- commandArgs(TRUE)

str <- read.delim(file=args[4],sep="\t",header=T,stringsAsFactors = F)

############################
### give filename  ####  change every time   ##
############################
filename <- args[1]
idcut <- as.numeric(args[2])
lencut <- as.numeric(args[3])
evacut <- as.numeric(args[5])
folder <- args[6]

### add "length" into filename dataframe
blastx <- read.delim(file=filename,sep="\t",header=F,stringsAsFactors = F)
colnames(blastx)[1:3]<- c("query","subject","identity")
colnames(blastx)[7:8]<- c("querystart","queryend")
colnames(blastx)[11:14]<- c("evalue","bitscore","slen","qlen")
blastx <- blastx[order(blastx$query,blastx$evalue,-blastx$bitscore),]
blastx <- blastx[!duplicated(blastx$query),]
blastx <- blastx[which(blastx$identity >= idcut & blastx$evalue <= evacut ),]
blastx$alignlength<-abs(blastx$queryend-blastx$querystart)+1
blastx$alignqueryratio <- 100*blastx$alignlength/blastx$qlen
blastx$alignsubjectratio <- 100*blastx$alignlength/blastx$slen
blastx <- blastx[which(blastx$alignqueryratio >= lencut & blastx$alignsubjectratio >= lencut ),]

blastx_classify <- merge(blastx,str,by.x="subject",by.y="SARG.Seq.ID",all.x=T)
write.table(blastx_classify,file=paste(folder,"/extracted_classification_",basename(filename),".txt",sep=''),sep="\t",quote=F,row.names=F,col.names=T)

