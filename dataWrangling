
#add sequence length to other metadata in R 

#read sequence length file from terminal
seqLength=read.table("~/Desktop/mar23_sequenceFiles/seqLength_noGap.txt")

#convert length entries to a new column 
seqLength_match=as.data.frame(matrix(seqLength$V1, ncol = 2, byrow = TRUE))
seqLength_match$V1 <- sub("\\|.*", "", seqLength_match$V1)

#rename columns 
colnames(seqLength_match)[1] <- "sampleID"
colnames(seqLength_match)[2] <- "sequenceLength"

#load metadata spreadsheet 
specimen_metadata=read.table("~/Desktop/mar23_sequenceFiles/specimen_metadata.txt",header=T,sep="\t")

#format the sample id to match the sequence length dataframe 
seqLength_match$sampleID <- sub(" ", "_", specimen_metadata$sampleID)
seqLength_match$sampleID <- sub("-", "", specimen_metadata$sampleID)
specimen_metadata$sampleID <- sub(" ", "_", specimen_metadata$sampleID)
specimen_metadata$sampleID <- sub("-", "", specimen_metadata$sampleID)

#merge 
specimenMetadata=merge(x = specimen_metadata, y = seqLength_match, by = "sampleID",all.x=TRUE,all.y=TRUE)

#export 
write.csv(specimenMetadata,file="specimenMetadata_withSeqlen.csv")
