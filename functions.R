library(Rsamtools)
library(data.table)
library(GenomicFeatures)
library(tidyr)
library(GenomicAlignments)
library(pryr)

# extracts the exons by transcript of a .gtf file
readTranscripts <- function(path){
  
  print("Reading gtf file...")
  old <- Sys.time()
  
  txdb <- makeTxDbFromGFF(path)
  transcripts <- exonsBy(txdb,"tx",use.names=T)
  print(paste0(length(transcripts)," transcripts read."))
  
  print("mem:")
  print(mem_used())
  print("time:")
  print(Sys.time()-old)
  
  return(transcripts)
}


# reads a fastq file and returns the headers as granges
fastqGRanges <- function(path){
  fastq <- readFastQHeadersAwk(path)
  
  print("Transforming to granges...")
  old <- Sys.time()
  
  granges <- GRanges(seqnames = Rle(fastq$transcript),
          ranges = IRanges(start=fastq$start,end=fastq$end),
          strand = Rle(strand(fastq$strand)),
          read_id = fastq$id,
          mate = fastq$mate,
          range_needs_check =fastq$range_needs_check
  )
  print(paste0(length(granges)," reads succesfully read..."))
  
  print("mem:")
  print(mem_used())
  print("time:")
  print(Sys.time()-old)
  
  return(granges)
}

#reads a fastq files and returns formatted headers
readFastQHeaders <- function(path){
  
  #read in file
  print("Reading fastq file...")
  old <- Sys.time()
  
  fastq <-  scan(path, what="", sep="\n")
  
  print("mem:")
  print(mem_used())
  print("time:")
  print(Sys.time()-old)
  
  
  #extract header rows
  print("Get headers...")
  old <- Sys.time()
  
  header <-data.table(header = fastq[grep("@",fastq)])
  
  print("mem:")
  print(mem_used())
  print("time:")
  print(Sys.time()-old)
  
  
  #format
  print("Formatting...")
  old <- Sys.time()
  
  header[,strand:="*"]
  header$header <- gsub("@","",header$header)
  header <- separate(header,col=header,into=c("id","info"),sep= "\\/")
  header <- separate(header,col=info,into=c("transcript","mate1","mate2"),sep=";")
  header <- melt(data=header,id.vars=c("id","transcript","strand"),variable.name="mate",value.name="range")
  header$mate <- as.integer(gsub("mate","",header$mate))
  header$range <- sapply(header$range,function(x){
    strsplit(x,":")[[1]][2]
  })
  header <- separate(header,col=range,into=c("start","end"),sep="-")
  header[,start:=as.integer(start)]
  header[,end:=as.integer(end)]
  
  #deal with na ends (not correct currently)
  header$range_needs_check <-F
  header[is.na(header$end)]$range_needs_check<-T 
  header[is.na(header$end)]$end<-76
  
  print("mem:")
  print(mem_used())
  print("time:")
  print(Sys.time()-old)
  
  
  return(header)
}

#reads a .bam file
readBamFile <- function(path){
  
  #read first mate
  print("reading first mate...")
  old <- Sys.time()
  param <- ScanBamParam( flag = scanBamFlag(isPaired = T,
                                            isProperPair = T,
                                            isUnmappedQuery = F,
                                            hasUnmappedMate = F,
                                            isSecondaryAlignment = F,
                                            isDuplicate = F,
                                            isSupplementaryAlignment = F,
                                            isFirstMateRead = T
  ))
  mate_1 <- readGAlignments(bam_path, param=param, use.names = T)
  
  print("mem:")
  print(mem_used())
  print("time:")
  print(Sys.time()-old)
  
  
  #read second mate
  print("reading second mate...")
  old <- Sys.time()
  param <- ScanBamParam( flag = scanBamFlag(isPaired = T,
                                            isProperPair = T,
                                            isUnmappedQuery = F,
                                            hasUnmappedMate = F,
                                            isSecondaryAlignment = F,
                                            isDuplicate = F,
                                            isSupplementaryAlignment = F,
                                            isSecondMateRead = T
  ))
  mate_2 <- readGAlignments(bam_path, param=param, use.names = T)
  
  print("mem:")
  print(mem_used())
  print("time:")
  print(Sys.time()-old)
  
  
  print("transforming to genomic ranges...")
  old <- Sys.time()
  #get the ranges
  bam <- list(mate_1=mate_1,mate_2=mate_2)
  
  bam_ranges_mate_1 <- granges(bam$mate_1)
  bam_ranges_mate_2 <- granges(bam$mate_2)
  bam_ranges_mate_1$mate <- 1
  bam_ranges_mate_2$mate <- 2
  
  print("mem:")
  print(mem_used())
  print("time:")
  print(Sys.time()-old)
  
  
  print("parsing junctions from cigar...")
  old <- Sys.time()
  #get junctions from cigar
  bam_ranges_mate_1$junctions <-  extractAlignmentRangesOnReference(cigar(bam$mate_1),pos=start(bam$mate_1))
  bam_ranges_mate_2$junctions <-  extractAlignmentRangesOnReference(cigar(bam$mate_2),pos=start(bam$mate_2))
  
  print("mem:")
  print(mem_used())
  print("time:")
  print(Sys.time()-old)
  
  
  # combine the mates into one granges object
  print("combining the two mates... ")
  old <- Sys.time()
  bam_ranges <- c(bam_ranges_mate_1,bam_ranges_mate_2)
  
  if(grepl("/", names(bam_ranges[1]), fixed = TRUE)){
    print("Parsing read names...")
    names(bam_ranges)<-unlist(strsplit(names(bam_ranges),"/",fixed=T))[ c(TRUE,FALSE) ]
  }
  bam_ranges$id <- names(bam_ranges)
  names(bam_ranges) <- paste(names(bam_ranges),bam_ranges$mate,sep="_")
  
  print("mem:")
  print(mem_used())
  print("time:")
  print(Sys.time()-old)
  
  print("finished reading the bam file! Number of reads:")
  print(length(bam_ranges))
  return(bam_ranges)
}

readFastQHeadersAwk <- function(path){
  
  awk.script <- "read_fastq.awk"
  
  print("Reading fastq file with awk...")
  old <- Sys.time()
  
  header <- fread(text=system2("awk",args=c("-f",awk.script,path),stdout=TRUE))
  
  print("mem:")
  print(mem_used())
  print("time:")
  print(Sys.time()-old)
  return(header)
}

readConfig <- function(path){
  config <- fread(path,header = F)
  path_variables <- as.list(config$V2)
  names(path_variables)<- gsub(":","",config$V1)
  print("Path variables: ")
  print(path_variables)
  return(path_variables)
}

