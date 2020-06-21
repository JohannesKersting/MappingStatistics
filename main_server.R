library(Rsamtools)
library(data.table)
library(GenomicFeatures)
library(GenomicRanges)
library(tidyr)
library(GenomicAlignments)
library(pryr)



#functions for io tasks
source("functions.R")
start_time <- Sys.time()

#get path variables from config file
config <- readConfig("config.txt")
bam_path <- config$bam_path
fastq_path <- config$fastq_path
gtf_path <- config$gtf_path
genomic_reads_data <- config$genomic_reads
output_folder <- config$output_folder

#path variables
#bam_path <- "/nfs/scratch/AS_benchmarking/Task1/error_0/contextmap2-output/subsample_01.bam"
#fastq_path <- "/nfs/scratch/AS_benchmarking/Task1/error_0/subsample_01_1.fastq"
#gtf_path <- "/nfs/scratch/AS_benchmarking/Task1/error_0/splicing_variants.gtf"
#genomic_reads_data <- "error_0_genomic_reads.RData"
#output_folder <- "error_0_contextmap2_subsample_01"

if(!file.exists(genomic_reads_data)){
  #reading input files
  print("Generating genomic_reads...")
  reads <- fastqGRanges(fastq_path)
  transcripts <- readTranscripts(gtf_path)
  
  
  
  #correct false ranges
  print("Correcting ranges...")
  old <- Sys.time()
  
  read_width <- width(reads[Position(function(x) x==FALSE, reads$range_needs_check)])
  print(paste0("Read width: ",read_width))
  
  transcript_ids_for_range_parsing_1 <- as.character(seqnames(reads[reads$mate==1&reads$range_needs_check]))
  transcript_ids_for_range_parsing_2 <- as.character(seqnames(reads[reads$mate==2&reads$range_needs_check]))
  
  transcipts_for_range_parsing_1 <- transcripts[transcript_ids_for_range_parsing_1]
  transcipts_for_range_parsing_2 <- transcripts[transcript_ids_for_range_parsing_2]
  
  width_per_transcript_1 <- sum(width(transcipts_for_range_parsing_1))
  width_per_transcript_2 <- sum(width(transcipts_for_range_parsing_2))
  
  corrected_ranges_1 <- IRanges(start=1,end=pmin(width_per_transcript_1,read_width))
  corrected_ranges_2 <- IRanges(start=pmax(1,width_per_transcript_2-read_width+1),end=width_per_transcript_2)
  
  ranges(reads[reads$mate==1&reads$range_needs_check])<-corrected_ranges_1
  ranges(reads[reads$mate==2&reads$range_needs_check])<-corrected_ranges_2
  
  print(paste0(sum(reads$range_needs_check)," ranges corrected."))
  
  print(Sys.time()-old)
  
  
  
  
  #transform transcriptomic coordinates into genomic coordinates
  print("Parsing read coordinates to genomic...")
  old <- Sys.time()
  
  mapping <- mapFromTranscripts(reads,transcripts,ignore.strand=F)
  
  
  print(paste0(length(mapping)," read coordinates parsed..."))
  print("Unparsed reads:")
  print(reads[-mapping$xHits])
  print(Sys.time()-old)
  print("mem:")
  print(mem_used())
  
  
  
  #get junctions
  print("Parsing junctions...")
  old <- Sys.time()
  
  mapping_list <- split(mapping,1:length(mapping))
  junctions <- intersect(mapping_list,transcripts[mapping$transcriptsHits],ignore.strand=T)
  
  print(Sys.time()-old)
  print("mem:")
  print(mem_used())
  
  #generate genomic reads
  print("Generate genomic reads...")
  old <- Sys.time()
  
  genomic_reads <- GRanges(seqnames = seqnames(mapping),
                           ranges = ranges(mapping),
                           strand = strand(mapping),
                           mate = reads$mate
  )
  names(genomic_reads)<-paste(reads$read_id,reads$mate,sep="_")
  genomic_reads$junctions <- ranges(junctions)
  
  print(Sys.time()-old)
  print("mem:")
  print(mem_used())
  
  #generate genomic reads
  print(paste0("Saving genomic reads as ", genomic_reads_data))
  old <- Sys.time()
  
  save(genomic_reads, file = genomic_reads_data)
  
  print(Sys.time()-old)
  print("mem:")
  print(mem_used())
}else{
  
  print(paste0("Genomic reads are already computed. Loading data ",genomic_reads_data))
  old <- Sys.time()
  load(genomic_reads_data)
  print(paste0(length(genomic_reads)," reads loaded!"))
  print(Sys.time()-old)
  print("mem:")
  print(mem_used())
  
}

print(paste0("Reading bam file ",bam_path))
bam <- readBamFile(bam_path)

print("Preparing lookup table...")
old <- Sys.time()

genomic_reads$index <- 1:length(genomic_reads)
lookup <- data.table(read_id=names(bam),short_id=bam$id,mate=bam$mate,bam_index=1:length(bam),read_index=genomic_reads[names(bam)]$index)
genomic_reads$index<-NULL

print(Sys.time()-old)
print("mem:")
print(mem_used())



print("Check chromosomes...")
old <- Sys.time()
#check if chromosom is correct
lookup[,correct_chrom:=(
  as.character(seqnames(genomic_reads[read_index]))==
    as.character(seqnames(bam[bam_index]))
)]
print(Sys.time()-old)

print("Check start pos...")
old <- Sys.time()
#check if start pos is correct
lookup[,correct_start:=(
  correct_chrom&(
    start(genomic_reads[read_index])==
      start( bam[bam_index])
  )
)]
print(Sys.time()-old)

print("Check end pos...")
old <- Sys.time()
#check if end pos is correct
lookup[,correct_end:=(
  correct_chrom&(
    end(genomic_reads[read_index])==
      end( bam[bam_index])
  )
)]
print(Sys.time()-old)

print("Check range...")
old <- Sys.time()
lookup[,correct_range:=(
  correct_chrom&(
    ranges(genomic_reads[read_index])==
      ranges( bam[bam_index])
  )
)]
print(Sys.time()-old)


#number of junctions new
print("Counting junctions...")
old <- Sys.time()
names(bam$junctions)<-1:length(bam)


genomic_reads_junctions <- as.data.table(unlist(genomic_reads$junctions))
bam_junctions <- as.data.table(unlist(bam$junctions))

genomic_reads_junctions[,names:=as.integer(names)]
bam_junctions[,names:=as.integer(names)]

genomic_reads_junctions[,id:=names(genomic_reads[names])]
bam_junctions[,id:=names(bam[names])]

lookup[,read_n_junc:=genomic_reads_junctions[,.N,by=names]$N[read_index]]
lookup[,bam_n_junc:=bam_junctions[,.N,by=names]$N[bam_index]]
print(Sys.time()-old)


#number of correct junctions
old <- Sys.time()
print("Checking number of correct junctions..")
merged <- merge(genomic_reads_junctions,bam_junctions,by=c("start","end","id"))
merged <- merged[lookup$correct_chrom[names.y]]
merged_per_read <- merged[,.N,by=names.y]
lookup$correct_n_junc<-0
lookup[merged_per_read$names.y]$correct_n_junc<-merged_per_read$N
lookup[,false_n_junc:=bam_n_junc-correct_n_junc]
print(Sys.time()-old)



print("Making summary...")
#count for statistics
old <- Sys.time()

counts <- data.table()
counts[,n_reads:= length(genomic_reads)]
counts[,n_proper_maps:=length(bam)]
counts[,n_correct_chrom:= sum(lookup$correct_chrom)]
counts[,n_correct_start:= sum(lookup$correct_start)]
counts[,n_correct_end:=sum(lookup$correct_end)]
counts[,n_correct_range:= sum(lookup$correct_range)]
counts[,n_all_junctions_correct:= sum(lookup$correct_n_junc==lookup$read_n_junc)]
counts[,n_junctions:= nrow(genomic_reads_junctions)]
counts[,n_correct_junctions:= sum(lookup$correct_n_junc)]
counts[,n_false_junctions:=sum(lookup$false_n_junc)]

print(Sys.time()-old)
print("Merging pairs...")
old <- Sys.time()

lookup_pairs <- merge(lookup[mate==1,-c("mate","read_id")],lookup[mate==2,-c("mate","read_id")],by="short_id",suffixes=c("_mate_1","_mate_2"))

print(Sys.time()-old)
print("Pair summary...")
old <- Sys.time()
counts[,n_pairs:= length(genomic_reads)/2]
counts[,n_proper_pairs:=nrow(lookup_pairs)]
counts[,n_correct_chrom_pairs:=sum(lookup_pairs$correct_chrom_mate_1&lookup_pairs$correct_chrom_mate_2)]
counts[,n_correct_range_pairs:=sum(lookup_pairs$correct_range_mate_1&lookup_pairs$correct_range_mate_2)]
counts[,n_all_junctions_correct_pairs:=sum((lookup_pairs$correct_n_junc_mate_1==lookup_pairs$read_n_junc_mate_1)&(lookup_pairs$correct_n_junc_mate_2==lookup_pairs$read_n_junc_mate_2))]
counts

print(Sys.time()-old)
old <- Sys.time()

if (!file.exists(output_folder)){
  dir.create(file.path(getwd(),output_folder))
}
print(paste0("Saving results to: ",file.path(getwd(),output_folder)))
save(bam, file = file.path(getwd(),output_folder,"bam_reads.RData"))
fwrite(counts, file.path(getwd(),output_folder,"summary.csv"))
fwrite(lookup, file.path(getwd(),output_folder,"full_info.csv"))
print(Sys.time()-old)


total_time <- Sys.time()-start_time
print("Total time:")
print(total_time)



