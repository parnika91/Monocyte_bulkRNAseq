## Script to run featureCounts (Rsubread) on mapped, SortedByCoordinates BAM files
library(Rsubread)

sampleTable <- read.csv("../Data/samplelist.csv", header = F)

dir_bam <- "../Results/bam"
#filenames <- file.path(dir_bam, paste0(sampleTable[,1], "Aligned.sortedByCoord.out.bam"))
filenames <- paste(dir_bam, list.files(dir_bam, pattern = ".bam"), sep = "/")

# files <- c()
# for(i in 1:length(filenames))
#   if(file.exists(filenames[i]))
#     files <- c(files, filenames[i])

gene_counts <- featureCounts(filenames,
              
              # annotation
              #annot.inbuilt="mm10",
              annot.ext="../Data/human.gtf",
              isGTFAnnotationFile=T,
              GTF.featureType="exon",
              GTF.attrType="gene_id",
              chrAliases=NULL,
              
              # level of summarization
              useMetaFeatures=TRUE,
              
              # overlap between reads and features
              allowMultiOverlap=FALSE,
              minOverlap=1,
              largestOverlap=FALSE,
              readExtension5=0,
              readExtension3=0,
              read2pos=NULL,
              
              # multi-mapping reads
              countMultiMappingReads=FALSE,
              fraction=FALSE,
              
              # read filtering
              minMQS=0,
              splitOnly=FALSE,
              nonSplitOnly=FALSE,
              primaryOnly=FALSE,
              ignoreDup=FALSE,
              
              # strandness
              strandSpecific=0,
              
              # exon-exon junctions
              juncCounts=FALSE,
              genome=NULL,
              
              # parameters specific to paired end reads
              isPairedEnd=TRUE,
              requireBothEndsMapped=FALSE,
              checkFragLength=FALSE,
              minFragLength=50,
              maxFragLength=600,
              countChimericFragments=TRUE,	
              autosort=TRUE,
              
              # miscellaneous
              nthreads=10,
              maxMOp=10)

saveRDS(gene_counts, "../Results/gene_counts/P1255_genecounts.rds")
write.csv2(gene_counts$counts, "../Results/gene_counts/P1255_genecounts.csv", quote = F)
