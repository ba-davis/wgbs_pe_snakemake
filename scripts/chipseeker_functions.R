

# Function to read and convert a single file to GRanges
read_and_convert_to_GRanges <- function(file_path) {
  df <- read.delim(file_path, header = TRUE)
  gr <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE,
                                  seqnames.field = "chr", start.field = "start",
                                  end.field = "end", ignore.strand = TRUE)
  return(gr)
}

new_annotate <- function(gr_obj, gr_name, inpath) {
  #my_comps <- names(gr_obj)
  
  # annotate via Chipseeker
  res <- annotatePeak(gr_obj,
    tssRegion = c(-3000, 200),
    TxDb = txdb,
    level = "transcript",
    assignGenomicAnnotation = TRUE,
    overlap="TSS",
    addFlankGeneInfo=TRUE
  )

  # anno pie
  pdf(paste0(inpath, "/", gr_name, ".anno_pie_chart.pdf"))
  plotAnnoPie(res)
  dev.off()

  # vennpie
  pdf(paste0(inpath, "/", gr_name, ".vennPie.pdf"))
  #par(mar=c(0.5,2,1,3))
  vennpie(res)
  dev.off()

  # convert to df
  foo <- as.data.frame(res)
  for(i in rownames(foo)) {
    if (grepl("Exon", foo[i,"annotation"], fixed = TRUE) || grepl("Intron", foo[i,"annotation"], fixed = TRUE) ) {
      #parse annotation column
      annotation = strsplit(foo[i,"annotation"], " ")
      #removes "," and "(" from string
      temp <- gsub('^.|.$', '', annotation[[1]][2])
      #get geneID and Transcript ID
      annotation <- strsplit(as.character(temp),"/")
      #check that annotation is consistent with geneId column
      if(annotation[[1]][2] != foo[i,"geneId"]) {
        flankTxIds <- strsplit(as.character(foo[i,"flank_txIds"]),";")
        flankGeneIds <- strsplit(as.character(foo[i,"flank_geneIds"]),";")
        flankGeneDistances <- strsplit(as.character(foo[i,"flank_gene_distances"]),";")
        #sanity check for identical number of flank genes
        temp = length(flankTxIds[[1]]) + length(flankGeneIds[[1]]) + length(flankGeneDistances[[1]])
        if(temp%%3 == 0) {
          if(length(flankGeneDistances[[1]]) == 1) {
            if(flankGeneDistances[[1]][1] == 0) {
              foo[i,"geneId"] = flankGeneIds[[1]][1]
              foo[i,"transcriptId"] = flankTxIds[[1]][1]
            }
          } else {
            #loop through list of flank genes
            for(j in 1:length(flankGeneDistances[[1]])) {
              #if peak is within a gene
              if(flankGeneDistances[[1]][j] == 0) {
                foo[i,"geneId"] = flankGeneIds[[1]][j]
                foo[i,"transcriptId"] = flankTxIds[[1]][j]
                break
              }
            }
          }
        }
      }
    }
  }
  return(foo)
}
