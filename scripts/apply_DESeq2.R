suppressPackageStartupMessages(library(DESeq2))

run_DESeq2 <- function(L, lfcThreshold = 0.0) {
  message("DESeq2")
  session_info <- sessionInfo()
  timing <- system.time({
    dds <- DESeqDataSetFromMatrix(countData = round(L$count), 
                                  colData = data.frame(condition = L$condt), 
                                  design = ~condition)
    dds <- DESeq(dds)
    res <- results(dds, contrast = c("condition", levels(factor(L$condt))[1], 
                                     levels(factor(L$condt))[2]), alpha = 0.05, lfcThreshold = lfcThreshold)
  })
  
  plotDispEsts(dds)
  plotMA(res)
  summary(res)
  
  list(session_info = session_info,
       timing = timing,
       res = res,
       df = data.frame(pval = res$pvalue,
                        padj = res$padj,
                        lfc = res$log2FoldChange,
                        row.names = rownames(res))
       )
}


run_DESeq2_multibatch <- function(L, lfcThreshold = 0.0) {
  message("DESeq2 Multibatch")
  session_info <- sessionInfo()
  timing <- system.time({
    dds <- DESeqDataSetFromMatrix(countData = round(L$count), 
                                  colData = data.frame(condition = L$condt, batch=L$batch), 
                                  design = ~condition+batch)
    dds <- DESeq(dds)
    res <- results(dds, contrast = c("condition", levels(factor(L$condt))[1], 
                                     levels(factor(L$condt))[2]), alpha = 0.05, lfcThreshold = lfcThreshold)
  })
  
  plotDispEsts(dds)
  plotMA(res)
  summary(res)
  
  list(session_info = session_info,
       timing = timing,
       res = res,
       df = data.frame(pval = res$pvalue,
                        padj = res$padj,
                        lfc = res$log2FoldChange,
                        row.names = rownames(res))
       )
}