suppressPackageStartupMessages(library(MAST))
suppressPackageStartupMessages(library(edgeR))

run_MASTcpm <- function(L) {
  message("MAST, CPM")
  session_info <- sessionInfo()
  timing <- system.time({
    stopifnot(all(names(L$condt) == colnames(L$count)))
    grp <- L$condt
    dge <- DGEList(counts = L$count)
    dge <- edgeR::calcNormFactors(dge)
    cpms <- edgeR::cpm(dge)
    sca <- FromMatrix(exprsArray = log2(cpms + 1), 
                      cData = data.frame(grp = grp))
    zlmdata <- zlm(~grp, sca)
    mast <- lrTest(zlmdata, "grp")
    lfcs <- getLogFC(zlmdata)
  })
  hist(mast[, "hurdle", "Pr(>Chisq)"], 50)
  
  list(session_info = session_info,
       timing = timing,
       res = mast,
       df = data.frame(pval = mast[, "hurdle", "Pr(>Chisq)"],
                       lfc = lfcs[, 'logFC'],
                       var_lfc = lfcs[, 'varLogFC'],
                       row.names = names(mast[, "hurdle", "Pr(>Chisq)"])))
}
