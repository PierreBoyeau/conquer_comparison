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

    summaryCond <- summary(zlmdata, doLRT=TRUE)
    summaryDt <- summaryCond$datatable
    fcHurdle <- merge(
                    summaryDt[contrast=='grp2' & component=='H',.(primerid, `Pr(>Chisq)`)],
                        #hurdle P values
                    summaryDt[contrast=='grp2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)],
                    by='primerid' #logFC coefficients
    )
    df = data.frame(pval = fcHurdle$"Pr(>Chisq)",
                lfc = fcHurdle$coef,
                row.names = fcHurdle$primerid)
    df <- df[order(as.numeric(row.names(df))),]

  })
  hist(mast[, "hurdle", "Pr(>Chisq)"], 50)
  
  list(session_info = session_info,
       timing = timing,
       res = mast,
       df = df)
}

run_MASTcpm_multibatch <- function(L) {
  message("MAST, CPM MULTI BATCHES")
  session_info <- sessionInfo()
  timing <- system.time({
    stopifnot(all(names(L$condt) == colnames(L$count)))
    grp <- L$condt
    dge <- DGEList(counts = L$count)
    dge <- edgeR::calcNormFactors(dge)
    cpms <- edgeR::cpm(dge)
    sca <- FromMatrix(exprsArray = log2(cpms + 1), 
                      cData = data.frame(grp = grp))
    zlmdata <- zlm(~L$condt + L$batch, sca)
    
    summaryCond <- summary(zlmdata, doLRT=TRUE)
    summaryDt <- summaryCond$datatable
    fcHurdle <- merge(
                    summaryDt[contrast=='grp2' & component=='H',.(primerid, `Pr(>Chisq)`)],
                        #hurdle P values
                    summaryDt[contrast=='grp2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)],
                    by='primerid' #logFC coefficients
    )
    df = data.frame(pval = fcHurdle$"Pr(>Chisq)",
                lfc = fcHurdle$coef,
                row.names = fcHurdle$primerid)
    df <- df[order(as.numeric(row.names(df))),]
    # mast <- lrTest(zlmdata, "grp")
    # lfcs <- getLogFC(zlmdata)
  })
  
  list(session_info = session_info,
       timing = timing,
       res = mast,
       df = df)
}
