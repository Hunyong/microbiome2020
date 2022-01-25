calculate.criterion <- function(pval, n.signal, n.gene) {
    real.P <- pval[1:n.signal]
    real.N <- pval[(n.signal + 1):n.gene]

    TP <- sum(real.P < 0.05)
    FN <- length(real.P) - TP

    TN <- sum(real.N > 0.05)
    FP <- length(real.N) - TN

    accuracy <- (TP + TN) / n.gene
    precision <- TP / (TP + FP)
    sensitivity <- TP / (TP + FN)
    FDR <- FP / (TP + FP)
    out <- list()
    out[["accuracy"]] <- accuracy
    out[["precision"]] <- accuracy
    out[["sensitivity"]] <- accuracy
    out[["FDR"]] <- FDR
    return(out)
}

#result$criterion <- apply(result$pval, 1, calculate.criterion, n.signal = n.signal, n.gene = n.gene)
