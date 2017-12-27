#' Kaplan-Merier function
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @import survival
#' @import logging
#' @param time, survival time with same length with data
#' @param status, survival status with same length with data, 1 = Dead, 0 = Alive
#' @param group, vector indicating the group information
#' @param pvalcorrection, fdr, bonferroni
#' @param min_samples, minimum number of samples for KM analysis
#' @return results
#' @export

drkaplanmeier <- function(time, status, group,
                           min_samples=5, pvalcorrection="fdr") {

    print("kaplan-meier analysis")

    ## start calculating
    Result <- matrix(NA, dim(group)[1], 2, byrow=TRUE)
    colnames(Result) <- c("KMp", "N")
    rownames(Result) <- rownames(group)

    foreach( i1 = c(1:nrow(group)) ) %dopar% {

        keeplink <- which(!is.na(time) & !is.na(group[i1,]) &
                              !is.na(status) & time>= 0)
        Time.survival <- as.numeric(as.vector(time[keeplink]))
        cen.status <- ifelse(as.vector(status[keeplink]) == 1, 1, 0)
        keep.group <- as.numeric(as.vector(group[i1, keeplink]))

        if (length(keep.group) >= min_samples) {
            Result[i1, "N"] <- length(keep.group)
            cutgroup <- keep.group
            if(length(unique(cutgroup)) > 1) {
                test.data1 <- list(time     = Time.survival,
                                   status   = cen.status,
                                   group    = as.factor(cutgroup))
                tempresult <- try( model1 <- survdiff(Surv(time, status) ~ group,
                                                data=test.data1, na.action=na.exclude))

                if(length(grep("error", tempresult, ignore.case = T)) == 0){
                    Result[i1, "KMp"] <- 1-pchisq(model1$chisq,
                                                  df=length(levels(factor(cutgroup)))-1)

                }######### end of if no error
            }######### end of if more than two groups
        }######### end of if each marker has at least 20 no NA samples
    }######### end of each group

    Result <- data.frame(Result, check.names = F)
    Result$KMpFDR <- p.adjust(Result$KMp, method = pvalcorrection)
    return(Result)
}
