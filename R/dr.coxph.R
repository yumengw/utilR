#' coxph function
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @import survival
#' @import logging
#' @param data, data matrix used for analysis, each row is a data vector, or categorical values
#' @param time, survival time with same length with data
#' @param status, survival status with same length with data
#' @param pvalcorrection, fdr, bonferroni
#' @param min_samples, minimum number of samples for KM analysis
#' @return results
#' @export

drcoxph <- function(job_dir, data, time, status, kmtop=0.5, outformat = "json",
                     min_samples=5, min_groups=2, pvalcorrection = "fdr") {

    print("survival analysis")
    proc_input = coxph_preprocessing(data, time, status)
    data = proc_input$data
    time = proc_input$time
    status = proc_input$status

    ## start calculating
    Result <- matrix(NA, dim(data)[1], 7, byrow=TRUE)
    rownames(Result) <- rownames(data)
    colnames(Result) <- c("N", "coef", "Exp(coef)", "Coxp", "CoxpFDR", "LRp", "LRpFDR")

    KMdata <- data.frame()
    keeplink <-  which(!is.na(time) & !is.na(status) & (time) > 0)
    KMdata <- data.frame(matrix(NA, length(keeplink), 2+nrow(data),
                                dimnames=list(colnames(data)[keeplink], c("time", "status", paste("grp::", rownames(data), sep="")))),
                                check.names=F)
    KMdata$time <- time[keeplink]
    KMdata$status <- status[keeplink]

    #foreach( i1 = c(1:nrow(data)) ) %dopar% {
    for (i1 in c(1:nrow(data)) ) {
        if (suppressWarnings(sum(is.na(as.numeric(t(data[i1,]))))) == ncol(data)) {
            dataused <- factor(as.vector(as.character(t(data[i1,])))) ## categorical variable
            used_patients <- colnames(data)
            keeplink <- which(!is.na(time) & !is.na(status) &
                                  (time) > 0 & !is.na(dataused))
            Time.survival <- as.numeric(as.vector(time[keeplink]))
            cen.status <- ifelse(as.vector(status[keeplink]) == 1, 1, 0)
            dataused <- dataused[keeplink]
            used_patients <- used_patients[keeplink]
            if (length(dataused) >= min_groups) {
                Result[i1, "N"] <- length(dataused)
                test.data1 <- list(time     = Time.survival,
                                   status   = cen.status,
                                   group    = dataused)
                ## cox
                tempresult <- try( model1 <- coxph(Surv(time, status) ~ group,
                                                   data=test.data1, na.action = na.exclude), silent=TRUE)
                if(length(grep("error", tempresult, ignore.case = T)) == 0){
                    Result[i1, c("coef", "Exp(coef)", "Coxp")] <-
                        summary(model1)$coefficients[1,c("coef", "exp(coef)", "Pr(>|z|)" )]
                }
                ## KM
                KMdata[used_patients, paste("grp::", rownames(data)[i1], sep="")] <- as.character(dataused)
                tempresult <- try( model1 <- survdiff(Surv(time, status) ~ group,
                                                      data=test.data1, na.action=na.exclude), silent=TRUE)
                if(length(grep("error", tempresult, ignore.case = T)) == 0){
                    Result[i1, "LRp"] <- 1-pchisq(model1$chisq,
                                                  df=length(levels(factor(dataused)))-1)
                }
            }######### end of if each marker has at least min_samples no NA samples
        }else{
            dataused <- as.vector(as.numeric(t(data[i1,]))) ## numerix variable
            used_patients <- colnames(data)
            keeplink <- which(!is.na(time) & !is.na(status) &
                                (time) > 0 & !is.na(dataused))
            Time.survival <- as.numeric(as.vector(time[keeplink]))
            cen.status <- ifelse(as.vector(status[keeplink]) == 1, 1, 0)
            dataused <- dataused[keeplink]
            used_patients <- used_patients[keeplink]
            if (length(dataused) >= min_samples) {
                Result[i1, "N"] <- length(dataused)
                test.data1 <- list(time     = Time.survival,
                                   status   = cen.status,
                                   group    = dataused)
                ## cox
                tempresult <- try( model1 <- coxph(Surv(time, status) ~ group,
                                                   data=test.data1, na.action = na.exclude), silent=TRUE)
                if(length(grep("error", tempresult, ignore.case = T)) == 0){
                    Result[i1, c("coef", "Exp(coef)", "Coxp")] <-
                        summary(model1)$coefficients[1,c("coef", "exp(coef)", "Pr(>|z|)" )]
                }
                ## KM
                risk.group <- ifelse(dataused > quantile(dataused, 1-kmtop), "High", "Medium")
                risk.group <-  ifelse(dataused <= quantile(dataused, kmtop), "Low", risk.group)
                keeper <- which(risk.group!="Medium")
                risk.group <- risk.group[keeper]
                Time.survival.tmp <- Time.survival[keeper]
                cen.status.tmp <- cen.status[keeper]
                dataused.tmp <- dataused[keeper]
                cutgroup.tmp <- risk.group
                used_patients <- used_patients[keeper]
                KMdata[used_patients, paste("grp::", rownames(data)[i1], sep="")] <- risk.group
                if (length(unique(cutgroup.tmp)) > 1) {
                    test.data1 <- list(time     = Time.survival.tmp,
                                       status   = cen.status.tmp,
                                       group    = as.factor(cutgroup.tmp))
                    tempresult <- try( model1 <- survdiff(Surv(time, status) ~ group,
                                                          data=test.data1, na.action=na.exclude))
                    if (length(grep("error", tempresult, ignore.case = T)) == 0){
                        Result[i1, "LRp"] <- 1-pchisq(model1$chisq, df=length(levels(factor(cutgroup.tmp)))-1)
                    }
                }######### has more than one group
            }######### end of if each marker has at least min_samples no NA samples
        }######## is numeric or categorical valirable
    }######### end of each gene
    Result <- data.frame(Result, check.names = F)
    Result$CoxpFDR <- p.adjust(Result$Coxp, method = pvalcorrection)
    Result$LRpFDR <- p.adjust(Result$LRp, method = pvalcorrection)
    ## used for output
    proc_input_reformat = as.data.frame(cbind(t(data), time, status))
    Result <- cbind(gene = rownames(Result), Result)
    ## return file list
    return(list(rawdata = writeResults(proc_input_reformat, job_dir, outformat),
                table = writeResults(Result, job_dir, outformat),
                KMplot = writeResults(KMdata, job_dir, outformat)))
}
