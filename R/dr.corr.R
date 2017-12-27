#' Correlation function
#' this function will calculate all correlations between data matrix 1 and 2 by row
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @import logging
#' @param datafiles contain file1, data matrix used for analysis, each column is a data vector
#' and file2, data matrix used for analysis, each column is a data vector
#' if file2 not given, calculate correlations within file1
#' @param method, pearson, spearman, kendall (default spearman)
#' @param alternative, two.sided, less, greater (default two.sided)
#' @param pvalcorrection, fdr, bonferroni
#' @param min_samples, minimum samples for correlation calculation
#' @return results
#' @export

drcorr <- function(datafiles, job_dir, outformat="json", method="spearman",
                   alternative="two.sided", min_samples = 5, pvalcorrection="fdr") {

    print("correlation analysis")
    fileSep <- .Platform$file.sep
    ## start calculating

    if (length(datafiles)==1) { #### within data1 correlation
        data1 = read.table(paste(job_dir, datafiles[1], sep=fileSep), sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F)
        if (length(grep('mRNA', colnames(data1))) > 0) {
            data1[data1<0] <- NA
            data1 <- log(data1 + 1e-6, 2)
        }
        Result <- data.frame(matrix(NA, dim(data1)[2]*(dim(data1)[2]-1)/2, 6, byrow=TRUE), check.names=F)
        all_data <- data1
        colname_tmp <- unlist(strsplit(colnames(data1)[1], " "))
        colname_1 <- paste(colname_tmp[length(colname_tmp)]," a",sep="")
        colname_2 <- paste(colname_tmp[length(colname_tmp)]," b",sep="")
        colnames(Result) <- c(colname_1, colname_2, "corr_coef", "N", "p_value", "FDR")
        Result_line_index <- 0
        foreach( i1 = c(1:(ncol(data1)-1)) ) %dopar% {
            foreach( i2 = c((i1+1):ncol(data1)) ) %dopar% {
                Result_line_index <- Result_line_index + 1
                keeplink <- which(!is.na(data1[,i1]) & !is.na(data1[,i2]))
                data1used <- as.numeric(as.vector(data1[keeplink, i1]))
                data2used <- as.numeric(as.vector(data1[keeplink, i2]))
                Result[Result_line_index, colname_1] <- colnames(data1)[i1]
                Result[Result_line_index, colname_2] <- colnames(data1)[i2]
                Result[Result_line_index, "N"] <- length(data1used)
                if (Result[Result_line_index, "N"] >= min_samples) {
                    tempresult <- try( tt <- cor.test(data1used, data2used,
                                                      method = method, alternative = alternative))
                    if(length(grep("error", tempresult, ignore.case = T)) == 0){
                        Result[Result_line_index, "corr_coef"] <- tempresult$estimate[[1]]
                        Result[Result_line_index, "p_value"] <- tempresult$p.value
                    }####### end if of no error
                }####### end if of more than min_samples data
            }####### end if of data2
        }####### end if of data1

    } else { #### corss data1, data2 correlation

        data1 = read.table(paste(job_dir, datafiles[1], sep=fileSep), sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F)
        data2 = read.table(paste(job_dir, datafiles[2], sep=fileSep), sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F)
        if (length(grep('mRNA', colnames(data1))) > 0) {
            data1[data1<0] <- NA
            data1 <- log(data1 + 1e-6, 2)
        }
        if (length(grep('mRNA', colnames(data2))) > 0) {
            data2[data2<0] <- NA
            data2 <- log(data2 + 1e-6, 2)
        }
        all_data <- cbind(data1, data2)
        Result <- data.frame(matrix(NA, dim(data1)[2] * dim(data2)[2], 6, byrow=TRUE), check.names=F)
        colname_tmp1 <- unlist(strsplit(colnames(data1)[1], " "))
        colname_tmp2 <- unlist(strsplit(colnames(data2)[1], " "))
        colname_1 <- colname_tmp1[length(colname_tmp1)]
        colname_2 <- colname_tmp2[length(colname_tmp2)]
        colnames(Result) <- c(colname_1, colname_2, "corr_coef", "N", "p_value", "FDR")

        Result_line_index <- 0
        foreach( i1 = c(1:ncol(data1)) ) %dopar% {
            foreach( i2 = c(1:ncol(data2)) ) %dopar% {
                Result_line_index <- Result_line_index + 1
                keeplink <- which(!is.na(data1[,i1]) & !is.na(data2[,i2]))
                data1used <- as.numeric(as.vector(data1[keeplink, i1]))
                data2used <- as.numeric(as.vector(data2[keeplink, i2]))
                Result[Result_line_index, colname_1] <- colnames(data1)[i1]
                Result[Result_line_index, colname_2] <- colnames(data2)[i2]
                Result[Result_line_index, "N"] <- length(data1used)
                if (Result[Result_line_index, "N"] >= min_samples) {
                    tempresult <- try( tt <- cor.test(data1used, data2used,
                                                    method = method, alternative = alternative))
                    if(length(grep("error", tempresult, ignore.case = T)) == 0){
                        Result[Result_line_index, "corr_coef"] <- tempresult$estimate[[1]]
                        Result[Result_line_index, "p_value"] <- tempresult$p.value
                    }####### end if of no error

                }####### end if of more than min_samples data
            }####### end if of data2
        }####### end if of data1
    } #data2 is not NA

    Result$FDR <- p.adjust(Result$p_value, method = pvalcorrection)
    ## write output
    colnames(Result)[colnames(Result)=="p_value"] <- paste("p_value(", method, ")", sep="")

    reformat_out_dataframe = structure(list(rownames=colnames(all_data),
                                            datalines=t(as.matrix(all_data))))
    outFilename <- paste(gsub("\\-", "", UUIDgenerate()), "json", sep=".")
    write(toJSON(reformat_out_dataframe, pretty = T),
          file=paste(job_dir, outFilename, sep=fileSep))

    return(list(table = writeResults(Result, job_dir, outformat),
                scatterplot = outFilename))
}
