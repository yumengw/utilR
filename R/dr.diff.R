#' Differential function
#' this function will calculate differentail of all data based on group
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @import logging
#' @import uuid
#' @param datafiles contain file1, data matrix used for analysis, each column is a data vector
#' and file2, group information for each patient
#' @param method, t-test, rank-test, (wilcox-test, kruskal, default test)
#' @param alternative, two.sided, less, greater (default two.sided)
#' @param pvalcorrection, fdr, bonferroni
#' @param min_samples, minimum samples for correlation calculation
#' @return results
#' @export

drdiff <- function(datafiles, job_dir, outformat="json", method="rank-test",
                   alternative="two.sided", min_samples = 5, pvalcorrection="fdr") {

    fileSep <- .Platform$file.sep
    group_file <- ''
    data_file <- ''
    for (eachfile in datafiles) {
        if (length(grep('_grp.tsv', eachfile)) > 0) group_file <- eachfile
        else data_file <- eachfile
    }

    data = read.table(paste(job_dir, data_file, sep=fileSep), sep="\t", header=T,
                      row.names=1, check.names=F, stringsAsFactors=F)
    group = read.table(paste(job_dir, group_file, sep=fileSep), sep="\t", header=T,
                       row.names=1, check.names=F, stringsAsFactors=F)

    if (length(grep('mRNA', colnames(data))) > 0) {
        data[data<0] <- NA
        data <- log(data + 1e-6, 2)
    }
    Result <- data.frame(matrix(NA, ncol(data), 6, byrow=TRUE), check.names=F)
    colname_tmp <- unlist(strsplit(colnames(data)[1], " "))
    colname_1 <- colname_tmp[length(colname_tmp)]
    colnames(Result) <- c(colname_1, "group", "N", "median", "p_value", "FDR")

    boxplot_data_all <- list()
    boxplot_data_outlier <- list()
    boxplot_data_header <- list()

    #foreach( i1 = c(1:(ncol(data))) ) %dopar% {
    for (i1 in c(1:ncol(data))) {

        group_count <- as.data.frame(table(group$group))
        group_count_used <- as.character(
                            group_count[group_count$Freq > min_samples, "Var1"])
        test.data <- data.frame(cbind(data[, i1], group$group), stringsAsFactors=F)
        colnames(test.data) <- c("value", "group")
        test.data <- test.data[!is.na(test.data$value) &
                                   test.data$group %in% group_count_used, ]

        test.data$value <- as.numeric(as.character(test.data$value))
        group_count <- as.data.frame(table(test.data$group))
        Result[i1,colname_1] <- colnames(data)[i1]
        Result[i1,"N"] <- paste(group_count$Freq, collapse="|")
        Result[i1,"group"] <- paste(group_count$Var1, collapse="|")
        median_tmp <- aggregate(test.data$value, list(test.data$group), median)
        Result[i1,"median"] <- paste(median_tmp$x, collapse="|")

        if (length(unique(test.data$group)) < 2)
            print("only one group cannot compare")

        if (method=="rank-test") {
            if (length(unique(test.data$group)) == 2) {
                tempresult <- try( tt <- wilcox.test(test.data$value ~ test.data$group))
                if (length(grep("error", tempresult, ignore.case = T)) == 0){
                    Result[i1,"p_value"] <- tempresult$p.value
                }
            }
            if (length(unique(test.data$group)) > 2) {
                tempresult <- try( tt <- kruskal.test(test.data$value ~ test.data$group))
                if (length(grep("error", tempresult, ignore.case = T)) == 0){
                    Result[i1,"p_value"] <- tempresult$p.value
                }
            }
        }

        if (method=="t-test") {
            if (length(unique(test.data$group)) != 2)
                print("Need to have two groups for t-test")
            else if (length(unique(test.data$group)) == 2) {
                tempresult <- try( tt <- t.test(test.data$value ~ test.data$group))
                if (length(grep("error", tempresult, ignore.case = T)) == 0){
                    Result[i1,"p_value"] <- tempresult$p.value
                }
            }
        }

        split_data <- split(test.data, test.data$group)
        unique_group <- unique(test.data$group)
        boxplot_data_each <- matrix(NA, length(unique_group), 5)
        boxplot_data_outlier_each <- list()
        for (i2 in c(1:length(unique_group))) {
            this_group_data <- split_data[[unique_group[i2]]]
            boxplot_data_each[i2, ] <- boxplot.stats(this_group_data$value)$stats
            boxplot_data_outlier_each[[i2]] <- boxplot(this_group_data$value, plot=FALSE)$out
        }
        boxplot_data_all[[i1]] <- boxplot_data_each
        boxplot_data_header[[i1]] <- t(matrix(unique_group))
        boxplot_data_outlier[[i1]] <- boxplot_data_outlier_each
    }

    Result$FDR <- p.adjust(Result$p_value, method=pvalcorrection)
    ## boxplot data
    reformat_out_dataframe = structure(list(
                                rownames=colnames(data),
                                colnames=boxplot_data_header,
                                datalines=boxplot_data_all))
    outFilename_data <- paste(gsub("\\-", "", UUIDgenerate()), outformat, sep=".")
    write(toJSON(reformat_out_dataframe, pretty = T),
            file=paste(job_dir, outFilename_data, sep=fileSep))
    ## boxplot outlier data
    reformat_out_dataframe = structure(list(
                                rownames=colnames(data),
                                colnames=boxplot_data_header,
                                datalines=boxplot_data_outlier))
    outFilename_outlier <- paste(gsub("\\-", "", UUIDgenerate()), outformat, sep=".")
    write(toJSON(reformat_out_dataframe, pretty = T),
          file=paste(job_dir, outFilename_outlier, sep=fileSep))

    return(list(table = writeResults(Result, job_dir, outformat),
               boxplot = outFilename_data,
               outlier = outFilename_outlier))
}
