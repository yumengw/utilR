#' mutWaterfall function
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @import survival
#' @import logging
#' @import maftools
#' @import uuid
#' @importFrom jsonlite toJSON
#' @param mafile, maf file name
#' @examples
#' drmutWaterfall(mafile, "/Users/yumengw/projects/drright/drright-modules/")
#' @return results
#' @export

drmutWaterfall <- function(mafile, job_dir, format="json") {
    fileSep <- .Platform$file.sep
    lmaf = read.maf(maf=paste(job_dir, mafile, sep=fileSep), useAll = TRUE, verbose = FALSE)
    top <- min(20, nrow(getGeneSummary(lmaf)))
    mat = oncoplot(maf = lmaf, top = top, removeNonMutated = T, writeMatrix = T, sortByMutation = T)
    sparse_mat = matrix(NA, nrow(mat)*ncol(mat), 3)
    mutTypes <- unique(as.vector(as.matrix(mat)))
    mutTypes <- mutTypes[!mutTypes==""] # remove empty
    row_count <- 1
    for (i in c(1:nrow(mat))) {
        for (j in c(1:ncol(mat))) {
            if (mat[i, j] != "") {
                sparse_mat[row_count, ] <- c(i, j, which(mutTypes==mat[i, j]))
                row_count <- row_count + 1
            }
        }
    }
    sparse_mat <- sparse_mat[!is.na(sparse_mat[,1]), ]

    reformat_out_dataframe = structure(list(
                                        colnames=colnames(mat),
                                        rownames=rownames(mat),
                                        muttypes=mutTypes,
                                        datalines=as.matrix(sparse_mat)))
    outFilename <- paste(gsub("\\-", "", UUIDgenerate()), format, sep=".")
    write(toJSON(reformat_out_dataframe, pretty = T),
          file=paste(job_dir, outFilename, sep=fileSep))
    ## return file list
    return(list(waterfall = outFilename))
}
