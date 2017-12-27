#' writeResults function
#' print results to uder folder
#' @param results, result data
#' @return nothing
#' @import uuid
#' @importFrom jsonlite toJSON
#' @examples
#' writeResults(data = results, format = "tsv", job_dir)
#' @export
writeResults <- function(out_dataframe, job_dir, format) {

    fileSep <- .Platform$file.sep
    outFilename <- paste(gsub("\\-", "", UUIDgenerate()), format, sep=".")
    switch(format,
           csv = {
               write.csv(out_dataframe, file=paste(job_dir, outFilename, sep=fileSep),
                         sep=",", quote=F, row.names=T, col.names=T)
           },
           rds = {
               saveRDS(out_dataframe, paste(job_dir, outFilename, sep=fileSep))
           },
           tsv = {
               write.table(out_dataframe, file=paste(job_dir, outFilename, sep=fileSep),
                           sep="\t", quote=F, row.names=T, col.names=T)
           },
           {
               reformat_out_dataframe = structure(list(
                                                    colnames=colnames(out_dataframe),
                                                    datalines=as.matrix(out_dataframe)))
               write(toJSON(reformat_out_dataframe, pretty = T),
                     file=paste(job_dir, outFilename, sep=fileSep))
           }
    )
    return(outFilename)
}
