#' preprocessing for coxph function
#' @param data, data matrix used for analysis, each row is a data vector, or categorical values
#' @param time, survival time with same length with data
#' @param status, survival status with same length with data
#' @return results
#' @export
coxph_preprocessing <- function(data, time, status) {

    usedpatient_id <- intersect(data$patient_id, time$patient_id)
    usedpatient_id <- intersect(usedpatient_id, status$patient_id)
    data <- data[match(usedpatient_id, data$patient_id), ]
    time <- time[match(usedpatient_id, time$patient_id), ]
    status <- status[match(usedpatient_id, status$patient_id), ]
    print(identical(data$patient_id, time$patient_id))
    print(identical(status$patient_id, time$patient_id))

    rownames(data) <- data$sample_id
    data$sample_id <- NULL
    data$patient_id <- NULL
    time$patient_id <- NULL
    status$patient_id <- NULL

    data <- t(data)
    time <- as.vector(as.numeric(t(time)))
    status <- as.vector(as.numeric(t(status)))
    return(list(data = data, time = time, status = status))
}
