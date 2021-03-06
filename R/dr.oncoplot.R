#' oncoplot function
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @import survival
#' @import logging
#' @import ComplexHeatmap
#' @import grid
#' @import uuid
#' @importFrom jsonlite toJSON
#' @param mafile, maf file name
#' @examples
#' drmutWaterfall(mafile, "/Users/yumengw/projects/drright/drright-modules/")
#' @return results
#' @export

droncoplot <- function(mafile, job_dir, format="json") {
    fileSep <- .Platform$file.sep
    mat = read.delim(paste(job_dir, mafile, sep=fileSep), sep="\t", header=T, check.names=F,
                     stringsAsFactors=F)
    mat[is.na(mat)] = ""
    rownames(mat) = mat[, 1]
    mat = mat[, -1]
    ######################## plot ########################

    mut_type = list(mut=c("3'UTR", "5'Flank", "5'UTR", "De_novo_Start_InFrame",
                          "De_novo_Start_OutOfFrame", "Frame_Shift_Del", "Frame_Shift_Ins",
                          "IGR", "In_Frame_Del", "In_Frame_Ins", "Intron", "lincRNA",
                          "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
                          "RNA", "Silent", "Splice_Site", "Start_Codon_Del",
                          "Start_Codon_SNP", "Stop_Codon_Del", "AMP", "HOMDEL"),
                    color=c("pink", "pink", "pink", "pink", "pink", "brown", "brown",
                            "pink", "black", "black", "pink", "yellow", "#3cba54",
                            "black", "black", "pink", "#CCCCCC", "black", "black",
                            "black", "purple", "#d62d20", "#437af8"))

    col <- as.vector(mut_type$color)
    names(col) <- mut_type$mut
    alter_fun = list(
        background = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
        },
        "3'UTR" = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["3'UTR"], col = NA))
        },
        "5'Flank" = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["5'Flank"], col = NA))
        },
        "5'UTR" = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["5'UTR"], col = NA))
        },
        "De_novo_Start_InFrame" = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["De_novo_Start_InFrame"], col = NA))
        },
        "De_novo_Start_OutOfFrame" = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["De_novo_Start_OutOfFrame"], col = NA))
        },
        "Frame_Shift_Del" = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
        },
        "Frame_Shift_Ins" = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))
        },
        "IGR" = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["IGR"], col = NA))
        },
        "In_Frame_Del" = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["In_Frame_Del"], col = NA))
        },
        "In_Frame_Ins" = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["In_Frame_Ins"], col = NA))
        },
        "Intron" = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Intron"], col = NA))
        },
        "lincRNA" = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["lincRNA"], col = NA))
        },
        "Missense_Mutation" = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Missense_Mutation"], col = NA))
        },
        "Nonsense_Mutation" = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
        },
        "Nonstop_Mutation" = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Nonstop_Mutation"], col = NA))
        },
        "RNA" = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["RNA"], col = NA))
        },
        "Silent" = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Silent"], col = NA))
        },
        "Splice_Site" = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Splice_Site"], col = NA))
        },
        "Start_Codon_Del" = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Start_Codon_Del"], col = NA))
        },
        "Start_Codon_SNP" = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Start_Codon_SNP"], col = NA))
        },
        "Stop_Codon_Del" = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Stop_Codon_Del"], col = NA))
        },
        "AMP" = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["AMP"], col = NA))
        },
        "HOMDEL" = function(x, y, w, h) {
            grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["HOMDEL"], col = NA))
        }
    )


    pdf(paste(job_dir, "CoMutplot.pdf", sep=fileSep), 7, 5)

    mutTypes <- unique(as.vector(as.matrix(mat)))
    mutTypes <- mutTypes[!mutTypes==""] # remove empty
    mutTypes <- unique(unlist(strsplit(mutTypes, ";")))

    ht = myoncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
                   alter_fun = alter_fun, col = col,
                   column_title = "",
                   remove_empty_columns = FALSE,
                   heatmap_legend_param = list(title = "Alternations", at = mutTypes,
                                               labels = mutTypes,
                                               nrow = 1, title_position = "leftcenter"))

    draw(ht[["ht"]], heatmap_legend_side = "bottom")
    dev.off()

    ######################################################
    mat <- mat[as.numeric(ht[["row_order"]]), as.numeric(ht[["col_order"]])]
    sparse_mat = matrix(NA, nrow(mat)*ncol(mat)*23, 3)

    row_count <- 1
    for (i in c(1:nrow(mat))) {
        for (j in c(1:ncol(mat))) {
            if (mat[i, j] != "") {
                tmp_variant <- unlist(strsplit(mat[i, j], ";"))
                for (k in c(1:length(tmp_variant))) {
                    sparse_mat[row_count, ] <- c(i, j, which(mutTypes==tmp_variant[k]))
                    row_count <- row_count + 1
                }
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

myoncoPrint = function(mat, get_type = function(x) x,
                     alter_fun = alter_fun_list, alter_fun_list = NULL, col,
                     row_order = oncoprint_row_order(),
                     column_order = oncoprint_column_order(),
                     show_pct = TRUE, pct_gp = row_names_gp, pct_digits = 0,
                     axis_gp = gpar(fontsize = 8),
                     show_row_barplot = TRUE,
                     row_barplot_width = unit(2, "cm"),
                     remove_empty_columns = FALSE,
                     heatmap_legend_param = list(title = "Alterations"),
                     top_annotation = HeatmapAnnotation(column_bar = anno_oncoprint_barplot(),
                                                        annotation_height = unit(2, "cm")),
                     top_annotation_height = top_annotation@size,
                     bottom_annotation = new("HeatmapAnnotation"),
                     bottom_annotation_height = bottom_annotation@size,
                     barplot_ignore = NULL,
                     row_title = character(0),
                     row_title_side = c("left", "right"),
                     row_title_gp = gpar(fontsize = 14),
                     row_title_rot = switch(row_title_side[1], "left" = 90, "right" = 270),
                     column_title = character(0),
                     column_title_side = c("top", "bottom"),
                     column_title_gp = gpar(fontsize = 14),
                     column_title_rot = 0,
                     show_row_names = TRUE,
                     row_names_gp = gpar(fontsize = 12),
                     show_column_names = FALSE,
                     column_names_gp = gpar(fontsize = 12),
                     split = NULL,
                     gap = unit(1, "mm"),
                     combined_name_fun = function(x) paste(x, collapse = "/"),
                     width = NULL,
                     ...) {

    if(length(names(list(...))) > 0) {
        if(any(names(list(...)) %in% c("show_column_barplot", "column_barplot_height"))) {
            stop("`show_column_barplot` and `column_barplot_height` is deprecated, please configure `top_annotation` directly.")
        }
    }

    if(!is.null(alter_fun_list)) {
        warning("`alter_fun_list` is deprecated, please `alter_fun` instead.")
    }

    # convert mat to mat_list
    if(inherits(mat, "data.frame")) {
        mat = as.matrix(mat)
    }
    if(inherits(mat, "matrix")) {
        all_type = unique(unlist(lapply(mat, get_type)))
        all_type = all_type[!is.na(all_type)]
        all_type = all_type[grepl("\\S", all_type)]

        mat_list = lapply(all_type, function(type) {
            m = sapply(mat, function(x) type %in% get_type(x))
            dim(m) = dim(mat)
            dimnames(m) = dimnames(mat)
            m
        })
        names(mat_list) = all_type
    } else if(inherits(mat, "list")) {
        mat_list = mat

        all_type = names(mat_list)
        mat_list = lapply(mat_list, function(x) {
            if(!is.matrix(x)) {
                stop("Expect a list of matrix (not data frames).")
            }
            oattr = attributes(x)
            x = as.logical(x)
            attributes(x) = oattr
            x
        })

        if(length(unique(sapply(mat_list, nrow))) > 1) {
            stop("All matrix in 'mat_list' should have same number of rows.")
        }

        if(length(unique(sapply(mat_list, ncol))) > 1) {
            stop("All matrix in 'mat_list' should have same number of columns.")
        }
    } else {
        stop("Incorrect type of 'mat'")
    }

    if(missing(alter_fun) && missing(col)) {
        if(length(mat_list) == 1) {
            af = function(x, y, w, h, v) {
                grid.rect(x, y, w, h, gp = gpar(fill = "#CCCCCC", col = NA))
                if(v[1]) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "red", col = NA))
            }
            col = "red"
        } else if(length(mat_list) == 2) {
            af = function(x, y, w, h, v) {
                grid.rect(x, y, w, h, gp = gpar(fill = "#CCCCCC", col = NA))
                if(v[1]) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "red", col = NA))
                if(v[2]) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = "blue", col = NA))
            }
            col = c("red", "blue")
        } else {
            stop("`alter_fun` should be specified.")
        }
        names(col) = names(mat_list)
    } else if(is.list(alter_fun)) {

        # validate the list first
        if(is.null(alter_fun$background)) alter_fun$background = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = "#CCCCCC", col = NA))
        sdf = setdiff(all_type, names(alter_fun))
        if(length(sdf) > 0) {
            stop(paste0("You should define shape function for: ", paste(sdf, collapse = ", ")))
        }

        alter_fun = alter_fun[unique(c("background", intersect(names(alter_fun), all_type)))]

        af = function(x, y, w, h, v) {
            if(!is.null(alter_fun$background)) alter_fun$background(x, y, w, h)

            alter_fun = alter_fun[names(alter_fun) != "background"]

            if(sum(v)) {
                for(nm in names(alter_fun)) {
                    if(v[nm]) alter_fun[[nm]](x, y, w, h)
                }
            }
        }
    } else {
        af = alter_fun
    }

    col = col[intersect(names(col), all_type)]

    # type as the third dimension
    arr = array(FALSE, dim = c(dim(mat_list[[1]]), length(all_type)), dimnames = c(dimnames(mat_list[[1]]), list(all_type)))
    for(i in seq_along(all_type)) {
        arr[, , i] = mat_list[[i]]
    }

    oncoprint_row_order = function() {
        order(rowSums(count_matrix), decreasing = TRUE)
    }

    oncoprint_column_order = function() {
        scoreCol = function(x) {
            score = 0
            for(i in 1:length(x)) {
                if(x[i]) {
                    score = score + 2^(length(x)-i*1/x[i])
                }
            }
            return(score)
        }
        scores = apply(count_matrix[row_order, ,drop = FALSE], 2, scoreCol)
        order(scores, decreasing=TRUE)
    }

    count_matrix = apply(arr, c(1, 2), sum)
    if(is.null(row_order)) row_order = seq_len(nrow(count_matrix))
    if(is.null(column_order)) column_order = seq_len(ncol(count_matrix))
    row_order = row_order
    if(is.character(column_order)) {
        column_order = structure(seq_len(dim(arr)[2]), names = dimnames(arr)[[2]])[column_order]
    }
    column_order = column_order
    names(column_order) = as.character(column_order)
    if(remove_empty_columns) {
        l = rowSums(apply(arr, c(2, 3), sum)) > 0
        arr = arr[, l, , drop = FALSE]
        column_order = structure(seq_len(sum(l)), names = which(l))[as.character(intersect(column_order, which(l)))]
    }

    # validate col
    sdf = setdiff(all_type, names(col))
    if(length(sdf) > 0) {
        stop(paste0("You should define colors for:", paste(sdf, collapse = ", ")))
    }

    # for each gene, percent of samples that have alterations
    pct = rowSums(apply(arr, 1:2, any)) / ncol(mat_list[[1]])
    pct = paste0(round(pct * 100, digits = pct_digits), "%")
    ha_pct = rowAnnotation(pct = row_anno_text(pct, just = "right", offset = unit(1, "npc"), gp = pct_gp), width = max_text_width(pct, gp = pct_gp))

    #####################################################################
    # row annotation which is a barplot
    anno_row_bar = function(index, k = NULL, N = NULL) {
        n = length(index)
        count = apply(arr, c(1, 3), sum)[index, , drop = FALSE]
        all_type = all_type[!(colnames(count) %in% barplot_ignore)]
        count = count[, setdiff(colnames(count), barplot_ignore), drop = FALSE]
        max_count = max(rowSums(count))
        pushViewport(viewport(xscale = c(0, max_count*1.1), yscale = c(0.5, n + 0.5)))
        for(i in seq_len(nrow(count))) {
            if(any(count[i, ] > 0)) {
                x = count[i, ]
                x = x[x > 0]
                x2 = cumsum(x)
                type = all_type[count[i, ] > 0]
                # row order is from top to end while coordinate of y is from bottom to top
                # so here we need to use n-i+1
                grid.rect(x2, n-i+1, width = x, height = 0.8, default.units = "native", just = "right", gp = gpar(col = NA, fill = col[type]))
            }
        }
        breaks = grid.pretty(c(0, max_count))
        if(k == 1) {
            grid.xaxis(at = breaks, label = breaks, main = FALSE, gp = axis_gp)
        }
        upViewport()
    }

    ha_row_bar = rowAnnotation(row_bar = anno_row_bar, width = row_barplot_width)

    ###################################################################
    # column annotation which is also a barplot
    anno_column_bar = function(index) {
        n = length(index)
        count = apply(arr, c(2, 3), sum)[index, , drop = FALSE]
        all_type = all_type[!(colnames(count) %in% barplot_ignore)]
        count = count[, setdiff(colnames(count), barplot_ignore), drop = FALSE]
        max_count = max(rowSums(count))
        pushViewport(viewport(yscale = c(0, max_count*1.1), xscale = c(0.5, n + 0.5)))
        for(i in seq_len(nrow(count))) {
            if(any(count[i, ] > 0)) {
                y = count[i, ]
                y = y[y > 0]
                y2 = cumsum(y)
                type = all_type[count[i, ] > 0]
                grid.rect(i, y2, height = y, width = 0.8, default.units = "native", just = "top", gp = gpar(col = NA, fill = col[type]))
            }
        }
        breaks = grid.pretty(c(0, max_count))
        grid.yaxis(at = breaks, label = breaks, gp = axis_gp)
        upViewport()
    }

    top_annotation = top_annotation

    #####################################################################
    # the main matrix
    pheudo = c(all_type, rep(NA, nrow(arr)*ncol(arr) - length(all_type)))
    dim(pheudo) = dim(arr)[1:2]
    dimnames(pheudo) = dimnames(arr)[1:2]

    if(length(list(...))) {
        if(any(names(list(...)) %in% c("rect_gp", "cluster_rows", "cluster_columns", "cell_fun"))) {
            stop("'rect_gp', 'cluster_rows', 'cluster_columns', 'cell_fun' are not allowed to use in `oncoPrint()`.")
        }
    }

    ht = Heatmap(pheudo, col = col, rect_gp = gpar(type = "none"),
                 cluster_rows = FALSE, cluster_columns = FALSE, row_order = row_order, column_order = column_order,
                 cell_fun = function(j, i, x, y, width, height, fill) {
                     z = arr[i, j, ]
                     names(z) = dimnames(arr)[[3]]
                     af(x, y, width, height, z)
                 },
                 top_annotation = top_annotation,
                 top_annotation_height = top_annotation_height,
                 bottom_annotation = bottom_annotation,
                 bottom_annotation_height = bottom_annotation_height,
                 row_title = row_title,
                 row_title_side = row_title_side,
                 row_title_gp = row_names_gp,
                 row_title_rot = row_title_rot,
                 column_title = column_title,
                 column_title_side = column_title_side,
                 column_title_gp = column_title_gp,
                 column_title_rot = column_title_rot,
                 show_row_names = show_row_names,
                 row_names_gp = row_names_gp,
                 show_column_names = show_column_names,
                 column_names_gp = column_names_gp,
                 heatmap_legend_param = heatmap_legend_param,
                 split = split,
                 gap = gap,
                 combined_name_fun = combined_name_fun,
                 width = width,
                 ...)

    ht@matrix_param$oncoprint = list()
    ht@matrix_param$oncoprint$arr = arr
    ht@matrix_param$oncoprint$barplot_ignore = barplot_ignore
    ht@matrix_param$oncoprint$all_type = all_type
    ht@matrix_param$oncoprint$axis_gp = axis_gp
    ht@matrix_param$oncoprint$col = col

    if(show_pct) {
        ht_list = ha_pct + ht
    } else {
        ht_list = ht
    }

    if(show_row_barplot) {
        ht_list = ht_list + ha_row_bar
    }
    return(list(ht = ht_list,
                row_order = row_order,
                col_order = as.numeric(names(column_order))))

}

# == title
# Unify a list of matrix
#
# == param
# -mat_list a list of matrix, all of them should have dimension names
# -default default values for the newly added rows and columns
#
# == details
# All matrix will be unified to have same row names and column names
#
# == value
# A list of matrix
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
unify_mat_list = function(mat_list, default = 0) {
    common_rn = unique(unlist(lapply(mat_list, rownames)))
    common_cn = unique(unlist(lapply(mat_list, colnames)))

    mat_list2 = lapply(seq_along(mat_list), function(i) {
        mat = matrix(default, nrow = length(common_rn), ncol = length(common_cn))
        dimnames(mat) = list(common_rn, common_cn)
        mat[rownames(mat_list[[i]]), colnames(mat_list[[i]])] = mat_list[[i]]
        mat
    })
    names(mat_list2) = names(mat_list)
    return(mat_list2)
}



# == title
# Column barplot annotation for oncoPrint
#
# == details
# This function is only used for column annotation
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
anno_oncoprint_barplot = function() {

    function(index) {
        object = get("object", envir = parent.frame(n = 5))
        arr = object@matrix_param$oncoprint$arr
        barplot_ignore = object@matrix_param$oncoprint$barplot_ignore
        all_type = object@matrix_param$oncoprint$all_type
        axis_gp = object@matrix_param$oncoprint$axis_gp
        col = object@matrix_param$oncoprint$col

        n = length(index)
        count = apply(arr, c(2, 3), sum)[index, , drop = FALSE]
        all_type = all_type[!(colnames(count) %in% barplot_ignore)]
        count = count[, setdiff(colnames(count), barplot_ignore), drop = FALSE]
        max_count = max(rowSums(count))
        pushViewport(viewport(yscale = c(0, max_count*1.1), xscale = c(0.5, n + 0.5)))
        for(i in seq_len(nrow(count))) {
            if(any(count[i, ] > 0)) {
                y = count[i, ]
                y = y[y > 0]
                y2 = cumsum(y)
                type = all_type[count[i, ] > 0]
                grid.rect(i, y2, height = y, width = 0.8, default.units = "native", just = "top", gp = gpar(col = NA, fill = col[type]))
            }
        }
        breaks = grid.pretty(c(0, max_count))
        grid.yaxis(at = breaks, label = breaks, gp = axis_gp)
        upViewport()
    }
}
