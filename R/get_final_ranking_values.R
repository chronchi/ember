#' Prepare the data for normalization
#'
#' @description
#' `prepare_data_normalization` checks if the data structure
#'     is right and if it has the appropriate slots
#'
#' @details
#' The function check if the object is a summarized experiment object
#' first, so the normalization procedures can be saved later on. It also
#' checks if the assay_to_use is available on the summarized experiment
#' object, without this it cannot proceed with the normalization. Lastly
#' it proceeds to subset the dataframe to only contain the necessary genes
#' to calculate the ranking obtained in the procedure described
#' in the paper.
#'
#' @param sum_exp Summarized experiment object. It should contain a "rank"
#'     slot that will be used to calculate the average rank of the stable
#'     genes.
#' @param assay_to_use A string. Which assay to use when calculating the
#'     average expression of the stable genes.
#' @param most_variable_genes A vector. Most variable genes determined by some
#'     procedure. Default is "small" and uses an internal list of 1000 genes
#'     to calculate the embedding. If using for calculating the regressed
#'     data, use the "big" option.
#' @param stable_genes A vector. A list with 44 genes that will be used
#'     to the qPCR-like normalization if NULL (default).
#' @param verbose. An integer. If verbose equals to 1, then it prints the
#'     number of stable genes, total genes in the dataframe and total number
#'     of samples
#' @examples
#' \dontrun{
#' prepare_data_normalization(
#'     tcga,
#'     "logFPKM_TMM",
#'     c("GAPDH"),
#'     c("ESR1", "GREB1")
#' )
#' }
#' @importFrom SummarizedExperiment assay
#' @importFrom singscore rankGenes
prepare_data_normalization <- function(
    sum_exp,
    assay_to_use,
    stable_genes = NULL,
    most_variable_genes = "small",
    verbose = 1
){

    sum_exp_classes <- c(
        "SummarizedExperiment",
        "RangedSummarizedExperiment"
    )

    if( !(is(sum_exp, "SummarizedExperiment")) ){
        stop(
            paste0(
                "Object is not from class(es) ",
                paste(sum_exp_classes, collapse = ", "),
                ". Please check your data structure."
            )
        )
    }

    if( !(assay_to_use %in% SummarizedExperiment::assayNames(sum_exp)) ){
        stop(
            paste0(
                "Assay not available, check if you specified the right ",
                "name or you included your assay in the data structure."
            )
        )
    }

    if (most_variable_genes == "small"){
        most_variable_genes <- most_variable_genes_small
    } else if (most_variable_genes == "big"){
        most_variable_genes <- most_variable_genes_big
    }

    if (is.null(stable_genes)) {
        stable_genes <- stable_genes_list
    }

    # subselect the necessary genes for downstream procedures
    sum_exp <- sum_exp[
        intersect(rownames(sum_exp), c(stable_genes, most_variable_genes)),
    ]

    if (verbose == 1){
        cat(paste0(
            "Total number of stable genes: ",
            length(intersect(rownames(sum_exp), stable_genes)), "\n",
            "Total number of genes: ", nrow(sum_exp), "\n",
            "Number of samples: ", ncol(sum_exp), "\n"
        ))
    }

    # calculate the ranking that will be used in the other functions
    SummarizedExperiment::assay(sum_exp, "rank") <- singscore::rankGenes(
        as.matrix(SummarizedExperiment::assay(sum_exp, assay_to_use))
    )

    sum_exp
}


#' Calculate the average ranking and average expression of stable genes
#'
#' @inheritParams prepare_data_normalization
#' @examples
#' \dontrun{
#' avg_ranking(tcga, "logFPKM_TMM", c("GAPDH"))
#' }
#' @importFrom SummarizedExperiment assay
#' @importFrom dplyr %>%
avg_ranking <- function(sum_exp, assay_to_use, stable_genes = NULL){

    if (is.null(stable_genes)) {
        stable_genes <- stable_genes_list
    }

    available_stable_genes <- intersect(rownames(sum_exp), stable_genes)

    rank_stable_genes <- SummarizedExperiment::assay(sum_exp, "rank")[available_stable_genes, ]
    expr_stable_genes <- SummarizedExperiment::assay(sum_exp, assay_to_use)[available_stable_genes, ]


    values <- list(
        avg_ranking = colMeans(rank_stable_genes)/nrow(sum_exp),
        avg_expression = colMeans(expr_stable_genes)
    )

    values$sd_ranking <- sd(values$avg_ranking)
    values$sd_expression <- sd(values$avg_expression)

    values
}

#' Calculate the average expression and ranking based on stable genes
#'
#' @inheritParams avg_ranking
#' @param avg_ranking_sds A list. Output from the function avg_ranking
#' @examples
#' \dontrun{
#' calculate_norm_ranks(tcga, "logFPKM_TMM", output_avg_ranking)
#' }
#' @importFrom SummarizedExperiment assay
calculate_norm_ranks <- function(sum_exp, assay_to_use, avg_ranking_sds){

    # for each patient calculate the division for each gene in each
    # element of average expression obtained from the avg_ranking function
    df <- sapply(
        1:ncol(sum_exp),
        function(i, sum_exp, avg_ranking_sds){
            expression_levels <- SummarizedExperiment::assay(sum_exp, assay_to_use)[,i]
            expression_levels/avg_ranking_sds$avg_expression[i]
        },
        sum_exp = sum_exp,
        avg_ranking_sds = avg_ranking_sds
    ) %>%
        `colnames<-`(colnames(sum_exp)) %>%
        `rownames<-`(rownames(sum_exp))

    # we now add the normalized average expression to the original summarized
    # experiment so it can be used to other analysis
    SummarizedExperiment::assay(sum_exp, "avg_expression") <- df

    # notice here we are converting the average normalized ranking back to
    # average rank, since we will divide the rank of the gene and not
    # the normalized rank. Also all this operations are sums, so we
    # can go back and forth here.
    df <- sapply(
        1:ncol(sum_exp),
        function(i, sum_exp, avg_ranking_sds){
            rank_genes <- SummarizedExperiment::assay(sum_exp, "rank")[,i]
            rank_genes/(avg_ranking_sds$avg_ranking[i]*nrow(sum_exp))
        },
        sum_exp = sum_exp,
        avg_ranking_sds = avg_ranking_sds
    ) %>%
        `colnames<-`(colnames(sum_exp)) %>%
        `rownames<-`(rownames(sum_exp))

    # we now add the normalized average ranking to the original summarized
    # experiment so it can be used to other analysis
    SummarizedExperiment::assay(sum_exp, "avg_ranking") <- df

    sum_exp$avg_ranking <- avg_ranking_sds$avg_ranking
    sum_exp$avg_expression <- avg_ranking_sds$avg_expression

    sum_exp

}

#' Wrapper for prepare_data_normalization avg_ranking and calculate_norm_ranks
#'
#' @inheritParams prepare_data_normalization
#' @returns A summarized experiment object containing the normalized
#'     data that can be applied later on to obtain the embedding using
#'     the loadings of the PCA.
#' @examples
#' \dontrun{
#' get_final_ranking_values(tcga, "logFPKM_TMM")
#' }
#' @importFrom SummarizedExperiment assay
#' @export
get_final_ranking_values <- function(
    sum_exp,
    assay_to_use,
    stable_genes = NULL,
    most_variable_genes = "small",
    verbose = 1
){

    sum_exp <- prepare_data_normalization(
        sum_exp = sum_exp,
        assay_to_use = assay_to_use,
        stable_genes = stable_genes,
        most_variable_genes = most_variable_genes,
        verbose = verbose
    )

    avg_ranking_sds <- avg_ranking(
        sum_exp = sum_exp,
        assay_to_use = assay_to_use,
        stable_genes = stable_genes
    )

    sum_exp <- calculate_norm_ranks(
        sum_exp = sum_exp,
        avg_ranking_sds = avg_ranking_sds,
        assay_to_use = assay_to_use
    )

    SummarizedExperiment::assay(sum_exp, "znorm") <- scale(SummarizedExperiment::assay(sum_exp, assay_to_use))

    sum_exp

}
