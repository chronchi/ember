#' Get regressed data
#'
#' @inheritParams prepare_data_normalization
#' @inheritParams get_pca_coordinates
#' @param remove_components A vector with the principal components
#'     to be removed. Default is PC1, PC2 and PC5
#'
#' @return A summarized experiment with the regressed data and the
#'     pca coordinates.
#' @export
get_regressed_data <- function(
    sum_exp,
    assay_to_use,
    which_pca = "big",
    stable_genes = NULL,
    remove_components = c("PC1", "PC2", "PC5")
){

    if (!(which_pca %in% c("small", "big"))){
        stop(
            paste0(
                "which_pca should be either 'big' or 'small', please ",
                "double check your function call."
            )
        )
    }

    # we first normalize the dataset and then get the pca coordinates
    # to later regress out the batch effects
    sum_exp <- get_final_ranking_values(
        sum_exp = sum_exp,
        assay_to_use = assay_to_use,
        stable_genes = stable_genes,
        most_variable_genes = which_pca
    )

    print("Normalization done.")

    # the steps for this function is to first get the pca coordinates
    # from the normalized dataset and then convert it back.
    pca_coordinates <- get_pca_coordinates(
        sum_exp = sum_exp,
        which_pca = which_pca
    )

    if (which_pca == "big"){
        pca_fit_to_use <- pca_fit_all
    } else if ( which_pca == "small") {
        pca_fit_to_use <- pca_fit
    }

    name_components <- colnames(pca_fit_to_use$loadings)
    keep_components <- setdiff(name_components, remove_components)

    genes_in_common <- intersect(
        rownames(sum_exp),
        rownames(pca_fit_to_use$loadings)
    )

    df_regressed <- as.matrix(
        SummarizedExperiment::colData(pca_coordinates)[, keep_components]
    ) %*% t(pca_fit_to_use$loadings[genes_in_common, keep_components]) %>%
        t

    sum_exp <- sum_exp[rownames(df_regressed), ]
    SummarizedExperiment::assay(sum_exp, "regressed") <- df_regressed

    sum_exp

}

#' Calculate scores for samples with regressed data
#'
#' @param sum_exp A summarized experiment object with the
#'     assay "regressed"
#' @param gene_sets A list with gene sets that will be used to
#'     calculate the scores
#'
#' @return A summarized experiment with the scores calculated
#'    as columns in the colData slot
#'
#' @importFrom SummarizedExperiment colData assay
#' @importFrom dplyr %>%
#' @export
get_scores_regressed <- function(
    sum_exp,
    gene_sets
){
    gene_sets_scores <- sapply(
        gene_sets,
        function(gene_set, df_pcs_regressed){

            gene_set <- intersect(
                gene_set,
                rownames(df_pcs_regressed)
            )
            colSums(df_pcs_regressed[gene_set, ])

        },
        df_pcs_regressed = assay(sum_exp, "regressed")
    ) %>%
        data.frame

    SummarizedExperiment::colData(sum_exp)[
        ,
        paste0("regressed_", names(gene_sets))
    ] <- gene_sets_scores

    sum_exp

}
