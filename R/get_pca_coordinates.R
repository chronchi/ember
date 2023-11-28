#' Get the PCA projection coordinates for a dataset
#'
#' @param sum_exp A summarized experiment object
#' @param assay_to_use A string stating the name of the slot in the
#'     Summarized Experiment object from which the PCs will be obtained.
#' @param which_pca A string stating which PCA object to use to get
#'     the coordinates. Two options available: "small" or "big". The
#'     "small" option uses the fit with only 1044 genes. The
#'     "big" option uses the one with all the 8000+ genes, used
#'     for the regression of the data.
#' @return A dataframe. The PC coordinates are returned for the new
#'     dataset
#' @examples
#' \dontrun{
#' get_pca_coordinates(
#'     scanb,
#'     which_pca = "small",
#'     assay_to_use = "avg_ranking"
#' )
#' }
#' @export
#' @importFrom SummarizedExperiment assay assayNames colData
#' @importFrom dplyr %>%
get_pca_coordinates <- function(
    sum_exp,
    which_pca = "small",
    assay_to_use = "avg_ranking"
){

    # first check if avg_ranking is there if it was
    # selected
    if (assay_to_use == "avg_ranking"){
        if( !("avg_ranking" %in% SummarizedExperiment::assayNames(sum_exp)) ){
            stop(
                paste0(
                    "avg_ranking not available in the object, make sure",
                    " the normalization was done."
                )
            )
        }
    }

    if (which_pca == "small"){
        pca_fit_for_coord <- pca_fit
    } else if (which_pca == "big"){
        pca_fit_for_coord <- pca_fit_all
    }

    loadings_pca <- pca_fit_for_coord$loadings

    # we now add 0 to average ranking for the genes that are not
    # available in the summarized experiment
    genes_for_pca <- rownames(loadings_pca)
    assay_matrix <- SummarizedExperiment::assay(sum_exp, assay_to_use) %>%
        as.matrix
    genes_not_available <- setdiff(genes_for_pca, rownames(assay_matrix))
    if (length(genes_not_available) > 0){
        assay_matrix <- rbind(
            assay_matrix,
            matrix(
                0,
                nrow = length(genes_not_available),
                ncol = ncol(sum_exp),
                dimnames = list(
                    genes_not_available,
                    colnames(assay_matrix)
                )
            )
        )
    }

    # calculate the pc coordinates
    assay_matrix <- assay_matrix[genes_for_pca, ]
    assay_matrix <- t(assay_matrix) %*% (loadings_pca %>% as.matrix)

    SummarizedExperiment::colData(sum_exp)[, colnames(assay_matrix)] <-
        assay_matrix

    sum_exp

}
