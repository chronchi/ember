#' Get gene sets formatted to apply the scores
#'
#' @details
#' This function is a wrapper on the function `msigdbr` from the package
#' msigdbr. The arguments are the ones that would be used for the mentioned
#' function
#'
#' @importFrom msigdbr msigdbr
#' @importFrom dplyr %>%
#' @export
get_gene_sets <- function(...){
    msigdbr::msigdbr(...) %>% sapply(
        .$gs_name %>% unique,
        function(x, gene_sets){
            gene_sets %>%
                dplyr::filter(gs_name == x) %>%
                dplyr::pull(gene_symbol)
        },
        gene_sets = .,
        USE.NAMES = TRUE,
        simplify = FALSE
    )
}
