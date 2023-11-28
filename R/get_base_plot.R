#' Plot PCA embedding of the 3 big cohorts to serve as a base plot
#'
#' @param x A string. PC for x-axis
#' @param y A string. PC for y-axis
#' @param color A string. Which column to use when coloring
#'     the dots.
#' @param size_dots An integer. Size of the points in the plot
#' @param alpha_val A number. Alpha value of the points
#' @param size_legend An integer. Size of the color legend
#' @param base_size An integer. Parameter for ggplot2::theme_bw
#' @return A plot of the TCGA, SCANB and METABRIC embedding
#' @import ggplot2
#' @importFrom dplyr %>%
#' @export
get_base_plot <- function(
    x = "PC3",
    y = "PC4",
    color = "pam50",
    size_dots = 2,
    alpha_val = 0.1,
    size_legend = 4,
    base_size = 10
){

    p <- df_pca_all_cohorts %>%
        dplyr::filter(
            pam50 %in% c("luma", "lumb", "basal", "normal", "her2")
        ) %>%
        ggplot2::ggplot(aes(x = !!sym(x), y = !!sym(y), color = !!sym(color))) +
        ggplot2::geom_point(size = size_dots, alpha = alpha_val) +
        ggplot2::scale_alpha(guide = 'none')

    if (color == "pam50"){
        p <- p + ggplot2::labs(
            color = "PAM50"
        ) +
        ggplot2::scale_color_manual(
           values = get_colors_pam50(df_pca_all_cohorts)
        )
    }

    p + ggplot2::guides(
        colour = ggplot2::guide_legend(
            override.aes = list(size = size_legend, alpha = 1)
        )
    ) +
        ggplot2::theme_bw(base_size = base_size)
}
