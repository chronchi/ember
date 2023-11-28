get_colors_pam50 <- function(df){

    # ggplot2 default palette scheme
    n <- 6
    hues <- seq(15, 375, length = n + 1)
    colors_pam50 <- hcl(h = hues, l = 65, c = 100)[1:n]
    names(colors_pam50) <- c(
        "basal",
        "her2",
        "lumb",
        "luma",
        "normal",
        "claudin-low"
    )

    # we only use the colors available in the dataframe, otherwise
    # it shows all categories when it is not actually necessary
    selected_colors <- intersect(
        names(colors_pam50),
        unique(df$pam50)
    )
    colors_pam50[selected_colors]
}
