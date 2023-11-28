# In this file we create the workflow for the PCA fits used in both
# the small and big models. It uses the rds files saved during the
# development of the tool

path_to_mol_land <- file.path(
    "~/Documents/manuscripts",
    "20220721_molecular_landscape",
    "results",
    "rds_files"
)

pca_fit <- readRDS(file.path(
    path_to_mol_land,
    "pca_merging/pca_fit.rds"
))

pca_fit_all <- readRDS(file.path(
    path_to_mol_land,
    "trying/pca_fit_all_genes_wo_outliers.rds"
))

stable_genes_list <- readRDS(file.path(
    path_to_mol_land,
    "pca_merging/stable_genes.rds"
))

most_variable_genes_small <- setdiff(
    rownames(pca_fit$loadings),
    stable_genes
)

most_variable_genes_big <- setdiff(
    rownames(pca_fit_all$loadings),
    stable_genes
)

# we also get the df_pca and save only the first four components
# of all the three cohorts as these are the ones that will be used
# when plotting only
df_pca_all_cohorts <- readRDS(file.path(
    path_to_mol_land,
    "pca_merging/df_pca_coordinates.rds"
)) %>%
    dplyr::select(PC1, PC2, PC3, PC4, pam50, cohort, er_status)

# xz 78mb <- this one is the slower to save
# bzip2 82.4mb
# gzip 81.4mb
usethis::use_data(
    pca_fit, pca_fit_all,
    stable_genes_list,
    most_variable_genes_big,
    most_variable_genes_small,
    df_pca_all_cohorts,
    internal = TRUE,
    overwrite = TRUE,
    compress = "xz"
)

