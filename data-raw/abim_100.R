library(dplyr)
library(SummarizedExperiment)

path_to_mol_land <- file.path(
    "~/Documents/manuscripts",
    "20220721_molecular_landscape",
    "data"
)

load(file.path(
    path_to_mol_land,
    "scanb_2022/ABiM.100.mymatrix.Rdata"
))

sheets <- "ABiM.100"
clin_data <- suppressWarnings({readxl::read_excel(
    file.path(
        "~/Documents/brisken-lab/BioRepo",
        "Data/20230125_scanb/",
        "Supplementary Data Table 1 - 2023-01-13.xlsx"
    ),
    sheet = sheets,
    progress = TRUE
)}) %>%
    dplyr::mutate(
        sample_name = GEX.assay,
        tumor_size = T.size,
        age = Age,
        pam50 = ProsignaFT.Subtype,
        node_status = dplyr::case_when(
            LN == 0 ~ "neg",
            LN == 1 ~ "pos",
            TRUE ~ NA
        )
    )

abim_100 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
        fpkm = ABiM.100.mymatrix,
        log2fpkm = log2(ABiM.100.mymatrix + 1)
    ),
    colData = clin_data %>%
        data.frame %>%
        `rownames<-`(clin_data$GEX.assay)
)

# load gene info
load(file.path(
    path_to_mol_land,
    "scanb_2022/Gene.ID.ann.Rdata"
))
gene_id_ann <- Gene.ID.ann

# convert ensembl genes to hugo IDs
rownames(abim_100) <- gene_id_ann[rownames(abim_100), "Gene.Name"]

# check what are the duplicated genes to see if we can safely drop them
dup_genes <- rownames(abim_100)[duplicated(rownames(abim_100))] %>%
    table

# Due to the lower amount of genes and the fact each gene has only one copy
# duplicated, we select the first one.
abim_100 <- abim_100[!duplicated(rownames(abim_100)), ]

usethis::use_data(abim_100, overwrite = TRUE, compress = "xz")

