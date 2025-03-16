#!user/bin/Rscript

#=======================================#
#|                                     |#
#|    scRNA-Seq results - analysis     |#
#|                                     |#
#=======================================#



#################################################################  Setting up, formatting and  ##################################################################
                                                                #  cleaning the relevant data  #


# Defining the required R packages for the full analysis
cran_packages <- c("tidyverse", "stringi", "BiocManager", "this.path", 
              "scales", "RCurl", "cowplot", "rebus", "ggsci",
              "progress", "metap", "doSNOW", "foreach", "scCustomize",
              "Matrix", "ggpubr", "R.utils", "devtools", "remotes", "RMTstat", "rstudioapi",
              "aplot", "circlize", "data.table", "devtools", "doParallel", "doRNG", "ggfun", "gghalves",
              "ggplotify", "ggridges", "irlba", "magrittr", "Matrix", "msigdbr", "pagoda2", "plyr", "pointr", 
              "RcppML", "reshape2", "reticulate", "rlang", "RMTstat", "RobustRankAggreg", "roxygen2", 
              "SeuratObject", "tidyselect", "tidytree", "VAM")

bioconductor_packages <- c("Seurat", "glmGamPoi", "multtest", "biomaRt", "AnnotationDbi",
              "EnsDb.Hsapiens.v86", "EnhancedVolcano", "graphite", "netgsa",
              "org.Hs.eg.db", "fgsea", "clusterProfiler", "SPIA", "PCAtools",
              "AUCell", "BiocParallel", "ComplexHeatmap", "decoupleR", "fgsea", "ggtree",
              "GSEABase", "GSVA", "Nebulosa", "scde", "singscore", "SummarizedExperiment",
              "UCell", "viper","sparseMatrixStats")

github_packages <- c("SeuratWrappers", "monocle3", "presto", "irGSEA")
github_install_names <- c("satijalab/seurat-wrappers", "cole-trapnell-lab/monocle3", "immunogenomics/presto", "chuiqin/irGSEA")                            


# This function will check if the required packages are installed and if not it installs them
# NOTE: by default the function will only check for CRAN packages, so if you want to install bioconductor
# or github packages you will have to feed those in separately
install_required_packages = function (CRAN_package_lst = NULL, BioC_package_lst = NULL, github_package_lst = NULL, github_install_param = NULL) {
    if (is.null(CRAN_package_lst) == FALSE) {
        for (cran_item in CRAN_package_lst) {
            if (is.element(cran_item, installed.packages()) == FALSE) {
                message("The package: ", cran_item, " is not installed. Installing package...")
                install.packages(cran_item)
            } else {
                message("The package: ", cran_item, " is installed. Load the package using the library function.")
            }
        }

    } else {
        message("No CRAN packages were requested.")
    }
      
    if (is.null(BioC_package_lst) == FALSE) {
        for (bioc_item in BioC_package_lst) {
            if (is.element(bioc_item, installed.packages()) == FALSE) {
                message("The package: ", bioc_item, " is not installed. Installing package...")
                BiocManager::install(bioc_item)
            } else {
                message("The package: ", bioc_item, " is installed. Load the package using the library function.")
            }
        }

    } else {
        message("\n", "No Bioconductor packages were requested.")

    }

    if (is.null(github_package_lst) == FALSE) {
        for (github_item in github_package_lst) {
            if (is.element(github_item, installed.packages()) == FALSE) {
                message("The package: ", github_item, " is not installed. Installing package...")
                for (install_param in github_install_param) {
                    devtools::install_github(install_param)
                }
            } else {
                message("The package: ", github_item, " is installed. Load the package using the library function.")
            }
        }

    } else {
        message("\n", "No github packages were requested.")
    }  
}
install_required_packages(CRAN_package_lst = cran_packages, BioC_package_lst = bioconductor_packages, github_package_lst = github_packages, github_install_param = github_install_names)


# As for a reproducible Seurat workflow
if (packageVersion("Seurat") != "4.3.0.1") {
    
    remotes::install_version("Seurat", version = "4.3.0.1")
}


# Used library packages
packages <- c("tidyverse", "stringi", "BiocManager", "Seurat", "biomaRt", "AnnotationDbi", "scales", "EnsDb.Hsapiens.v86", "this.path",
                "RCurl", "cowplot", "rebus", "ggsci", "EnhancedVolcano", "progress", "metap", "doSNOW", "foreach", "scCustomize", "graphite",
                    "org.Hs.eg.db", "fgsea", "clusterProfiler", "ggpubr", "SeuratWrappers", "multtest", "PCAtools", "RMTstat", "irGSEA")

lapply(packages, library, character.only = TRUE)


# Setting wd + path and listing the read output tables
local_dir <- this.path::this.dir()
setwd(local_dir)
bam_files <- paste0(getwd(), "/", "Re-seq_bam")


# Define the script version for output files, for easier version control, and make a dedicated folder in the WD
define_version_create_folder = function() {
    #define the current script version
    #NOTE: the version numbers must match the defined pattern. I  updated the function to accept sub-versions below two digits.
    if (nchar(stringr::str_extract(rstudioapi::getActiveDocumentContext()$path, "_V.*")) == 7) {
        script_version <- stringr::str_extract(rstudioapi::getActiveDocumentContext()$path, "_V" %R% rebus::one_or_more(ASCII_ALNUM) %R% one_or_more(PUNCT) %R%
                                      one_or_more(ASCII_ALNUM))
    } else {
        script_version <- stringr::str_extract(rstudioapi::getActiveDocumentContext()$path, "_V" %R% rebus::one_or_more(ASCII_ALNUM) %R% one_or_more(PUNCT) %R%
                                      one_or_more(ASCII_ALNUM) %R% one_or_more(PUNCT) %R% one_or_more(ASCII_ALNUM))
    }
    
    #check if the current version already has a folder or not and if not it creates one
    if (dir.exists(paths = paste0("Analysis", script_version, "/", "Plots")) == TRUE) {
        message("The current script directory is already present. No new directory will be created.")
    } else {
        message("The current script vesrsion has no directorey yet.", "\n", "Creating one now using the path:", paste0(getwd(), "Analysis_", script_version, "/", "Plots"))
        dir.create(path = paste0("Analysis", script_version, "/", "Plots"), recursive = TRUE)
        dir.create(path = paste0("Analysis", script_version, "/", "GSEA_results"), recursive = TRUE)
        dir.create(path = paste0("Analysis", script_version, "/", "GSEA_results", "/", "Plots"), recursive = TRUE)
        
    }
    assign("script_version", script_version, envir = .GlobalEnv)
    assign("analysis_dir", paste0("Analysis", script_version), envir = .GlobalEnv)
}
define_version_create_folder()


# This function will read and compile all .bam files into a single df, and names the columns according to the sample names
# note:the function automatically creates the GeneCountTable object in global environment, no assignment needed
read_RPG.tabs = function (path = getwd(), filenames, test_set = FALSE) {
    RPG.tab_lst <- list.files(path = path, 
                              pattern = filenames, 
                              all.files = TRUE,
                              recursive = TRUE) #super neat function which goes into sub-folders to find the pattern designated file
    
    if (test_set == TRUE) {
        test_lst <- list()
        for (e in 1:5) {
            test_lst[[e]] <- read.table(paste0(path, "/", RPG.tab_lst[e]),
                                        stringsAsFactors = FALSE, 
                                        header = FALSE)
        } 
    } else {
        message("No test_set was requested, compiling GeneCountTable")
    }
    
    #this function will allow to visualize a progression bar
    progressBar <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                    total = length(RPG.tab_lst),
                                    complete = "=",
                                    incomplete = "-",
                                    current = ">",
                                    clear = FALSE,
                                    width = 100)
    
    
        for (e in seq_along(RPG.tab_lst)) {
        progressBar$tick()
            
        tmp <- read.table(paste0(path, "/", RPG.tab_lst[e]),
                          stringsAsFactors = FALSE,
                          header = FALSE)
        tmp <- tmp[grep("^ENS", tmp[, 1]), c(1:2)]
        
        if (e == 1) {
            GeneCountTable <- tmp
        } else {
            GeneCountTable <- cbind(GeneCountTable, tmp[, 2])
        }
        
    }
    
        
    samples <- stringr::str_remove(RPG.tab_lst, "^lane1")
    samples <- stringr::str_remove(samples, "_.*")
    read_df_names <- c("geneID", samples)
    colnames(GeneCountTable) <- (read_df_names)
    
    if (test_set == TRUE)  {
        assign("test_lst", test_lst, envir = .GlobalEnv)
        assign("GeneCountTable", GeneCountTable, envir = .GlobalEnv)
        message("test_set was requested, compiling GeneCountTable and test_set")
    } else {
        assign("GeneCountTable", GeneCountTable, envir = .GlobalEnv)
    }
    
}
read_RPG.tabs(path = bam_files, filenames = "ReadsPerGene.out.tab", test_set = FALSE)


# Writing the compiled GeneCountTable for future use, so it can be loaded directly
write_csv(GeneCountTable, paste0(analysis_dir, "/", "Unstranded_GeneCountTable_re-seq", script_version, ".csv"))


# From this point on we can work with the clean gene count table
if (!exists("GeneCountTable", where = .GlobalEnv) && file.exists(paste0(analysis_dir, "/", "Unstranded_GeneCountTable_re-seq", script_version, ".csv"))) {
    GeneCountTable <- read.csv(file = paste0(analysis_dir, "/", "Unstranded_GeneCountTable_re-seq", script_version, ".csv"))
}
head(GeneCountTable)


# Saving the geneIDs into a separate string for gene name conversion
geneIDs <- GeneCountTable$geneID
rownames(GeneCountTable) <- geneIDs
head(GeneCountTable)


# Load the unified gene list if available, otherwise build it from scratch
if (file.exists(paste0("Additional_data_files", "/", "Unified_gene_names_table", script_version, ".csv")) && !exists("unified_gene_names_table", where = .GlobalEnv)) {
    unif_gene_names <- read.csv(file = paste0("Additional_data_files", "/", "Unified_gene_names_table", script_version, ".csv"))
    message("The unified gene names table was found and loaded as: ", deparse(substitute(unified_gene_names_table)))
} else {
    
    # Retrieving gene names using the ensembl IDs
    biomaRt::listEnsembl() #lists the available biomart keynames, here I need "gene"
    ensembl <- biomaRt::useEnsembl(biomart = "genes") #creates an ensembl object needed for the database connection
    datasets <- biomaRt::listDatasets(ensembl) #lists the organism based datasest from which we need to chose one for the db connection
    ensembl_connection <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl") #this connects the database and the dataset for the search query
    attributes <- biomaRt::listAttributes(ensembl_connection) #needed for the query, describes the desired info like gene names
    filters <- biomaRt::listFilters(ensembl_connection) # needed for the query, defined by upon which data we search for(?) like here the ensemblIds
    
    gene_names <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                        filters = "ensembl_gene_id",
                        values = geneIDs,
                        mart = ensembl_connection)
    
    
    
    # If you run into missing IDs you can check which ones are missing
    missingIDs <- geneIDs[geneIDs %in% gene_names$ensembl_gene_id == FALSE]
    missingNames <- gene_names[stringi::stri_isempty(gene_names$external_gene_name) == TRUE, 1]
    unif_missing <- c(missingIDs, missingNames)
    
    # An alternative method for gene id mapping if the direct query does not work properly
    # like shorter name list than id list...
    AnnotationDbi::keytypes(EnsDb.Hsapiens.v86)
    AnnotationDbi::columns(EnsDb.Hsapiens.v86)
    
    gene_names_2 <- AnnotationDbi::mapIds(x = EnsDb.Hsapiens.v86,
                           keys = unif_missing,
                           column = "SYMBOL",
                           keytype = "GENEID")
    gene_names_2 <- data.frame(names(gene_names_2), gene_names_2)
    rownames(gene_names_2) <- NULL
    colnames(gene_names_2) <- colnames(gene_names)
    clean_gene_names <- gene_names[stringi::stri_isempty(gene_names$external_gene_name) == FALSE, ]
    unif_gene_names <- rbind(clean_gene_names, gene_names_2)
    unif_gene_names <- dplyr::arrange(unif_gene_names, ensembl_gene_id)
    
    
    # Removing unnecessary variables
    rm(gene_names, gene_names_2, clean_gene_names, missingIDs, missingNames, unif_missing)
    
    
    # Save the unified gene names for later
    write_csv(x = unif_gene_names, file = paste0("Additional_data_files", "/", "Unified_gene_names_table", script_version, ".csv"))
    
}




## Loading the metadata.csv (only do this if you don't have the extended version yet)


# This if statement will check if the processes/extended metadata is present or not. If yes it will not run the following code chunk
if (file.exists(paste0("Additional_data_files", "/", "Re-seq_sample_extended_fixed_updated_metadata_with_markers_V3.csv"))) {
    message("The processed metadata file is already present in the working directory.",
    "\n", "No further processing will be done.") 
} else if (file.exists(paste0("Additional_data_files", "/", "Re-seq_sample_extended_fixed_updated_metadata_with_markers_V3.csv"))) {
    message("The extended metadata file is already present in the working directory.", 
    "\n", "No further processing will be done.")
} else {
    message("The process metadata is not present in the working directory.",
    "\n", "The metadata will be processed and saved now.")

    # The metadata will be loaded and processed by the following code chunk
    metad <- read.csv("Index_table_reseq.csv", sep = ";")
    glimpse(metad)
    metad <- metad[, 2:3]
    colnames(metad) <- c("sampleID", "RespGroup")
    metad$sampleID <- str_replace_all(metad$sampleID, "\\(A\\)", "A")
    patientID <- str_remove(metad$sampleID, "P" %R% one_or_more(ASCII_ALNUM))
    metad <- mutate(metad, patientID = patientID)
    metadata <- metad[order(metad$patientID), ]
    rm(metad)


    # Sub-setting sample codes for additional entries to the metadata df
    age <- if_else(str_detect(metadata$sampleID, "^A"),
            true = str_sub(metadata$sampleID, 4, 5),
            false = str_sub(metadata$sampleID, 3, 4))

    collection <- if_else(str_detect(metadata$sampleID, "^A"),
                          true = str_sub(metadata$sampleID, 6, 10),
                          false = str_sub(metadata$sampleID, 5, 9))

    cycle <- if_else(str_detect(metadata$sampleID, "^A"),
                          true = str_sub(metadata$sampleID, 11, 11),
                          false = str_sub(metadata$sampleID, 10, 10))

    plate <- if_else(str_detect(metadata$sampleID, "^A"),
                     true = str_sub(metadata$sampleID, 12, 13),
                     false = str_sub(metadata$sampleID, 11, 12))

    well <- if_else(str_detect(metadata$sampleID, "^A"),
                    true = str_sub(metadata$sampleID, 14, 15),
                    false = str_sub(metadata$sampleID, 13, 14))

    treatment_type <- if_else(str_detect(metadata$sampleID, "^A"),
                              true = "AL",
                              false = "L")


    # Adding additional information to the metadata and formatting it for a seurat object (has to have the samples as rownames, it has to be a dataframe and not a tibble)
    metadata <- mutate(metadata, Age = age,
                       Collection = collection,
                       Cycle = cycle,
                       Treatment = treatment_type,
                       Plate = plate,
                       Well = well)
    rm(age, collection, cycle, plate, well, treatment_type)
    write_csv(metadata, "Re-seq_sample_extended_metadata.csv")
    # IMPORTANT NOTE: in this case the noted treatment cycles and the real cycles and treatments
    # do not correspond, so the metadata has to be corrected accordingly. This was done separately
    # and the corrected data will be loaded at the next section

}

#################################################################   End Section  ##################################################################






#################################################################   Creating the Seurat object and  ##################################################################
                                                                #   prepping for the data analysis  #


# NOTE (for this I need a count table and a prepped metadata dataframe)




## Continuing with the full data


# Load the metadata table if available and not already loaded
read_metadata <- function(path = getwd(), filename, sep) {
    if (exists("meta", envir = .GlobalEnv) == TRUE) {
        stop("The metadata is already loaded as the following object: meta")
    } else if (exists("metad", envir = .GlobalEnv) == TRUE) {
        stop("The metadata is already loaded as the following object: metad")
    } else {
        meta <- read.csv(paste0(path, "/", filename), sep = sep)
    }
    assign("meta", meta, envir = .GlobalEnv)
    message("The metadata is loaded as the following object: meta")
}
read_metadata(filename = paste0("Additional_data_files", "/", "Re-seq_sample_extended_fixed_updated_metadata_with_markers_V3.csv"), sep = ";")
glimpse(meta)
head(meta)


# NOTE:these steps seem to be necessary for the seurat object to utilize the metadata properly
rownames(meta) <- meta$sampleID
meta <- mutate(meta, sampleID = NULL)
head(meta)

GCT <- mutate(GeneCountTable, geneID = NULL)
CTC_reseq.obj <- CreateSeuratObject(counts = GCT,
                                    project = "CTC_reseq",
                                    assay = "RNA_seq",
                                    meta.data = meta,
                                    min.cells = 0)
glimpse(CTC_reseq.obj@meta.data)


# Save the created new Seurat object
saveRDS(object = CTC_reseq.obj, file = paste0(getwd(), "/", analysis_dir, "/", "Original_CTC_Seurat_object", script_version, ".rds"))


# Following the creation of the seurat object remove some not needed objects to save memory
rm(meta, GCT)

#################################################################   End Section  ##################################################################






#################################################################   Quality control  ##################################################################
#                    #


# Load the original CTC Seurat object if not loaded yet
if (exists("CTC_reseq.obj", where = .GlobalEnv)) {
  message("The CTC_reseq.obj is already present in the .GlobalEnv")
} else {
  message("Loading the CTC_reseq.obj")
  CTC_reseq.obj <- readRDS(file = paste0(getwd(), "/", analysis_dir, "/", "Original_CTC_Seurat_object", script_version, ".rds"))
  
}


# Load mitochondrial geneIDs and names (got it from ENSEMBL) for the QC
MT_genes <- read_csv(paste0("Additional_data_files", "/", "ENSEMBL_MT_genes.csv"))
mito.genes <- rownames(CTC_reseq.obj)[rownames(CTC_reseq.obj) %in% MT_genes$`Gene stable ID`]
print(mito.genes)
rm(mito.genes)


# Alternatively one can load it from the scCustomize package
#ensembl_mito_id


# Add mitochondrial gene percentage and genes per nCount percentages (log10 transformed for better visibility) to our
# seurat.obj
CTC_reseq.obj[["mt.percent"]] <- Seurat::PercentageFeatureSet(object = CTC_reseq.obj, features = MT_genes$`Gene stable ID`)
CTC_reseq.obj[["nFeature.per.nCount"]] <- log10(CTC_reseq.obj$nFeature_RNA_seq) / log10(CTC_reseq.obj$nCount_RNA_seq)


# At this point I will separate the metadata from the seurat object
CTC_meta.data <- CTC_reseq.obj@meta.data
CTC_meta.data_clean <- CTC_meta.data[!is.na(CTC_meta.data$mt.percent), ]
CTC_meta.data_clean <- mutate(CTC_meta.data_clean, Cell_names = rownames(CTC_meta.data_clean))
rm(CTC_meta.data)


# This function simply runs a block of code, plotting and saving the QC metrics one by one (the only purpose of this is conciseness)
generate_pre_qc_plots = function(make_plots = TRUE) {
  if (make_plots == TRUE) {
    message("Pre-QC plots were requested. Running the plotting code block.")
    
    
    # Visualize the cell numbers / resistance group and save it
    cn_p <- ggplot2::ggplot(CTC_meta.data_clean,
                            aes(x = RespGroup, fill = RespGroup)) +
      ggplot2::geom_bar(color = "black") +
      ggplot2::geom_text(aes(label = after_stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + #thx Bing
      ggsci::scale_fill_tron() +
      ggplot2::labs(x = "Responder groups", y = "Cell numbers", fill = "Response groups", color = "Response groups") +
      #ggplot2::ggtitle("Cell numbers after sequencing") +
      ggplot2::guides(fill = guide_legend(title = "Response groups")) +
      ggplot2::theme_classic() +
      ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
    print(cn_p)
    
    ggplot2::ggsave(filename = paste0("Cell numbers", script_version, ".png"), plot = cn_p,
                    device = "png", path = paste0(analysis_dir, "/", "Plots/"),
                    width = 3000, height = 1800, units = "px", dpi = 320)
    rm(cn_p)
    
    
    # Visualizing the mt.percent across cells using density plot and save it
    mt_pd <- ggplot2::ggplot(CTC_meta.data_clean,
                             aes(x = mt.percent, fill = RespGroup)) +
      ggplot2::geom_density(alpha = 0.2) +
      ggsci::scale_fill_tron() +
      ggplot2::labs(x = "MT gene expression percentage", y = "MT gene percent distribution") +
      ggplot2::theme_classic() +
      ggplot2::geom_vline(xintercept = 75, color = "red") +
      #ggplot2::ggtitle("Mitochondrial gene expression") +
      ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
    print(mt_pd)
    
    ggplot2::ggsave(filename = paste0("Mitochondrial gene expression", script_version, ".png"), plot = mt_pd,
                    device = "png", path = paste0(analysis_dir, "/", "Plots/"),
                    width = 3000, height = 1800, units = "px", dpi = 320)
    rm(mt_pd)
    
    
    # Visualizing the nCount numbers/treatment groups and save it
    # The scales package helps to remove scientific notation form the axises
    nCount_p <- ggplot2::ggplot(CTC_meta.data_clean,
                                aes(x = nCount_RNA_seq, fill = RespGroup)) +
      ggplot2::geom_density(alpha = 0.2) +
      ggsci::scale_fill_tron() +
      ggplot2::labs(x = "nCounts", y = expression("log"[10]* " Count density"), fill = "Response group") +
      ggplot2::scale_x_log10(labels = scales::comma) + 
      ggplot2::theme_classic() +
      ggplot2::geom_vline(xintercept = 4000, color = "red") +
      #ggplot2::ggtitle("mRNA count distribution") +
      ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
    print(nCount_p)
    
    ggplot2::ggsave(filename = paste0("mRNA count distribution", script_version, ".png"), plot = nCount_p,
                    device = "png", path = paste0(analysis_dir, "/", "Plots/"),
                    width = 3000, height = 1800, units = "px", dpi = 320)
    rm(nCount_p)
    
    
    # Visualizing the nFeature numbers/treatment groups and save it
    nFeature_p <- ggplot2::ggplot(CTC_meta.data_clean,
                                  aes(x = nFeature_RNA_seq, fill = RespGroup)) +
      ggplot2::geom_density(alpha = 0.2) +
      ggsci::scale_fill_tron() +
      ggplot2::labs(x = "nFeatures", y = expression("log"[10]* " Feature density"), fill = "Response group") +
      ggplot2::scale_x_log10(labels = scales::comma) + 
      ggplot2::theme_classic() +
      #ggplot2::ggtitle("Gene count distribution") +
      ggplot2::geom_vline(xintercept = 1000, color = "red") +
      ggplot2::geom_vline(xintercept = 7800, color = "orange") +
      ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
    print(nFeature_p)
    
    ggplot2::ggsave(filename = paste0("Gene count distribution", script_version, ".png"), plot = nFeature_p,
                    device = "png", path = paste0(analysis_dir, "/", "Plots/"),
                    width = 3000, height = 1800, units = "px", dpi = 320)
    rm(nFeature_p)
    
    
    # Cell complexity plotting by checking the nFeature/nCount ratio (the higher the nFeature/nCount
    # the more complex the cells are) and save it
    nFpernC_p <-ggplot2::ggplot(CTC_meta.data_clean,
                                aes(x = nFeature.per.nCount, fill = RespGroup)) +
      ggplot2::geom_density(alpha = 0.2) +
      ggsci::scale_fill_tron() +
      ggplot2::labs(x = "nFeatures/nCounts", y = expression("log" [10]* " density"), fill = "Response group") +
      ggplot2::theme_classic() +
      #ggplot2::ggtitle("Cell complexity distribution (nFeature/nCount)") +
      ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
    print(nFpernC_p)
    
    ggplot2::ggsave(filename = paste0("Cell complexity distribution", script_version, ".png"), plot = nFpernC_p,
                    device = "png", path = paste0(analysis_dir, "/", "Plots/"),
                    width = 3000, height = 1800, units = "px", dpi = 320)
    rm(nFpernC_p)
    
    
    # nCount vs nGene (a high ratio is better) and save it
    nCvnF_p <- ggplot2::ggplot(data = CTC_meta.data_clean,
                               aes(x = nCount_RNA_seq, y = nFeature_RNA_seq, color = mt.percent)) +
      ggplot2::geom_point(size = 3, alpha = 0.5) + 
      ggplot2::scale_color_gradient(low = "grey90", high = "red", limits = c(0, 100)) +
      ggsci::scale_fill_tron() +
      ggplot2::labs(x = "mRNA counts", y = "Gene counts", color = "MT gene expression \n percentage") +
      #stat_smooth(method = lm) +
      ggplot2::scale_x_log10(labels = scales::comma) +
      ggplot2::scale_y_log10(labels = scales::comma) +
      ggplot2::geom_hline(yintercept = 1000, linetype = "dashed") +
      ggplot2::geom_hline(yintercept = 7800, linetype = "dashed") +
      ggplot2::geom_vline(xintercept = 4000, linetype = "dashed") +
      ggplot2::theme_classic() +
      #facet_wrap(facets = ~RespGroup)
      #ggplot2::ggtitle("nCount vs nFeature - QC plot") +
      ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
    print(nCvnF_p)
    
    ggplot2::ggsave(filename = paste0("nCount vs nFeature - QC plot", script_version, ".png"), plot = nCvnF_p,
                    device = "png", path = paste0(analysis_dir, "/", "Plots/"),
                    width = 3800, height = 1900, units = "px", dpi = 320)
    rm(nCvnF_p)
    
    
    # End of run message
    message("The requested pre-QC plots were printed and saved into the: ", paste0(analysis_dir, "/", "Plots/"), " directory.")
    
  } else {
    message("Pre-QC plots were not requested. Skipping the plotting code block.")
  }
}
generate_pre_qc_plots()


# Clean objects
rm(CTC_meta.data_clean)

#################################################################   End Section  ##################################################################






#################################################################   Data filtering  ##################################################################
#                   #




## Filtering setup


# I set the lower cutoff for nFeature to 1000 genes and the upper cutoff to 7800.
# Later analysis showed that the population above this cutoff clusters together,
# probably the sign of multiple cells in one well, so the cutoff seems justified.




## Filtering


# Filter the Seurat object according to the above decided parameters
filt_CTC.obj <- subset(CTC_reseq.obj, 
                       subset = nCount_RNA_seq >= 4000 & nFeature_RNA_seq >= 1000 &
                         nFeature_RNA_seq <= 7800 & mt.percent < 75)


# Save the filtered Seurat object for later use
saveRDS(object = filt_CTC.obj, file = paste0(getwd(), "/", analysis_dir, "/", "Filtered_CTC_Seurat_object", script_version, ".rds"))


# Load the filtered Seurat object if not loaded yet
if (exists("filt_CTC.obj", where = .GlobalEnv)) {
  message("The filt_CTC.obj is already present in the .GlobalEnv")
} else {
  message("Loading the filt_CTC.obj")
  filt_CTC.obj <- readRDS(file = paste0(getwd(), "/", analysis_dir, "/", "Filtered_CTC_Seurat_object", script_version, ".rds"))
  
}




## Re-affirming the QC metrics on the trimmed, filtered data


# This function simply runs a block of code, plotting and saving the filtered data one by one (the only purpose of this is conciseness)
generate_filtered_data_plots = function(print_plots = TRUE, save_plots = TRUE) {
  # Defining the additional color palette
  cluster_colors <- ggsci::pal_tron(palette = "legacy")(7)
  
  # Checking the new cell numbers
  cn_p <- ggplot2::ggplot(filt_CTC.obj@meta.data,
                          aes(x = RespGroup, fill = RespGroup)) +
    ggplot2::geom_bar(color = "black") +
    ggplot2::geom_text(aes(label = after_stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) +
    ggsci::scale_fill_tron() +
    ggplot2::labs(x = "Response groups", y = "Cell numbers", fill = "Response groups") +
    #ggplot2::ggtitle("Filtered cell numbers after sequencing - Response groups") +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
  
  # Checking the new cell numbers
  dcn_p <- ggplot2::ggplot(filt_CTC.obj@meta.data,
                           aes(x = detailedResp, fill = detailedResp)) +
    ggplot2::geom_bar(color = "black") +
    ggplot2::geom_text(aes(label = after_stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) +
    ggplot2::scale_fill_manual(values = cluster_colors[c(5, 1, 2, 3)]) +
    ggplot2::labs(x = "Response groups", y = "Cell numbers", fill = "Response groups") +
    #ggplot2::ggtitle("Filtered cell numbers after sequencing - Detailed response") +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
  
  # Visualizing the mt.percent across cells using density plot and save it
  mt_pd <- ggplot2::ggplot(filt_CTC.obj@meta.data,
                           aes(x = mt.percent, fill = RespGroup)) +
    ggplot2::geom_density(alpha = 0.2) +
    ggsci::scale_fill_tron() +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Cell percentage", y = "MT gene percent distribution", fill = "Response groups") +
    #ggplot2::ggtitle("Filtered mitochondrial gene expression") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
  
  # Visualizing the nCount numbers/treatment groups and save it
  nCount_p <- ggplot2::ggplot(filt_CTC.obj@meta.data,
                              aes(x = nCount_RNA_seq, fill = RespGroup)) +
    ggplot2::geom_density(alpha = 0.2) +
    ggsci::scale_fill_tron() +
    ggplot2::labs(x = "nCounts", y = expression("log"[10]* " Count density"), fill = "Response group") +
    ggplot2::scale_x_log10(labels = scales::comma) + 
    ggplot2::theme_classic() +
    #geom_vline(xintercept = 2000, color = "red") +
    #ggplot2::ggtitle("Filtered mRNA count distribution") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
  
  # Visualizing the nFeature numbers/treatment groups and save it
  nFeature_p <- ggplot2::ggplot(filt_CTC.obj@meta.data,
                                aes(x = nFeature_RNA_seq, fill = RespGroup)) +
    ggplot2::geom_density(alpha = 0.2) +
    ggsci::scale_fill_tron() +
    ggplot2::labs(x = "nFeatures", y = expression("log"[10]* " Feature density"), fill = "Response group") +
    ggplot2::scale_x_log10(labels = scales::comma) + 
    ggplot2::theme_classic() +
    #ggplot2::ggtitle("Filtered gene count distribution") +
    #geom_vline(xintercept = 1000, color = "red") +
    #geom_vline(xintercept = 6500, color = "orange") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
  
  # Cell complexity plotting by checking the nFeature/nCount ratio (the higher the nFeature/nCount
  # the more complex the cells are) and save it
  nFpernC_p <-ggplot2::ggplot(filt_CTC.obj@meta.data,
                              aes(x = nFeature.per.nCount, fill = RespGroup)) +
    ggplot2::geom_density(alpha = 0.2) +
    ggsci::scale_fill_tron() +
    ggplot2::labs(x = "nFeatures/nCounts", y = expression("log"[10]* " density"), fill = "Response group") +
    ggplot2::theme_classic() +
    #ggplot2::ggtitle("Cell complexity distribution (nFeture/nCount)") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
  
  # nCount vs nGene (a high ratio is better) and save it
  nCvnF_p <- ggplot2::ggplot(data = filt_CTC.obj@meta.data,
                             aes(x = nCount_RNA_seq, y = nFeature_RNA_seq, color = mt.percent)) +
    ggplot2::geom_point(size = 3, alpha = 0.5) + 
    ggplot2::scale_color_gradient(low = "grey90", high = "red", limits = c(0, 100)) +
    #stat_smooth(method = lm) +
    ggplot2::labs(x = "mRNA counts", y = "Gene counts", color = "MT gene expression \n percentage") +
    ggplot2::scale_x_log10(labels = scales::comma) +
    ggplot2::scale_y_log10(labels = scales::comma) +
    #geom_hline(yintercept = 1000, linetype = "dashed") +
    #geom_hline(yintercept = 6500, linetype = "dashed") +
    #geom_vline(xintercept = 4000, linetype = "dashed") +
    ggplot2::theme_classic() +
    #facet_wrap(facets = ~RespGroup)
    #ggplot2::ggtitle("Filtered nCount vs nFeature - QC plot") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
  
  
  # This if statement controls if the plots are printed
  if (print_plots == TRUE) {
    
    message("Filtered data plots were requested. Running the plotting code block.")
    
    
    
    print(cn_p)
    print(dcn_p)
    print(mt_pd)
    print(nCount_p)
    print(nFeature_p)
    print(nFpernC_p)
    print(nCvnF_p)
    
    
    
    
    
    # End of run message
    message("The requested filtered data plots were printed as requested.")
    
  } else {
    
    message("The plots were not printed as it was not requested")
    
  }
  
  
  # This if statement control if the plots will be saved
  if (save_plots == TRUE) {
    
    ggplot2::ggsave(filename = paste0("Filtered cell numbers", script_version, ".png"), plot = cn_p,
                    device = "png", path = paste0(analysis_dir, "/", "Plots/"),
                    width = 3000, height = 1800, units = "px", dpi = 320)
    
    ggplot2::ggsave(filename = paste0("Filtered cell numbers - detailed response", script_version, ".png"), plot = dcn_p,
                    device = "png", path = paste0(analysis_dir, "/", "Plots/"),
                    width = 3000, height = 1800, units = "px", dpi = 320)
    
    ggplot2::ggsave(filename = paste0("Filtered mitochondrial gene expression", script_version, ".png"), plot = mt_pd,
                    device = "png", path = paste0(analysis_dir, "/", "Plots/"),
                    width = 3000, height = 1800, units = "px", dpi = 320)
    
    ggplot2::ggsave(filename = paste0("Filtered mRNA count distribution", script_version, ".png"), plot = nCount_p,
                    device = "png", path = paste0(analysis_dir, "/", "Plots/"),
                    width = 3000, height = 1800, units = "px", dpi = 320)
    
    ggplot2::ggsave(filename = paste0("Filtered gene count distribution", script_version, ".png"), plot = nFeature_p,
                    device = "png", path = paste0(analysis_dir, "/", "Plots/"),
                    width = 3000, height = 1800, units = "px", dpi = 320)
    
    ggplot2::ggsave(filename = paste0("Filtered cell complexity distribution", script_version, ".png"), plot = nFpernC_p,
                    device = "png", path = paste0(analysis_dir, "/", "Plots/"), units = "px",
                    width = 3000, height = 1800, dpi = 320)
    
    ggplot2::ggsave(filename = paste0("Filtered nCount vs nFeature - QC plot", script_version, ".png"), plot = nCvnF_p,
                    device = "png", path = paste0(analysis_dir, "/", "Plots/"),
                    width = 3800, height = 1900, units = "px", dpi = 320)
    
    message("The plots were saved to the folder: ", paste0(analysis_dir, "/", "Plots/"))
    
    
  } else {
    
    message("WARNING: the plots are not saved as a save was not requested")
    
  }
  
}
generate_filtered_data_plots(print_plots = TRUE, save_plots = TRUE)


# Remove unnecessary variables
rm(GeneCountTable)


#################################################################   End Section  ##################################################################






#################################################################   Normalization and variance-  ##################################################################
                                                                #          stabilization         #




## First things first, we want to look at factors commonly influencing clustering behavior like cell cycle status or mt gene expression
## and see if they have enough effect on cell clustering or not in order to determine if want to regress them out or not


# Cell cycle
# First we check for cell cycle effects for this we do a rough nCount normalization with the
# NormalizeData() function (will divide the nCounts with the cell number and does a log10 transform)
ccPhase_CTC.obj <- Seurat::NormalizeData(filt_CTC.obj)


# This if else statement will check if the S and G2M cell cycle markers are present in the working directory and if yes are they loaded in the global Env.
# If present but not loaded it will load them, if not present and not loaded it will attempt to make them 
if ((file.exists(paste0("Additional_data_files", "/", "s_phase_IDs.csv")) && file.exists(paste0("Additional_data_files", "/", "g2m_phase_IDs.csv"))) && (!exists(x = "s.IDs", where = .GlobalEnv) & !exists(x = "g2m.IDs", where = .GlobalEnv))) {
  
  s.IDs <- read.csv(file = paste0("Additional_data_files", "/", "s_phase_IDs.csv"))
  g2m.IDs <- read.csv(file = paste0("Additional_data_files", "/", "g2m_phase_IDs.csv"))
  message("The cell cycle marker tables for S and G2M phases are loaded as the following variables: ", deparse(substitute(s.IDs)), " and " ,deparse(substitute(g2m.IDs)))
  
} else if (exists(x = "s.IDs", where = .GlobalEnv) & exists(x = "g2m.IDs", where = .GlobalEnv)) {
  
  message("The cell cycle marker tables for S and G2M phases are already loaded in the .GlobalEnv.")
  
  
} else {
  
  message("The cell cycle marker tables for S and G2M are not found in the working directory. Running the code block to manually create them:")
  
  # Load cell cycle markers
  # Cell cycle markers are available as part of the seurat package in gene name format
  s.genes <- cc.genes.updated.2019$s.genes
  g2m.genes <- cc.genes.updated.2019$g2m.genes
  
  
  # Retrieving gene names using the ensembl IDs
  biomaRt::listEnsembl() #lists the available biomart keynames, here I need "gene"
  ensembl <- biomaRt::useEnsembl(biomart = "genes") #creates an ensembl object needed for the database connection
  datasets <- biomaRt::listDatasets(ensembl) #lists the organism based datasest from which we need to chose one for the db connection
  ensembl_connection <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl") #this connects the database and the dataset for the search query
  attributes <- biomaRt::listAttributes(ensembl_connection) #needed for the query, describes the desired info like gene names
  filters <- biomaRt::listFilters(ensembl_connection) # needed for the query, defined by upon which data we search for(?) like here the ensemblIds
  
  
  # In this case I did not convert my ensemblIDs to gene names as the conversion missed
  # some IDs, so for the cell cycle identification and matching I have to convert 
  # the gene names into ensembl IDs like I tried before
  s.IDs <- getBM(attributes = c("ensembl_gene_id","external_gene_name"),
                 filters = "external_gene_name",
                 values = s.genes,
                 mart = ensembl_connection)
  g2m.IDs <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                   filters = "external_gene_name",
                   values = g2m.genes,
                   mart = ensembl_connection)
  
  
  # Write the resulting dataframes into the working directory
  write_csv(x = s.IDs, file = "s_phase_IDs.csv")
  write_csv(x = g2m.IDs, file = "g2m_phase_IDs.csv")
  
}


# Clean the workspace a bit by removing unnecessary objects
rm(ensembl, datasets, ensembl_connection, attributes, filters)


# Now I can look for the cell cycle (cc) genes in the seurat object and assign
# a cc score
ccPhase_CTC.obj <- Seurat::CellCycleScoring(ccPhase_CTC.obj,
                                            s.features = s.IDs$ensembl_gene_id,
                                            g2m.features = g2m.IDs$ensembl_gene_id,
                                            set.ident = FALSE)

# For the next step I need to reduce variation, for this I have to find the feature
# causing high variation (usually the very highly expressed genes) and scale the data
# seurat has built in functions for this
ccPhase_CTC.obj <- Seurat::FindVariableFeatures(ccPhase_CTC.obj,
                                                selection.method = "vst",
                                                nfeatures = 3000,
                                                verbose = TRUE)
ccPhase_CTC.obj <- Seurat::ScaleData(ccPhase_CTC.obj)


# Now I can do the PCA and plot the results
ccPhase_CTC.obj <- Seurat::RunPCA(ccPhase_CTC.obj, approx = FALSE, verbose = FALSE)


# Testing the significance of PCs with the Jackstraw method and plotting the PCs
# NOTE: this is not possible on an objects transformed by SCTransform!
test_object <- Seurat::JackStraw(object = ccPhase_CTC.obj, reduction = "pca", assay = "RNA_seq",
                                 dims = 30)
Seurat::ScoreJackStraw(test_object, dims = 1:30, do.plot = TRUE)
rm(test_object)


# Plotting and saving the result using the cell cycle phase overlay
CellCyc_split_p <- Seurat::DimPlot(ccPhase_CTC.obj,
                                   reduction = "pca",
                                   group.by = "Phase",
                                   split.by = "Phase",
                                  pt.size = 3,
                                  cols = ggsci::pal_tron(palette = "legacy", alpha = 0.5)(7))
CellCyc_split_p2 <- CellCyc_split_p +
  #ggsci::scale_fill_tron() +
  #ggsci::scale_color_tron() +
  ggplot2::labs(x = "PC 1", y = "PC 2", color = "Cell cycle phase", fill = "Cell cycle phase") +
  #ggplot2::ggtitle("Cell cycle phases (split)") +
  ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
print(CellCyc_split_p2)

ggplot2::ggsave(filename = paste0("Cell cycle phase clustering (split)", script_version, ".png"), CellCyc_split_p2,
                device = "png", path = paste0(analysis_dir, "/", "Plots/"),
                width = 3000, height = 1800, units = "px")
rm(CellCyc_split_p, CellCyc_split_p2)


CellCyc_p <- Seurat::DimPlot(ccPhase_CTC.obj, 
                             reduction = "pca",
                             group.by = "Phase",
                             pt.size = 3,
                             cols = ggsci::pal_tron(palette = "legacy", alpha = 0.5)(7))
CellCyc_p2 <- CellCyc_p +
  #ggsci::scale_fill_tron() +
  #ggsci::scale_color_tron() +
  ggplot2::labs(x = "PC 1", y = "PC 2", color = "Cell cycle phase", fill = "Cell cycle phase") +
  #ggplot2::ggtitle("Cell cycle phases") +
  ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 18))
print(CellCyc_p2)

ggplot2::ggsave(filename = paste0("Cell cycle phase clustering", script_version, ".png"), CellCyc_p2,
                device = "png", path = paste0(analysis_dir, "/", "Plots/"),
                width = 3000, height = 1800, units = "px")
rm(CellCyc_p, CellCyc_p2)


# Removing the ccPhase and other unnecessary objects
rm(ccPhase_CTC.obj, CTC_reseq.obj)

# In this particular case both the cell cycle and mitotic gene expression data scatters together, and there are no visible clusters
# based on the cell cycle phase or mitotic gene expression so they do not have to be regressed out

#################################################################   End Section  ##################################################################






#################################################################   Final normalization and regressing out  ##################################################################
                                                                #        sources of unwanted variation      #




## Normalization and clustering


# As for the final more accurate method of normalization SCTransform, which is estimating the variance of the raw filtered data,
# and identifying the most variable genes.

# You might load the already prepared SCT transformed seurat object if available in the working directory
if (file.exists(paste0(getwd(), "/", analysis_dir, "/", "SCTransformed_fully_prepared_CTC_Seurat_object", script_version, ".rds")) && !exists(x = "sct_CTC.obj", where = .GlobalEnv)) {
  sct_CTC.obj <- readRDS(file = paste0(getwd(), "/", analysis_dir, "/", "SCTransformed_fully_prepared_CTC_Seurat_object", script_version, ".rds"))
  message("The standard SCTransformed Seurat object is present in the working directory and will be loaded as: ", deparse(substitute(sct_CTC.obj)))
  
} else {
  
  # If not present, load the filtered Surat object for further processing
  if (!exists(x = "filt_CTC.obj", where = .GlobalEnv)) {
    filt_CTC.obj <- readRDS(file = paste0(getwd(), "/", analysis_dir, "/", "Filtered_CTC_Seurat_object", script_version, ".rds"))
    message("The filtered CTC Seurat object was loaded as: ", deparse(substitute(filt_CTC.obj)))
  }
  
  
  # Run SCTransform 
  sct_CTC.obj <- SCTransform(filt_CTC.obj, assay = "RNA_seq", method = "glmGamPoi")
  
  
  # Save the generated sCTransformed seurat object
  saveRDS(object = sct_CTC.obj, file = paste0(getwd(), "/", analysis_dir, "/", "SCTransformed_CTC_Seurat_object", script_version, ".rds"))
  
  
  # Just as Jonathan, I will also repeat the cell cycle and mt scoring so I can overlay this info later onto the clustering
  sct_CTC.obj <- CellCycleScoring(sct_CTC.obj, 
                                  s.features = s.IDs$ensembl_gene_id, 
                                  g2m.features = g2m.IDs$ensembl_gene_id, 
                                  set.ident = FALSE)
  
  
  # Check how many PCs are responsible for the most variation. This is important
  # for the downstream steps. One way to check this is to plot the PCs with an ElbowPlot (plotting and saving the result)
  sct_CTC.obj <- RunPCA(sct_CTC.obj, assay = "SCT", approx = FALSE, verbose = FALSE, npcs = 50) #approx = FALSE to run with normal vsd method
  PC_weight_p <- ElbowPlot(sct_CTC.obj, ndims = 50, reduction = "pca")
  PC_weight_p2 <- PC_weight_p +
    labs(x = "PCs") +
    ggtitle("PC weights") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 18))
  print(PC_weight_p2)
  
  ggplot2::ggsave(filename = paste0("Principal component weights", script_version, ".png"), PC_weight_p2,
                  device = "png", path = paste0(analysis_dir, "/", "Plots/"),
                  width = 1500, height = 1000, units = "px", dpi = 320)
  rm(PC_weight_p, PC_weight_p2)
  
  
  # Running dimensionality reduction
  sct_CTC.obj <- RunUMAP(sct_CTC.obj, dims = 1:20, reduction = "pca", verbose = FALSE)
  
  
  # Finding neighbors and clustering
  sct_CTC.obj <- FindNeighbors(sct_CTC.obj, reduction = "pca", dims = 1:20)
  sct_CTC.obj <- FindClusters(sct_CTC.obj, resolution = c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2))
  
  
  # Set various identities to check the clustering, and to determine the optimal clustering resolution
  Idents(object = sct_CTC.obj) <- "SCT_snn_res.0.2"
  Idents(sct_CTC.obj) <- "SCT_snn_res.0.4"
  Idents(sct_CTC.obj) <- "SCT_snn_res.0.6"
  Idents(sct_CTC.obj) <- "SCT_snn_res.0.8"
  Idents(sct_CTC.obj) <- "SCT_snn_res.1"
  Idents(sct_CTC.obj) <- "SCT_snn_res.1.2"
  Idents(sct_CTC.obj) <- "SCT_snn_res.1.4"
  Idents(sct_CTC.obj) <- "SCT_snn_res.1.6"
  Idents(sct_CTC.obj) <- "SCT_snn_res.1.8"
  Idents(sct_CTC.obj) <- "SCT_snn_res.2"
  
  
  # Visualize the various clustering resolutions
  UMAP_p <- DimPlot(sct_CTC.obj, reduction = "umap")
  print(UMAP_p)
  rm(UMAP_p)
  
  PCA_p <- PCAPlot(sct_CTC.obj, reduction = "pca")
  print(PCA_p)
  rm(PCA_p)
  
  
  # Re-run FindCluster with the preferred seeing
  sct_CTC.obj <- FindClusters(sct_CTC.obj, resolution = 0.6)
  
  
  # Saving the fully prepared CTC seurat object
  saveRDS(object = sct_CTC.obj, file = paste0(getwd(), "/", analysis_dir, "/", "SCTransformed_fully_prepared_CTC_Seurat_object", script_version, ".rds"))
  
  
  # Closing message
  message("The block ran succesfully. \n",
          "The newly prepared SCT transformed Seurat object and the SCT transformed, fully prepared Seurat object were saved into the version folder ",
          script_version, " as follows: \n",
          "SCTransformed_CTC_Seurat_object", script_version, ".rds \n",
          "SCTransformed_fully_prepared_CTC_Seurat_object", script_version, ".rds")
  
}


# Cluster visualization (PCA, UMAP and t-SNE plots)
cluster_colors <- pal_tron(palette = "legacy", alpha = 0.5)(7)
cluster_colors_2 <- pal_futurama(palette = "planetexpress")(12)
cluster_colors_3 <- pal_startrek(palette = "uniform")(7)
cluster_colors_4 <- pal_rickandmorty(palette = "schwifty")(12)


# This function will print and save the dimensionality reduction plots
plot_DR_cluster = function (sct_object, print_plots = TRUE, save_plots = TRUE) {
  
  cluster_colors <- pal_tron(palette = "legacy", alpha = 0.5)(7)
  cluster_colors_2 <- pal_futurama(palette = "planetexpress")(12)
  cluster_colors_3 <- pal_startrek(palette = "uniform")(7)
  cluster_colors_4 <- pal_rickandmorty(palette = "schwifty")(12)
  
  # Create the PCA plot
  PCA_p <- DimPlot(sct_object, 
    reduction = "pca",
    group.by = "seurat_clusters", 
    pt.size = 5,
    cols = ggsci::pal_tron(palette = "legacy", alpha = 0.5)(7)[c(6, 7)])
  PCA_p2 <- PCA_p +
    #scale_fill_manual(values = cluster_colors[c(6, 7)]) +
    #scale_color_manual(values = cluster_colors[c(6, 7)]) +
    labs(x = "PC 1", y = "PC 2", color = "Clusters", fill = "Clusters") +
    #ggtitle("PCA plot") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
  
  
  # Create the UMAP plot
  UMAP_p <- DimPlot(sct_object, 
    reduction = "umap", 
    group.by = "seurat_clusters",
    pt.size = 5,
    cols = ggsci::pal_tron(palette = "legacy", alpha = 0.5)(7)[c(6, 7)])
  UMAP_p2 <- UMAP_p +
    #scale_fill_manual(values = cluster_colors[c(6, 7)]) +
    #scale_color_manual(values = cluster_colors[c(6, 7)]) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Clusters", fill = "Clusters") +
    #ggtitle("UMAP plot") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 18))
  

  # This if statement controls if the plost will be saved
  if (save_plots == TRUE) {
    
    ggsave(filename = paste0("PCA_clustering_gb", script_version, ".png"), PCA_p2,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    
    ggsave(filename = paste0("UMAP_clustering_gb", script_version, ".png"), UMAP_p2,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    
    message("The plots were saved to the folder: ", paste0(analysis_dir, "/", "Plots/"))
    
  } else {
    
    message("WARNING: the plots are not saved as a save was not requested")
    
  }
  
  
  # This if statement controls if the plots are printed
  if (print_plots == TRUE){
    
    print(PCA_p2)
    print(UMAP_p2)
    
    message("The plots were printed as requested.")
    
  } else {
    
    message("The plots were not printed as print was not requested.")
    
  }
  
  
  
}
plot_DR_cluster(sct_object = sct_CTC.obj, print_plots = TRUE, save_plots = TRUE)


# Removing unnecessary objects
rm(g2m.IDs, s.IDs, MT_genes, QC_PCA, QC_PCA_par, filt_CTC.obj, sct_variance)

#################################################################   End Section  ##################################################################






#################################################################   Plotting different metrics over the UMAP  ##################################################################
                                                                #         clusters to if they correlate       #


# This function simply runs a block of code, plotting and saving additional data overplayed on the UMAP (the only purpose of this is conciseness)
plot_UMAP_ovelrays = function (print_plots = TRUE, save_plots = TRUE) {
  
  # Define the used color plalettes
  grad_col <- pal_locuszoom("default", alpha = 0.5)(7)
  cluster_colors <- pal_tron(palette = "legacy")(7)
  cluster_colors_2 <- pal_futurama(palette = "planetexpress")(12)
  cluster_colors_3 <- pal_startrek(palette = "uniform")(7)
  cluster_colors_4 <- pal_rickandmorty(palette = "schwifty")(12)
  
  # Overlaying nFeature
  UMAP_nF_p <- FeaturePlot(sct_CTC.obj,
                           reduction = "umap",
                           features = "nFeature_SCT",
                           pt.size = 5)
  UMAP_nF_p2 <- UMAP_nF_p +
    #scale_color_gradient(low = "lightgrey", high = cluster_colors[7]) +
    scale_color_gradient(low = "lightgrey", high = grad_col[5]) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Clusters", fill = "Clusters") +
    #ggtitle("UMAP plot - nFeatures") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
  
  # Overlaying nCount
  UMAP_nC_p <- FeaturePlot(sct_CTC.obj,
                           reduction = "umap",
                           features = "nCount_SCT",
                           pt.size = 5)
  UMAP_nC_p2 <- UMAP_nC_p +
    scale_color_gradient(low = "lightgrey", high = grad_col[5]) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Clusters", fill = "Clusters") +
    #ggtitle("UMAP plot - nCounts") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
  
  # Overlay mt.percent
  UMAP_mt.percent <- FeaturePlot(sct_CTC.obj,
                                 reduction = "umap",
                                 features = "mt.percent",
                                 pt.size = 5)
  UMAP_mt.percent2 <- UMAP_mt.percent +
    scale_color_gradient(low = "lightgrey", high = grad_col[5]) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Clusters", fill = "Clusters") +
    #ggtitle("UMAP plot - MT gene expression") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
  
  # Overlay cell cycle
  UMAP_CellCyc <- DimPlot(sct_CTC.obj,
                          reduction = "umap",
                          group.by = "Phase",
                          pt.size = 5,
                          cols = ggsci::pal_startrek(palette = "uniform", alpha = 0.5)(7)[c(4, 5, 7)])
  UMAP_CellCyc2 <- UMAP_CellCyc +
    #scale_fill_manual(values = cluster_colors_3[c(4, 5, 7)]) +
    #scale_color_manual(values = cluster_colors_3[c(4, 5, 7)]) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Cell cycle phase", fill = "Cell cycle phase") +
    #ggtitle("UMAP plot - Cell cycle phases") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
  
  # Overlay treatment cycle
  UMAP_TreatCyc <- DimPlot(sct_CTC.obj,
                           reduction = "umap",
                           group.by = "RealCycle",
                           pt.size = 5,
                           cols = ggsci::pal_startrek(palette = "uniform", alpha = 0.5)(7)[c(4, 5, 7)])
  UMAP_TreatCyc2 <- UMAP_TreatCyc +
    #scale_fill_manual(values = cluster_colors_3[c(4, 5, 7)]) +
    #scale_color_manual(values = cluster_colors_3[c(4, 5, 7)]) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Treatment cycles", fill = "Treatment cycles") +
    #ggtitle("UMAP plot - Treatment cycles") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
  
  # Overlay treatment type
  UMAP_TreatType <- DimPlot(sct_CTC.obj,
                            reduction = "umap",
                            group.by = "Treatment",
                            pt.size = 5,
                            cols = ggsci::pal_startrek(palette = "uniform", alpha = 0.5)(7)[c(4, 5, 7)])
  UMAP_TreatType2 <- UMAP_TreatType +
    #scale_fill_manual(values = cluster_colors_3[c(4, 5, 7)]) +
    #scale_color_manual(values = cluster_colors_3[c(4, 5, 7)]) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Treatment", fill = "Treatment") +
    #ggtitle("UMAP plot - Treatment types") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
  
  # Overlay RespGroup
  UMAP_RespG <- DimPlot(sct_CTC.obj,
                        reduction = "umap",
                        group.by = "RespGroup",
                        pt.size = 5,
                        cols = ggsci::pal_tron(palette = "legacy", alpha = 0.5)(7)[c(1, 2)])
  UMAP_RespG2 <- UMAP_RespG +
    #scale_fill_manual(values = cluster_colors[c(1, 2)]) +
    #scale_color_manual(values = cluster_colors[c(1, 2)]) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Response group", fill = "Response group") +
    #ggtitle("UMAP plot - Response groups") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
  
  # Overlay Detailed response
  UMAP_DetRespG <- DimPlot(sct_CTC.obj,
                           reduction = "umap",
                           group.by = "detailedResp",
                           pt.size = 5,
                           cols = ggsci::pal_tron(palette = "legacy", alpha = 0.5)(7)[c(5, 1, 2, 3)])
  UMAP_DetRespG2 <- UMAP_DetRespG +
    #scale_fill_manual(values = cluster_colors[c(5, 1, 2, 3)]) +
    #scale_color_manual(values = cluster_colors[c(5, 1, 2, 3)]) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Response profile", fill = "Response profile") +
    #ggtitle("UMAP plot - Detailed treatment response") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
  
  # Overlay Luthetium-177 doses 
  UMAP_Lu <- DimPlot(sct_CTC.obj,
                     reduction = "umap",
                     group.by = "Lu177_GBq",
                     pt.size = 5,
                     cols = ggsci::pal_rickandmorty(palette = "schwifty", alpha = 0.5)(12)[c(4, 5, 6, 7, 8, 9)])
  UMAP_Lu2 <- UMAP_Lu +
    #scale_fill_manual(values = cluster_colors_4[c(4, 5, 6, 7, 8, 9)]) +
    #scale_color_manual(values = cluster_colors_4[c(4, 5, 6, 7, 8, 9)]) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Lu-177 doses (GBq)", fill = "Lutetium-177 doses (GBq)") +
    #ggtitle("UMAP plot - Lu-177 treatment doses") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
  
  # Overlay Actinium-225 doses
  UMAP_Ac <- DimPlot(sct_CTC.obj,
                     reduction = "umap",
                     group.by = "Ac225_MBq",
                     pt.size = 5,
                     cols = ggsci::pal_rickandmorty(palette = "schwifty", alpha = 0.5)(12)[c(4, 5, 6, 7)])
  UMAP_Ac2 <- UMAP_Ac +
    #scale_fill_manual(values = cluster_colors_4[c(4, 5, 6, 7)]) +
    #scale_color_manual(values = cluster_colors_4[c(4, 5, 6, 7)]) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Ac-225 doses (MBq)", fill = "Ac-225 doses (MBq)") +
    #ggtitle("UMAP plot - Ac-225 treatment doses") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
  
  # Overlay marker status
  UMAP_MarkerStat <- DimPlot(sct_CTC.obj,
                             reduction = "umap",
                             group.by = "Marker_status",
                             pt.size = 5,
                             cols = ggsci::pal_startrek(palette = "uniform", alpha = 0.5)(7)[c(4, 5, 7)])
  UMAP_MarkerStat2 <- UMAP_MarkerStat +
    #scale_fill_manual(values = cluster_colors_3[c(4, 5, 7)]) +
    #scale_color_manual(values = cluster_colors_3[c(4, 5, 7)]) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Marker status", fill = "Marker status") +
    #ggtitle("UMAP plot - Surface marker distribution") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
  
  
  
  # This if statement controls if th plots are printed
  if (print_plots == TRUE) {
    
    message("Overlay plots were requested. Running the plotting code block.")
    
    
    print(UMAP_nF_p2)
    print(UMAP_nC_p2)
    print(UMAP_mt.percent2)
    print(UMAP_CellCyc2)
    print(UMAP_TreatCyc2)
    print(UMAP_TreatType2)
    print(UMAP_RespG2)
    print(UMAP_DetRespG2)
    print(UMAP_Lu2)
    print(UMAP_Ac2)
    print(UMAP_MarkerStat2)
    
  } else {
    
    message("The overlay plots were printed as requested.")
    
  }
  
  
  # This if statement controls if the plots are saved
  if (save_plots == TRUE) {
    
    ggsave(filename = paste0("UMAP_nFeature", script_version, ".png"), UMAP_nF_p2,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    
    ggsave(filename = paste0("UMAP_nCount", script_version, ".png"), UMAP_nC_p2,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    
    ggsave(filename = paste0("UMAP_mt.percent", script_version, ".png"), UMAP_mt.percent2,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    
    ggsave(filename = paste0("UMAP_Cell_cycle", script_version, ".png"), UMAP_CellCyc2,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    
    ggsave(filename = paste0("UMAP_Treatment_cycles", script_version, ".png"), UMAP_TreatCyc2,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    
    ggsave(filename = paste0("UMAP_Treatment_type", script_version, ".png"), UMAP_TreatType2,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    
    ggsave(filename = paste0("UMAP_Response_group", script_version, ".png"), UMAP_RespG2,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    
    ggsave(filename = paste0("UMAP_Detailed_treatment_response", script_version, ".png"), UMAP_DetRespG2,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    
    ggsave(filename = paste0("UMAP_Lu-177_treatment_doses", script_version, ".png"), UMAP_Lu2,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    
    ggsave(filename = paste0("UMAP_Ac-225_treatment_doses", script_version, ".png"), UMAP_Ac2,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    
    ggsave(filename = paste0("UMAP_Marker_Status", script_version, ".png"), UMAP_MarkerStat2,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    rm(UMAP_MarkerStat, UMAP_MarkerStat2)
    
    
    message("The plots were saved to the folder: ", paste0(analysis_dir, "/", "Plots/"))
    
  } else {
    
    message("WARNING: the plots are not saved as a save was not requested")
    
  }
  
}
plot_UMAP_ovelrays(print_plots = TRUE, save_plots = TRUE)




## Visualizing the overlay data in a more understandable manner, using bar charts


# This function will calculate cluster composition stats and plots and saves them as requested (the only purpose of this is conciseness)
calculate_and_plot_cluster_composition = function (make_plots = TRUE, return_statistics = TRUE, print_statistics = TRUE) {
  
  # Define the color palette for plotting
  cluster_colors <- ggsci::pal_tron(palette = "legacy")(7)
  
  sct_CTC_meta_df <- sct_CTC.obj@meta.data
  sct_CTC_meta_df$TreatmentStatus <- ifelse(sct_CTC_meta_df$RealCycle == 0, "Non-treated", "Treated")
  sct_CTC_meta_df$patients <- str_remove(sct_CTC_meta_df$patientID, "." %R% END) 
  
  
  # Visualizing data from the variables side
  # Visualizing the treatment cycle data
  # Creating the plotting df
  Treatment_cycles_summary_df <- data.frame(Cycles = c(0, 1, 2),
                                            Total_cell_numbers = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 0)),
                                                                   nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 1)),
                                                                   nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 2))),
                                            Cell_numbers_c0 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 0, seurat_clusters == 0)),
                                                                nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 1, seurat_clusters == 0)),
                                                                nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 2, seurat_clusters == 0))),
                                            Cell_numbers_c1 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 0, seurat_clusters == 1)),
                                                                nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 1, seurat_clusters == 1)),
                                                                nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 2, seurat_clusters == 1))),
                                            Cell_percent_c0 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 0, seurat_clusters == 0)) /
                                                                  nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 0)) * 100,
                                                                nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 1, seurat_clusters == 0)) /
                                                                  nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 1)) * 100,
                                                                nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 2, seurat_clusters == 0)) /
                                                                  nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 2)) * 100),
                                            Cell_percent_c1 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 0, seurat_clusters == 1)) /
                                                                  nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 0)) * 100,
                                                                nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 1, seurat_clusters == 1)) /
                                                                  nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 1)) * 100,
                                                                nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 2, seurat_clusters == 1)) /
                                                                  nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 2)) * 100))
  
  
  
  # Transforming the plotting df to a long format
  Treatment_cycles_summary_df_long <- pivot_longer(Treatment_cycles_summary_df, cols = c(Cell_percent_c0, Cell_percent_c1), names_to = "Clusters", values_to = "Percentages")
  
  
  # Visualizing the treatment type data
  # Creating the plotting df
  Treatment_summary_df <- data.frame(Treatment_type = c("AL", "L", "NT"),
                                     Total_cell_numbers = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "AL")),
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "L")),
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "NT"))),
                                     Cell_numbers_c0 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "AL", seurat_clusters == 0)),
                                                         nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "L", seurat_clusters == 0)),
                                                         nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "NT", seurat_clusters == 0))),
                                     Cell_numbers_c1 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "AL", seurat_clusters == 1)),
                                                         nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "L", seurat_clusters == 1)),
                                                         nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "NT", seurat_clusters == 1))),
                                     Cell_percent_c0 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "AL", seurat_clusters == 0)) /
                                                           nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "AL")) * 100,
                                                         nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "L", seurat_clusters == 0)) /
                                                           nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "L")) * 100,
                                                         nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "NT", seurat_clusters == 0)) /
                                                           nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "NT")) * 100),
                                     Cell_percent_c1 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "AL", seurat_clusters == 1)) /
                                                           nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "AL")) * 100,
                                                         nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "L", seurat_clusters == 1)) /
                                                           nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "L")) * 100,
                                                         nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "NT", seurat_clusters == 1)) /
                                                           nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "NT")) * 100))
  
  
  # Transforming the plotting df to a long format
  Treatment_summary_df_long <- pivot_longer(Treatment_summary_df, cols = c(Cell_percent_c0, Cell_percent_c1), names_to = "Clusters", values_to = "Percentages")
  
  
  # Visualizing the cell cycle data
  # Creating the plotting df
  Cell_cycle_summary_df <- data.frame(CC_phase = c("G1", "S", "G2M"),
                                      Total_cell_numbers = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "G1")),
                                                             nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "S")),
                                                             nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "G2M"))),
                                      Cell_numbers_c0 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "G1", seurat_clusters == 0)),
                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "S", seurat_clusters == 0)),
                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "G2M", seurat_clusters == 0))),
                                      Cell_numbers_c1 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "G1", seurat_clusters == 1)),
                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "S", seurat_clusters == 1)),
                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "G2M", seurat_clusters == 1))),
                                      Cell_percent_c0 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "G1", seurat_clusters == 0)) /
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "G1")) * 100,
                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "S", seurat_clusters == 0)) /
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "S")) * 100,
                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "G2M", seurat_clusters == 0)) /
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "G2M")) * 100),
                                      Cell_percent_c1 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "G1", seurat_clusters == 1)) /
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "G1")) * 100,
                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "S", seurat_clusters == 1)) /
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "S")) * 100,
                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "G2M", seurat_clusters == 1)) /
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "G2M")) * 100))
  
  
  # Transforming the plotting df to a long format
  Cell_cycle_summary_df_long <- pivot_longer(Cell_cycle_summary_df, cols = c(Cell_percent_c0, Cell_percent_c1), names_to = "Clusters", values_to = "Percentages")
  print(Cell_cycle_summary_df_long)
  
  
  # Visualizing the response group data
  # Creating the plotting df
  Response_summary_df <- data.frame(response_group = c("Responder", "Nonresponder"),
                                    Total_cell_numbers = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Responder")),
                                                           nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Nonresponder"))),
                                    Cell_numbers_c0 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Responder", seurat_clusters == 0)),
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Nonresponder", seurat_clusters == 0))),
                                    Cell_numbers_c1 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Responder", seurat_clusters == 1)),
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Nonresponder", seurat_clusters == 1))),
                                    Cell_percent_c0 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Responder", seurat_clusters == 0)) /
                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Responder")) * 100,
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Nonresponder", seurat_clusters == 0)) /
                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Nonresponder")) * 100),
                                    Cell_percent_c1 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Responder", seurat_clusters == 1)) /
                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Responder")) * 100,
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Nonresponder", seurat_clusters == 1)) /
                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Nonresponder")) * 100))
  
  
  # Transforming the plotting df to a long format
  Response_summary_df_long <- pivot_longer(Response_summary_df, cols = c(Cell_percent_c0, Cell_percent_c1), names_to = "Clusters", values_to = "Percentages")
  
  
  
  # Visualizing the response group data
  # Creating the plotting df
  Treatment_status_df <- data.frame(treatment_status = c("Non-treated", "Treated"),
                                    Total_cell_numbers = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, TreatmentStatus == "Non-treated")),
                                                           nrow(dplyr::filter(.data = sct_CTC_meta_df, TreatmentStatus == "Treated"))),
                                    Cell_numbers_c0 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, TreatmentStatus == "Non-treated", seurat_clusters == 0)),
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, TreatmentStatus == "Treated", seurat_clusters == 0))),
                                    Cell_numbers_c1 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, TreatmentStatus == "Non-treated", seurat_clusters == 1)),
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, TreatmentStatus == "Treated", seurat_clusters == 1))),
                                    Cell_percent_c0 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, TreatmentStatus == "Non-treated", seurat_clusters == 0)) /
                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, TreatmentStatus == "Non-treated")) * 100,
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, TreatmentStatus == "Treated", seurat_clusters == 0)) /
                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, TreatmentStatus == "Treated")) * 100),
                                    Cell_percent_c1 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, TreatmentStatus == "Non-treated", seurat_clusters == 1)) /
                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, TreatmentStatus == "Non-treated")) * 100,
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, TreatmentStatus == "Treated", seurat_clusters == 1)) /
                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, TreatmentStatus == "Treated")) * 100))
  
  
  # Transforming the plotting df to a long format
  Treatment_status_df_long <- pivot_longer(Treatment_status_df, cols = c(Cell_percent_c0, Cell_percent_c1), names_to = "Clusters", values_to = "Percentages")
  
  
  # Visualizing the response group data
  # Creating the plotting df
  Detailed_response_summary_df <- data.frame(response_group = c("NT", "PD", "PR", "SD"),
                                             Total_cell_numbers = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "NT")),
                                                                    nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "PD")),
                                                                    nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "PR")),
                                                                    nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "SD"))),
                                             Cell_numbers_c0 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "NT", seurat_clusters == 0)),
                                                                 nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "PD", seurat_clusters == 0)),
                                                                 nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "PR", seurat_clusters == 0)),
                                                                 nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "SD", seurat_clusters == 0))),
                                             Cell_numbers_c1 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "NT", seurat_clusters == 1)),
                                                                 nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "PD", seurat_clusters == 1)),
                                                                 nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "PR", seurat_clusters == 1)),
                                                                 nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "SD", seurat_clusters == 1))),
                                             Cell_percent_c0 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "NT", seurat_clusters == 0)) /
                                                                   nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "NT")) * 100,
                                                                 nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "PD", seurat_clusters == 0)) /
                                                                   nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "PD")) * 100,
                                                                 nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "PR", seurat_clusters == 0)) /
                                                                   nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "PR")) * 100,
                                                                 nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "SD", seurat_clusters == 0)) /
                                                                   nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "SD")) * 100),
                                             Cell_percent_c1 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "NT", seurat_clusters == 1)) /
                                                                   nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "NT")) * 100,
                                                                 nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "PD", seurat_clusters == 1)) /
                                                                   nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "PD")) * 100,
                                                                 nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "PR", seurat_clusters == 1)) /
                                                                   nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "PR")) * 100,
                                                                 nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "SD", seurat_clusters == 1)) /
                                                                   nrow(dplyr::filter(.data = sct_CTC_meta_df, detailedResp == "SD")) * 100))
  
  
  # Transforming the plotting df to a long format
  Detailed_response_summary_df_long <- pivot_longer(Detailed_response_summary_df, cols = c(Cell_percent_c0, Cell_percent_c1), names_to = "Clusters", values_to = "Percentages")
  
  # Visualizing the surface marker distribution data
  # Creating the plotting df form the markers side
  Surface_marker_summary_df <- data.frame(marker_status = c("EpCAM+", "EpCAM+PSMA+", "PSMA+"),
                                          Total_cell_numbers = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker_status == "EpCAM+")),
                                                                 nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker_status == "EpCAM+PSMA+")),
                                                                 nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker_status == "PSMA+"))),
                                          Cell_numbers_c0 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker_status == "EpCAM+", seurat_clusters == 0)),
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker_status == "EpCAM+PSMA+", seurat_clusters == 0)),
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker_status == "PSMA+", seurat_clusters == 0))),
                                          Cell_numbers_c1 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker_status == "EpCAM+", seurat_clusters == 1)),
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker_status == "EpCAM+PSMA+", seurat_clusters == 1)),
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker_status == "PSMA+", seurat_clusters == 1))),
                                          Cell_percent_c0 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker_status == "EpCAM+", seurat_clusters == 0)) /
                                                                nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker_status == "EpCAM+")) * 100,
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker_status == "EpCAM+PSMA+", seurat_clusters == 0)) /
                                                                nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker_status == "EpCAM+PSMA+")) * 100,
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker_status == "PSMA+", seurat_clusters == 0)) /
                                                                nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker_status == "PSMA+")) * 100),
                                          Cell_percent_c1 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker_status == "EpCAM+", seurat_clusters == 1)) /
                                                                nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker_status == "EpCAM+")) * 100,
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker_status == "EpCAM+PSMA+", seurat_clusters == 1)) /
                                                                nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker_status == "EpCAM+PSMA+")) * 100,
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker_status == "PSMA+", seurat_clusters == 1)) /
                                                                nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker_status == "PSMA+")) * 100))
  
  
  # Transform the plotting df to the long format for plotting
  Surface_marker_summary_df_long <- pivot_longer(Surface_marker_summary_df, cols = c(Cell_percent_c0, Cell_percent_c1), names_to = "Clusters", values_to = "Percentages")
  
  
  # Visualizing data from the cluster side
  # Visualizing the treatment cycle data
  # Creating the plotting df
  Treatment_cycles_summary_df_clust <- data.frame(Clusters = c("Cluster 0", "Cluster 1"),
                                                  Total_cell_numbers = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")), 
                                                                         nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1"))),
                                                  Cycle_0_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", RealCycle == "0")),
                                                                nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", RealCycle == "0"))),
                                                  Cycle_1_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", RealCycle == "1")),
                                                                nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", RealCycle == "1"))),
                                                  Cycle_2_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", RealCycle == "2")),
                                                                nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", RealCycle == "2"))),
                                                  Cycle_0_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", RealCycle == "0")) /
                                                                  nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                                nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", RealCycle == "0")) /
                                                                  nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100),
                                                  Cycle_1_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", RealCycle == "1")) /
                                                                  nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                                nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", RealCycle == "1")) /
                                                                  nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100),
                                                  Cycle_2_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", RealCycle == "2")) /
                                                                  nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                                nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", RealCycle == "2")) /
                                                                  nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100))
  
  
  # Transform the plotting df to a long format for plotting
  Treatment_cycles_summary_df_clust_long <- pivot_longer(Treatment_cycles_summary_df_clust, cols = c(Cycle_0_p, Cycle_1_p, Cycle_2_p), names_to = "Cycles", values_to = "Percentages")
  
  
  # Visualizing the treatment type data
  # Creating the plotting df
  Treatment_summary_df_clust <- data.frame(Clusters = c("Cluster 0", "Cluster 1"),
                                           Total_cell_numbers = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")), 
                                                                  nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1"))),
                                           AL_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Treatment == "AL")),
                                                    nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Treatment == "AL"))),
                                           L_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Treatment == "L")),
                                                   nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Treatment == "L"))),
                                           NT_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Treatment == "NT")),
                                                    nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Treatment == "NT"))),
                                           AL_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Treatment == "AL")) /
                                                      nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                    nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Treatment == "AL")) /
                                                      nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100),
                                           L_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Treatment == "L")) /
                                                     nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                   nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Treatment == "L")) /
                                                     nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100),
                                           NT_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Treatment == "NT")) /
                                                      nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                    nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Treatment == "NT")) /
                                                      nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100))
  
  
  # Transforming the plotting df to a long format
  Treatment_summary_df_clust_long <- pivot_longer(Treatment_summary_df_clust, cols = c(AL_p, L_p, NT_p), names_to = "Treatment_type", values_to = "Percentages")
  
  
  # Visualizing the cell cycle data
  # Creating the plotting df
  Cell_cycle_summary_df_clust <- data.frame(Clusters = c("Cluster 0", "Cluster 1"),
                                            Total_cell_numbers = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")), 
                                                                   nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1"))),
                                            S_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Phase == "S")),
                                                    nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Phase == "S"))),
                                            G1_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Phase == "G1")),
                                                     nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Phase == "G1"))),
                                            G2M_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Phase == "G2M")),
                                                      nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Phase == "G2M"))),
                                            S_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Phase == "S")) /
                                                      nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                    nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Phase == "S")) /
                                                      nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100),
                                            G1_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Phase == "G1")) /
                                                       nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                     nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Phase == "G1")) /
                                                       nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100),
                                            G2M_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Phase == "G2M")) /
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                      nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Phase == "G2M")) /
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100))
  
  
  # Transforming the plotting df to a long format
  Cell_cycle_summary_df_clust_long <- pivot_longer(Cell_cycle_summary_df_clust, cols = c(S_p, G1_p, G2M_p), names_to = "CC_phase", values_to = "Percentages")
  
  
  # Visualizing the response group data
  # Creating the plotting df
  Response_summary_df_clust <- data.frame(Clusters = c("Cluster 0", "Cluster 1"),
                                          Total_cell_numbers = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")), 
                                                                 nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1"))),
                                          Resp_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", RespGroup == "Responder")),
                                                     nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", RespGroup == "Responder"))),
                                          Nonresp_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", RespGroup == "Nonresponder")),
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", RespGroup == "Nonresponder"))),
                                          Resp_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", RespGroup == "Responder")) /
                                                       nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                     nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", RespGroup == "Responder")) /
                                                       nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100),
                                          Nonresp_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", RespGroup == "Nonresponder")) /
                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", RespGroup == "Nonresponder")) /
                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100))
  
  
  # Transforming the plotting df to a long format
  Response_summary_df_clust_long <- pivot_longer(Response_summary_df_clust, cols = c(Resp_p, Nonresp_p), names_to = "Response_group", values_to = "Percentages")
  
  
  # Visualizing the response group data
  # Creating the plotting df
  Treatment_status_df_clust <- data.frame(Clusters = c("Cluster 0", "Cluster 1"),
                                          Total_cell_numbers = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")), 
                                                                 nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1"))),
                                          NTreat_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", TreatmentStatus == "Non-treated")),
                                                       nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", TreatmentStatus == "Non-treated"))),
                                          Treat_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", TreatmentStatus == "Treated")),
                                                      nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", TreatmentStatus == "Treated"))),
                                          NTreat_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", TreatmentStatus == "Non-treated")) /
                                                         nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                       nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", TreatmentStatus == "Non-treated")) /
                                                         nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100),
                                          Treat_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", TreatmentStatus == "Treated")) /
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                      nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", TreatmentStatus == "Treated")) /
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100))
  
  
  # Transforming the plotting df to a long format
  Treatment_status_df_clust_long <- pivot_longer(Treatment_status_df_clust, cols = c(NTreat_p, Treat_p), names_to = "Treatment_status", values_to = "Percentages")
  
  
  # Calculating the detailed response stat
  # Creating the plotting df
  Detailed_response_summary_df_clust <- data.frame(Clusters = c("Cluster 0", "Cluster 1"),
                                                   Total_cell_numbers = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")), 
                                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1"))),
                                                   NT_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", detailedResp == "NT")),
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", detailedResp == "NT"))),
                                                   PD_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", detailedResp == "PD")),
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", detailedResp == "PD"))),
                                                   PR_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", detailedResp == "PR")),
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", detailedResp == "PR"))),
                                                   SD_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", detailedResp == "SD")),
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", detailedResp == "SD"))),
                                                   NT_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", detailedResp == "NT")) /
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", detailedResp == "NT")) /
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100),
                                                   PD_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", detailedResp == "PD")) /
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", detailedResp == "PD")) /
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100),
                                                   PR_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", detailedResp == "PR")) /
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", detailedResp == "PR")) /
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100),
                                                   SD_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", detailedResp == "SD")) /
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", detailedResp == "SD")) /
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100))
  
  
  # Transforming the plotting df to a long format
  Detailed_response_summary_df_clust_long <- pivot_longer(Detailed_response_summary_df_clust, cols = c(NT_p, PD_p, PR_p, SD_p),
                                                          names_to = "Response_group", values_to = "Percentages")
  
  # Visualizing the marker status data
  # Creating the plotting df form the cluster side
  Surface_marker_summary_df_clust <- data.frame(Clusters = c("Cluster 0", "Cluster 1"),
                                                Total_cell_numbers = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")), 
                                                                       nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1"))),
                                                EpCAM_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Marker_status == "EpCAM+")),
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Marker_status == "EpCAM+"))),
                                                EpCAM_PSMA_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Marker_status == "EpCAM+PSMA+")),
                                                                 nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Marker_status == "EpCAM+PSMA+"))),
                                                PSMA_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Marker_status == "PSMA+")),
                                                           nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Marker_status == "PSMA+"))),
                                                EpCAM_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Marker_status == "EpCAM+")) /
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Marker_status == "EpCAM+")) /
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100),
                                                EpCAM_PSMA_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Marker_status == "EpCAM+PSMA+")) /
                                                                   nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                                 nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Marker_status == "EpCAM+PSMA+")) /
                                                                   nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100),
                                                PSMA_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Marker_status == "PSMA+")) /
                                                             nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                           nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Marker_status == "PSMA+")) /
                                                             nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100))
  
  
  # Transform the plotting df to a long format for plotting
  Surface_marker_summary_df_clust_long <- pivot_longer(Surface_marker_summary_df_clust, cols = c(EpCAM_p, EpCAM_PSMA_p, PSMA_p), names_to = "Markers", values_to = "Percentages")
  
  
  # This if statement control if the calculated stat dfs will be assigned to the global env and printed to the console
  if (return_statistics == TRUE) {
    variable_side_cluster_stats <- list(
      Treatment_cycles_summary_df,
      Treatment_summary_df,
      Cell_cycle_summary_df,
      Response_summary_df,
      Treatment_status_df,
      Detailed_response_summary_df,
      Surface_marker_summary_df
    )
    
    names(variable_side_cluster_stats) <- c("Treatment_cycles_summary_df",
                                            "Treatment_summary_df",
                                            "Cell_cycle_summary_df",
                                            "Response_summary_df",
                                            "Treatment_status_df",
                                            "Detailed_response_summary_df",
                                            "Surface_marker_summary_df")
    
    
    variable_side_cluster_stats_plotting_dfs <- list(
      Treatment_cycles_summary_df_long,
      Treatment_summary_df_long,
      Cell_cycle_summary_df_long,
      Response_summary_df_long,
      Treatment_status_df_long,
      Detailed_response_summary_df_long,
      Surface_marker_summary_df_long
    )
    
    names(variable_side_cluster_stats) <- c("Treatment_cycles_summary_df_long",
                                            "Treatment_summary_df_long",
                                            "Cell_cycle_summary_df_long",
                                            "Response_summary_df_long",
                                            "Treatment_status_df_long",
                                            "Detailed_response_summary_df_long",
                                            "Surface_marker_summary_df_long")
    
    
    cluster_side_cluster_stats <- list(
      Treatment_cycles_summary_df_clust,
      Treatment_summary_df_clust,
      Cell_cycle_summary_df_clust,
      Response_summary_df_clust,
      Treatment_status_df_clust,
      Detailed_response_summary_df_clust,
      Surface_marker_summary_df_clust
    )
    
    names(variable_side_cluster_stats) <- c("Treatment_cycles_summary_df_clust",
                                            "Treatment_summary_df_clust",
                                            "Cell_cycle_summary_df_clust",
                                            "Response_summary_df_clust",
                                            "Treatment_status_df_clust",
                                            "Detailed_response_summary_df_clust",
                                            "Surface_marker_summary_df_clust")
    
    
    cluster_side_cluster_stats_plotting_dfs <- list(
      Treatment_cycles_summary_df_clust_long,
      Treatment_summary_df_clust_long,
      Cell_cycle_summary_df_clust_long,
      Response_summary_df_clust_long,
      Treatment_status_df_clust_long,
      Detailed_response_summary_df_clust_long,
      Surface_marker_summary_df_clust_long
    )
    
    names(variable_side_cluster_stats) <- c("Treatment_cycles_summary_df_clust_long",
                                            "Treatment_summary_df_clust_long",
                                            "Cell_cycle_summary_df_clust_long",
                                            "Response_summary_df_clust_long",
                                            "Treatment_status_df_clust_long",
                                            "Detailed_response_summary_df_clust_long",
                                            "Surface_marker_summary_df_clust_long")
    
    stat_lists <- list(variable_side_cluster_stats, variable_side_cluster_stats_plotting_dfs, 
                       cluster_side_cluster_stats, cluster_side_cluster_stats_plotting_dfs)
    
    names(stat_lists) <- c("variable_side_cluster_stats", "variable_side_cluster_stats_plotting_dfs", 
                           "cluster_side_cluster_stats", "cluster_side_cluster_stats_plotting_dfs")
    
    
    
    assign(x = "stat_lists",
           value = stat_lists,
           envir = .GlobalEnv)
    
    
    # This if statement control if the calculated stat dfs will be printed to the consol
    if (print_statistics == TRUE) {
      for (index in seq_along(stat_lists)) {
        for (element in stat_lists[index]) {
          print(element)
        }
      }
      
    } else {
      
      message("The stat table print was not requested. \n")
      
    }
    
    
    # Return statistics message
    message("Return statistcs was requested, the stat tables were assigned to a new list: stat_lists")
    
  } else {
    
    message("No return statisctics was requested, teherefore no stats will be assigned or printed. \n")
    
  }
  
  
  # This if statement controls if the plots will be printed and saved
  if (make_plots == TRUE) {
    
    # Plotting the treatment cycle distribution
    treatment_cyc_bar_p <- ggplot(Treatment_cycles_summary_df_long,
                                  aes(x = Cycles, y = Percentages, fill = Clusters)) +
      geom_col(color = "black", position = "dodge") +
      geom_text(aes(label = round(Percentages, digits = 2)), vjust = -0.5, position = position_dodge(width = 0.9)) +
      scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
      #scale_x_discrete(labels = c("Double positive","EpCAM positive","PSMA positive")) +
      scale_fill_manual(values = cluster_colors[c(6, 7)], label = c("Cluster 0", "Cluster 1")) +
      labs(fill = "CTC cluster", x = expression("Treatment cycle"), y = expression("Cell number percentages")) +
      #ggtitle("CTC clusters - cell distribution based on treatment cycle") +
      #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
            strip.text = element_blank(), text = element_text(size = 24))
    print(treatment_cyc_bar_p)
    
    ggsave(filename = paste0("Treatment cycles in clusters", script_version, ".png"), treatment_cyc_bar_p,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    rm(treatment_cyc_bar_p)
    
    
    # Plotting the treatment type distribution
    treatment_type_bar_p <- ggplot(Treatment_summary_df_long,
                                   aes(x = factor(Treatment_type, levels = c("NT", "AL", "L")), y = Percentages, fill = Clusters)) +
      geom_col(color = "black", position = "dodge") +
      geom_text(aes(label = round(Percentages, digits = 2)), vjust = -0.5, position = position_dodge(width = 0.9)) +
      scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
      #scale_x_discrete(labels = c("Double positive","EpCAM positive","PSMA positive")) +
      scale_fill_manual(values = cluster_colors[c(6, 7)], label = c("Cluster 0", "Cluster 1")) +
      labs(fill = "CTC cluster", x = expression("Treatment type"), y = expression("Cell number percentages")) +
      #ggtitle("CTC clusters - cell distribution based on treatment types") +
      #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
            strip.text = element_blank(), text = element_text(size = 24))
    print(treatment_type_bar_p)
    
    ggsave(filename = paste0("Treatment types over clusters", script_version, ".png"), treatment_type_bar_p,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    rm(treatment_type_bar_p)
    
    
    # Plotting the cell cycle distribution
    cell_cyc_bar_p <- ggplot(Cell_cycle_summary_df_long,
                             aes(x = factor(CC_phase, levels = c("G1", "S", "G2M")), y = Percentages, fill = Clusters)) +
      geom_col(color = "black", position = "dodge") +
      geom_text(aes(label = round(Percentages, digits = 2)), vjust = -0.5, position = position_dodge(width = 0.9)) +
      scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
      #scale_x_discrete(labels = c("Double positive","EpCAM positive","PSMA positive")) +
      scale_fill_manual(values = cluster_colors[c(6, 7)], label = c("Cluster 0", "Cluster 1")) +
      labs(fill = "CTC cluster", x = expression("Cell cycle phase"), y = expression("Cell number percentages")) +
      #ggtitle("CTC clusters - cell distribution based on cell cycle phase") +
      #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
            strip.text = element_blank(), text = element_text(size = 24))
    print(cell_cyc_bar_p)
    
    ggsave(filename = paste0("Cell cycle over clusters", script_version, ".png"), cell_cyc_bar_p,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    rm(cell_cyc_bar_p)
    
    
    # Plotting the treatment response distribution from response view
    response_group_bar_p <- ggplot(Response_summary_df_long,
                                   aes(x = factor(response_group, levels = c("Responder", "Nonresponder")), y = Percentages, fill = Clusters)) +
      geom_col(color = "black", position = "dodge") +
      geom_text(aes(label = round(Percentages, digits = 2)), vjust = -0.5, position = position_dodge(width = 0.9)) +
      scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
      #scale_x_discrete(labels = c("Double positive","EpCAM positive","PSMA positive")) +
      scale_fill_manual(values = cluster_colors[c(6, 7)], label = c("Cluster 0", "Cluster 1")) +
      labs(fill = "CTC cluster", x = expression("Response group"), y = expression("Cell number percentages")) +
      #ggtitle("CTC clusters - cell distribution based on treatment response") +
      #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
            strip.text = element_blank(), text = element_text(size = 24))
    print(response_group_bar_p)
    
    ggsave(filename = paste0("Treatment response over clusters", script_version, ".png"), response_group_bar_p,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    rm(response_group_bar_p)
    
    
    # Plotting the treatment status distribution from response view
    treatment_group_bar_p <- ggplot(Treatment_status_df_long,
                                    aes(x = factor(treatment_status, levels = c("Non-treated", "Treated")), y = Percentages, fill = Clusters)) +
      geom_col(color = "black", position = "dodge") +
      geom_text(aes(label = round(Percentages, digits = 2)), vjust = -0.5, position = position_dodge(width = 0.9)) +
      scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
      #scale_x_discrete(labels = c("Double positive","EpCAM positive","PSMA positive")) +
      scale_fill_manual(values = cluster_colors[c(6, 7)], label = c("Cluster 0", "Cluster 1")) +
      labs(fill = "CTC cluster", x = expression("Treatment status"), y = expression("Cell number percentages")) +
      #ggtitle("CTC clusters - cell distribution based on treatment status") +
      #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
            strip.text = element_blank(), text = element_text(size = 24))
    print(treatment_group_bar_p)
    
    ggsave(filename = paste0("Treatment status over clusters", script_version, ".png"), treatment_group_bar_p,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    rm(treatment_group_bar_p)
    
    
    # Plotting the detailed treatment response distribution from response view
    detailed_response_bar_p <- ggplot(Detailed_response_summary_df_long,
                                      aes(x = factor(response_group, levels = c("NT", "PD", "PR", "SD")), y = Percentages, fill = Clusters)) +
      geom_col(color = "black", position = "dodge") +
      geom_text(aes(label = round(Percentages, digits = 2)), vjust = -0.5, position = position_dodge(width = 0.9)) +
      scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
      #scale_x_discrete(labels = c("Double positive","EpCAM positive","PSMA positive")) +
      scale_fill_manual(values = cluster_colors[c(6, 7)], label = c("Cluster 0", "Cluster 1")) +
      labs(fill = "CTC cluster", x = expression("Detailed treatment response"), y = expression("Cell number percentages")) +
      #ggtitle("CTC clusters - cell distribution based on detailed treatment response") +
      #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
            strip.text = element_blank(), text = element_text(size = 24))
    print(detailed_response_bar_p)
    
    ggsave(filename = paste0("Detailed treatment response over clusters", script_version, ".png"), detailed_response_bar_p,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    rm(detailed_response_bar_p)
    
    
    # Plotting the surface marker distribution
    surface_marker_bar_p <- ggplot(Surface_marker_summary_df_long,
                                   aes(x = factor(marker_status, levels = c("EpCAM+", "EpCAM+PSMA+", "PSMA+")), y = Percentages, fill = Clusters)) +
      geom_col(color = "black", position = "dodge") +
      geom_text(aes(label = round(Percentages, digits = 2)), vjust = -0.5, position = position_dodge(width = 0.9)) +
      scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
      #scale_x_discrete(labels = c("Double positive","EpCAM positive","PSMA positive")) +
      scale_fill_manual(values = cluster_colors[c(6, 7)], label = c("Cluster 0", "Cluster 1")) +
      labs(fill = "CTC cluster", x = expression("Surface marker status"), y = expression("Cell number percentages")) +
      #ggtitle("CTC clusters - cell distribution based on the surface marker status") +
      #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
            strip.text = element_blank(), text = element_text(size = 24))
    print(surface_marker_bar_p)
    
    ggsave(filename = paste0("Surface marker over clusters", script_version, ".png"), surface_marker_bar_p,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    rm(surface_marker_bar_p)
    
    
    # Plotting the treatment cycle distribution from the cluster view
    treatment_cyc_bar_p <- ggplot(Treatment_cycles_summary_df_clust_long,
                                  aes(x = Clusters, y = Percentages, fill = as.factor(Cycles))) +
      geom_col(color = "black", position = "dodge") +
      geom_text(aes(label = round(Percentages, digits = 2)), vjust = -0.5, position = position_dodge(width = 0.9)) +
      scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
      scale_x_discrete(labels = c("Cluster 0", "Cluster 1", "Cluster 2")) +
      scale_fill_manual(values = cluster_colors[c(4, 6, 7)], label = c("Cycle 0", "Cycle 1", "Cycle 2")) +
      labs(fill = "Treatment cycles", x = expression("CTC clusters"), y = expression("Cell number percentages")) +
      #ggtitle("CTC clusters - cell distribution based on treatment cycle") +
      #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
            strip.text = element_blank(), text = element_text(size = 24))
    print(treatment_cyc_bar_p)
    
    ggsave(filename = paste0("Cluster_composition-Treatment_cycles", script_version, ".png"), treatment_cyc_bar_p,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    rm(treatment_cyc_bar_p)
    
    
    # Plotting the treatment type distribution from the cluster view
    treatment_type_bar_p <- ggplot(Treatment_summary_df_clust_long,
                                   aes(x = Clusters, y = Percentages, fill = factor(Treatment_type, levels = c("NT_p", "AL_p", "L_p")))) +
      geom_col(color = "black", position = "dodge") +
      geom_text(aes(label = round(Percentages, digits = 2)), vjust = -0.5, position = position_dodge(width = 0.9)) +
      scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
      scale_x_discrete(labels = c("Cluster 0", "Cluster 1", "Cluster 2")) +
      scale_fill_manual(values = cluster_colors[c(4, 6, 7)], label = c("NT", "AL", "L")) +
      labs(fill = "Treatment types", x = expression("CTC clusters"), y = expression("Cell number percentages")) +
      #ggtitle("CTC clusters - cell distribution based on treatment types") +
      #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
            strip.text = element_blank(), text = element_text(size = 24))
    print(treatment_type_bar_p)
    
    ggsave(filename = paste0("Cluster_composition-Treatment_types", script_version, ".png"), treatment_type_bar_p,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    rm(treatment_type_bar_p)
    
    
    # Plotting the cell cycle distribution from the cluster view
    cell_cyc_bar_p <- ggplot(Cell_cycle_summary_df_clust_long,
                             aes(x = Clusters, y = Percentages, fill = factor(CC_phase, levels = c("G1_p", "S_p", "G2M_p")))) +
      geom_col(color = "black", position = "dodge") +
      geom_text(aes(label = round(Percentages, digits = 2)), vjust = -0.5, position = position_dodge(width = 0.9)) +
      scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
      scale_x_discrete(labels = c("Cluster 0", "Cluster 1", "Cluster 2")) +
      scale_fill_manual(values = cluster_colors[c(4, 6, 7)], label = c("G1", "S", "G2/M")) +
      labs(fill = "Cell cycle \nphase", x = expression("CTC cluster"), y = expression("Cell number percentages")) +
      #ggtitle("CTC clusters - cell distribution based on cell cycle phase") +
      #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
            strip.text = element_blank(), text = element_text(size = 24))
    print(cell_cyc_bar_p)
    
    ggsave(filename = paste0("Cluster_composition-Cell_cycle_phases", script_version, ".png"), cell_cyc_bar_p,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    rm(cell_cyc_bar_p)
    
    
    # Plotting the treatment response distribution from the cluster view
    response_group_bar_p <- ggplot(Response_summary_df_clust_long,
                                   aes(x = Clusters, y = Percentages, fill = factor(Response_group, levels = c("Resp_p", "Nonresp_p")))) +
      geom_col(color = "black", position = "dodge") +
      geom_text(aes(label = round(Percentages, digits = 2)), vjust = -0.5, position = position_dodge(width = 0.9)) +
      scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
      scale_x_discrete(labels = c("Cluster 0", "Cluster 1")) +
      scale_fill_manual(values = cluster_colors[c(2, 1)], label = c("Responder", "Nonresponder")) +
      labs(fill = "Response group", x = expression("CTC cluster"), y = expression("Cell number percentages")) +
      #ggtitle("CTC clusters - cell distribution based on treatment response") +
      #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
            strip.text = element_blank(), text = element_text(size = 24))
    print(response_group_bar_p)
    
    ggsave(filename = paste0("Cluster_composition-Response_groups", script_version, ".png"), response_group_bar_p,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    rm(response_group_bar_p)
    
    
    # Plotting the treatment status distribution from the cluster view
    treatment_status_bar_p <- ggplot(Treatment_status_df_clust_long,
                                     aes(x = Clusters, y = Percentages, fill = factor(Treatment_status, levels = c("NTreat_p", "Treat_p")))) +
      geom_col(color = "black", position = "dodge") +
      geom_text(aes(label = round(Percentages, digits = 2)), vjust = -0.5, position = position_dodge(width = 0.9)) +
      scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
      scale_x_discrete(labels = c("Cluster 0", "Cluster 1")) +
      scale_fill_manual(values = cluster_colors[c(2, 1)], label = c("Non-treated", "Treated")) +
      labs(fill = "Treatment status", x = expression("CTC cluster"), y = expression("Cell number percentages")) +
      #ggtitle("CTC clusters - cell distribution based on treatment response") +
      #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
            strip.text = element_blank(), text = element_text(size = 24))
    print(treatment_status_bar_p)
    
    ggsave(filename = paste0("Cluster_composition-Treatment_status
                             ", script_version, ".png"), treatment_status_bar_p,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    rm(treatment_status_bar_p)
    
    
    # Plotting the detailed treatment response distribution from the cluster view
    detailed_response_group_bar_p <- ggplot(Detailed_response_summary_df_clust_long,
                                            aes(x = Clusters, y = Percentages, fill = factor(Response_group, levels = c("NT_p", "PD_p", "PR_p", "SD_p")))) +
      geom_col(color = "black", position = "dodge") +
      geom_text(aes(label = round(Percentages, digits = 2)), vjust = -0.5, position = position_dodge(width = 0.9)) +
      scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
      scale_x_discrete(labels = c("Cluster 0", "Cluster 1")) +
      scale_fill_manual(values = cluster_colors[c(5, 1, 2, 3)], label = c("NT", "PD", "PR", "SD")) +
      labs(fill = "Detailed treatment \nresponse", x = expression("CTC cluster"), y = expression("Cell number percentages")) +
      #ggtitle("CTC clusters - cell distribution based on detailed treatment response") +
      #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
            strip.text = element_blank(), text = element_text(size = 24))
    print(detailed_response_group_bar_p)
    
    ggsave(filename = paste0("Cluster_composition-Detailed_treatment_response", script_version, ".png"), detailed_response_group_bar_p,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    rm(detailed_response_group_bar_p)
    
    
    # Plotting the surface marker distribution from the cluster view
    surface_m_bar_p <- ggplot(Surface_marker_summary_df_clust_long,
                              aes(x = Clusters, y = Percentages, fill = as.factor(Markers))) +
      geom_col(color = "black", position = "dodge") +
      geom_text(aes(label = round(Percentages, digits = 2)), vjust = -0.5, position = position_dodge(width = 0.9)) +
      scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
      scale_x_discrete(labels = c("Cluster 0", "Cluster 1", "Cluster 2")) +
      scale_fill_manual(values = cluster_colors[c(4, 6, 7)], label = c("EpCAM+", "EpCAM+PSMA+", "PSMA+")) +
      labs(fill = "Treatment cycles", x = expression("CTC clusters"), y = expression("Cell number percentages")) +
      #ggtitle("CTC clusters - cell distribution based on surface markers") +
      #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
            strip.text = element_blank(), text = element_text(size = 24))
    print(surface_m_bar_p)
    
    ggsave(filename = paste0("Cluster_composition-Surface_markers", script_version, ".png"), surface_m_bar_p,
           device = "png", path = paste0(analysis_dir, "/", "Plots/"),
           width = 3000, height = 1800, units = "px", dpi = 320)
    rm(surface_m_bar_p)
    
    
    # Plotting message
    message("The requested plots have been printed and saved into the: ", paste0(analysis_dir, "/", "Plots/"), "folder. \n")
    
  } else {
    
    # No plotting message
    message("No plots were requested, therefore no plots will be printed or saved.")
    
  }
  
  
}
calculate_and_plot_cluster_composition(make_plots = TRUE, return_statistics = TRUE, print_statistics = TRUE)

#################################################################   End Section  ##################################################################






#################################################################  Condition influence on clustering  ##################################################################
                                                                #       (Statistical testing)         #




## Calculating the association between various conditions and clustering using the Chi square test


# This function will simply calculate if the association between the culstering and various conditions is significant or not (the purpose of this function is
# to make this code block more concise)

# This is the Chi-square test of independence 
# NOTE: it takes the previously created stats list as the main data input and looks at the same conditions (Treatment cycles, treatment types, cell cycle phase
# overall and detailed treatment response)
calculating_association_stats = function (stat_list, return_contingency_tables = TRUE, return_statisctics = TRUE, print_statistics = TRUE) {
  
  ## Creating contingency tables
  # NOTE: cols should be conditions and rows should be outcomes
  
  
  # Treatment cycles
  treatment_cyc_ct <- matrix(ncol = 3, nrow = 2)
  colnames(treatment_cyc_ct) <- c("Cycle_0", "Cycle_1", "Cycle_2")
  rownames(treatment_cyc_ct) <- c("Cluster_0", "Cluster_1")
  treatment_cyc_ct[1, ] <- stat_list[[1]][[1]]$Cell_numbers_c0
  treatment_cyc_ct[2, ] <- stat_list[[1]][[1]]$Cell_numbers_c1
  
  
  # Treatment type
  treatment_type_ct <- matrix(ncol = 2, nrow = 2)
  colnames(treatment_type_ct) <- c("AL", "L")
  rownames(treatment_type_ct) <- c("Cluster_0", "Cluster_1")
  treatment_type_ct[1, ] <- stat_list[[1]][[2]]$Cell_numbers_c0[c(1, 2)]
  treatment_type_ct[2, ] <- stat_list[[1]][[2]]$Cell_numbers_c1[c(1, 2)]
  
  
  # Cell cycle phases
  cell_cyc_ct <- matrix(ncol = 3, nrow = 2)
  colnames(cell_cyc_ct) <- c("G1", "S", "G2M")
  rownames(cell_cyc_ct) <- c("Cluster_0", "Cluster_1")
  cell_cyc_ct[1, ] <- stat_list[[1]][[3]]$Cell_numbers_c0
  cell_cyc_ct[2, ] <- stat_list[[1]][[3]]$Cell_numbers_c1
  
  
  # Overall treatment response
  overall_resp_ct <- matrix(ncol = 2, nrow = 2)
  colnames(overall_resp_ct) <- c("Responder", "Nonresponder")
  rownames(overall_resp_ct) <- c("Cluster_0", "Cluster_1")
  overall_resp_ct[1, ] <- stat_list[[1]][[4]]$Cell_numbers_c0
  overall_resp_ct[2, ] <- stat_list[[1]][[4]]$Cell_numbers_c1
  
  
  # Overall treatment status
  treatment_status_ct <- matrix(ncol = 2, nrow = 2)
  colnames(treatment_status_ct) <- c("Non-treated", "Treated")
  rownames(treatment_status_ct) <- c("Cluster_0", "Cluster_1")
  treatment_status_ct[1, ] <- stat_list[[1]][[5]]$Cell_numbers_c0
  treatment_status_ct[2, ] <- stat_list[[1]][[5]]$Cell_numbers_c1
  
  
  # Detailed treatment response
  detailed_resp_ct <- matrix(ncol = 3, nrow = 2)
  colnames(detailed_resp_ct) <- c("PD", "PR", "SD")
  rownames(detailed_resp_ct) <- c("Cluster_0", "Cluster_1")
  detailed_resp_ct[1, ] <- stat_list[[1]][[6]]$Cell_numbers_c0[c(2:4)]
  detailed_resp_ct[2, ] <- stat_list[[1]][[6]]$Cell_numbers_c1[c(2:4)]
  
  
  # Surface marker distribution
  surface_marker_ct <- matrix(ncol = 3, nrow = 2)
  colnames(surface_marker_ct) <- c("EPCAM+PSMA-", "EPCAM+PSMA+", "EPCAM-PSMA+")
  rownames(surface_marker_ct) <- c("Cluster_0", "Cluster_1")
  surface_marker_ct[1, ] <- stat_list[[1]][[7]]$Cell_numbers_c0
  surface_marker_ct[2, ] <- stat_list[[1]][[7]]$Cell_numbers_c1
  
  
  
  # Compile the tables together into a list
  contingency_lst <- list (treatment_cyc_ct, treatment_type_ct, cell_cyc_ct,
                           overall_resp_ct, treatment_status_ct, detailed_resp_ct,
                           surface_marker_ct)
  names(contingency_lst) <- c("treatment_cyc_ct", "treatment_type_ct", "cell_cyc_ct",
                              "overall_resp_ct", "treatment_status_ct", "detailed_resp_ct",
                              "surface_marker_ct")
  
  
  ## Calculate the association with the Chi-square test based on the
  ## dimensions of the contingency tables
  association_stat_res <- vector(mode = "list", length = length(contingency_lst))
  names(association_stat_res) <- c("treatment_cyc_ct", "treatment_type_ct", "cell_cyc_ct",
                                   "overall_resp_ct", "treatment_status_ct", "detailed_resp_ct",
                                   "surface_marker_ct")
  for (i in seq_along(contingency_lst)) {
    
    association_stat_res[[i]] <- chisq.test(contingency_lst[[i]])
    
  }
  
  
  ## This if statement controls if a list of contingency tables is returned
  if (return_contingency_tables == TRUE) {
    
    assign("contingency_tables_lst", contingency_lst, envir = .GlobalEnv)
    message("The contingency tables are assigned to a new variable: contingency_tables_lst")
    
  } else {
    
    message("The return of the contingency tables was not requested.")
    
  }
  
  
  ## This if statement controls if a new list of statistics is returned
  if (return_statisctics == TRUE) {
    
    assign("association_stat_res", association_stat_res, envir = .GlobalEnv)
    message("The stat results are assinged to a new variable: association_stat_res")
    
  } else {
    
    message("The return of the stat results was not requested.")
    
  }
  
  
  ## This if statement controls if the stats and contingency tables are printed
  if (print_statistics == TRUE) {
    
    print(contingency_lst)
    print(association_stat_res)
    
  } else {
    
    message("As requested, the contingency tables and stat results will not be printed.")
    
  }
  
}
calculating_association_stats(stat_list = stat_lists, return_contingency_tables = TRUE, return_statisctics = TRUE, print_statistics = TRUE)


# Save the results (use capture.output of print to get the nice format)
save_association_stat_results = function (condition_infulance_stat_lst) {
  
  result_to_save <- c()
  for (i in seq_along(condition_infulance_stat_lst)) {
    
    result_to_save <- capture.output(print(condition_infulance_stat_lst[[i]]))
    write_lines(x = result_to_save, file = paste0(analysis_dir, "/", "Association_test_sesult", "_",
                                                  names(condition_infulance_stat_lst[i]), script_version, ".txt"))
    
  }
  
}
save_association_stat_results(association_stat_res)


# Remove unnecessary variables
rm(stat_lists, association_stat_res)

#################################################################   End Section  ##################################################################






#################################################################   Differential gene expression analysis  ##################################################################
                                                                #                                          #


# Finding differentially expressed (DE) genes using the FindMarkers function to compare two
# clusters (if the Log2FC is positive, then it is upregulated in ident.1 vs ident.2 if it is negative then it is downregulated ident.1 vs ident.2)
DE_Cluster_0v1 <- FindMarkers(sct_CTC.obj, ident.1 = "0", ident.2 = "1", assay = "SCT")
head(DE_Cluster_0v1)


# This is the parallelized version of the gene name matching function
# NOTE: the original geneName_matching_par does not work with the results generated by FindAllMarkers, as the rownames are altered
# if by .1, .2 and so on if DEGs are matching between groups. Therefore I updated the function with a selectable DE_lst_ID column.

geneName_matching_par_V2 = function (DE_lst, gene_lst, ID_col, name_col, n_cores = 2, DE_lst_IDs_are_rows = TRUE, DE_lst_IDs = NULL) {
  # initializing the progress bar
  progressBar <- progress::progress_bar$new(format = "Progress = (:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                            total = nrow(DE_lst),
                                            complete = "=",
                                            incomplete = "-",
                                            current = ">",
                                            clear = FALSE,
                                            width = 100)
  prog_range <- seq_len(nrow(DE_lst))
  progress = function(n) {
    progressBar$tick(tokens = list(Progress = prog_range[n]))
  }
  
  # This section is froeach specific and is required in order to visualize the progress bar
  # with the foreach function
  opts <- list(progress = progress)
  
  # This part defines the number of cores which should be used and registers them as sockets
  cores <- n_cores
  clust <- parallel::makeCluster(cores)
  doSNOW::registerDoSNOW(cl = clust)
  
  # Here I define the forach package and am trying to import hte dopar keyword from it
  # as otherwise it does not seem to work in VS Code
  #' @importFrom foreach %dopar%
  
  # Here I set up the foreach loop (which should be able to do the parallellized run)
  position <- c()
  out <- c()
  
  if (DE_lst_IDs_are_rows == TRUE) {
    
    geneNames <- foreach::foreach(i = seq_len(nrow(DE_lst)), .options.snow = opts, .combine = c) %dopar% {
      position <- grep(rownames(DE_lst)[i], gene_lst[, ID_col])
      out[i] <- gene_lst[, name_col][position]
    }
    
  } else {
    
    geneNames <- foreach::foreach(i = seq_len(nrow(DE_lst)), .options.snow = opts, .combine = c) %dopar% {
      position <- grep(DE_lst[, DE_lst_IDs][i], gene_lst[, ID_col])
      out[i] <- gene_lst[, name_col][position]
    }    
    
  }
  
  
  # Closing the parallel sockets
  parallel::stopCluster(cl = clust)
  
  # Closing the function with assigning the new object and naming it
  out_lst <- cbind(DE_lst, geneNames)
  original_name <- deparse(substitute(DE_lst))
  output_name <- paste0(original_name, "_", "geneNames")
  assign(output_name, out_lst, envir = .GlobalEnv)
  message("The gene name matched DE_cluster_0v1 list was assigned to a new object named:", "\n", output_name)
}
geneName_matching_par_V2(DE_lst = DE_Cluster_0v1, gene_lst = unif_gene_names, ID_col = "ensembl_gene_id", name_col = "external_gene_name", n_cores = 5)

head(DE_Cluster_0v1_geneNames)


# Save the gene name matched DEG list
write_csv(x = DE_Cluster_0v1_geneNames, file = paste0(analysis_dir, "/", "DE_Cluster_0v1_geneNames", script_version, ".csv"))


# Checking if there is differential expression between the resp groups (not according to the clustering)
Idents(sct_CTC.obj) <- "RespGroup"
DE_RespvNonresp <- FindMarkers(sct_CTC.obj, ident.1 = "Responder", ident.2 = "Nonresponder", assay = "SCT")

head(DE_RespvNonresp)
summary(DE_RespvNonresp$p_val_adj)


# Gene name matching with treatment cycle DEG markers
geneName_matching_par_V2(DE_lst = DE_RespvNonresp, gene_lst = unif_gene_names, ID_col = "ensembl_gene_id", name_col = "external_gene_name", n_cores = 5)
head(DE_RespvNonresp_geneNames)


# Saving the list of RespvsNonresp DEG markers
write_csv(x = DE_RespvNonresp_geneNames, file = paste0(analysis_dir, "/", "DE_RespvNonresps_geneNames", script_version, ".csv"))

# Remove unnecessary objects
rm(DE_RespvNonresp, CON_Cluster_0v1, CON_Cluster_0v1_geneNames)




## Assigning ENTREZ IDs to the DE gene list for future utility and eas of use


# This if Statment will check if the final version of the UnifiedEntrezIds df is present in the main working directory. If yes,
# the file will be loaded, if not the processing will be ran.
if (exists("unifiedEntrezIDs", where = .GlobalEnv)) {
  
  message("The unifiedEntrezIDs variable is already present in the .GlobalEnv.")  
  
} else if (file.exists(paste0("Additional_data_files", "/", "UnifiedEntrezIDs_final.csv")) == TRUE) {
  
  unifiedEntrezIDs <- read.csv(file = paste0("Additional_data_files", "/", "UnifiedEntrezIDs_final.csv"))
  message("The final, fully prepared version of the ENTREZ IDs list is present in the working directory, \n",
          "and will be loaded as: ", "unifiedEntrezIDs")
  
} else {
  
  message("The final, fully prepared version of the ENTREZ IDs list was not found in the working directory, \n",
          "Running the following code block to prpare it. \n")
  
  
  # This if statement will check if the entrezID table is already present in the working directory and if yes it will load it
  if (file.exists(paste0("Additional_data_files", "/", "Symbols_ENSEMBL_ENTREZ_ID_table.csv")) & !exists(x = "entrezIDs", where = .GlobalEnv)) {
    
    entrezIDs <- read.csv(file = paste0("Additional_data_files", "/", "Symbols_ENSEMBL_ENTREZ_ID_table.csv"))
    message("The geneNames, entrezIDs and ensemblIDs containing table has been loaded as a variable: ", deparse(substitute(entrezIDs)))
    
  } else {
    
    # Retrieving ENTREZ IDs using the ensembl IDs
    listEnsembl() #lists the available biomart keynames, here I need "gene"
    ensembl <- useEnsembl(biomart = "genes")
    datasets <- listDatasets(ensembl) #lists the organism based datasest from which we need to chose one for the db connection
    ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl") #creates a mart object used for the queries 
    attributes <- listAttributes(ensembl) #needed for the query, describes the desired info like gene names
    filters <- listFilters(ensembl) # needed for the query, defined by upon which data we search for(?) like here the ensemblIds
    
    
    # Writing the database query using the previous settings
    entrezIDs <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name", "uniprot_gn_symbol"),
                       filters = "ensembl_gene_id",
                       values = geneIDs,
                       mart = ensembl)
    
    # NOTE: as with the gene names, some IDs are missing, therefore I will try to runt he missing IDs through another database and see if they map onto
    # ENTREZ IDs or not
    # Additionally not every ENSEMBL ID returns an NA, so I will have to check if some were dropped and add them back
    
    
    # Save the created table into the working directory
    write_csv(x = entrezIDs, file = paste0(analysis_dir, "/", "Symbols_ENSEMBL_ENTREZ_ID_table.csv"))
    
    message("The entrezID variable was created and saved into the working directory.")
    
  }
  
  
  
  # Check for the dropped ENSEMBL IDs
  dropped_IDs <- geneIDs[!geneIDs %in% unique(entrezIDs$ensembl_gene_id)]
  
  
  # Add back these ENSEMBL IDs with an NA to ENTREZ IDS (and other IDs if present)
  dropped_IDs_df <- data.frame(ensembl_gene_id = dropped_IDs,
                               entrezgene_id = NA,
                               external_gene_name = NA,
                               uniprot_gn_symbol = NA)
  
  
  # I will bind the dropped ENSEMBL IDs to the main entrezIDs list
  entrezIDs <- rbind(entrezIDs, dropped_IDs_df)
  
  
  # NOTE: during retrieval, a non matching ID/gene name will be left empty, which makes later processing more difficult
  # To amend this, this function will fill the empty positions found int he ENTREZ ID dataframe with NAs for easier handling
  fill_empty_positions = function(dataframe) {
    # Initialize a progress bar
    progressBar <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                              total = nrow(dataframe),
                                              complete = "=",
                                              incomplete = "-",
                                              current = ">",
                                              clear = FALSE,
                                              width = 100)
    
    # This nested loop will look at each position and determines if there is an NA (then it does nothing)
    # or if there a missing value, then it dilles it with NA
    for(i in seq_len(nrow(dataframe))) {
      for(col in 2:4) {
        
        if (!is.na(dataframe[i, col]) && stringi::stri_isempty(dataframe[i, col])) {
          dataframe[i, col] <- NA
          
        }
      } 
      # Start the progress bar
      progressBar$tick()
    }
    
    return(dataframe)
    message("The missing positions have been filled with NAs.")
  }
  entrezIDs <- fill_empty_positions(entrezIDs)
  
  # NOTE: As the functions runs with a nested for loop it is time and resource intensive.
  # I tried parallelizing it, but gained no performance boost.
  
  
  # As the previous step take a considerable amount of time, I will save the output so it can be loaded if needed
  write_csv(x = entrezIDs, file = paste0(analysis_dir, "/", "ENSEMBL_ID_list_with_ENTREZ_IDs_and_gene_names", script_version, ".csv"))
  
  
  # Load the saved entrezIDs df if not loaded yet
  if (!exists(x = "entrezIDs", where = .GlobalEnv)) {
    entrezIDs <- read.csv(file = paste0("Additional_data_files", "/", "ENSEMBL_ID_list_with_ENTREZ_IDs_and_gene_names", script_version, ".csv"))
    message("The ENTREZ ID list is loaded as the following object: ", deparse(substitute(entrezIDs)))
  }
  
  
  # First I will sub-select the missing ENTREZ IDs 
  missing_entrezIDs <- entrezIDs[is.na(entrezIDs$entrezgene_id), ]
  
  
  # Next I wil use the AnnotationDbi package and the org.Hs.eg.db database to select the proper keys and columns for the search
  columns(org.Hs.eg.db)
  keytypes(org.Hs.eg.db)
  
  
  # Now I will make the call based on the missing IDs df 
  entrezIDs_p2 <- AnnotationDbi::select(org.Hs.eg.db, keys = missing_entrezIDs$ensembl_gene_id,
                                        columns = c("ENSEMBL", "ENTREZID", "SYMBOL", "UNIPROT"),
                                        keytype = "ENSEMBL")
  
  
  # This function will merge the two parts to form a full ID list
  merge_ID_lists = function (ID_lst_1, ID_lst_2, new_obj_name = "unifiedEntrezIDs") {
    
    # Fist assign the non NA elements of the first list (main list) to a new temporary list
    tmpLstP1 <- dplyr::filter(ID_lst_1, entrezgene_id != "NA")
    
    # In order to bind the second list with the first one they need to have the same colnames
    # so I assign the second list to a temporary list and change it's colnames to match the colnames
    # of the first list
    tmpLstP2 <- ID_lst_2
    colnames(tmpLstP2) <- colnames(tmpLstP1)
    
    # Unify the two lists by binding them torgwther by rows
    unifiedEntrezIDs <- rbind(tmpLstP1, tmpLstP2)
    
    # Name the new object and assign it into the global environment
    objectName <- new_obj_name
    assign(objectName, unifiedEntrezIDs, envir = .GlobalEnv)
    message("The unified ID lists were assigned to the .GlobalEnv ans a new object: ", objectName)
  }
  merge_ID_lists(entrezIDs, entrezIDs_p2)
  
  
  # Inspect the new unified ENTREZ ID list and see how many NAs it still has
  head(unifiedEntrezIDs)
  summary(is.na(unifiedEntrezIDs$entrezgene_id))
  
  
  # NOTE: during mapping, it happens that multiple ENTREZ IDs are mapped to the same ENSEMBL ID
  # so in order to check if these multi mappings are also covering different genes I will write a function
  # which sub-selects the multi mapped IDs and returns a unique names list where the names are ENSEMBL IDs
  # while the elements are gene names (if available). This will help me determine if one ENSEMBL ID maps
  # onto multiple genes or not, and if not, I can randomly drop the multiple ENTREZ IDs mapped to a single
  # ENSEMBL ID, as they will be covering the same gene
  check_multimapped_IDs = function(multi_ID_dataframe) {
    # Let's count the number of times each Ensembl ID appears in `Ensembl` column
    multi_mapped <- dplyr::count(multi_ID_dataframe, ensembl_gene_id, name = "entrez_id_count")
    
    # Arrange by the genes with the highest number of Entrez IDs mapped
    multi_mapped <- dplyr::arrange(multi_mapped, dplyr::desc(entrez_id_count))
    
    # Filter out the single mapped entries, leaving only the multi-mapped ones
    multi_mapped <- dplyr::filter(multi_mapped, entrez_id_count > 1)
    
    # Initialize the progress bar
    progressBar <- progress::progress_bar$new(format = "(:spin) [:bar] [:percent] [Elapsed time: :elapsedfull] || Estimated time remaining: :eta",
                                              total = nrow(multi_mapped),
                                              complete = "=",
                                              incomplete = "-",
                                              current = ">",
                                              clear = FALSE,
                                              width = 100)
    
    #loop over the dataframe and grep each NESEMBL ID and the associated unique external gene names
    duplicated_genes <- list()
    for(element in multi_mapped$ensembl_gene_id) {
      # Start the progress bar
      progressBar$tick()
      
      duplicated_genes[[element]] <- unique(multi_ID_dataframe$external_gene_name[grep(element, multi_ID_dataframe$ensembl_gene_id)])
    }
    return(duplicated_genes)
  }
  multimapped_genes <- check_multimapped_IDs(unifiedEntrezIDs)
  
  
  # NOTE: after inspecting the resulting named list, it seems to me that there is no multi mapping on the gene name
  # level, therefore a random drop of entrezIDs (basically just calling unique on the ENSEMBL IDs) should be safe to do
  unifiedEntrezIDs <- unifiedEntrezIDs[!duplicated(unifiedEntrezIDs$ensembl_gene_id), ]
  
  
  # Save the final unified entrezID dataframe for further use
  write_csv(x = unifiedEntrezIDs, file = paste0(analysis_dir, "/", "UnifiedEntrezIDs_final.csv"))
  
}


# I will combine the gene name and the ENTREZ ID DE lists by feeding the gene name list into the 
# ENTREZ ID matching function for the cluster DEG list
entrezID_matching_par_V2 = function (DE_lst, entrezID_lst, ENSEMBL_ID_col, ENTREZ_ID_col, n_cores = 2, new_object_name = NA,
                                     DE_lst_IDs_are_rows = TRUE, DE_lst_IDs = NULL) {
  # Initializing the progress bar
  progressBar <- progress::progress_bar$new(format = "Progress = (:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                            total = nrow(DE_lst),
                                            complete = "=",
                                            incomplete = "-",
                                            current = ">",
                                            clear = FALSE,
                                            width = 100)
  
  # This part is necessary for the proper functioning of the progress bar together with the foreach function
  prog_range <- seq_len(nrow(DE_lst))
  progress = function(n) {
    progressBar$tick(tokens = list(Progress = prog_range[n]))
  }
  
  # This section is foreach specific and is required in order to visualize the progress bar
  # with the foreach function
  opts <- list(progress = progress)
  
  # This part defines the number of cores which should be used and registers them as sockets
  cores <- n_cores
  clust <- parallel::makeCluster(cores)
  doSNOW::registerDoSNOW(cl = clust)
  
  # Here I define the foreach package and am trying to import hte dopar keyword from it
  # as otherwise it does not seem to work in VS Code
  #' @importFrom foreach %dopar%
  
  # Here I set up the foreach loop (which should be able to do the parallelized run)
  position <- c()
  out <- c()
  
  if (DE_lst_IDs_are_rows == TRUE) {
    
    entrezIDs <- foreach::foreach(i = seq_len(nrow(DE_lst)), .options.snow = opts, .combine = c) %dopar% {
      position <- grep(rownames(DE_lst)[i], entrezID_lst[, ENSEMBL_ID_col])
      out[i] <- entrezID_lst[, ENTREZ_ID_col][position]
    }
    
  } else {
    
    entrezIDs <- foreach::foreach(i = seq_len(nrow(DE_lst)), .options.snow = opts, .combine = c) %dopar% {
      position <- grep(DE_lst[, DE_lst_IDs][i], entrezID_lst[, ENSEMBL_ID_col])
      out[i] <- entrezID_lst[, ENTREZ_ID_col][position]
    }    
    
  }
  
  # Closing the parallel sockets
  parallel::stopCluster(cl = clust)
  
  # Closing the function with assigning the new object and naming it
  out_lst <- cbind(DE_lst, entrezIDs)
  
  # This new if statment allows the user to specify a new object anme, but could just go with the standard nameing
  # the function uses
  if(is.na(new_object_name)) {
    original_name <- deparse(substitute(DE_lst))
    output_name <- paste0(original_name, "_", "entrezIDs")
  } else {
    output_name <- new_object_name
  }
  
  assign(output_name, out_lst, envir = .GlobalEnv)
  message("The gene name matched dataframe was assigned to a new object named:", "\n", output_name)
}
entrezID_matching_par_V2(DE_Cluster_0v1_geneNames, unifiedEntrezIDs, "ensembl_gene_id", "entrezgene_id", n_cores = 5, new_object_name = "DE_Cluster_0v1_geneNames_entrezIDs")
head(DE_Cluster_0v1_geneNames_entrezIDs)




## Further processing the DE gene list (continuing with the geneName + ENTREZ ID list)
## NOTE: as the Fisher's test showed, the two DE lists are not significantly different
## therefore I will continue only with the cluster DEG list


# Attaching some additional info to the cluster DE gene list and saving the whole list
DE_Cluster_0v1_geneNames_entrezIDs$ExprChange <- ifelse(DE_Cluster_0v1_geneNames_entrezIDs$avg_log2FC >= 0 & DE_Cluster_0v1_geneNames_entrezIDs$p_val_adj <= 0.05, "Upregulated",
                                                        ifelse(DE_Cluster_0v1_geneNames_entrezIDs$avg_log2FC <= 0 & DE_Cluster_0v1_geneNames_entrezIDs$p_val_adj <= 0.05, "Downregulated", "Not significant"))

# NOTE: the ranking here will be solely based on the p_adj value (as the sign() function will turn all log2FC into a negative, 0 or positive 1)
# if you don't want that you can omit the sign() function
DE_Cluster_0v1_geneNames_entrezIDs$Ranking <- DE_Cluster_0v1_geneNames_entrezIDs$avg_log2FC * (-log10(DE_Cluster_0v1_geneNames_entrezIDs$p_val_adj))


# Save the resulting cluster DEG list
write_csv(x = DE_Cluster_0v1_geneNames_entrezIDs, file = paste0(analysis_dir, "/", "DEG_clusters_0v1", script_version, ".csv"))




## As for further analysis, I will be mainly using the gene names and ENTREZ IDs, I will remove genes from the DEG list which have no ENTREZ ID or gene name


# Load the extended DEG list if not loaded yet
if (!exists(x = "DE_Cluster_0v1_geneNames_entrezIDs", where = .GlobalEnv)) {
  DE_Cluster_0v1_geneNames_entrezIDs <- read.csv(file = paste0(analysis_dir, "/", "DEG_clusters_0v1", script_version, ".csv"))
  message("The extended DEG list was loaded under the name: ", deparse(substitute(DE_Cluster_0v1_geneNames_entrezIDs)))
}


# Filter the original list and assign it to a new, clean variable by entrezIDs
DEG_cluster_0v1_clean <- dplyr::filter(DE_Cluster_0v1_geneNames_entrezIDs, entrezIDs != "NA")


# Filter the original list and assign it to a new, clean variable by geneNames
DEG_cluster_0v1_symbols_clean <- dplyr::filter(DE_Cluster_0v1_geneNames_entrezIDs, geneNames != "NA")


# Save the new, clean DEG lists for future use
write_csv(x = DEG_cluster_0v1_clean, file = paste0(analysis_dir, "/", "DEG_clusters_0v1_clean", script_version, ".csv"))

write_csv(x = DEG_cluster_0v1_symbols_clean, file = paste0(analysis_dir, "/", "DEG_clusters_0v1_symbols_clean", script_version, ".csv"))


# Remove unnecessary objects
rm(DE_Cluster_0v1, DE_Cluster_0v1_entrezIDs, DE_Cluster_0v1_geneNames, DE_Cluster_0v1_geneNames_entrezIDs)




## Visualizing DE genes using a volcano plot


# Load the DEG list if not loaded yet
if (!exists(x = "DEG_cluster_0v1_symbols_clean", where = .GlobalEnv)) {
  DEG_cluster_0v1_symbols_clean <- read.csv(file = paste0(analysis_dir, "/", "DEG_clusters_0v1_symbols_clean", script_version, ".csv"))
  message("The clean DEG lists were loaded under the names: ", deparse(substitute(DE_Cluster_0v1_geneNames_entrezIDs)), " and ", deparse(substitute(DE_Cluster_0v1_geneNames_entrezIDs)))
}


# Annotated visualization with ggplot2 for geneNames
DEG_cluster_0v1_symbols_clean$Significance <- ifelse(DEG_cluster_0v1_symbols_clean$avg_log2FC > 1 & DEG_cluster_0v1_symbols_clean$p_val_adj < 0.05, "Upregulated",
                                                     ifelse(DEG_cluster_0v1_symbols_clean$avg_log2FC < -1 & DEG_cluster_0v1_symbols_clean$p_val_adj < 0.05, "Downregulated", "Not significant")) #labeling genes for coloring
upreg <- DEG_cluster_0v1_symbols_clean[DEG_cluster_0v1_symbols_clean$Significance == "Upregulated", ] #sub-setting the upregulated genes
topupreg <- head(dplyr::arrange(upreg, desc(Ranking))$geneNames, 15) #selecting the top upreg genes
downreg <- DEG_cluster_0v1_symbols_clean[DEG_cluster_0v1_symbols_clean$Significance == "Downregulated", ] #sub-setting the downregulated genes
topdownreg <- head(dplyr::arrange(downreg, Ranking)$geneNames, 15) #selecting the top downreg genes

DEG_cluster_0v1_symbols_clean$DElabel <- ifelse(DEG_cluster_0v1_symbols_clean$geneNames %in% topupreg | DEG_cluster_0v1_symbols_clean$geneNames %in% topdownreg,
                                                DEG_cluster_0v1_symbols_clean$geneNames, NA) #marking the top DEGs in the dataframe

volcano_p <- ggplot(data = DEG_cluster_0v1_symbols_clean , 
                    aes(x = avg_log2FC, y = -log10(p_val_adj), col = Significance, label = DElabel)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_point(data = subset(DEG_cluster_0v1_symbols_clean, !is.na(DElabel)),
             shape = 21, color = "black", stroke = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = -1, linetype = "dashed") +
  scale_color_manual(values = cluster_colors[c(2, 5, 3)],
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  scale_x_continuous(limits = c(-22, 22), breaks = seq(-22, 22, 2)) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 2)) +
  labs(color = "Expression status", x = expression("log"[2] * "FC"), y = expression("-log"[10] * "p-value")) +
  #ggtitle("Differentially expressed genes - geneNames") +
  geom_text_repel(show.legend = FALSE, max.overlaps = Inf, color = "black", size = 5) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
print(volcano_p)

ggsave(filename = paste0("DEGs_Volcano_plot_symbols", script_version, ".png"), volcano_p,
       device = "png", path = paste0(analysis_dir, "/", "Plots/"),
       width = 4800, height = 2400, units = "px", dpi = 320)
rm(volcano_p)

# Construct an upreg-downreg df of gene names for saving
topupreg <- head(dplyr::arrange(upreg, p_val_adj), 50) #selecting the top upreg genes
topdownreg <- head(dplyr::arrange(downreg, Ranking), 50) #selecting the top downreg genes
top_upreg_and_downreg_DEGs_symbols <- data.frame(top_upregulated = topupreg, top_downregulated = topdownreg)


# Save the upreg-downreg df of gene names
write_csv(x = top_upreg_and_downreg_DEGs_symbols, file = paste0(analysis_dir, "/", "top_upreg_and_downreg_DEGs_symbols", script_version, ".csv"))


# Remove unnecessary variables
rm(upreg, downreg, top50upreg, topupreg, topdownreg, top_upreg_and_downreg_DEGs, top_upreg_and_downreg_DEGs_symbols)

#################################################################   End Section  ##################################################################







#################################################################               irGSEA                  ##################################################################
                                                                #   (Gene Set Enrichment Analysis)    #




## A wrapper function for the irGSEA.score function which will run it on multiple entries
## NOTE: only the options I would change are part of the wrapper arguments
run_serial_irGSEA = function(seurat_obj, seurat_assay = "SCT", seurat_slot = "data", category_lst, subcategory_lst, object_names, n_cores = 5) {
  #load the required package
  require(package = irGSEA)
  
  #define the output list
  output_list <- vector(mode = "list", length = length(category_lst))
  
  #run the irGSEA.score function on all elements of the category/subcategory lists
  for (item in seq_along(category_lst)) {
    output_list[[item]] <- irGSEA::irGSEA.score(object = seurat_obj, assay = seurat_assay, 
                                                slot = seurat_slot, seeds = 123, ncores = n_cores,
                                                min.cells = 1, min.feature = 0,
                                                custom = FALSE, geneset = NULL, msigdb = TRUE, 
                                                species = "Homo sapiens", category = category_lst[[item]],  
                                                subcategory = subcategory_lst[[item]], geneid = "ensembl",
                                                method = c("AUCell", "UCell", "singscore", 
                                                           "ssgsea", "JASMINE", "viper"),
                                                aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                                                kcdf = 'Gaussian')
    
  }
  
  #assign the proper names to the irGSEA list objects
  names(output_list) <- object_names
  
  #return the output list
  return(output_list)
  
}




## A wrapper function to save the seurat objects generated by irGSEA.score
## NOTE: it was tailored for interactive use, therefore it uses the readline function instead of the scan function to read user input
save_irGSEA_list = function(irGSEA_list, directory, filename) {
  
  #this function will allow to visualize a progression bar
  progressBar <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                            total = length(irGSEA_list),
                                            complete = "=",
                                            incomplete = "-",
                                            current = ">",
                                            clear = FALSE,
                                            width = 100)
  
  #An if statement to determine if the save directory exists and if not what should be the course of action
  if (!dir.exists(directory)) {
    
    message("The target directory does not exist. Would you like to create one? \n",
            "1: Yes \n",
            "2: No")
    
    input <- readline(prompt = "Enter your choice (1 or 2): ")
  
    if (!is.null(input) && input == "1") {
      dir.create(path = directory)
      message("The directory was created, saving files...")
      
    } else if (!is.null(input) && input == "2") {
      stop("No target directory was created.")
      
    } else {
      stop("Wrong input. It has to be '1' or '2'.")
      
    }
    
    
  }
  
  #A for loop responsible for saving the input list elements
  for (element in seq_along(irGSEA_list)) {
    progressBar$tick()
    
    saveRDS(object = irGSEA_list[[element]], file = paste0(directory, "/", filename, "_", names(irGSEA_list)[element], script_version, ".rds"))
    
  }
  
}




## A wrapper function to save the seurat objects generated by irGSEA.score
## NOTE: it was tailored for interactive use, therefore it uses the readline function instead of the scan function to read user input
load_irGSEA_list = function(directory, pattern = NULL) {
  #Define a vetor to hold the names of the input files
  input_files <- NULL
  
  #Define the input file list
  if (is.null(pattern)) {
    input_files <- list.files(path = directory, include.dirs = FALSE)
    
  } else {
    input_files <- list.files(path = directory, pattern = pattern, include.dirs = FALSE)
    
  }
  
  #this function will allow to visualize a progression bar
  progressBar <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                            total = length(input_files),
                                            complete = "=",
                                            incomplete = "-",
                                            current = ">",
                                            clear = FALSE,
                                            width = 100)
  
  #Define the output list
  output_list <- vector(mode = "list", length = length(input_files))
  
  
  for (e in seq_along(input_files)) {
    progressBar$tick()
    
    output_list[[e]] <- readRDS(file = paste0(directory, "/", input_files[e]))
    
  }
  
  #Extract the filenames for the list
  #NOTE: the (?<=) operator is a special operator called look_behind, this will define the part of the string
  #after which the search will start, while the second part [^_V]+ will give me the end boundary
  filenames <- stringr::str_extract(string = input_files, pattern = "(?<=GSEA_)[^_V]+")
  
  #Name the list elements 
  names(output_list) <- filenames
  
  #Return the output list
  return(output_list)
  
  
}




## A wrapper function to save the seurat objects generated by irGSEA.hub
## NOTE: only the options I would change are part of the wrapper arguments
run_serial_hub_gene_search = function(irGSEA_list, seurat_assay = "SCT", seurat_slot = "data", irGSEA_methods = c("AUCell", "UCell", "singscore", "ssgsea", "JASMINE", "viper"), irGSEA_pathways, n_cores = 5,
                                      correlation_type = "rank", correlation_color = c("#0073c2","white","#efc000")) {
  #load the required package
  require(package = irGSEA)
  
  #define the output list
  output_list <- vector(mode = "list", length = length(irGSEA_list))
  
  #this function will allow to visualize a progression bar
  progressBar <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                            total = length(irGSEA_list),
                                            complete = "=",
                                            incomplete = "-",
                                            current = ">",
                                            clear = FALSE,
                                            width = 100)
  
  #run the irGSEA.score function on all elements of the category/subcategory lists
  for (item in seq_along(irGSEA_list)) {
    progressBar$tick()
    
    output_list[[item]] <- irGSEA::irGSEA.hub(irGSEA_list[[item]],
                                              assay = seurat_assay,
                                              method = irGSEA_methods,
                                              show.geneset = irGSEA_pathways[[item]],
                                              ncores = 5,
                                              type = correlation_type,
                                              correlation.color = correlation_color)
    
  }
  
  #assign the proper names to the irGSEA list objects
  names(output_list) <- names(irGSEA_list)
  
  #return the output list
  return(output_list)
  
}



## A wrapper function for the irGSEA.integrate function which will run it on multiple entries
## NOTE: only the options I would change are part of the wrapper arguments
run_serial_integration <- function(irGSEA_list, methods = c("AUCell", "UCell", "singscore", "ssgsea", "JASMINE", "viper")) {
  # Load the required package
  require(irGSEA)
  
  # Define the output list
  output_list <- vector(mode = "list", length = length(irGSEA_list))
  
  # Run the irGSEA.integrate function on all elements of the input list
  for (element in seq_along(irGSEA_list)) {
    # Use tryCatch to handle any errors
    result <- tryCatch(
      {
        # Attempt to run the integration
        irGSEA::irGSEA.integrate(object = irGSEA_list[[element]], method = methods)
      },
      error = function(e) {
        # Handle the error
        message(paste("Error in element", element, ":", e$message))
        return(NULL)  # Return NULL in case of error
      }
    )
    # Assign the result to the output list
    output_list[[element]] <- result
  }
  
  # Assign the proper names to the output list objects
  names(output_list) <- names(irGSEA_list)
  
  # Return the output list
  return(output_list)
}




## A wrapper function for the various irGSEA visualization functions which will run it on multiple entries
run_serial_visualization = function(integrated_irGSEA_list, plotting_function = "bubble", n_show = 50, cluster_color = NULL, direction_color = NULL,
                                    significance_color = NULL, bar_method = NULL, plot_method = "RRA", method_color = NULL, clust_rows = TRUE) {
  
  #load the required package
  require(package = irGSEA)
  
  #define the output list
  output_list <- vector(mode = "list", length = length(integrated_irGSEA_list))
  
  #this function will allow to visualize a progression bar
  progressBar <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                            total = length(integrated_irGSEA_list),
                                            complete = "=",
                                            incomplete = "-",
                                            current = ">",
                                            clear = FALSE,
                                            width = 100)
  
  #run the various plotting function on all elements of the input list
  for (element in seq_along(integrated_irGSEA_list)) {
    progressBar$tick()
    
    if (!is.character(plotting_function)) {
      stop(message = "The 'plotting_function' argument must be of type character.")
      
    } else {
      
      if (plotting_function == "bubble") {
        
        output_list[[element]] <- irGSEA::irGSEA.bubble(object = integrated_irGSEA_list[[element]], top = n_show,
                                                        cluster.color = cluster_color, direction.color = direction_color,
                                                        significance.color = significance_color, cluster_rows = clust_rows,
                                                        method = plot_method)
        
      } else if (plotting_function == "heatmap") {
        
        output_list[[element]] <- irGSEA::irGSEA.heatmap(object = integrated_irGSEA_list[[element]], top = n_show,
                                                         cluster.color = cluster_color, direction.color = direction_color,
                                                         significance.color = significance_color, cluster_rows = clust_rows,
                                                         method = plot_method, heatmap.width = unit(10, "npc"), heatmap.heigh = unit(10, "npc"),
                                                         rowname.fointsize = 15)
        
      } else if (plotting_function == "barplot") {
        
        output_list[[element]] <- irGSEA::irGSEA.barplot(object = integrated_irGSEA_list[[element]], significance.color = significance_color,
                                                         color.cluster = cluster_color, method = bar_method, color.method = method_color)
        
      } else {
        
        stop(message = "The 'plotting_function' argument must be either 'bubble', 'heatmap' or 'barplot'.")
        
      }
      
    }
    
  }
  
  #assign the proper names to the output list objects
  names(output_list) <- names(integrated_irGSEA_list)
  
  #return the output list
  return(output_list)
  
}




## A function to save the irGSEA plots
save_irGSEA_plots = function(irGSEA_plot_list, directory, filename) {
  #this function will allow to visualize a progression bar
  progressBar <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                            total = length(irGSEA_plot_list),
                                            complete = "=",
                                            incomplete = "-",
                                            current = ">",
                                            clear = FALSE,
                                            width = 100)
  
  #An if statement to determine if the save directory exists and if not what should be the course of action
  if (!dir.exists(directory)) {
    
    message("The target directory does not exist. Would you like to create one? \n",
            "1: Yes \n",
            "2: No")
    
    input <- readline(prompt = "Enter your choice (1 or 2): ")
    
    if (!is.null(input) && input == "1") {
      dir.create(path = directory)
      message("The directory was created, saving files...")
      
    } else if (!is.null(input) && input == "2") {
      stop("No target directory was created.")
      
    } else {
      stop("Wrong input. It has to be '1' or '2'.")
      
    }
    
    
  }
  
  #A for loop responsible for saving the input list elements
  for (element in seq_along(irGSEA_plot_list)) {
    progressBar$tick()
    
    ggsave(filename = paste0(filename, "_", names(irGSEA_plot_list)[element], script_version, ".png"), device = "png", path = directory, plot = irGSEA_plot_list[[element]],
      width = 5000, height = 3000, units = "px", dpi = 320)
    
  }
  
  
}




## Load the irGSEA results if not present
save_dir <- paste0(getwd(), "/", analysis_dir, "/", "GSEA_results")
if (!exists("irGSEA_sct_lst", where = .GlobalEnv)) {
  irGSEA_sct_lst <- load_irGSEA_list(directory = save_dir, pattern = ".rds")
  
}




## Run the irGSEA analysis using various gene sets


# Print the full list of gene sets available inside the MSigDB
print(msigdbr::msigdbr_collections(), n = 23)


# Prep the necessary lists for the analysis
gene_cat_lst <- list("C2", "C2", "C2", "C2", "C3", "C4", "C4", "C5", "C5", "C5", "C6", "H")
gene_subcat_lst <- list("CP:BIOCARTA", "CP:KEGG", "CP:REACTOME", "CP:WIKIPATHWAYS", "TFT:GTRD", "CGN", "CM", "GO:BP", "GO:MF", "GO:CC", NULL,  NULL)
seurat_names <- c("BIOCARTA", "KEGG", "REACTOME", "WIKIPATHWAYS", "TFT:GTRD", "CGN", "CM", "GO:BP", "GO:MF", "GO:CC", "ONCOGENIC_SIGN", "HALLMARK")


# Start the analysis run by calling the irGSEA.score wrapper
irGSEA_sct_lst <- run_serial_irGSEA(seurat_obj = sct_CTC.obj, seurat_assay = "SCT", seurat_slot = "data", category_lst = gene_cat_lst,
                            subcategory_lst = gene_subcat_lst, object_names = seurat_names, n_cores = 5)


# Save the resulting irGSEA list
save_irGSEA_list(irGSEA_list = irGSEA_sct_lst, directory = save_dir, filename = "irGSEA")


# Run the hub gene analysis to find the hub genes for the interesting pathways
interesing_pathways <- list(c("HALLMARK-APOPTOSIS", "HALLMARK-PROTEIN-SECRETION", "HALLMARK-MTORC1-SIGNALING", "HALLMARK-DNA-REPAIR", "HALLMARK-APICAL-JUNCTION", "HALLMARK-COAGULATION"),
                            c("BIOCARTA-RAB-PATHWAY", "BIOCARTA-SUMO-PATHWAY", "BIOCARTA-IGF1R-PATHWAY", "BIOCARTA-NPC-PATHWAY", "BIOCARTA-CTBP1-PATHWAY", "BIOCARTA-CHEMICAL-PATHWAY"),
                            c("GOBP-POSITIVE-REGULATION-OF-PLASMINOGEN-ACTIVATION", "GOBP-POSITIVE-REGULATION-OF-RNA-POLYMERASE-II-TRANSCRIPTION-PREINITIATION-COMPLEX-ASSEMBLY",
                              "GOBP-POSITIVE-REGULATION-OF-PROTEIN-K63-LINKED-UBIQUITINATION"),
                            c("REACTOME-ER-TO-GOLGI-ANTEROGRADE-TRANSPORT", "REACTOME-TRANSPORT-TO-THE-GOLGI-AND-SUBSEQUENT-MODIFICATION", "REACTOME-UPTAKE-AND-FUNCTION-OF-DIPHTHERIA-TOXIN", "REACTOME-GRB7-EVENTS-IN-ERBB2-SIGNALING",
                              "REACTOME-SUMOYLATION-OF-DNA-REPLICATION-PROTEINS"))
hub_genes <- run_serial_hub_gene_search(irGSEA_list = irGSEA_sct_lst[c(7, 1, 4, 10)], irGSEA_pathways = interesing_pathways)
hub_genes_exp <- run_serial_hub_gene_search(irGSEA_list = irGSEA_sct_lst[c(7, 1, 4, 10)], irGSEA_pathways = interesing_pathways, correlation_type = "expression")


# Integrate the produced irGSEA.score data for visualization
irGSEA_sct_lst_integrated <- run_serial_integration(irGSEA_list = irGSEA_sct_lst)


# Plot the integrated irGSEA results
cluster_colors <- pal_tron(palette = "legacy")(7)

irGSEA_dotplots <- run_serial_visualization(integrated_irGSEA_list = irGSEA_sct_lst_integrated, plotting_function = "heatmap", cluster_color = cluster_colors[c(6, 5)],
                                            n_show = 10, clust_rows = TRUE, plot_method = "RRA")


# Save the generated irGSEA plots
plot_save_dir <- paste0(save_dir, "/", "Plots")
save_irGSEA_plots(irGSEA_plot_list = irGSEA_dotplots, directory = plot_save_dir, filename = "irGSEA_heatmap")


#################################################################   End Section  ##################################################################






#################################################################   Plotting DEGs based on enrichment data    ##################################################################
                                                                #                and literature               #




## Looking at EMT and CSC markers based on the BIORAD SAB target list H96 and other interesting markers


# Load the aforementioned lists
CSC_list <- read.csv(paste0(getwd(), "/", "Additional_data_files", "/", "CSC_SAB_target_list_H96.csv"), sep = ";")
EMT_list <- read.csv(paste0(getwd(), "/", "Additional_data_files", "/", "EMT_SAB_target_list_H96.csv"), sep = ";")
Krat_panel <- c("TP53", "CHEK2", "BRCA1", "BRCA2", "MSH2", "MSH6", "NBN", "PMS1")
sign_hub_genes <- c("COPB2", "SNX2", "AURKA", "POLR3C", "CCNO", "RAB4A", "RAB5A", "RAB9A", "RAB11A", "SUMO2", "SUMO1", "UBA2", "CDCA8", "IGF1R", "PIK3CA", "RANGAP1", "CASP7",
                    "BID", "HPN", "ENO1", "TP53", "CAND1", "UBE2V1", "UBE2V2", "UBE2N", "BIRC2", "TRAPPC2L", "TXNRD1", "ERBB2", "ERBB3")
surf_mark <- c("EPCAM", "FOLH1")



# Filter the DEG list to select the proper genes 
CSC_DEGs <- DEG_cluster_0v1_symbols_clean[DEG_cluster_0v1_symbols_clean$geneNames %in% CSC_list$CSC_SAB_target_list_H96, ]
EMT_DEGs <- DEG_cluster_0v1_symbols_clean[DEG_cluster_0v1_symbols_clean$geneNames %in% EMT_list$EMT_SAB_target_list_H96, ]
hub_DEGs <- DEG_cluster_0v1_symbols_clean[DEG_cluster_0v1_symbols_clean$geneNames %in% sign_hub_genes, ]
Krat_DEGs <- DEG_cluster_0v1_symbols_clean[DEG_cluster_0v1_symbols_clean$geneNames %in% Krat_panel, ]
Surf_mark_DEGs <- DEG_cluster_0v1_symbols_clean[DEG_cluster_0v1_symbols_clean$geneNames %in% surf_mark, ]



# Visualize the genes based on their DEG expression levels
ggplot(data = filter(CSC_DEGs, ExprChange != "Not significant"), aes(avg_log2FC, geneNames, fill = avg_log2FC, colour = avg_log2FC)) +
  geom_bar(stat = "identity", color = "black") +
  scale_x_continuous(limits = c(-5, 10), breaks = seq(-5, 10, 2.5)) +
  scale_fill_gradient(low = "#6EE2FFFF", high = "#F7C530FF", limits = c(-5, 10))+
  labs(fill = expression(""), x = expression("log"[2] * "FC"), y = expression("Gene symbols")) +
  ggtitle("Cancer stem cell marker expression in Cluster 0 vs Cluster 1") +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), text = element_text(size = 18))


ggplot(data = filter(EMT_DEGs, ExprChange != "Not significant"), aes(avg_log2FC, geneNames, fill = avg_log2FC, colour = avg_log2FC)) +
  geom_bar(stat = "identity", color = "black") +
  scale_x_continuous(limits = c(-5, 10), breaks = seq(-5, 10, 2.5)) +
  scale_fill_gradient(low = "#6EE2FFFF", high = "#F7C530FF", limits = c(-5, 10))+
  labs(fill = expression(""), x = expression("log"[2] * "FC"), y = expression("Gene symbols")) +
  ggtitle("EMT marker expression in Cluster 0 vs Cluster 1") +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), text = element_text(size = 18))


ggplot(data = filter(Krat_DEGs, ExprChange != "Not significant"), aes(avg_log2FC, geneNames, fill = avg_log2FC, colour = avg_log2FC)) +
  geom_bar(stat = "identity", color = "black") +
  scale_x_continuous(limits = c(-5, 10), breaks = seq(-5, 10, 2.5)) +
  scale_fill_gradient(low = "#6EE2FFFF", high = "#F7C530FF", limits = c(-5, 10))+
  labs(fill = expression(""), x = expression("log"[2] * "FC"), y = expression("Gene symbols")) +
  ggtitle("Krat. marker expression in Cluster 0 vs Cluster 1") +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), text = element_text(size = 18))


ggplot(data = filter(Surf_mark_DEGs, ExprChange != "Not significant"), aes(avg_log2FC, geneNames, fill = avg_log2FC, colour = avg_log2FC)) +
  geom_bar(stat = "identity", color = "black") +
  scale_x_continuous(limits = c(-5, 10), breaks = seq(-5, 10, 2.5)) +
  scale_fill_gradient(low = "#6EE2FFFF", high = "#F7C530FF", limits = c(-5, 10))+
  labs(fill = expression(""), x = expression("log"[2] * "FC"), y = expression("Gene symbols")) +
  ggtitle("EpCAM/PSMA expression in Cluster 0 vs Cluster 1") +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), text = element_text(size = 18))


hub_DEGs_p <- ggplot(data = filter(hub_DEGs, ExprChange != "Not significant"), aes(avg_log2FC, geneNames, fill = avg_log2FC, colour = avg_log2FC)) +
  geom_bar(stat = "identity", color = "black") +
  scale_x_continuous(limits = c(-5, 10), breaks = seq(-5, 10, 2.5)) +
  scale_fill_gradient(low = "#6EE2FFFF", high = "#F7C530FF", limits = c(-5, 10))+
  labs(fill = expression(""), x = expression("log"[2] * "FC"), y = expression("Gene symbols")) +
  #ggtitle("Hub-gene expression in Cluster 0 vs Cluster 1") +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), text = element_text(size = 24))
print(hub_DEGs_p)


ggsave(filename = paste0("Hub-gene_expression_in_Cluster_0_vs_Cluster_1", script_version, ".png"), hub_DEGs_p,
       device = "png", path = paste0(analysis_dir, "/", "GSEA_results/", "Plots/"),
       width = 2000, height = 2500, units = "px", dpi = 320)
rm(hub_DEGs_p)

#################################################################   End Section  ##################################################################

##############################################################################################################################################################################
################################################################   This is the end of the analysis for now  ##################################################################
##############################################################################################################################################################################

