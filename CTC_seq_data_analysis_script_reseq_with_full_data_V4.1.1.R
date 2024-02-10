#!user/bin/Rscript

#=================================#
#|                               |#
#|    scS results - analysis     |#
#|                               |#
#=================================#



#################################################################  Setting up, formatting and  ##################################################################
                                                                #  cleaning the relevant data  #


# Defining the required R packages for the full analysis
cran_packages <- c("tidyverse", "stringi", "BiocManager",
              "scales", "RCurl", "cowplot", "rebus", "ggsci",
              "progress", "metap", "doSNOW", "foreach", "scCustomize",
              "Matrix", "ggpubr", "R.utils", "devtools", "remotes", "RMTstat")

bioconductor_packages <- c("Seurat", "glmGamPoi", "multtest", "biomaRt", "AnnotationDbi",
              "EnsDb.Hsapiens.v86", "EnhancedVolcano", "graphite", "netgsa",
              "org.Hs.eg.db", "fgsea", "clusterProfiler", "SPIA", "PCAtools")

github_packages <- c("SeuratWrappers", "monocle3", "presto")
github_install_names <- c("satijalab/seurat-wrappers", "cole-trapnell-lab/monocle3", "immunogenomics/presto")                            


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
packages <- c("tidyverse", "stringi", "BiocManager", "Seurat", "biomaRt", "AnnotationDbi", "scales", "EnsDb.Hsapiens.v86",
                "RCurl", "cowplot", "rebus", "ggsci", "EnhancedVolcano", "progress", "metap", "doSNOW", "foreach", "scCustomize", "graphite",
                    "org.Hs.eg.db", "fgsea", "clusterProfiler", "ggpubr", "SeuratWrappers", "multtest", "PCAtools", "RMTstat")

lapply(packages, library, character.only = TRUE)


# Setting wd + path and listing the read output tables
setwd("C:/Work/DKFZ/Alpha_project_data_anal/Reseq")
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
    if (dir.exists(paths = paste0(script_version, "/", "Plots")) == TRUE) {
        message("The current script directory is already present. No new directory will be created.")
    } else {
        message("The current script vesrsion has no directorey yet.", "\n", "Creating one now using the path:", paste0(getwd(), script_version, "/", "Plots"))
        dir.create(path = paste0(script_version, "/", "Plots"), recursive = TRUE)
        dir.create(path = paste0(script_version, "/", "GSEA_results"), recursive = TRUE)
        dir.create(path = paste0(script_version, "/", "GSEA_results", "/", "Plots"), recursive = TRUE)
        
    }
    assign("script_version", script_version, envir = .GlobalEnv)
}
define_version_create_folder()


# Set the new wd to be the script version specific folder for ease of use (and ease of version control)
#setwd(paste0("C:/Work/DKFZ/Alpha_project_data_anal/Reseq", "/", script_version))


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
    
    #this function and sequence replaces the for (i in ...) and allows parallel processing. By default the function returns a list
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
write_csv(GeneCountTable, paste0(script_version, "/", "Unstranded_GeneCountTable_re-seq", script_version, ".csv"))


# From this point on we can work with the clean gene count table
GeneCountTable <- read.csv(file = paste0(script_version, "/", "Unstranded_GeneCountTable_re-seq", script_version, ".csv"))
head(GeneCountTable)


# Saving the geneIDs into a separate string for gene name conversion
geneIDs <- GeneCountTable$geneID
rownames(GeneCountTable) <- geneIDs
head(GeneCountTable)


# Load the unified gene list if available, otherwie build it from scratch
if (file.exists("Unified_gene_names_table.csv") && !exists("unified_gene_names_table", where = .GlobalEnv)) {
    unif_gene_names <- read.csv(file = "Unified_gene_names_table.csv")
    message("The unified gene names table was found and loaded as: ", deparse(substitute(unified_gene_names_table)))
} else {
    
    # Retrieving gene names using the ensembl IDs
    biomaRt::listEnsembl() #lists the available biomart keynames, here I need "gene"
    biomaRt::ensembl <- useEnsembl(biomart = "genes") #creates an ensembl object needed for the database connection
    datasets <- listDatasets(ensembl) #lists the organism based datasest from which we need to chose one for the db connection
    ensembl_connection <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl") #this connects the database and the dataset for the search query
    attributes <- listAttributes(ensembl_connection) #needed for the query, describes the desired info like gene names
    filters <- listFilters(ensembl_connection) # needed for the query, defined by upon which data we search for(?) like here the ensemblIds
    
    gene_names <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                        filters = "ensembl_gene_id",
                        values = geneIDs,
                        mart = ensembl_connection)
    
    
    
    # If you run into missing IDs you can check which ones are missing
    missingIDs <- geneIDs[geneIDs %in% gene_names$ensembl_gene_id == FALSE]
    missingNames <- gene_names[stringi::stri_isempty(gene_names$external_gene_name) == TRUE, 1]
    unif_missing <- c(missingIDs, missingNames)
    
    # An alternative method for gene id mapping if the direct query does not work properly
    # like shorter name list than id list...
    keytypes(EnsDb.Hsapiens.v86)
    columns(EnsDb.Hsapiens.v86)
    
    gene_names_2 <- mapIds(x = EnsDb.Hsapiens.v86,
                           keys = unif_missing,
                           column = "SYMBOL",
                           keytype = "GENEID")
    gene_names_2 <- data.frame(names(gene_names_2), gene_names_2)
    rownames(gene_names_2) <- NULL
    colnames(gene_names_2) <- colnames(gene_names)
    clean_gene_names <- gene_names[stringi::stri_isempty(gene_names$external_gene_name) == FALSE, ]
    unif_gene_names <- rbind(clean_gene_names, gene_names_2)
    unif_gene_names <- arrange(unif_gene_names, ensembl_gene_id)
    
    
    # Removing unnecessary variables
    rm(gene_names, gene_names_2, clean_gene_names, missingIDs, missingNames, unif_missing)
    
    
    # Save the unified gene names for later
    write_csv(x = unif_gene_names, file = "Unified_gene_names_table.csv")
    
}


# Re-ordering the GeneCountTable (GCT) by rows to match the re-ordered unif_gene_names df (by ensembl IDs)
# Removing the geneID column to allow ordering by sample (column) name and ordering the counts
#oGCT <- arrange(GeneCountTable, geneID)
#oGCT <- as.data.frame(oGCT)

#rownames(oGCT) <- unif_gene_names$ensembl_gene_id #optional, to track which row is which gene without messing with the ordering




## Loading the metadata.csv (only do this if you don't have the extended version yet)


# This if statement will check if the processes/extended metadata is present or not. If yes it will not run the following code chunk
if (file.exists("Re-seq_sample_extended_fixed_updated_metadata_with_markers_V2.csv")) {
    message("The processed metadata file is already present in the working directory.",
    "\n", "No further processing will be done.") 
} else if (file.exists("Re-seq_sample_extended_fixed_updated_metadata_with_markers_V2.csv")) {
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
read_metadata(filename = "Re-seq_sample_extended_fixed_updated_metadata_with_markers_V2.csv", sep = ";")
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
saveRDS(object = CTC_reseq.obj, file = paste0(getwd(), "/", script_version, "/", "Original_CTC_Seurat_object", script_version, ".rds"))


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
    CTC_reseq.obj <- readRDS(file = paste0(getwd(), "/", script_version, "/", "Original_CTC_Seurat_object", script_version, ".rds"))

}


# Load mitochondrial geneIDs and names (got it from ENSEMBL) for the QC
MT_genes <- read_csv("ENSEMBL_MT_genes.csv")
mito.genes <- rownames(CTC_reseq.obj)[rownames(CTC_reseq.obj) %in% MT_genes$`Gene stable ID`]
print(mito.genes)
rm(mito.genes)


# Alternatively one can load it from the scCustomize package
#ensembl_mito_id


# Add mitochondrial gene percentage and genes per nCount percentages (log10 transformed for better visibility) to our
# seurat.obj
CTC_reseq.obj[["mt.percent"]] <- PercentageFeatureSet(CTC_reseq.obj, features = MT_genes$`Gene stable ID`)
CTC_reseq.obj[["nFeature.per.nCount"]] <- log10(CTC_reseq.obj$nFeature_RNA_seq) / log10(CTC_reseq.obj$nCount_RNA_seq)


# At this point I will separate the metadata from the seurat object, so I won't screw up the object
# while trying out different things
CTC_meta.data <- CTC_reseq.obj@meta.data
CTC_meta.data_clean <- CTC_meta.data[!is.na(CTC_meta.data$mt.percent), ]
CTC_meta.data_clean <- mutate(CTC_meta.data_clean, Cell_names = rownames(CTC_meta.data_clean))
rm(CTC_meta.data)


# This function simply runs a block of code, plotting and saving the QC metrics one by one (the only purpose of this is conciseness)
generate_pre_qc_plots = function(make_plots = TRUE) {
    if (make_plots == TRUE) {
        message("Pre-QC plots were requested. Running the plotting code block.")
        
        
        # Visualize the cell numbers / resistance group and save it
        cn_p <- ggplot(CTC_meta.data_clean,
                       aes(x = RespGroup, fill = RespGroup)) +
            geom_bar(color = "black") +
            geom_text(aes(label = after_stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + #thx Bing
            scale_fill_tron() +
            labs(x = "Responder groups", y = "Cell numbers", fill = "Response groups", color = "Response groups") +
            ggtitle("Cell numbers after sequencing") +
            guides(fill = guide_legend(title = "Response groups")) +
            theme_classic() +
            theme(plot.title = element_text(hjust = 0.5))
        print(cn_p)
        
        ggsave(filename = paste0("Cell numbers", script_version, ".png"), plot = cn_p,
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        rm(cn_p)
        
        
        # Visualizing the mt.percent across cells using density plot and save it
        mt_pd <- ggplot(CTC_meta.data_clean,
                        aes(x = mt.percent, fill = RespGroup)) +
            geom_density(alpha = 0.2) +
            scale_fill_tron() +
            labs(x = "MT gene expression percentage", y = "MT gene percent distribution") +
            theme_classic() +
            ggtitle("Mitochondrial gene expression") +
            theme(plot.title = element_text(hjust = 0.5))
        print(mt_pd)
        
        ggsave(filename = paste0("Mitochondrial gene expression", script_version, ".png"), plot = mt_pd,
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        rm(mt_pd)
        
        
        # Visualizing the nCount numbers/treatment groups and save it
        # The scales package helps to remove scientific notation form the axises
        nCount_p <- ggplot(CTC_meta.data_clean,
                           aes(x = nCount_RNA_seq, fill = RespGroup)) +
            geom_density(alpha = 0.2) +
            scale_fill_tron() +
            labs(x = "nCounts", y = expression("log"[10]* " Count density"), fill = "Response group") +
            scale_x_log10(labels = scales::comma) + 
            theme_classic() +
            geom_vline(xintercept = 4000, color = "red") +
            ggtitle("mRNA count distribution") +
            theme(plot.title = element_text(hjust = 0.5))
        print(nCount_p)
        
        ggsave(filename = paste0("mRNA count distribution", script_version, ".png"), plot = nCount_p,
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        rm(nCount_p)
        
        
        # Visualizing the nFeature numbers/treatment groups and save it
        nFeature_p <- ggplot(CTC_meta.data_clean,
                             aes(x = nFeature_RNA_seq, fill = RespGroup)) +
            geom_density(alpha = 0.2) +
            scale_fill_tron() +
            labs(x = "nFeatures", y = expression("log"[10]* " Feature density"), fill = "Response group") +
            scale_x_log10(labels = scales::comma) + 
            theme_classic() +
            ggtitle("Gene count distribution") +
            geom_vline(xintercept = 1000, color = "red") +
            geom_vline(xintercept = 7000, color = "orange") +
            theme(plot.title = element_text(hjust = 0.5))
        print(nFeature_p)
        
        ggsave(filename = paste0("Gene count distribution", script_version, ".png"), plot = nFeature_p,
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        rm(nFeature_p)
        
        
        # Cell complexity plotting by checking the nFeature/nCount ratio (the higher the nFeature/nCount
        # the more complex the cells are) and save it
        nFpernC_p <-ggplot(CTC_meta.data_clean,
                           aes(x = nFeature.per.nCount, fill = RespGroup)) +
            geom_density(alpha = 0.2) +
            scale_fill_tron() +
            labs(x = "nFeatures/nCounts", y = expression("log" [10]* " density"), fill = "Response group") +
            theme_classic() +
            ggtitle("Cell complexity distribution (nFeature/nCount)") +
            theme(plot.title = element_text(hjust = 0.5))
        print(nFpernC_p)
        
        ggsave(filename = paste0("Cell complexity distribution", script_version, ".png"), plot = nFpernC_p,
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        rm(nFpernC_p)
        
        
        # nCount vs nGene (a high ratio is better) and save it
        nCvnF_p <- ggplot(data = CTC_meta.data_clean,
                          aes(x = nCount_RNA_seq, y = nFeature_RNA_seq, color = mt.percent)) +
            geom_point() + 
            scale_color_gradient(low = "grey90", high = "red") +
            scale_fill_tron() +
            labs(x = "mRNA counts", y = "Gene counts", color = "MT gene expression percentage") +
            #stat_smooth(method = lm) +
            scale_x_log10(labels = scales::comma) +
            scale_y_log10(labels = scales::comma) +
            geom_hline(yintercept = 1000, linetype = "dashed") +
            geom_hline(yintercept = 6500, linetype = "dashed") +
            geom_vline(xintercept = 4000, linetype = "dashed") +
            theme_classic() +
            #facet_wrap(facets = ~RespGroup)
            ggtitle("nCount vs nFeature - QC plot") +
            theme(plot.title = element_text(hjust = 0.5))
        print(nCvnF_p)
        
        ggsave(filename = paste0("nCount vs nFeature - QC plot", script_version, ".png"), plot = nCvnF_p,
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        rm(nCvnF_p)
        
        
        # End of run message
        message("The requested pre-QC plots were printed and saved into the: ", paste0(script_version, "/", "Plots/"), " directory.")
        
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


# I set the lower cutoff here to 1000 genes (analysis tutorial, 250 genes, other sources 200 genes
# Jonathan lower 500) and the upper cutoff to 7000 (Jonathan upper 5000) there seems to be no consensus
# regarding this metric. I ultimately would not draw an upper limit as the population with
# high gen/mRNA counts seem well integrated to the corresponding population. Regardless the advice is
# that one should look the 3 parameter (nCount, nGene, MTratio) together and set the 
# thresholds accordingly (that is why I love the UvGplot)



## Filtering and trimming


# Filter the trimmed Seurat object according to the above decided parameters
filt_CTC.obj <- subset(CTC_reseq.obj, 
                           subset = nCount_RNA_seq >= 4000 & nFeature_RNA_seq >= 1000 &
                                    nFeature_RNA_seq <= 6500)


# Save the filtered Seurat object for later use
saveRDS(object = filt_CTC.obj, file = paste0(getwd(), "/", script_version, "/", "Filtered_CTC_Seurat_object", script_version, ".rds"))


# Load the filtered Seurat object if not loaded yet
if (exists("filt_CTC.obj", where = .GlobalEnv)) {
    message("The filt_CTC.obj is already present in the .GlobalEnv")
} else {
    message("Loading the filt_CTC.obj")
    filt_CTC.obj <- readRDS(file = paste0(getwd(), "/", script_version, "/", "Filtered_CTC_Seurat_object", script_version, ".rds"))

}




## Another advised step is to remove genes with 0 counts across all cells


#There are two ways to do this:


# 1st

# For this first extract the counts from the seurat object using the GetAssayData function
#counts <- GetAssayData(filt_CTC.obj, slot = "counts")
#head(counts)

# Next, create a boolean matrix, describing if the count is bigger than 0 or not
#nonzero <- counts > 0
#head(nonzero)

# Sums all TRUE values and returns TRUE if more than 3 cells express this gene
# (15% of the nonresponder group and 5% of the responder group )
#keep_genes <- Matrix::rowSums(nonzero) >= 3
#summary(keep_genes)

# Only keeping those genes expressed in more than 3 cells
#filtered_counts <- counts[keep_genes, ]
#glimpse(filtered_counts)

# Re-create the seurat object containing the filtered gene
#filt_CTC.obj <- CreateSeuratObject(counts = filtered_counts,
#                                    project = "CTC_reseq",
#                                    assay = "RNA_seq",
#                                    meta.data = filt_CTC.obj@meta.data,
#                                    min.cells = 0)


# 2nd

# Creating a new seurat object from the filt_CTC.obj, using the min.cell argument with the
# desired number
#counts <- GetAssayData(filt_CTC.obj, slot = "counts")
#filt_CTC.obj <- CreateSeuratObject(counts = counts,
#                                    project = "CTC_reseq",
#                                    assay = "RNA_seq",
#                                    meta.data = filt_CTC.obj@meta.data,
#                                    min.cells = 3)

#rm(counts, CTC_meta.data_clean)




## Re-affirming the QC metrics on the trimmed, filtered data


# This function simply runs a block of code, plotting and saving the filtered data one by one (the only purpose of this is conciseness)
generate_filtered_data_plots = function(print_plots = TRUE, save_plots = TRUE) {
    # Defining the additional color palette
    cluster_colors <- pal_tron(palette = "legacy")(7)
    
    # Checking the new cell numbers
    cn_p <- ggplot(filt_CTC.obj@meta.data,
                   aes(x = RespGroup, fill = RespGroup)) +
        geom_bar(color = "black") +
        geom_text(aes(label = after_stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) +
        scale_fill_tron() +
        labs(x = "Response groups", y = "Cell numbers", fill = "Response groups") +
        ggtitle("Filtered cell numbers after sequencing - Response groups") +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5))
    
    # Checking the new cell numbers
    dcn_p <- ggplot(filt_CTC.obj@meta.data,
                   aes(x = detailedResp, fill = detailedResp)) +
        geom_bar(color = "black") +
        geom_text(aes(label = after_stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) +
        scale_fill_manual(values = cluster_colors[c(5, 1, 2, 3)]) +
        labs(x = "Response groups", y = "Cell numbers", fill = "Response groups") +
        ggtitle("Filtered cell numbers after sequencing - Detailed response") +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5))
    
    # Visualizing the mt.percent across cells using density plot and save it
    mt_pd <- ggplot(filt_CTC.obj@meta.data,
                    aes(x = mt.percent, fill = RespGroup)) +
        geom_density(alpha = 0.2) +
        scale_fill_tron() +
        theme_classic() +
        labs(x = "Cell percentage", y = "MT gene percent distribution", fill = "Response groups") +
        ggtitle("Filtered mitochondrial gene expression") +
        theme(plot.title = element_text(hjust = 0.5))
    
    # Visualizing the nCount numbers/treatment groups and save it
    nCount_p <- ggplot(filt_CTC.obj@meta.data,
                       aes(x = nCount_RNA_seq, fill = RespGroup)) +
        geom_density(alpha = 0.2) +
        scale_fill_tron() +
        labs(x = "nCounts", y = expression("log"[10]* " Count density"), fill = "Response group") +
        scale_x_log10(labels = scales::comma) + 
        theme_classic() +
        #geom_vline(xintercept = 2000, color = "red") +
        ggtitle("Filtered mRNA count distribution") +
        theme(plot.title = element_text(hjust = 0.5))
    
    # Visualizing the nFeature numbers/treatment groups and save it
    nFeature_p <- ggplot(filt_CTC.obj@meta.data,
                         aes(x = nFeature_RNA_seq, fill = RespGroup)) +
        geom_density(alpha = 0.2) +
        scale_fill_tron() +
        labs(x = "nFeatures", y = expression("log"[10]* " Feature density"), fill = "Response group") +
        scale_x_log10(labels = scales::comma) + 
        theme_classic() +
        ggtitle("Filtered gene count distribution") +
        #geom_vline(xintercept = 1000, color = "red") +
        #geom_vline(xintercept = 6500, color = "orange") +
        theme(plot.title = element_text(hjust = 0.5))
    
    # Cell complexity plotting by checking the nFeature/nCount ratio (the higher the nFeature/nCount
    # the more complex the cells are) and save it
    nFpernC_p <-ggplot(filt_CTC.obj@meta.data,
                       aes(x = nFeature.per.nCount, fill = RespGroup)) +
        geom_density(alpha = 0.2) +
        scale_fill_tron() +
        labs(x = "nFeatures/nCounts", y = expression("log"[10]* " density"), fill = "Response group") +
        theme_classic() +
        ggtitle("Cell complexity distribution (nFeture/nCount)") +
        theme(plot.title = element_text(hjust = 0.5))
    
    # nCount vs nGene (a high ratio is better) and save it
    nCvnF_p <- ggplot(data = filt_CTC.obj@meta.data,
                      aes(x = nCount_RNA_seq, y = nFeature_RNA_seq, color = mt.percent)) +
        geom_point() + 
        scale_color_gradient(low = "grey90", high = "red") +
        #stat_smooth(method = lm) +
        labs(x = "mRNA counts", y = "Gene counts", color = "MT gene expression percentage") +
        scale_x_log10(labels = scales::comma) +
        scale_y_log10(labels = scales::comma) +
        #geom_hline(yintercept = 1000, linetype = "dashed") +
        #geom_hline(yintercept = 6500, linetype = "dashed") +
        #geom_vline(xintercept = 4000, linetype = "dashed") +
        theme_classic() +
        #facet_wrap(facets = ~RespGroup)
        ggtitle("Filtered nCount vs nFeature - QC plot") +
        theme(plot.title = element_text(hjust = 0.5))
    
    
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
        
        ggsave(filename = paste0("Filtered cell numbers", script_version, ".png"), plot = cn_p,
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        
        ggsave(filename = paste0("Filtered cell numbers - detailed response", script_version, ".png"), plot = dcn_p,
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        
        ggsave(filename = paste0("Filtered mitochondrial gene expression", script_version, ".png"), plot = mt_pd,
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        
        ggsave(filename = paste0("Filtered mRNA count distribution", script_version, ".png"), plot = nCount_p,
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        
        ggsave(filename = paste0("Filtered gene count distribution", script_version, ".png"), plot = nFeature_p,
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        
        ggsave(filename = paste0("Filtered cell complexity distribution", script_version, ".png"), plot = nFpernC_p,
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        
        ggsave(filename = paste0("Filtered nCount vs nFeature - QC plot", script_version, ".png"), plot = nCvnF_p,
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        
        message("The plots were saved to the folder: ", paste0(script_version, "/", "Plots/"))
        
        
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
ccPhase_CTC.obj <- NormalizeData(filt_CTC.obj)


# This if else statement will check if the S and G2M cell cycle markers are present in the working directory and if yes are they loaded in the global Env.
# If present but not loaded it will load them, if not present and not loaded it will attempt to make them 
if ((file.exists("s_phase_IDs.csv") & file.exists("g2m_phase_IDs.csv")) && (!exists(x = "s.IDs", where = .GlobalEnv) & !exists(x = "g2m.IDs", where = .GlobalEnv))) {
    
    s.IDs <- read.csv(file = "s_phase_IDs.csv")
    g2m.IDs <- read.csv(file = "g2m_phase_IDs.csv")
    message("The cell cycle marker tables for S and G2M phases are loaded as the following variables: ", deparse(substitute(s.IDs)), " and " ,deparse(substitute(g2m.IDs)))
    
} else if (exists(x = "s.IDs", where = .GlobalEnv) & exists(x = "g2m.IDs", where = .GlobalEnv)) {
    
    message("The cell cycle marker tables for S and G2M phases are already loaded in the .GlobalEnv.")
    
    
} else {
    
    message("The cell cycle marker tables for S and G2M are not found in the working directory. Running the code block to manually create them:")
    
    # Load cell cycle markers
    # Cell cycle markers are available as part of the seurat package in gene name format
    s.genes <- cc.genes.updated.2019$s.genes
    g2m.genes <- cc.genes.updated.2019$g2m.genes
    
    
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
ccPhase_CTC.obj <- CellCycleScoring(ccPhase_CTC.obj,
                                    s.features = s.IDs$ensembl_gene_id,
                                    g2m.features = g2m.IDs$ensembl_gene_id,
                                    set.ident = FALSE)

# For the next step I need to reduce variation, for this I have to find the feature
# causing high variation (usually the very highly expressed genes) and scale the data
# seurat has built in functions for this
ccPhase_CTC.obj <- FindVariableFeatures(ccPhase_CTC.obj,
                                        selection.method = "vst",
                                        nfeatures = 3000,
                                        verbose = TRUE)
ccPhase_CTC.obj <- ScaleData(ccPhase_CTC.obj)


# Now I can do the PCA and plot the results
ccPhase_CTC.obj <- RunPCA(ccPhase_CTC.obj, approx = FALSE, verbose = FALSE)


# Testing the significance of PCs with the Jackstraw method and plotting the PCs
# NOTE: this is not possible on an objects transformed by SCTransform!
test_object <- JackStraw(object = ccPhase_CTC.obj, reduction = "pca", assay = "RNA_seq",
                        dims = 30)
ScoreJackStraw(test_object, dims = 1:30, do.plot = TRUE)
rm(test_object)


# Plotting and saving the result using the cell cycle phase overlay
CellCyc_split_p <- DimPlot(ccPhase_CTC.obj,
                           reduction = "pca",
                           group.by = "Phase",
                           split.by = "Phase")
CellCyc_split_p2 <- CellCyc_split_p +
                        scale_fill_tron() +
                        scale_color_tron() +
                        labs(x = "PC 1", y = "PC 2", color = "Cell cycle phase", fill = "Cell cycle phase") +
                        ggtitle("Cell cycle phases (split)")
print(CellCyc_split_p2)

ggsave(filename = paste0("Cell cycle phase clustering (split)", script_version, ".png"), CellCyc_split_p2,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 1500, height = 1000, units = "px")
rm(CellCyc_split_p, CellCyc_split_p2)


CellCyc_p <- DimPlot(ccPhase_CTC.obj, 
                     reduction = "pca",
                     group.by = "Phase")
CellCyc_p2 <- CellCyc_p +
                scale_fill_tron() +
                scale_color_tron() +
                labs(x = "PC 1", y = "PC 2", color = "Cell cycle phase", fill = "Cell cycle phase") +
                ggtitle("Cell cycle phases")
print(CellCyc_p2)

ggsave(filename = paste0("Cell cycle phase clustering", script_version, ".png"), CellCyc_p2,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 1500, height = 1000, units = "px")
rm(CellCyc_p, CellCyc_p2)


# Removing the ccPhase and other unnecessary objects
rm(ccPhase_CTC.obj, CTC_reseq.obj)


# Mitochondrial gene expression
# First we take a look at the mitotic gene expression percentage and check the quartile values
#summary(ccPhase_CTC.obj@meta.data$mt.percent)


# Then, we turn the mt.percent into a categorical factor vector, based on the quartile values
#ccPhase_tCTC.obj@meta.data$mtFr <- cut(ccPhase_tCTC.obj@meta.data$mt.percent,
#                                        breaks = c(-Inf, 7.8936, 14.2639, 26.8550, Inf),
#                                        labels = c("Low", "Medium", "Medium high", "High"))

# Run a PCA (might not be necessary...)
#ccPhase_tCTC.obj <- RunPCA(ccPhase_tCTC.obj, approx = FALSE, verbose = FALSE)


# Plot the result using the mt.fr overlay
#DimPlot(ccPhase_tCTC.obj,
#        reduction = "pca",
#        group.by = "mtFr",
#        split.by = "mtFr")


# Overlay the mt.percent using the FeaturePlot
#mt.per_p <- FeaturePlot(ccPhase_CTC.obj,
#                        features = "mt.percent")
#print(mt.per_p)
#rm(mt.per_p)

# In this particular case both the cell cycle and mitotic gene expression data scatters together, and there are no visible clusters
# based on the cell cycle phase or mitotic gene expression so they do not have to be regressed out

#################################################################   End Section  ##################################################################






#################################################################   Final normalization and regressing out  ##################################################################
                                                                #        sources of unwanted variation      #




## Normalization and clustering without integration


# As for the final more accurate method of normalization SCTransform, which is estimating the variance of the raw filtered data,
# and identifying the most variable genes. First i want to try and look at the data without integration.


# You might load the already prepared SCT transformed seurat object if available in the working directory
# NOTE: V3.9 is used as the standard as this script was forked for the Monocle3 analysis, and general refactoring
if (file.exists(paste0(getwd(), "/", "SCTransformed_fully_prepared_Seurat_object_V3.9.rds")) && !exists(x = "sct_CTC.obj", where = .GlobalEnv)) {
    sct_CTC.obj <- readRDS(file = paste0(getwd(), "/", "SCTransformed_fully_prepared_Seurat_object_V3.9.rds"))
    message("The standard SCTransformed Seurat object is present in the working directory and will be loaded as: ", deparse(substitute(sct_CTC.obj)))

} else {

    # If not present, load the filtered Surat object for further processing
    if (!exists(x = "filt_CTC.obj", where = .GlobalEnv)) {
        filt_CTC.obj <- readRDS(file = paste0(getwd(), "/", script_version, "/", "Filtered_CTC_Seurat_object", script_version, ".rds"))
        message("The filtered CTC Seurat object was loaded as: ", deparse(substitute(filt_CTC.obj)))
    }


    # Try both the newer and older SCTransform version to see which one works better
    sct_CTC.obj <- SCTransform(filt_CTC.obj, assay = "RNA_seq", method = "glmGamPoi")

    #sct_CTC.obj <- SCTransform(filt_CTC.obj, assay = "RNA_seq")
    # NOTE: The older version seems to work better, as it gives more defined clusters


    # Save the generated sCTransformed seurat object
    saveRDS(object = sct_CTC.obj, file = paste0(getwd(), "/", script_version, "/", "SCTransformed_CTC_Seurat_object", script_version, ".rds"))


    # Just as Jonathan, I will also repeat the cell cycle and mt scoring so I can overlay this info later onto the clustering
    sct_CTC.obj <- CellCycleScoring(sct_CTC.obj, 
                                    s.features = s.IDs$ensembl_gene_id, 
                                    g2m.features = g2m.IDs$ensembl_gene_id, 
                                    set.ident = FALSE)

    #summary(sct_CTC.obj@meta.data$mt.percent)
    #sct_CTC.obj@meta.data$mtFr <- cut(sct_CTC.obj@meta.data$mt.percent,
    #                                       breaks = c(-Inf, 9.1067, 15.4310, 26.4454, Inf),
    #                                       labels = c("Low", "Medium", "Medium high", "High"))


    # NOTE: there seems to be some genes missing here after SCTransform, so I'm unsure how reliable this is


    # An important thing to check is which PCs are responsible for the most variation. This is important
    # for the downstream steps. One way to check this is to plot the PCs with an ElbowPlot (plotting and saving the result)
    sct_CTC.obj <- RunPCA(sct_CTC.obj, assay = "SCT", approx = FALSE, verbose = FALSE, npcs = 50) #approx = FALSE to run with normal vsd method
    PC_weight_p <- ElbowPlot(sct_CTC.obj, ndims = 50, reduction = "pca")
    PC_weight_p2 <- PC_weight_p +
        labs(x = "PCs") +
        ggtitle("PC weights") +
        theme(plot.title = element_text(hjust = 0.5))
    print(PC_weight_p2)

    ggsave(filename = paste0("Principal component weights", script_version, ".png"), PC_weight_p2,
        device = "png", path = paste0(script_version, "/", "Plots/"),
        width = 1500, height = 1000, units = "px")
    rm(PC_weight_p, PC_weight_p2)


    # Another way to determine how many PCs one wants to use is by using the PCAtools package
    # NOTE: for the PCAtools to work one cannot use the seurat object itself, so some prep work is required


    # For this check to work we will need to re-do the PCA with the PCATools. In order to reproduce the same PCA
    # as we got with seurat we have to set the same seed as Seurat uses for the PCATools package and 
    # pull the variance data from the SCTransformed object (from the SCT slot)
    set.seed(42)
    sct_variance <- sct_CTC.obj@assays$SCT@scale.data[VariableFeatures(sct_CTC.obj), ]


    # Now that we set the seed and pulled the variance data we can run the pca
    # NOTE: in order to reproduce the Seurat PCA properly wee need to specifying the approximation method Seurat uses (IrlbaParam()) 
    # the number of PC's to retain (default is 50) and import the metadata for each cell (stored in slot @meta.data)
    QC_PCA <- PCAtools::pca(sct_variance, rank = 50, BSPARAM = BiocSingular::IrlbaParam(), metadata = sct_CTC.obj@meta.data)


    # Following the pca one can determine the number of PCs to keep using different methods:
    # 1st - finding an elbow point (similar to what one would visually do)
    PCAtools::findElbowPoint(QC_PCA$variance)

    # 2nd - using a parallel pca which will give the number of PCs to keep as $n
    QC_PCA_par <- PCAtools::parallelPCA(sct_variance, max.rank = 50, BSPARAM = BiocSingular::IrlbaParam(), threshold = 0.05)
    QC_PCA_par$n

    # 3rd - using the following methods
    PCAtools::chooseMarchenkoPastur(sct_variance, var.explained =  QC_PCA$sdev^2, noise = 4)
    PCAtools::chooseGavishDonoho(sct_variance, var.explained =  QC_PCA$sdev^2, noise = 4)

    # IMPORTANT NOTE: the various methods usually do not agree on how many PCs one should keep, so you probably should try a couple and
    # see if the UMAP and clustering makes sense based on what you set or not...
    
    
    # Unset seed
    set.seed(NULL)
    

    # Running dimensionality reduction (trying both UMAP and tSNE :) )
    sct_CTC.obj <- RunUMAP(sct_CTC.obj, dims = 1:20, reduction = "pca", verbose = FALSE)
    sct_CTC.obj <- RunTSNE(sct_CTC.obj, dims = 1:20, reduction = "pca", perplexity = 5)


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

    TSNE_p <- DimPlot(sct_CTC.obj, reduction = "tsne")
    print(TSNE_p)
    rm(TSNE_p)

    PCA_p <- PCAPlot(sct_CTC.obj, reduction = "pca")
    print(PCA_p)
    rm(PCA_p)


    # Re-run FindCluster with the preferred seeing
    sct_CTC.obj <- FindClusters(sct_CTC.obj, resolution = 0.6)


    # Saving the fully prepared CTC seurat object
    saveRDS(object = sct_CTC.obj, file = paste0(getwd(), "/", script_version, "/", "SCTransformed_fully_prepared_CTC_Seurat_object", script_version, ".rds"))

    
    # Closing message
    message("The block ran succesfully. \n",
            "The newly prepared SCT transformed Seurat object and the SCT transformed, fully prepared Seurat object were saved into the version folder ",
            script_version, " as follows: \n",
            "SCTransformed_CTC_Seurat_object", script_version, ".rds \n",
            "SCTransformed_fully_prepared_CTC_Seurat_object", script_version, ".rds")
    
}


# Cluster visualization (PCA, UMAP and t-SNE plots)
cluster_colors <- pal_tron(palette = "legacy")(7)
cluster_colors_2 <- pal_futurama(palette = "planetexpress")(12)
cluster_colors_3 <- pal_startrek(palette = "uniform")(7)
cluster_colors_4 <- pal_rickandmorty(palette = "schwifty")(12)


# This function will print and save the dimensionality reduction plots
plot_DR_cluster = function (print_plots = TRUE, save_plots = TRUE) {
    
    cluster_colors <- pal_tron(palette = "legacy")(7)
    cluster_colors_2 <- pal_futurama(palette = "planetexpress")(12)
    cluster_colors_3 <- pal_startrek(palette = "uniform")(7)
    cluster_colors_4 <- pal_rickandmorty(palette = "schwifty")(12)
    
    # Create the PCA plot
    PCA_p <- DimPlot(sct_CTC.obj, reduction = "pca", group.by = "seurat_clusters")
    PCA_p2 <- PCA_p +
        scale_fill_manual(values = cluster_colors[c(6, 7)]) +
        scale_color_manual(values = cluster_colors[c(6, 7)]) +
        labs(x = "PC 1", y = "PC 2", color = "Clusters", fill = "Clusters") +
        ggtitle("PCA plot")
    
    
    # Create the UMAP plot
    UMAP_p <- DimPlot(sct_CTC.obj, reduction = "umap", group.by = "seurat_clusters")
    UMAP_p2 <- UMAP_p +
        scale_fill_manual(values = cluster_colors[c(6, 7)]) +
        scale_color_manual(values = cluster_colors[c(6, 7)]) +
        labs(x = "UMAP 1", y = "UMAP 2", color = "Clusters", fill = "Clusters") +
        ggtitle("UMAP plot")
    
    
    # Create the TSNE plot
    TSNE_p <- DimPlot(sct_CTC.obj, reduction = "tsne", group.by = "seurat_clusters")
    TSNE_p2 <- TSNE_p +
        scale_fill_manual(values = cluster_colors[c(6, 7)]) +
        scale_color_manual(values = cluster_colors[c(6, 7)]) +
        labs(x = "tSNE 1", y = "tSNE 2", color = "Clusters", fill = "Clusters") +
        ggtitle("tSNE plot")
    
    
    
    # This if statement controls if the plost will be saved
    if (save_plots == TRUE) {
        
        ggsave(filename = paste0("PCA_clustering_gb", script_version, ".png"), PCA_p2,
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        
        ggsave(filename = paste0("UMAP_clustering_gb", script_version, ".png"), UMAP_p2,
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        
        ggsave(filename = paste0("TSNE_clustering_gb", script_version, ".png"), TSNE_p2,
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        
        message("The plots were saved to the folder: ", paste0(script_version, "/", "Plots/"))
        
    } else {
        
        message("WARNING: the plots are not saved as a save was not requested")
        
    }
    
    
    # This if statement controls if the plots are printed
    if (print_plots == TRUE){
        
        print(PCA_p2)
        print(UMAP_p2)
        print(TSNE_p2)
        
        message("The plots were printed as requested.")
        
    } else {
        
        message("The plots were not printed as print was not requested.")
        
    }
    
    
    
}
plot_DR_cluster(print_plots = TRUE, save_plots = TRUE)


# Removing unnecessary objects
rm(g2m.IDs, s.IDs, MT_genes, QC_PCA, QC_PCA_par, filt_CTC.obj, sct_variance)

#################################################################   End Section  ##################################################################






#################################################################   Plotting different metrics over the UMAP  ##################################################################
                                                                #         clusters to if they correlate       #


# Defining a color for gradient plots
grad_col <- pal_locuszoom("default")(7)


# This function simply runs a block of code, plotting and saving additional data overplayed on the UMAP (the only purpose of this is conciseness)
plot_UMAP_ovelrays = function (print_plots = TRUE, save_plots = TRUE) {
    
    # Define the used color plalettes
    grad_col <- pal_locuszoom("default")(7)
    cluster_colors <- pal_tron(palette = "legacy")(7)
    cluster_colors_2 <- pal_futurama(palette = "planetexpress")(12)
    cluster_colors_3 <- pal_startrek(palette = "uniform")(7)
    cluster_colors_4 <- pal_rickandmorty(palette = "schwifty")(12)
    
    # Overlaying nFeature
    UMAP_nF_p <- FeaturePlot(sct_CTC.obj,
                             reduction = "umap",
                             features = "nFeature_SCT")
    UMAP_nF_p2 <- UMAP_nF_p +
        #scale_color_gradient(low = "lightgrey", high = cluster_colors[7]) +
        scale_color_gradient(low = "lightgrey", high = grad_col[5]) +
        labs(x = "UMAP 1", y = "UMAP 2", color = "Clusters", fill = "Clusters") +
        ggtitle("UMAP plot - nFeatures")
    
    # Overlaying nCount
    UMAP_nC_p <- FeaturePlot(sct_CTC.obj,
                             reduction = "umap",
                             features = "nCount_SCT")
    UMAP_nC_p2 <- UMAP_nC_p +
        scale_color_gradient(low = "lightgrey", high = grad_col[5]) +
        labs(x = "UMAP 1", y = "UMAP 2", color = "Clusters", fill = "Clusters") +
        ggtitle("UMAP plot - nCounts")
    
    # Overlay mt.percent
    UMAP_mt.percent <- FeaturePlot(sct_CTC.obj,
                                   reduction = "umap",
                                   features = "mt.percent")
    UMAP_mt.percent2 <- UMAP_mt.percent +
        scale_color_gradient(low = "lightgrey", high = grad_col[5]) +
        labs(x = "UMAP 1", y = "UMAP 2", color = "Clusters", fill = "Clusters") +
        ggtitle("UMAP plot - MT gene expression")
    
    # Overlay cell cycle
    UMAP_CellCyc <- DimPlot(sct_CTC.obj,
                            reduction = "umap",
                            group.by = "Phase")
    UMAP_CellCyc2 <- UMAP_CellCyc +
        scale_fill_manual(values = cluster_colors_3[c(4, 5, 7)]) +
        scale_color_manual(values = cluster_colors_3[c(4, 5, 7)]) +
        labs(x = "UMAP 1", y = "UMAP 2", color = "Cell cycle phase", fill = "Cell cycle phase") +
        ggtitle("UMAP plot - Cell cycle phases")
    
    # Overlay treatment cycle
    UMAP_TreatCyc <- DimPlot(sct_CTC.obj,
                             reduction = "umap",
                             group.by = "RealCycle")
    UMAP_TreatCyc2 <- UMAP_TreatCyc +
        scale_fill_manual(values = cluster_colors_3[c(4, 5, 7)]) +
        scale_color_manual(values = cluster_colors_3[c(4, 5, 7)]) +
        labs(x = "UMAP 1", y = "UMAP 2", color = "Treatment cycles", fill = "Treatment cycles") +
        ggtitle("UMAP plot - Treatment cycles")
    
    # Overlay treatment type
    UMAP_TreatType <- DimPlot(sct_CTC.obj,
                              reduction = "umap",
                              group.by = "Treatment")
    UMAP_TreatType2 <- UMAP_TreatType +
        scale_fill_manual(values = cluster_colors_3[c(4, 5, 7)]) +
        scale_color_manual(values = cluster_colors_3[c(4, 5, 7)]) +
        labs(x = "UMAP 1", y = "UMAP 2", color = "Treatment", fill = "Treatment") +
        ggtitle("UMAP plot - Treatment types")
    
    # Overlay RespGroup
    UMAP_RespG <- DimPlot(sct_CTC.obj,
                          reduction = "umap",
                          group.by = "RespGroup")
    UMAP_RespG2 <- UMAP_RespG +
        scale_fill_manual(values = cluster_colors[c(1, 2)]) +
        scale_color_manual(values = cluster_colors[c(1, 2)]) +
        labs(x = "UMAP 1", y = "UMAP 2", color = "Response group", fill = "Response group") +
        ggtitle("UMAP plot - Response groups")
    
    # Overlay Detailed response
    UMAP_DetRespG <- DimPlot(sct_CTC.obj,
                             reduction = "umap",
                             group.by = "detailedResp")
    UMAP_DetRespG2 <- UMAP_DetRespG +
        scale_fill_manual(values = cluster_colors[c(5, 1, 2, 3)]) +
        scale_color_manual(values = cluster_colors[c(5, 1, 2, 3)]) +
        labs(x = "UMAP 1", y = "UMAP 2", color = "Response profile", fill = "Response profile") +
        ggtitle("UMAP plot - Detailed treatment response")
    
    # Overlay Luthetium-177 doses
    UMAP_Lu <- DimPlot(sct_CTC.obj,
                       reduction = "umap",
                       group.by = "Lu177_GBq")
    UMAP_Lu2 <- UMAP_Lu +
        scale_fill_manual(values = cluster_colors_4[c(4, 5, 6, 7, 8, 9)]) +
        scale_color_manual(values = cluster_colors_4[c(4, 5, 6, 7, 8, 9)]) +
        labs(x = "UMAP 1", y = "UMAP 2", color = "Lu-177 doses (GBq)", fill = "Lutetium-177 doses (GBq)") +
        ggtitle("UMAP plot - Lu-177 treatment doses")
    
    # Overlay Actinium-225 doses
    UMAP_Ac <- DimPlot(sct_CTC.obj,
                       reduction = "umap",
                       group.by = "Ac225_MBq")
    UMAP_Ac2 <- UMAP_Ac +
        scale_fill_manual(values = cluster_colors_4[c(4, 5, 6, 7)]) +
        scale_color_manual(values = cluster_colors_4[c(4, 5, 6, 7)]) +
        labs(x = "UMAP 1", y = "UMAP 2", color = "Ac-225 doses (MBq)", fill = "Ac-225 doses (MBq)") +
        ggtitle("UMAP plot - Ac-225 treatment doses")
    
    # Overlay marker status
    UMAP_MarkerStat <- DimPlot(sct_CTC.obj,
                               reduction = "umap",
                               group.by = "Marker.status")
    UMAP_MarkerStat2 <- UMAP_MarkerStat +
        scale_fill_manual(values = cluster_colors_3[c(4, 5, 7)]) +
        scale_color_manual(values = cluster_colors_3[c(4, 5, 7)]) +
        labs(x = "UMAP 1", y = "UMAP 2", color = "Marker status", fill = "Marker status") +
        ggtitle("UMAP plot - Surface marker distribution")
    
    
    
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
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        
        ggsave(filename = paste0("UMAP_nCount", script_version, ".png"), UMAP_nC_p2,
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        
        ggsave(filename = paste0("UMAP_mt.percent", script_version, ".png"), UMAP_mt.percent2,
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        
        ggsave(filename = paste0("UMAP_Cell_cycle", script_version, ".png"), UMAP_CellCyc2,
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        
        ggsave(filename = paste0("UMAP_Treatment_cycles", script_version, ".png"), UMAP_TreatCyc2,
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        
        ggsave(filename = paste0("UMAP_Treatment_type", script_version, ".png"), UMAP_TreatType2,
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        
        ggsave(filename = paste0("UMAP_Response_group", script_version, ".png"), UMAP_RespG2,
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        
        ggsave(filename = paste0("UMAP_Detailed_treatment_response", script_version, ".png"), UMAP_DetRespG2,
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        
        ggsave(filename = paste0("UMAP_Lu-177_treatment_doses", script_version, ".png"), UMAP_Lu2,
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        
        ggsave(filename = paste0("UMAP_Ac-225_treatment_doses", script_version, ".png"), UMAP_Ac2,
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        
        ggsave(filename = paste0("UMAP_Marker_Status", script_version, ".png"), UMAP_MarkerStat2,
               device = "png", path = paste0(script_version, "/", "Plots/"),
               width = 1500, height = 1000, units = "px")
        rm(UMAP_MarkerStat, UMAP_MarkerStat2)
        
        
        message("The plots were saved to the folder: ", paste0(script_version, "/", "Plots/"))
        
    } else {
        
        message("WARNING: the plots are not saved as a save was not requested")
        
    }
    
}
plot_UMAP_ovelrays(print_plots = TRUE, save_plots = TRUE)




## Visualizing the overlay data in a more understandable manner, using bar charts


# This function will calculate cluster composition stats and plots and saves them as requested (the only purpose of this is conciseness)
calculate_and_plot_cluster_composition = function (make_plots = TRUE, return_statistics = TRUE, print_statistics = TRUE) {
        
        # Preparing the main dataframe for the visualization
        #head(sct_CTC.obj@meta.data)
        
        sct_CTC_meta_df <- sct_CTC.obj@meta.data
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
                                                Total_cell_numbers = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "EpCAM+")),
                                                                       nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "EpCAM+PSMA+")),
                                                                       nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "PSMA+"))),
                                                Cell_numbers_c0 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "EpCAM+", seurat_clusters == 0)),
                                                                    nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "EpCAM+PSMA+", seurat_clusters == 0)),
                                                                    nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "PSMA+", seurat_clusters == 0))),
                                                Cell_numbers_c1 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "EpCAM+", seurat_clusters == 1)),
                                                                    nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "EpCAM+PSMA+", seurat_clusters == 1)),
                                                                    nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "PSMA+", seurat_clusters == 1))),
                                                Cell_percent_c0 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "EpCAM+", seurat_clusters == 0)) /
                                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "EpCAM+")) * 100,
                                                                    nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "EpCAM+PSMA+", seurat_clusters == 0)) /
                                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "EpCAM+PSMA+")) * 100,
                                                                    nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "PSMA+", seurat_clusters == 0)) /
                                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "PSMA+")) * 100),
                                                Cell_percent_c1 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "EpCAM+", seurat_clusters == 1)) /
                                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "EpCAM+")) * 100,
                                                                    nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "EpCAM+PSMA+", seurat_clusters == 1)) /
                                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "EpCAM+PSMA+")) * 100,
                                                                    nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "PSMA+", seurat_clusters == 1)) /
                                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "PSMA+")) * 100))
        
        
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
                                                      EpCAM_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Marker.status == "EpCAM+")),
                                                                  nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Marker.status == "EpCAM+"))),
                                                      EpCAM_PSMA_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Marker.status == "EpCAM+PSMA+")),
                                                                       nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Marker.status == "EpCAM+PSMA+"))),
                                                      PSMA_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Marker.status == "PSMA+")),
                                                                 nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Marker.status == "PSMA+"))),
                                                      EpCAM_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Marker.status == "EpCAM+")) /
                                                                      nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                                  nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Marker.status == "EpCAM+")) /
                                                                      nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100),
                                                      EpCAM_PSMA_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Marker.status == "EpCAM+PSMA+")) /
                                                                           nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                                       nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Marker.status == "EpCAM+PSMA+")) /
                                                                           nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100),
                                                      PSMA_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Marker.status == "PSMA+")) /
                                                                     nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                                 nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Marker.status == "PSMA+")) /
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
                Detailed_response_summary_df,
                Surface_marker_summary_df
            )
            
            names(variable_side_cluster_stats) <- c("Treatment_cycles_summary_df",
                                                    "Treatment_summary_df",
                                                    "Cell_cycle_summary_df",
                                                    "Response_summary_df",
                                                    "Detailed_response_summary_df",
                                                    "Surface_marker_summary_df")
            
            
            variable_side_cluster_stats_plotting_dfs <- list(
                Treatment_cycles_summary_df_long,
                Treatment_summary_df_long,
                Cell_cycle_summary_df_long,
                Response_summary_df_long,
                Detailed_response_summary_df_long,
                Surface_marker_summary_df_long
            )
            
            names(variable_side_cluster_stats) <- c("Treatment_cycles_summary_df_long",
                                                    "Treatment_summary_df_long",
                                                    "Cell_cycle_summary_df_long",
                                                    "Response_summary_df_long",
                                                    "Detailed_response_summary_df_long",
                                                    "Surface_marker_summary_df_long")
            
            
            cluster_side_cluster_stats <- list(
                Treatment_cycles_summary_df_clust,
                Treatment_summary_df_clust,
                Cell_cycle_summary_df_clust,
                Response_summary_df_clust,
                Detailed_response_summary_df_clust,
                Surface_marker_summary_df_clust
            )
            
            names(variable_side_cluster_stats) <- c("Treatment_cycles_summary_df_clust",
                                                    "Treatment_summary_df_clust",
                                                    "Cell_cycle_summary_df_clust",
                                                    "Response_summary_df_clust",
                                                    "Detailed_response_summary_df_clust",
                                                    "Surface_marker_summary_df_clust")
            
            
            cluster_side_cluster_stats_plotting_dfs <- list(
                Treatment_cycles_summary_df_clust_long,
                Treatment_summary_df_clust_long,
                Cell_cycle_summary_df_clust_long,
                Response_summary_df_clust_long,
                Detailed_response_summary_df_clust_long,
                Surface_marker_summary_df_clust_long
            )
            
            names(variable_side_cluster_stats) <- c("Treatment_cycles_summary_df_clust_long",
                                                    "Treatment_summary_df_clust_long",
                                                    "Cell_cycle_summary_df_clust_long",
                                                    "Response_summary_df_clust_long",
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
            
            
            #print(Treatment_cycles_summary_df)
            #print(Treatment_cycles_summary_df_long)
            #print(Treatment_summary_df)
            #print(Treatment_summary_df_long)
            #print(Cell_cycle_summary_df)
            #print(Response_summary_df)
            #print(Response_summary_df_long)
            #print(Surface_marker_summary_df)
            #print(Surface_marker_summary_df_long)
            #print(Treatment_cycles_summary_df_clust)
            #print(Treatment_cycles_summary_df_clust_long)
            #print(Treatment_summary_df_clust)
            #print(Treatment_summary_df_clust_long)
            #print(Cell_cycle_summary_df_clust)
            #print(Cell_cycle_summary_df_clust_long)
            #print(Response_summary_df_clust)
            #print(Response_summary_df_clust_long)
            #print(Surface_marker_summary_df_clust)
            #print(Surface_marker_summary_df_clust_long)
            
            
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
                geom_text(aes(label = as.integer(Percentages)), vjust = -0.5, position = position_dodge(width = 0.9)) +
                scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
                #scale_x_discrete(labels = c("Double positive","EpCAM positive","PSMA positive")) +
                scale_fill_manual(values = cluster_colors[c(6, 7)], label = c("Cluster 0", "Cluster 1")) +
                labs(fill = "CTC cluster", x = expression("Treatment cycle"), y = expression("Cell number percentages")) +
                ggtitle("CTC clusters - cell distribution based on treatment cycle") +
                #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
                theme_classic() +
                theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
                      strip.text = element_blank())
            print(treatment_cyc_bar_p)
            
            ggsave(filename = paste0("Treatment cycles in clusters", script_version, ".png"), treatment_cyc_bar_p,
                   device = "png", path = paste0(script_version, "/", "Plots/"),
                   width = 1500, height = 1000, units = "px")
            rm(treatment_cyc_bar_p)
            
            
            # Plotting the treatment type distribution
            treatment_type_bar_p <- ggplot(Treatment_summary_df_long,
                                           aes(x = factor(Treatment_type, levels = c("NT", "AL", "L")), y = Percentages, fill = Clusters)) +
                geom_col(color = "black", position = "dodge") +
                geom_text(aes(label = as.integer(Percentages)), vjust = -0.5, position = position_dodge(width = 0.9)) +
                scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
                #scale_x_discrete(labels = c("Double positive","EpCAM positive","PSMA positive")) +
                scale_fill_manual(values = cluster_colors[c(6, 7)], label = c("Cluster 0", "Cluster 1")) +
                labs(fill = "CTC cluster", x = expression("Treatment type"), y = expression("Cell number percentages")) +
                ggtitle("CTC clusters - cell distribution based on treatment types") +
                #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
                theme_classic() +
                theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
                      strip.text = element_blank())
            print(treatment_type_bar_p)
            
            ggsave(filename = paste0("Treatment types over clusters", script_version, ".png"), treatment_type_bar_p,
                   device = "png", path = paste0(script_version, "/", "Plots/"),
                   width = 1500, height = 1000, units = "px")
            rm(treatment_type_bar_p)
            
            
            # Plotting the cell cycle distribution
            cell_cyc_bar_p <- ggplot(Cell_cycle_summary_df_long,
                                     aes(x = factor(CC_phase, levels = c("G1", "S", "G2M")), y = Percentages, fill = Clusters)) +
                geom_col(color = "black", position = "dodge") +
                geom_text(aes(label = as.integer(Percentages)), vjust = -0.5, position = position_dodge(width = 0.9)) +
                scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
                #scale_x_discrete(labels = c("Double positive","EpCAM positive","PSMA positive")) +
                scale_fill_manual(values = cluster_colors[c(6, 7)], label = c("Cluster 0", "Cluster 1")) +
                labs(fill = "CTC cluster", x = expression("Cell cycle phase"), y = expression("Cell number percentages")) +
                ggtitle("CTC clusters - cell distribution based on cell cycle phase") +
                #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
                theme_classic() +
                theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
                      strip.text = element_blank())
            print(cell_cyc_bar_p)
            
            ggsave(filename = paste0("Cell cycle over clusters", script_version, ".png"), cell_cyc_bar_p,
                   device = "png", path = paste0(script_version, "/", "Plots/"),
                   width = 1500, height = 1000, units = "px")
            rm(cell_cyc_bar_p)
           
            
            # Plotting the treatment response distribution from response view
            response_group_bar_p <- ggplot(Response_summary_df_long,
                                           aes(x = factor(response_group, levels = c("Responder", "Nonresponder")), y = Percentages, fill = Clusters)) +
                geom_col(color = "black", position = "dodge") +
                geom_text(aes(label = as.integer(Percentages)), vjust = -0.5, position = position_dodge(width = 0.9)) +
                scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
                #scale_x_discrete(labels = c("Double positive","EpCAM positive","PSMA positive")) +
                scale_fill_manual(values = cluster_colors[c(6, 7)], label = c("Cluster 0", "Cluster 1")) +
                labs(fill = "CTC cluster", x = expression("Response group"), y = expression("Cell number percentages")) +
                ggtitle("CTC clusters - cell distribution based on treatment response") +
                #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
                theme_classic() +
                theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
                      strip.text = element_blank())
            print(response_group_bar_p)
            
            ggsave(filename = paste0("Treatment response over clusters", script_version, ".png"), response_group_bar_p,
                   device = "png", path = paste0(script_version, "/", "Plots/"),
                   width = 1500, height = 1000, units = "px")
            rm(response_group_bar_p)
            
            
            # Plotting the surface marker distribution
            surface_marker_bar_p <- ggplot(Surface_marker_summary_df_long,
                                           aes(x = factor(marker_status, levels = c("EpCAM+", "EpCAM+PSMA+", "PSMA+")), y = Percentages, fill = Clusters)) +
                geom_col(color = "black", position = "dodge") +
                geom_text(aes(label = as.integer(Percentages)), vjust = -0.5, position = position_dodge(width = 0.9)) +
                scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
                #scale_x_discrete(labels = c("Double positive","EpCAM positive","PSMA positive")) +
                scale_fill_manual(values = cluster_colors[c(6, 7)], label = c("Cluster 0", "Cluster 1")) +
                labs(fill = "CTC cluster", x = expression("Surface marker status"), y = expression("Cell number percentages")) +
                ggtitle("CTC clusters - cell distribution based on the surface marker status") +
                #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
                theme_classic() +
                theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
                      strip.text = element_blank())
            print(surface_marker_bar_p)
            
            ggsave(filename = paste0("Surface marker over clusters", script_version, ".png"), surface_marker_bar_p,
                   device = "png", path = paste0(script_version, "/", "Plots/"),
                   width = 1500, height = 1000, units = "px")
            rm(surface_marker_bar_p)
            
            
            # Plotting the treatment cycle distribution from the cluster view
            treatment_cyc_bar_p <- ggplot(Treatment_cycles_summary_df_clust_long,
                                          aes(x = Clusters, y = Percentages, fill = as.factor(Cycles))) +
                geom_col(color = "black", position = "dodge") +
                geom_text(aes(label = as.integer(Percentages)), vjust = -0.5, position = position_dodge(width = 0.9)) +
                scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
                scale_x_discrete(labels = c("Cluster 0", "Cluster 1", "Cluster 2")) +
                scale_fill_manual(values = cluster_colors[c(4, 6, 7)], label = c("Cycle 0", "Cycle 1", "Cycle 2")) +
                labs(fill = "Treatment cycles", x = expression("CTC clusters"), y = expression("Cell number percentages")) +
                ggtitle("CTC clusters - cell distribution based on treatment cycle") +
                #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
                theme_classic() +
                theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
                      strip.text = element_blank())
            print(treatment_cyc_bar_p)
            
            ggsave(filename = paste0("Cluster_composition-Treatment_cycles", script_version, ".png"), treatment_cyc_bar_p,
                   device = "png", path = paste0(script_version, "/", "Plots/"),
                   width = 1500, height = 1000, units = "px")
            rm(treatment_cyc_bar_p)
            
            
            # Plotting the treatment type distribution from the cluster view
            treatment_type_bar_p <- ggplot(Treatment_summary_df_clust_long,
                                           aes(x = Clusters, y = Percentages, fill = factor(Treatment_type, levels = c("NT_p", "AL_p", "L_p")))) +
                geom_col(color = "black", position = "dodge") +
                geom_text(aes(label = as.integer(Percentages)), vjust = -0.5, position = position_dodge(width = 0.9)) +
                scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
                scale_x_discrete(labels = c("Cluster 0", "Cluster 1", "Cluster 2")) +
                scale_fill_manual(values = cluster_colors[c(4, 6, 7)], label = c("NT", "AL", "L")) +
                labs(fill = "Treatment types", x = expression("CTC clusters"), y = expression("Cell number percentages")) +
                ggtitle("CTC clusters - cell distribution based on treatment types") +
                #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
                theme_classic() +
                theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
                      strip.text = element_blank())
            print(treatment_type_bar_p)
            
            ggsave(filename = paste0("Cluster_composition-Treatment_types", script_version, ".png"), treatment_type_bar_p,
                   device = "png", path = paste0(script_version, "/", "Plots/"),
                   width = 1500, height = 1000, units = "px")
            rm(treatment_type_bar_p)
            
            
            # Plotting the cell cycle distribution from the cluster view
            cell_cyc_bar_p <- ggplot(Cell_cycle_summary_df_clust_long,
                                     aes(x = Clusters, y = Percentages, fill = factor(CC_phase, levels = c("G1_p", "S_p", "G2M_p")))) +
                geom_col(color = "black", position = "dodge") +
                geom_text(aes(label = as.integer(Percentages)), vjust = -0.5, position = position_dodge(width = 0.9)) +
                scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
                scale_x_discrete(labels = c("Cluster 0", "Cluster 1", "Cluster 2")) +
                scale_fill_manual(values = cluster_colors[c(4, 6, 7)], label = c("G1", "S", "G2/M")) +
                labs(fill = "Cell cycle phase", x = expression("CTC cluster"), y = expression("Cell number percentages")) +
                ggtitle("CTC clusters - cell distribution based on cell cycle phase") +
                #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
                theme_classic() +
                theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
                      strip.text = element_blank())
            print(cell_cyc_bar_p)
            
            ggsave(filename = paste0("Cluster_composition-Cell_cycle_phases", script_version, ".png"), cell_cyc_bar_p,
                   device = "png", path = paste0(script_version, "/", "Plots/"),
                   width = 1500, height = 1000, units = "px")
            rm(cell_cyc_bar_p)
            
            
            # Plotting the treatment response distribution from the cluster view
            response_group_bar_p <- ggplot(Response_summary_df_clust_long,
                                           aes(x = Clusters, y = Percentages, fill = factor(Response_group, levels = c("Resp_p", "Nonresp_p")))) +
                geom_col(color = "black", position = "dodge") +
                geom_text(aes(label = as.integer(Percentages)), vjust = -0.5, position = position_dodge(width = 0.9)) +
                scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
                scale_x_discrete(labels = c("Cluster 0", "Cluster 1")) +
                scale_fill_manual(values = cluster_colors[c(2, 1)], label = c("Responder", "Nonresponder")) +
                labs(fill = "Response group", x = expression("CTC cluster"), y = expression("Cell number percentages")) +
                ggtitle("CTC clusters - cell distribution based on treatment response") +
                #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
                theme_classic() +
                theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
                      strip.text = element_blank())
            print(response_group_bar_p)
            
            ggsave(filename = paste0("Cluster_composition-Response_groups", script_version, ".png"), response_group_bar_p,
                   device = "png", path = paste0(script_version, "/", "Plots/"),
                   width = 1500, height = 1000, units = "px")
            rm(response_group_bar_p)
            
            
            # Plotting the detailed treatment response distribution from the cluster view
            detailed_response_group_bar_p <- ggplot(Detailed_response_summary_df_clust_long,
                                           aes(x = Clusters, y = Percentages, fill = factor(Response_group, levels = c("NT_p", "PD_p", "PR_p", "SD_p")))) +
                geom_col(color = "black", position = "dodge") +
                geom_text(aes(label = as.integer(Percentages)), vjust = -0.5, position = position_dodge(width = 0.9)) +
                scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
                scale_x_discrete(labels = c("Cluster 0", "Cluster 1")) +
                scale_fill_manual(values = cluster_colors[c(5, 1, 2, 3)], label = c("NT", "PD", "PR", "SD")) +
                labs(fill = "Detailed treatment response", x = expression("CTC cluster"), y = expression("Cell number percentages")) +
                ggtitle("CTC clusters - cell distribution based on detailed treatment response") +
                #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
                theme_classic() +
                theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
                      strip.text = element_blank())
            print(detailed_response_group_bar_p)
            
            ggsave(filename = paste0("Cluster_composition-Detailed_treatment_response", script_version, ".png"), detailed_response_group_bar_p,
                   device = "png", path = paste0(script_version, "/", "Plots/"),
                   width = 1500, height = 1000, units = "px")
            rm(response_group_bar_p)
            
            
            # Plotting the surface marker distribution from the cluster view
            surface_m_bar_p <- ggplot(Surface_marker_summary_df_clust_long,
                                      aes(x = Clusters, y = Percentages, fill = as.factor(Markers))) +
                geom_col(color = "black", position = "dodge") +
                geom_text(aes(label = as.integer(Percentages)), vjust = -0.5, position = position_dodge(width = 0.9)) +
                scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
                scale_x_discrete(labels = c("Cluster 0", "Cluster 1", "Cluster 2")) +
                scale_fill_manual(values = cluster_colors[c(4, 6, 7)], label = c("EpCAM+", "EpCAM+PSMA+", "PSMA+")) +
                labs(fill = "Treatment cycles", x = expression("CTC clusters"), y = expression("Cell number percentages")) +
                ggtitle("CTC clusters - cell distribution based on surface markers") +
                #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
                theme_classic() +
                theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
                      strip.text = element_blank())
            print(surface_m_bar_p)
            
            ggsave(filename = paste0("Cluster_composition-Surface_markers", script_version, ".png"), surface_m_bar_p,
                   device = "png", path = paste0(script_version, "/", "Plots/"),
                   width = 1500, height = 1000, units = "px")
            rm(surface_m_bar_p)
            
            
            # Plotting message
            message("The requested plots have been printed and saved into the: ", paste0(script_version, "/", "Plots/"), "folder. \n")
            
        } else {
            
            # No plotting message
            message("No plots were requested, therefore no plots will be printed or saved.")
            
        }
        
        
}
calculate_and_plot_cluster_composition(make_plots = TRUE, return_statistics = TRUE, print_statistics = TRUE)

#################################################################   End Section  ##################################################################






#################################################################  Condition influence on clustering  ##################################################################
                                                                #       (Statistical testing)         #




## Calculating the association between various conditions and clustering using Fisher's exact test an Chi square test


# This function will simply calculate if the association between the culstering and various conditions is significant or not (the purpose of this function is
# to make this code block more concise)

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
    
    
    # Detailed treatment response
    detailed_resp_ct <- matrix(ncol = 3, nrow = 2)
    colnames(detailed_resp_ct) <- c("PD", "PR", "SD")
    rownames(detailed_resp_ct) <- c("Cluster_0", "Cluster_1")
    detailed_resp_ct[1, ] <- stat_list[[1]][[5]]$Cell_numbers_c0[c(2:4)]
    detailed_resp_ct[2, ] <- stat_list[[1]][[5]]$Cell_numbers_c1[c(2:4)]
    
    
    
    # Compile the tables together into a list
    contingency_lst <- list (treatment_cyc_ct, treatment_type_ct, cell_cyc_ct,
                             overall_resp_ct, detailed_resp_ct)
    names(contingency_lst) <- c("treatment_cyc_ct", "treatment_type_ct", "cell_cyc_ct",
                                "overall_resp_ct", "detailed_resp_ct")
    
    
    ## Calculate the association with Fisher's exact test or Chi-square test based on the
    ## dimensions of the contingency tables
    association_stat_res <- vector(mode = "list", length = length(contingency_lst))
    names(association_stat_res) <- c("treatment_cyc_ct", "treatment_type_ct", "cell_cyc_ct",
                                     "overall_resp_ct", "detailed_resp_ct")
    for (i in seq_along(contingency_lst)) {
        
        if (ncol(contingency_lst[[i]]) > 2) {
            
            association_stat_res[[i]] <- chisq.test(contingency_lst[[i]])
            
        } else {
            
            association_stat_res[[i]] <- fisher.test(contingency_lst[[i]])
            
        }
        
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
        write_lines(x = result_to_save, file = paste0(script_version, "/", "Association_test_sesult", "_",
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
DE_Cluster_0v1 <- FindMarkers(sct_CTC.obj, ident.1 = "0", ident.2 = "1")
head(DE_Cluster_0v1)


# This is the parallelized version of the gene name matching function
geneName_matching_par = function (DE_lst, gene_lst, ID_col, name_col, n_cores = 2) {
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
    geneNames <- foreach::foreach(i = seq_len(nrow(DE_lst)), .options.snow = opts, .combine = c) %dopar% {
        position <- grep(rownames(DE_lst)[i], gene_lst[, ID_col])
        out[i] <- gene_lst[, name_col][position]
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
geneName_matching_par(DE_Cluster_0v1, unif_gene_names, "ensembl_gene_id", "external_gene_name")

head(DE_Cluster_0v1_geneNames)


# NOTE: this version is less time intensive (runs under 5 minutes) as it runs the process on multiple cores in parallel, however
# it relies on on a foreach "loop" which is more complicated then a simple for loop


# Save the gene name matched DEG list
write_csv(x = DE_Cluster_0v1_geneNames, file = paste0(script_version, "/", "DE_Cluster_0v1_geneNames", script_version, ".csv"))


# Checking if there is differential expression between the resp groups (not according to the clustering)
Idents(sct_CTC.obj) <- "RespGroup"
DE_RespvNonresp <- FindMarkers(sct_CTC.obj, ident.1 = "Responder", ident.2 = "Nonresponder")

head(DE_RespvNonresp)
summary(DE_RespvNonresp$p_val_adj)


# Checking if there is differential expression between treatment cycles
Idents(sct_CTC.obj) <- "RealCycle"
DE_TreatmentCycles <- FindAllMarkers(sct_CTC.obj, assay = "SCT")

head(DE_TreatmentCycles)
summary(DE_TreatmentCycles$p_val_adj)

print(dplyr::filter(DE_TreatmentCycles, cluster == "0", p_val_adj <= 0.05))
print(dplyr::filter(DE_TreatmentCycles, cluster == "1", p_val_adj <= 0.05))
print(dplyr::filter(DE_TreatmentCycles, cluster == "2", p_val_adj <= 0.05))


# NOTE: I checked each treatment cycle separately to see where the difference come form, and it seems that there is no significant difference
# between Cycle 0 - 2 and Cycle 1 - 2, only between Cycle 0 - 1, however here the p_val_adj values are different.
# To try to deal with this I will run a simple FindMarkers with Cycle 0 and Cycle 1 as identities.


# Checking if there is differential expression between treatment cycles 0 and 1
DE_TreatmentCycles_0v1 <- Seurat::FindMarkers(sct_CTC.obj, ident.1 = "0", ident.2 = "1")


# Gene name matching with treatment cycle DEG markers

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

geneName_matching_par_V2(DE_lst = DE_TreatmentCycles_0v1, gene_lst = unif_gene_names, ID_col = "ensembl_gene_id", name_col = "external_gene_name")

head(DE_TreatmentCycles_0v1_geneNames)


# Saving the list of treatment cycle DEG markers
write_csv(x = DE_TreatmentCycles_0v1_geneNames, file = paste0(script_version, "/", "DE_TreatmentCycles_geneNames", script_version, ".csv"))


# Looking for conserved markers between the clusters (the function needs a grouping variable, which I set to the same across all cells)
Idents(sct_CTC.obj) <- "seurat_clusters"
sct_CTC.obj[["grouping.var"]] <- "RNA"
CON_Cluster_0v1 <- FindConservedMarkers(sct_CTC.obj, ident.1 = "0", ident.2 = "1", assay = "SCT", grouping.var = "grouping.var")
head(CON_Cluster_0v1)


# Gene name matching with the conserved markers
geneName_matching_par_V2(CON_Cluster_0v1, DE_lst_IDs_are_rows = TRUE, gene_lst =  unif_gene_names,
                         ID_col = "ensembl_gene_id", name_col = "external_gene_name")

head(CON_Cluster_0v1_geneNames)


# Saving the list of consereved markers
write_csv(x = CON_Cluster_0v1_geneNames, file = paste0(script_version, "/", "CON_Cluster_0v1_geneNames", script_version, ".csv"))


# Remove unnecessary objects
rm(DE_RespvNonresp, CON_Cluster_0v1, CON_Cluster_0v1_geneNames)




## Assigning ENTREZ IDs to the DE gene list


# NOTE: this step will be needed later anyway for the SPIA analysis. For now I want to see how the DE list and GSEA will look like
# if I use ENTREZ IDs instead of gene names (in this case, only annotated genes should appear, no pseudo genes)


# This if Statment will check if the final version of the UnifiedEntrezIds df is present in the main working directory. If yes,
# the file will be loaded, if not the processing will be ran.
if (exists("unifiedEntrezIDs", where = .GlobalEnv)) {
    
  message("The unifiedEntrezIDs variable is already present in the .GlobalEnv.")  
    
} else if (file.exists("UnifiedEntrezIDs_final.csv") == TRUE) {
    
    unifiedEntrezIDs <- read.csv(file = "UnifiedEntrezIDs_final.csv")
    message("The final, fully prepared version of the ENTREZ IDs list is present in the working directory, \n",
            "and will be loaded as: ", "unifiedEntrezIDs")
    
} else {
    
    message("The final, fully prepared version of the ENTREZ IDs list was not found in the working directory, \n",
            "Running the following code block to prpare it. \n")
    
    
    # This if statement will check if the entrezID table is already present in the working directory and if yes it will load it
    if (file.exists("Symbols_ENSEMBL_ENTREZ_ID_table.csv") & !exists(x = "entrezIDs", where = .GlobalEnv)) {
        
        entrezIDs <- read.csv(file = "Symbols_ENSEMBL_ENTREZ_ID_table.csv")
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
        write_csv(x = entrezIDs, file = "Symbols_ENSEMBL_ENTREZ_ID_table.csv")
        
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
    write_csv(x = entrezIDs, file = paste0(script_version, "/", "ENSEMBL_ID_list_with_ENTREZ_IDs_and_gene_names", script_version, ".csv"))
    
    
    # Load the saved entrezIDs df if not loaded yet
    if (!exists(x = "entrezIDs", where = .GlobalEnv)) {
        entrezIDs <- read.csv(file = paste0(script_version, "/", "ENSEMBL_ID_list_with_ENTREZ_IDs_and_gene_names", script_version, ".csv"))
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
    write_csv(x = unifiedEntrezIDs, file = "UnifiedEntrezIDs_final.csv")
    
}


# Attaching ENTREZ IDs to the DE gene lists using the following functions
# This is the original, well trusted for loop based implementation
entrezID_matching = function (DE_lst, entrezID_lst, ENSEMBL_ID_col, ENTREZ_ID_col, new_object_name = NA) {
    progressBar <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                    total = nrow(DE_lst),
                                    complete = "=",
                                    incomplete = "-",
                                    current = ">",
                                    clear = FALSE,
                                    width = 100)
    position <- c()
    for(i in seq_along(rownames(DE_lst))) {
        progressBar$tick()
        
        position <-  grep(rownames(DE_lst)[i], entrezID_lst[, ENSEMBL_ID_col])
        DE_lst$entrezID[i] <- entrezID_lst[, ENTREZ_ID_col][position]
        
    }
    
    # This new if statment allows the user to specify a new object anme, but could just go with the standard nameing
    # the function uses
    if(is.na(new_object_name)) {
        original_name <- deparse(substitute(DE_lst))
        output_name <- paste0(original_name, "_", "entrezIDs")
    } else {
        output_name <- new_object_name
    }

    assign(output_name, DE_lst, envir = .GlobalEnv)
    message("The entrezID matched ", DE_lst, " list was assigned to a new object named:", "/n", output_name)
}
entrezID_matching(DE_Cluster_0v1, unifiedEntrezIDs, "ensembl_gene_id", "entrezgene_id")

# NOTE: this for loop is highly time and resource intensive and will take over 5 minutes to run.
# in exchange however it will insert a positionally accurate name for the DE genes in a well understood manner


# This is the parallelized version of the matching function and the upgraded V2, which should work with the FindAllMarkers output
entrezID_matching_par = function (DE_lst, entrezID_lst, ENSEMBL_ID_col, ENTREZ_ID_col, n_cores = 2, new_object_name = NA) {
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
    entrezIDs <- foreach::foreach(i = seq_len(nrow(DE_lst)), .options.snow = opts, .combine = c) %dopar% {
        position <- grep(rownames(DE_lst)[i], entrezID_lst[, ENSEMBL_ID_col])
        out[i] <- entrezID_lst[, ENTREZ_ID_col][position]
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

# NOTE: this version is less time intensive (runs under 5 minutes) as it runs the process on multiple cores in parallel, however
# it relies on on a foreach "loop" which is more complicated then a simple for loop


# Attaching entrezIDs to the cluster DEG list
entrezID_matching_par_V2(DE_Cluster_0v1, unifiedEntrezIDs, "ensembl_gene_id", "entrezgene_id")

head(DE_Cluster_0v1_entrezIDs)


# Attaching entrezIDs to the treatment cycle DEG list
entrezID_matching_par_V2(DE_lst = DE_TreatmentCycles_0v1, entrezID_lst = unifiedEntrezIDs, ENSEMBL_ID_col = "ensembl_gene_id",
                         ENTREZ_ID_col = "entrezgene_id")

head(DE_TreatmentCycles_entrezIDs)


# I will combine the gene name and the ENTREZ ID DE lists by feeding the gene name list into the 
# ENTREZ ID matching function for the cluster DEG list
entrezID_matching_par_V2(DE_Cluster_0v1_geneNames, unifiedEntrezIDs, "ensembl_gene_id", "entrezgene_id")

head(DE_Cluster_0v1_geneNames_entrezIDs)


# I will combine the gene name and the ENTREZ ID DE lists by feeding the gene name list into the 
# ENTREZ ID matching function for the treatment cycle DEG list
entrezID_matching_par_V2(DE_lst = DE_TreatmentCycles_0v1_geneNames, entrezID_lst = unifiedEntrezIDs, ENSEMBL_ID_col = "ensembl_gene_id",
                         ENTREZ_ID_col = "entrezgene_id")

head(DE_TreatmentCycles_0v1_geneNames_entrezIDs)


# Save the raw generated ENSEMBL, ENTREZ ID and gene name (SYMBOL) containing table (no additional modifications)
write_csv(x = DE_Cluster_0v1_geneNames_entrezIDs, file = paste0(script_version, "/", "DEG_clusters_0v1_raw", script_version, ".csv"))


# Save the raw generated ENSEMBL, ENTREZ ID and gene name (SYMBOL) containing table (no additional modifications)
write_csv(x = DE_TreatmentCycles_0v1_geneNames_entrezIDs, file = paste0(script_version, "/", "DEG_treatment_cycles_raw", script_version, ".csv"))


# Remove unnecessary variables
rm(ensembl, datasets, attributes, filters, entrezIDs, dropped_IDs, dropped_IDs_df, missing_entrezIDs, entrezIDs_p2)





## Check if the two DE lists are significantly different using Fisher's exact test
## NOTE: the H0 here is that the two lists are significantly different from each other


# Construct the empty contingency list
DE_lst_contingency_t <- matrix(nrow = 2, ncol = 2)
colnames(DE_lst_contingency_t) <- c("in_cluster_lst", "not_in_cluster_lst")
rownames(DE_lst_contingency_t) <- c("in_cycle_lst", "not_in_cycle_lst")
print(DE_lst_contingency_t)


# Calculated the number of detected genes not found in either of the lists
not_in_either <- unif_gene_names$ensembl_gene_id[!unif_gene_names$ensembl_gene_id %in% rownames(DE_Cluster_0v1_geneNames_entrezIDs)]
not_in_either <- not_in_either[!not_in_either %in% rownames(DE_TreatmentCycles_0v1_geneNames_entrezIDs)]


# Fill the list with the appropriate counts based on similarities and differences
DE_lst_contingency_t[1, 1] <- as.numeric(summary(rownames(DE_Cluster_0v1_geneNames_entrezIDs) %in% rownames(DE_TreatmentCycles_0v1_geneNames_entrezIDs))[3])
DE_lst_contingency_t[1, 2] <- as.numeric(summary(rownames(DE_Cluster_0v1_geneNames_entrezIDs) %in% rownames(DE_TreatmentCycles_0v1_geneNames_entrezIDs))[2])
DE_lst_contingency_t[2, 1] <- as.numeric(summary(rownames(DE_TreatmentCycles_0v1_geneNames_entrezIDs) %in% rownames(DE_Cluster_0v1_geneNames_entrezIDs))[2])
DE_lst_contingency_t[2, 2] <- length(not_in_either)
print(DE_lst_contingency_t)


# Conduct Fisher's exact test and capture the nice table result for saving
DE_lists_similarity_test <- capture.output(print(fisher.test(DE_lst_contingency_t)))


# Save the results as .txt
write_lines(DE_lists_similarity_test, file = paste0(script_version, "/", "TreatmentCycel_vs_Cluster_DEG_similarity_test", script_version, ".txt"))


# Remove unnecessary variables
rm(DE_TreatmentCycles, DE_TreatmentCycles_0v1, DE_TreatmentCycles_0v1_geneNames,
   DE_TreatmentCycles_0v1_entrezIDs, DE_TreatmentCycles_0v1_geneNames_entrezIDs,
   DE_lst_contingency_t, DE_lists_similarity_test)




## Further processing the DE gene list (continuing with the geneName + ENTREZ ID list)
## NOTE: as the Fisher's test showed, the two DE lists are not significantly different
## therefore I will continue only with the cluster DEG list


# Attaching some additional info to the cluster DE gene list and saving the whole list
DE_Cluster_0v1_geneNames_entrezIDs$ExprChange <- ifelse(DE_Cluster_0v1_geneNames_entrezIDs$avg_log2FC >= 0 & DE_Cluster_0v1_geneNames_entrezIDs$p_val_adj <= 0.05, "Upregulated",
                                    ifelse(DE_Cluster_0v1_geneNames_entrezIDs$avg_log2FC <= 0 & DE_Cluster_0v1_geneNames_entrezIDs$p_val_adj <= 0.05, "Downregulated", "Not significant"))

# NOTE: the ranking here will be solely based on the p_adj value (as the sign() function will turn all log2FC into a negative, 0 or positive 1)
# if you don't want that (which we want to do here) you can omit the sign() function
DE_Cluster_0v1_geneNames_entrezIDs$Ranking <- DE_Cluster_0v1_geneNames_entrezIDs$avg_log2FC * (-log10(DE_Cluster_0v1_geneNames_entrezIDs$p_val_adj))


# Save the resulting cluster DEG list
write_csv(x = DE_Cluster_0v1_geneNames_entrezIDs, file = paste0(script_version, "/", "DEG_clusters_0v1", script_version, ".csv"))




## As for further analysis, I will be mainly using the gene names and ENTREZ IDs, I will remove genes from the DEG list which have no ENTREZ ID or gene name


# Load the extended DEG list if not loaded yet
if (!exists(x = "DE_Cluster_0v1_geneNames_entrezIDs", where = .GlobalEnv)) {
    DE_Cluster_0v1_geneNames_entrezIDs <- read.csv(file = paste0(script_version, "/", "DEG_clusters_0v1", script_version, ".csv"))
    message("The extended DEG list was loaded under the name: ", deparse(substitute(DE_Cluster_0v1_geneNames_entrezIDs)))
}


# Filter the original list and assign it to a new, clean variable by entrezIDs
DEG_cluster_0v1_clean <- dplyr::filter(DE_Cluster_0v1_geneNames_entrezIDs, entrezIDs != "NA")


# Filter the original list and assign it to a new, clean variable by geneNames
DEG_cluster_0v1_symbols_clean <- dplyr::filter(DE_Cluster_0v1_geneNames_entrezIDs, geneNames != "NA")


# Save the new, clean DEG lists for future use
write_csv(x = DEG_cluster_0v1_clean, file = paste0(script_version, "/", "DEG_clusters_0v1_clean", script_version, ".csv"))

write_csv(x = DEG_cluster_0v1_symbols_clean, file = paste0(script_version, "/", "DEG_clusters_0v1_symbols_clean", script_version, ".csv"))


# Remove unnecessary objects
rm(DE_Cluster_0v1, DE_Cluster_0v1_entrezIDs, DE_Cluster_0v1_geneNames, DE_Cluster_0v1_geneNames_entrezIDs)




## Visualizing DE genes using a volcano plot


# Load the DEG list if not loaded yet
if (!exists(x = "DEG_cluster_0v1_clean", where = .GlobalEnv) || !exists(x = "DEG_cluster_0v1_symbols_clean", where = .GlobalEnv)) {
    DEG_cluster_0v1_clean <- read.csv(file = paste0(script_version, "/", "DEG_clusters_0v1_clean", script_version, ".csv"))
    DEG_cluster_0v1_symbols_clean <- read.csv(file = paste0(script_version, "/", "DEG_clusters_0v1_symbols_clean", script_version, ".csv"))
    message("The clean DEG lists were loaded under the names: ", deparse(substitute(DE_Cluster_0v1_geneNames_entrezIDs)), " and ", deparse(substitute(DE_Cluster_0v1_geneNames_entrezIDs)))
}


# Annotated visualization with ggplot2 for entrezIDs
DEG_cluster_0v1_clean$Significance <- ifelse(DEG_cluster_0v1_clean$avg_log2FC > 0 & DEG_cluster_0v1_clean$p_val_adj < 0.05, "Upregulated",
                                    ifelse(DEG_cluster_0v1_clean$avg_log2FC < 0 & DEG_cluster_0v1_clean$p_val_adj < 0.05, "Downregulated", "Not significant")) #labeling genes for coloring
upreg <- DEG_cluster_0v1_clean[DEG_cluster_0v1_clean$Significance == "Upregulated", ] #sub-setting the upregulated genes
topupreg <- head(upreg$geneName[order(upreg$p_val_adj)], 10) #selecting the top upreg genes
downreg <- DEG_cluster_0v1_clean[DEG_cluster_0v1_clean$Significance == "Downregulated", ] #sub-setting the downregulated genes
topdownreg <- head(downreg$geneName[order(downreg$p_val_adj)], 10) #selecting the top downreg genes


top50upreg <- head(DEG_cluster_0v1_clean$geneName[order(DEG_cluster_0v1_clean$p_val_adj)], 50) #choosing the top 50 differentially expressed genes
DEG_cluster_0v1_clean$DElabel <- ifelse(DEG_cluster_0v1_clean$geneNames %in% topupreg | DEG_cluster_0v1_clean$geneNames %in% topdownreg,
                                        DEG_cluster_0v1_clean$geneNames, NA) #marking the top DEGs in the dataframe

volcano_p <- ggplot(data = DEG_cluster_0v1_clean, 
                aes(x = avg_log2FC, y = -log10(p_val_adj), col = Significance, label = DElabel)) +
            geom_point() +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
            scale_color_manual(values = cluster_colors[c(2, 5, 3)],
                            labels = c("Downregulated", "Not significant", "Upregulated")) +
            scale_x_continuous(limits = c(-15, 20), breaks = seq(-15, 20, 2)) +
            scale_y_continuous(limits = c(0, 14), breaks = seq(0, 14, 2)) +
            labs(color = "Expression status", x = expression("log"[2] * "FC"), y = expression("-log"[10] * "p-value")) +
            ggtitle("Differentially expressed genes - entrezIDs") +
            geom_text_repel(show.legend = FALSE, max.overlaps = Inf) +
            theme_classic() +
            theme(plot.title = element_text(hjust = 0.5))
print(volcano_p)

ggsave(filename = paste0("DEGs_Volcano_plot_entrezIDs", script_version, ".png"), volcano_p,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 3500, height = 3000, units = "px")
rm(volcano_p)


# Construct an upreg-downreg df of entrezIDs for saving
top_upreg_and_downreg_DEGs <- data.frame(top_upregulated = topupreg, top_downregulated = topdownreg)


# Save the upreg-downreg df of entrezIDs
write_csv(x = top_upreg_and_downreg_DEGs, file = paste0(script_version, "/", "top_upreg_and_downreg_DEGs", script_version, ".csv"))


# Annotated visualization with ggplot2 for geneNames
DEG_cluster_0v1_symbols_clean $Significance <- ifelse(DEG_cluster_0v1_symbols_clean $avg_log2FC > 0 & DEG_cluster_0v1_symbols_clean $p_val_adj < 0.05, "Upregulated",
                                    ifelse(DEG_cluster_0v1_symbols_clean $avg_log2FC < 0 & DEG_cluster_0v1_symbols_clean $p_val_adj < 0.05, "Downregulated", "Not significant")) #labeling genes for coloring
upreg <- DEG_cluster_0v1_symbols_clean [DEG_cluster_0v1_symbols_clean $Significance == "Upregulated", ] #sub-setting the upregulated genes
topupreg <- head(upreg$geneName[order(upreg$p_val_adj)], 10) #selecting the top upreg genes
downreg <- DEG_cluster_0v1_symbols_clean [DEG_cluster_0v1_symbols_clean $Significance == "Downregulated", ] #sub-setting the downregulated genes
topdownreg <- head(downreg$geneName[order(downreg$p_val_adj)], 10) #selecting the top downreg genes


top50upreg <- head(DEG_cluster_0v1_symbols_clean $geneName[order(DEG_cluster_0v1_symbols_clean $p_val_adj)], 50) #choosing the top 50 differentially expressed genes
DEG_cluster_0v1_symbols_clean $DElabel <- ifelse(DEG_cluster_0v1_symbols_clean $geneNames %in% topupreg | DEG_cluster_0v1_symbols_clean $geneNames %in% topdownreg,
                                        DEG_cluster_0v1_symbols_clean $geneNames, NA) #marking the top DEGs in the dataframe

volcano_p <- ggplot(data = DEG_cluster_0v1_symbols_clean , 
                aes(x = avg_log2FC, y = -log10(p_val_adj), col = Significance, label = DElabel)) +
            geom_point() +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
            scale_color_manual(values = cluster_colors[c(2, 5, 3)],
                            labels = c("Downregulated", "Not significant", "Upregulated")) +
            scale_x_continuous(limits = c(-15, 20), breaks = seq(-15, 20, 2)) +
            scale_y_continuous(limits = c(0, 14), breaks = seq(0, 14, 2)) +
            labs(color = "Expression status", x = expression("log"[2] * "FC"), y = expression("-log"[10] * "p-value")) +
            ggtitle("Differentially expressed genes - geneNames") +
            geom_text_repel(show.legend = FALSE, max.overlaps = Inf) +
            theme_classic() +
            theme(plot.title = element_text(hjust = 0.5))
print(volcano_p)

ggsave(filename = paste0("DEGs_Volcano_plot_symbols", script_version, ".png"), volcano_p,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 3500, height = 3000, units = "px")
rm(volcano_p)

# Construct an upreg-downreg df of gene names for saving
top_upreg_and_downreg_DEGs_symbols <- data.frame(top_upregulated = topupreg, top_downregulated = topdownreg)


# Save the upreg-downreg df of gene names
write_csv(x = top_upreg_and_downreg_DEGs_symbols, file = paste0(script_version, "/", "top_upreg_and_downreg_DEGs_symbols", script_version, ".csv"))


# Remove unnecessary variables
rm(upreg, downreg, top50upreg, topupreg, topdownreg, top_upreg_and_downreg_DEGs, top_upreg_and_downreg_DEGs_symbols)

#################################################################   End Section  ##################################################################






#################################################################               GSEA                  ##################################################################
                                                                #   (Gene Set Enrichment Analysis)    #


# NOTE: for the gene set enrichment analysis I will use the clusterProfiler r package




## Creating and processing a new DE list  of genes for the GSEA


# Create the new DEG list
# Note: this will include all genes in the set, not just the DE genes.
Seurat::Idents(sct_CTC.obj) <- "seurat_clusters"
GSEA_DE_lst <- Seurat::FindMarkers(sct_CTC.obj, 
                           ident.1 = "0", 
                           ident.2 = "1", 
                           logfc.threshold = -Inf, 
                           min.pct = -Inf, 
                           min.diff.pct = -Inf, 
                           verbose = TRUE)


# Attaching gene names to the DE genes using the following functions
# This is the parallelized version of the matching function (defined before)
geneName_matching_par_V2(GSEA_DE_lst, unif_gene_names, "ensembl_gene_id", "external_gene_name") 

# Attaching the ENTREZ IDs to the named DEG list
entrezID_matching_par_V2(GSEA_DE_lst_geneNames, unifiedEntrezIDs, "ensembl_gene_id", "entrezgene_id")

head(GSEA_DE_lst_geneNames_entrezIDs)




## Create a ranked gene list for the GSEA based on the log2FC and the P-adjusted values


# Attach the gene names and ENTREZ IDs to the DEG list
# Note: the ranking here will be solely based on the p_adj value (as the sign() function will turn all log2FC into a negative, 0 or positive 1)
# if you don't want that you can omit the sign() function
#GSEA_DE_lst_geneNames$Ranking <- sign(GSEA_DE_lst_geneNames$avg_log2FC) * (-log10(GSEA_DE_lst_geneNames$p_val_adj))

GSEA_DE_lst_geneNames_entrezIDs$Log_and_p_Ranking <- GSEA_DE_lst_geneNames_entrezIDs$avg_log2FC * (-log10(GSEA_DE_lst_geneNames_entrezIDs$p_val_adj))


# Adding an additional column with the expression change significance
GSEA_DE_lst_geneNames_entrezIDs$Significance <- ifelse(GSEA_DE_lst_geneNames_entrezIDs$avg_log2FC > 0 & GSEA_DE_lst_geneNames_entrezIDs$p_val_adj < 0.05, "Upregulated",
                                                ifelse(GSEA_DE_lst_geneNames_entrezIDs$avg_log2FC < 0 & GSEA_DE_lst_geneNames_entrezIDs$p_val_adj < 0.05, "Downregulated", "Not significant")) #labeling genes for coloring


# Save the generated GSEA_DEG list for further use
write_csv(x = GSEA_DE_lst_geneNames_entrezIDs, file = paste0(script_version, "/", "GSEA_DEG_list", script_version, ".csv"))




## DEG list QC


# Load the GSEA DEG list if not present yet
if (!exists("GSEA_DE_lst_geneNames_entrezIDs", where = .GlobalEnv)) {
    GSEA_DE_lst_geneNames_entrezIDs <- read.csv(file = paste0(script_version, "/", "GSEA_DEG_list", script_version, ".csv"))
    message("The GSEA DEG list was loaded as the following object: ", deparse(substitute(GSEA_DE_lst_geneNames_entrezIDs)))
}


# In order to run run GSEA we should remove any NAs (not converted IDs) and duplicated entries
# This will be done based on either the ENTREZ IDs or the geneNames 
summary(is.na(GSEA_DE_lst_geneNames_entrezIDs$entrezIDs))
summary(is.na(GSEA_DE_lst_geneNames_entrezIDs$geneNames))
GSEA_DE_lst_geneNames_entrezIDs[!is.na(GSEA_DE_lst_geneNames_entrezIDs$entrezIDs), ]
GSEA_DE_lst_geneNames_entrezIDs[!is.na(GSEA_DE_lst_geneNames_entrezIDs$geneNames), ]


# Remove the NAs from the GSEA DEG lists
GSEA_DE_lst_geneNames_entrezIDs_symbols <- GSEA_DE_lst_geneNames_entrezIDs[!is.na(GSEA_DE_lst_geneNames_entrezIDs$geneNames), ]

GSEA_DE_lst_geneNames_entrezIDs <- GSEA_DE_lst_geneNames_entrezIDs[!is.na(GSEA_DE_lst_geneNames_entrezIDs$entrezIDs), ]
GSEA_DE_lst_geneNames_entrezIDs <- GSEA_DE_lst_geneNames_entrezIDs[!is.na(GSEA_DE_lst_geneNames_entrezIDs$geneNames), ]


# Next, check if there are any duplicated entries
summary(duplicated(GSEA_DE_lst_geneNames_entrezIDs$entrezIDs))
summary(duplicated(GSEA_DE_lst_geneNames_entrezIDs$geneNames))
GSEA_DE_lst_geneNames_entrezIDs[duplicated(GSEA_DE_lst_geneNames_entrezIDs$entrezIDs), ]
GSEA_DE_lst_geneNames_entrezIDs[duplicated(GSEA_DE_lst_geneNames_entrezIDs$geneNames), ]

GSEA_DE_lst_geneNames_entrezIDs_symbols[duplicated(GSEA_DE_lst_geneNames_entrezIDs$geneNames), ]

# Subset the duplicated entries
dup_entrezIDs <- GSEA_DE_lst_geneNames_entrezIDs$entrezIDs[duplicated(GSEA_DE_lst_geneNames_entrezIDs$entrezIDs)]
dup_geneName <- GSEA_DE_lst_geneNames_entrezIDs$geneNames[duplicated(GSEA_DE_lst_geneNames_entrezIDs$geneNames)]


# Inspect the duplicated entries
dplyr::filter(GSEA_DE_lst_geneNames_entrezIDs, entrezIDs %in% dup_entrezIDs)
dplyr::filter(GSEA_DE_lst_geneNames_entrezIDs, geneNames %in% dup_geneName)

dplyr::filter(GSEA_DE_lst_geneNames_entrezIDs_symbols, geneNames %in% dup_geneName)

# Fix the duplication issue
# NOTE: in the case of the duplicated entrezIDs, I double checked the genes and one of the ID is assigned wrong
# therefore I will fix the aforementioned ID
GSEA_DE_lst_geneNames_entrezIDs["ENSG00000205571", "entrezIDs"] <- "6607"


# NOTE: in the case of the duplicated gene name, it is caused by a readthrough to the neighboring gene.
# I will fix this  by amending the gene name with one found on NCBI
GSEA_DE_lst_geneNames_entrezIDs_symbols["ENSG00000168255", "geneNames"] <- "POLR2J3-UPK3BL2"


# Assign the clean GSEA DEG lists to a new variable
GSEA_DEG_list_clean <- GSEA_DE_lst_geneNames_entrezIDs

GSEA_DEG_list_symbols_clean <- GSEA_DE_lst_geneNames_entrezIDs_symbols

# Save the clean GSEA DEG lists for further use
write_csv(x = GSEA_DEG_list_clean, file = paste0(script_version, "/", "GSEA_DEG_list_clean", script_version, ".csv"))

write_csv(x = GSEA_DEG_list_symbols_clean, file = paste0(script_version, "/", "GSEA_DEG_list_symbols_clean", script_version, ".csv"))

# Remove unnecessary variables
rm(GSEA_DE_lst, GSEA_DE_lst_geneNames, GSEA_DE_lst_geneNames_entrezIDs, GSEA_DEG_list_clean_symbols)



## Create a raked vector, which is the input format for the GSEA from the clusterProfiler package


# Load the clean GSEA DEG list if not present yet
if (!exists("GSEA_DEG_list_clean", where = .GlobalEnv) || !exists("GSEA_DEG_list_symbols_clean", where = .GlobalEnv)) {
    GSEA_DEG_list_clean <- read.csv(file = paste0(script_version, "/", "GSEA_DEG_list_clean", script_version, ".csv"))
    GSEA_DEG_list_symbols_clean <- read.csv(file = paste0(script_version, "/", "GSEA_DEG_list_symbols_clean", script_version, ".csv"))
    message("The clean GSEA DEG lists were loaded as the following variables: ", deparse(substitute(GSEA_DEG_list_clean)), " and ", deparse(substitute(GSEA_DEG_list_symbols_clean)))
}


# Create the ranked list using the ENTREZ IDs
GSEA_DEG_ranked_list <- GSEA_DEG_list_clean$Log_and_p_Ranking
names(GSEA_DEG_ranked_list) <- GSEA_DEG_list_clean$entrezIDs
head(GSEA_DEG_ranked_list)


# Create the ranked list using the geneNames
GSEA_DEG_symbols_ranked_list <- GSEA_DEG_list_symbols_clean$Log_and_p_Ranking
names(GSEA_DEG_symbols_ranked_list) <- GSEA_DEG_list_symbols_clean$geneNames
head(GSEA_DEG_symbols_ranked_list)



## Run some basic QC on the ranked list to sort out infinite elements and to see if the distribution is correct


# Check if there are any infinite values present or not
summary(GSEA_DEG_ranked_list)

summary(GSEA_DEG_symbols_ranked_list)


# Fix any remaining issues like NAs or Inf values
GSEA_DEG_ranked_list <- GSEA_DEG_ranked_list[!is.na(GSEA_DEG_ranked_list)]

GSEA_DEG_symbols_ranked_list <- GSEA_DEG_symbols_ranked_list[!is.na(GSEA_DEG_symbols_ranked_list)]


# First order the ranked list, then plot it to check how the values are distributed, and how they look like
GSEA_DEG_list_ord <- sort(GSEA_DEG_ranked_list, decreasing = TRUE)
plot(GSEA_DEG_list_ord)

GSEA_DEG_symbols_list_ord <- sort(GSEA_DEG_symbols_ranked_list, decreasing = TRUE)
plot(GSEA_DEG_symbols_list_ord)


# This function will plot the top DEG genes in the GSEA DEG list using ggplot2
# NOTE: the function will take two input df names (as a string). The 1st position is for 
# the df with ENTREZ IDs and the second is for the df with gene SYMBOLS
plot_top_GSEA_DEGs = function(DEG_lists = c("GSEA_DEG_list_clean", "GSEA_DEG_list_symbols_clean"),
                              num_of_top_genes = 50, make_plots = TRUE) {
    
    if (make_plots == TRUE) {
        
        # Get the plot dfs
        ID_plot_df <- get(DEG_lists[1])
        SYMBOL_plot_df <- get(DEG_lists[2])
        
        
        # Make a DEG_ranked_df for a ggplot representation of the top 50 upreg and downreg genes
        GSEA_DEG_ranked_ID_df <- ID_plot_df[, c(6, 7, 8)]
        GSEA_DEG_ranked_ID_df_ord <- dplyr::arrange(GSEA_DEG_ranked_ID_df, desc(Log_and_p_Ranking))
        
        GSEA_DEG_ranked_SYMBOL_df <- ID_plot_df[, c(6, 7, 8)]
        GSEA_DEG_ranked_SYMBOL_df_ord <- dplyr::arrange(GSEA_DEG_ranked_SYMBOL_df, desc(Log_and_p_Ranking))
        
        
        # Plot cluster 0 top upregulated
        p1_ID <- ggplot(GSEA_DEG_ranked_ID_df_ord[1:num_of_top_genes, ],
                        aes(x = entrezIDs, y = Log_and_p_Ranking)) +
                        geom_point() +
                        ggtitle(paste0("Cluster 0 top ", num_of_top_genes, "upregulated genes - ENTREZ IDs")) +
                        theme_classic() +
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        print(p1_ID)
        
        p1_SYMBOL <- ggplot(GSEA_DEG_ranked_SYMBOL_df_ord[1:num_of_top_genes, ],
                            aes(x = geneNames, y = Log_and_p_Ranking)) +
                            geom_point() +
                            ggtitle(paste0("Cluster 0 top ", num_of_top_genes, "upregulated genes - gene SYMBOLS")) +
                            theme_classic() +
                            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        print(p1_SYMBOL)
        
        # Plot cluster 1 top upregulated
        p2_ID <- ggplot(GSEA_DEG_ranked_ID_df_ord[nrow(GSEA_DEG_ranked_ID_df_ord):(nrow(GSEA_DEG_ranked_ID_df_ord) - num_of_top_genes), ],
                        aes(x = entrezIDs, y = Log_and_p_Ranking)) +
                        geom_point() +
                        ggtitle(paste0("Cluster 0 top ", num_of_top_genes, "upregulated genes - ENTREZ IDs")) +
                        theme_classic() +
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        print(p2_ID)
        
        p2_SYMBOL <- ggplot(GSEA_DEG_ranked_SYMBOL_df_ord[nrow(GSEA_DEG_ranked_SYMBOL_df_ord):(nrow(GSEA_DEG_ranked_SYMBOL_df_ord) - num_of_top_genes), ],
                            aes(x = entrezIDs, y = Log_and_p_Ranking)) +
                            geom_point() +
                            ggtitle(paste0("Cluster 0 top ", num_of_top_genes, "upregulated genes - gene SYMBOLS")) +
                            theme_classic() +
                            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        print(p2_SYMBOL)
        
        # End message
        message("The requested GSEA top ", num_of_top_genes, " DEGs were successfully plotted.")
        
    } else {
        
        message("No GSEA DEG plots were requested.")
        
    }
    
}

plot_top_GSEA_DEGs(DEG_lists = c("GSEA_DEG_list_clean", "GSEA_DEG_list_symbols_clean"),
                   num_of_top_genes = 50, make_plots = TRUE)




## Prepare the reference .gmt and other necessary files for the GSEA run and run the analysis

# NOTE: this GSEA run will be done based on the entrezIDs. It can be done based on the gene names as well with an adjusted
# named, ranked list.


# Define the path to the .gmt (reference sets) files
gmt_files <- paste0(getwd(), "/", "GSEA_entrez_ref_sets")
gmt_files_symbols <- paste0(getwd(), "/", "GSEA_symbols_ref_sets")


# List the reference pathway files
pathways <- list.files(gmt_files)
pathways_symbols <- list.files(gmt_files_symbols)


# To bulk read in .gmt files for the clusterProfiler package
read_gmt = function(path) {
    # Listing the files to load
    file_names <- list.files(path = path)
    
    # Creates a progress bar for the file loading
    progressBar <- progress::progress_bar$new(format = "Reading gmt - (:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                    total = length(file_names),
                                    complete = "=",
                                    incomplete = "-",
                                    current = ">",
                                    clear = FALSE,
                                    width = 100)
    
    # This loop, will ensure loading and processing of the gmt files in a sequential order
    tmp_gmt <- list()
    name_to_assign <- c()
    for (e in seq_along(file_names)) {
        # Initiate progress bar
        progressBar$tick()
        message("\n" ,"loading .gmt file: ", e)
        
        # Creating the name for the new pathway object
        name_to_assign <- file_names[e]
        
        #loading a gmt file into a temporary object
        tmp_gmt <- clusterProfiler::read.gmt(paste0(path, "/", file_names[e]))
        
        #assigning the trimmed files into new objects in the .GlobalEnv
        assign(name_to_assign, tmp_gmt, envir = .GlobalEnv)
        message("A reference pathway (.gmt) was assigned to a new object:", name_to_assign)
        
    }
    
}

read_gmt(gmt_files)
read_gmt(gmt_files_symbols)


# This function is a re-worked wrapper function, which will do multiple GSEA runs based on the set of reference
# gene sets provided, assigns the result into a new object, saves the result as a .csv file, prints the GSEA plots
# and saves them as .png files and assigns the results either as indivisual files, or as a list
run_serial_GSEA = function(DE_gene_list,
                              random_seed = 1234,
                              verbose = TRUE,
                              min_size = 5,
                              max_size = 500,
                              p_val_cutoff = 0.05,
                              p_val_adjust_method = "BH",
                              GSEA_method = "fgsea",
                              nPermSimple = 10000,
                              ref_gene_lst,
                              GSEA_csv_file_name,
                              GSEA_result_object_name,
                              GSEA_objects_to_list = TRUE,
                              GSEA_result_list_name,
                              plot_title = paste0("GSEA", pathway_names[e], "set"),
                              plot_file_name,
                              plot_print_requested = TRUE,
                              save_folder = "GSEA_results") {
    
    
    # An if statement to check if a custom save folder was requested. If yes, it creates the folder and subfolder and saves the results there.
    if (save_folder == "GSEA_results") {
        
        # do nothing
        
    } else if (file.exists(paste0(script_version, "/", save_folder)) == TRUE) {
        
        message("The requested custom save folder ", save_folder, "already exists. Files will be saved there")
        
    } else {
        
        dir.create(path = paste0(script_version, "/", save_folder), recursive = TRUE)
        
        dir.create(path = paste0(script_version, "/", save_folder, "/", "Plots"), recursive = TRUE)
        
        message("A custom save folder under the name: ", save_folder, "was requested and created.", "\n", 
                "the GSEA results will be saved into the new folder.")
        
    }
    
    
    # Extracting parts of the reference pathway names to be used for result names
    pathway_names <- stringr::str_extract(ref_gene_lst, rebus::one_or_more(ASCII_ALNUM) %R% "." %R% rebus::one_or_more(ASCII_ALNUM) %R% "." %R% rebus::one_or_more(ASCII_ALNUM))
    
    
    # Creating the progress bar
    progressBar <- progress::progress_bar$new(format = "Running GSEA - (:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                              total = length(ref_gene_lst),
                                              complete = "=",
                                              incomplete = "-",
                                              current = ">",
                                              clear = FALSE,
                                              width = 100)
    
    
    # Creating an empty multi-element list to store the run results in
    GSEA_result_list <- vector(mode = "list", length = length(ref_gene_lst))
    names(GSEA_result_list) <- ref_gene_lst
    
    
    # This for loop should run the GSEA function over the reference gene/entrezID lists one by one
    for (e in seq_along(ref_gene_lst)) {
        # Initiating the progress bar
        progressBar$tick()
        
        
        # Print starting message
        message("\n" ,"Initiating GSEA on set ", ref_gene_lst[e])
        
        
        # Running GSEA using the clusterProfiler package
        # the function uses permutation testing to adjust the p-values. as this is based on a random number, for the sake of
        # reproducibility, one should set the seed
        set.seed(random_seed)
        
        
        # Run GSEA on selected dataset and assign the result into the GSEA_result variable
        try(GSEA_result <- clusterProfiler::GSEA(geneList = DE_gene_list,
                                             minGSSize = min_size,
                                             maxGSSize = max_size,
                                             pvalueCutoff = p_val_cutoff,
                                             pAdjustMethod = p_val_adjust_method,
                                             TERM2GENE = get(ref_gene_lst[e]),
                                             by = GSEA_method,
                                             verbose = verbose,
                                             nPermSimple = nPermSimple))
        
        # Convert the GSEA result into a df so it can be saved as a .csv
        GSEA_result_df <- as.data.frame(GSEA_result)
        
        # Save the GSEA result as a .csv
        write.csv(GSEA_result_df, file = paste0(script_version, "/", save_folder, "/", ref_gene_lst[e], "_", GSEA_csv_file_name, script_version, ".csv"))
        message("A GSEA run was succesfully carried out, using the reference set: ", ref_gene_lst[e], "\n",
                "The result was saved as a .csv file:", GSEA_csv_file_name, "\n",
                "The result is assigned to the object:", GSEA_result_object_name, "\n")
        
        
        # This if statment will control if the GSEA result will be saved into a list or as an individual object
        GSEA_result_list[[e]] <- GSEA_result
    }
            
    
        # Initiating progress bar for the plotting loop
        progressBar <- progress::progress_bar$new(format = "Plotting GSEA - (:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                                  total = length(GSEA_result_list),
                                                  complete = "=",
                                                  incomplete = "-",
                                                  current = ">",
                                                  clear = FALSE,
                                                  width = 100)
        
        
        # Initiating the plotting with a message
        message("Plotting the GSEA results and saving the resulting plots:\n")
        
        
        # This loop will print and save the GSEA plot as .png files
        for (i in seq_along(GSEA_result_list)) {
            # Initiating the progress bar
            progressBar$tick()
            
            for (e in seq_len(nrow(GSEA_result_list[[i]]@result))) {
                
                # Call the gseaplot function to plot the gsea result
                gsea_res_p <- clusterProfiler::gseaplot(x = GSEA_result_list[[i]],
                                                        geneSetID = GSEA_result_list[[i]]@result[e, 1],
                                                        title = paste0(plot_title, "- ", GSEA_result_list[[i]]@result[e, 1]))
                
                # An if statement to see if plot printing was requested or not
                if (plot_print_requested == TRUE) {
                    # Visualize the resulting plot
                    print(gsea_res_p)
                }
                
                # Saving the resulting plot as a .png
                ggplot2::ggsave(filename = paste0(plot_file_name, "-", GSEA_result_list[[i]]@result[e, 1], script_version, ".png"), gsea_res_p,
                                device = "png", path = paste0(script_version, "/", save_folder, "/", "Plots", "/"),
                                width = 3500, height = 3000, units = "px")   
                
            }
                
        }
        message("The resulting GSEA plots were printed and saved into the folder: ", paste0(script_version, "/", save_folder, "/", "Plots"))
        
    
    
    if (GSEA_objects_to_list == TRUE) {
        
        # Assign the result as an object in the .GlobalEnv
        assign(GSEA_result_list_name, GSEA_result_list, envir = .GlobalEnv)
        
    } else {
        
        # Assign the result as an object in the .GlobalEnv
        assign(GSEA_result_object_name, GSEA_result, envir = .GlobalEnv)
        
    }
    
    
}


# Call the serial GSEA run on the ENTREZ ID reference list
run_serial_GSEA(DE_gene_list = GSEA_DEG_list_ord,
                   ref_gene_lst = pathways,
                   GSEA_csv_file_name = "GSEA_logP_entrez_result",
                   GSEA_result_object_name = "GSEA_logP_entrez_object",
                   GSEA_objects_to_list = TRUE,
                   GSEA_result_list_name = "GSEA_entrez_list",
                   plot_file_name = "GSEA_logP_entrez_plot",
                   plot_print_requested = FALSE,
                   save_folder = "GSEA_logP_entrez_result")


# Call the serial GSEA run on the gene SYMBOL refrence list
run_serial_GSEA(DE_gene_list = GSEA_DEG_symbols_list_ord,
                ref_gene_lst = pathways_symbols,
                GSEA_csv_file_name = "GSEA_logP_symbols_result",
                GSEA_result_object_name = "GSEA_logP_symbols_object",
                GSEA_objects_to_list = TRUE,
                GSEA_result_list_name = "GSEA_symbols_list",
                plot_file_name = "GSEA_logP_symbols_plot",
                plot_print_requested = FALSE,
                save_folder = "GSEA_logP_symbols_result")


# NOTE: compare the output results just in case. However when I did the comparison the results were the same, despite the
# difference in input lengths between ENTREZ IDs and gene SYMBOLS. That is probably because the missing gene SYMBOL - ENTREZ IDs
# are the ones covering the obscure pseudo genes, which are not well represented in the reference list


# Remove unnecessary variables
rm(list = pathways)
rm(list = pathways_symbols)

#################################################################   End Section  ##################################################################






#################################################################        GSEA result plotting         ##################################################################
                                                                #   (Gene Set Enrichment Analysis)    #




## Test plots for now, to  see which one should we go with (I chose the 4 best ones, I think)


# Dot plot
enrichplot::dotplot(GSEA_entrez_list[[16]])


# UpSet plot
enrichplot::upsetplot(GSEA_entrez_list[[16]])


# Tree plot
hallmark_set <- enrichplot::pairwise_termsim(GSEA_entrez_list[[16]])
enrichplot::treeplot(hallmark_set)


# Ridge plot
enrichplot::ridgeplot(GSEA_entrez_list[[16]])


#################################################################   End Section  ##################################################################






#################################################################                  SPIA                    ##################################################################
                                                                #   (Signaling Pathway Impact Analysis)    #




## Data preparation
## NOTE: SPIA takes two vectors as an input. One has to contain the DEG log2FC values (all the genes above a user defined p_val_adj) 
## as a named vector (names = ENTREZ IDs, values = avg_log2FC) (significant genes). The second input vector contains the ENTREZ IDs
## of all genes detected during sequencing (background genes).




## Preparing the significant genes vector with a p_val_adj <= 0.05 cutoff


# First subset the DEG list by the selected p_val_adj
significant_DEG_df <- dplyr::filter(DEG_cluster_0v1_clean, p_val_adj <= 0.05)


# NOTE: it is very important to remove duplicated ENTREZ IDs otherwise SPIA will throw an error:
# de must be a vector of log2 fold changes. The names of de should be included in the refference array!

summary(duplicated(significant_DEG_df$entrezIDs))
significant_DEG_df[grep(significant_DEG_df$entrezIDs[duplicated(significant_DEG_df$entrezIDs)], significant_DEG_df$entrezIDs), ]


# Remove/rename any potential duplicates
significant_DEG_df["ENSG00000205571", ]$entrezIDs <- "6607"


# Assign the avg_log2FC values to a new vector and name the vector with the ENTREZ IDs (this will be out SPIA input)
significant_DEG_lst <- significant_DEG_df$avg_log2FC
#significant_DEG_lst <- round(x = significant_DEG_lst, digits = 2)
names(significant_DEG_lst) <- as.vector(significant_DEG_df$entrezIDs)


# Preparing the background genes vector
background_gene_lst <- unifiedEntrezIDs$entrezgene_id




## Running SPIA
SPIA_result <- SPIA::spia(de = significant_DEG_lst, all = background_gene_lst, organism = "hsa",
                          nB = 5000, plots = TRUE, verbose = TRUE, combine = "fisher")


# Save the result of the SPIA run
write_csv(x = SPIA_result, file = paste0(script_version, "/", "SPIA_result_table", script_version, ".csv"))


# Plot the SPIA result
spia_p <- SPIA::plotP(x = SPIA_result, threshold = 0.05)

























#################################################################   End Section  ##################################################################






#################################################################   Overlaying interesting genes    ##################################################################
                                                                #           (DE or not)             #



# Surface markers (split form the main file for better customization)
surf_mark_df <- data.frame(EpCAM = sct_CTC.obj[["SCT"]]@data[c("ENSG00000119888"), ],
                           PSMA = sct_CTC.obj[["SCT"]]@data[c("ENSG00000086205"), ],
                           Cluster = sct_CTC.obj$seurat_clusters,
                           RespGroup = sct_CTC.obj$RespGroup,
                           Cycles = sct_CTC.obj$RealCycle)

surf_mark_df_long <- pivot_longer(data = surf_mark_df, cols = c(EpCAM, PSMA), names_to = "Marker", values_to = "ExprLvl")

# EpCAM responder
EpCAM_expr_over_cyc_resp_p <- ggplot(data = dplyr::filter(surf_mark_df_long, Marker == "EpCAM", RespGroup == "Responder"),
                                aes(x = factor(Cycles), y = ExprLvl, fill = factor(Cycles))) +
    #geom_violin(adjust = 1, trim = TRUE, scale = "width") +
    geom_boxplot(show.legend = FALSE, position = position_dodge(width = 1, preserve = "single")) +
    geom_dotplot(show.legend = FALSE, binaxis = "y", stackdir = "center", binwidth = 0.05, alpha = 0.5, position = position_dodge(width = 1, preserve = "total"),
                 dotsize = 2) +
    scale_fill_tron() +
    #stat_compare_means(label.y = 10, label.x = 1.2, method = "kruskal.test") +
    geom_pwc(method = "wilcox_test", label = "p = {p}", ref.group = "0") + #this geom is for pairwise comparison and is amazing!
    labs(x = expression("Treatment cycles"), y = expression("Log"[e] * " expression")) +
    ggtitle("EpCAM expression level changes in the responder group over treatment cycles") +
    #geom_jitter() +
    #facet_wrap(vars(Marker)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
print(EpCAM_expr_over_cyc_resp_p)

ggsave(filename = paste0("EpCAM_expression_over_cycles_Responders", script_version, ".png"), EpCAM_expr_over_cyc_resp_p,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 3500, height = 3000, units = "px")
rm(EpCAM_expr_over_cyc_resp_p)


# EpCAM nonresponder
EpCAM_expr_over_cyc_nresp_p <- ggplot(data = dplyr::filter(surf_mark_df_long, Marker == "EpCAM", RespGroup == "Nonresponder"),
                                      aes(x = factor(Cycles), y = ExprLvl, fill = factor(Cycles))) +
    #geom_violin(adjust = 1, trim = TRUE, scale = "width") +
    geom_boxplot(show.legend = FALSE, position = position_dodge(width = 1, preserve = "single")) +
    geom_dotplot(show.legend = FALSE, binaxis = "y", stackdir = "center", binwidth = 0.05, alpha = 0.5, position = position_dodge(width = 1, preserve = "total"),
                 dotsize = 2) +
    scale_fill_tron() +
    #stat_compare_means(label.y = 10, label.x = 1.2, method = "kruskal.test") +
    geom_pwc(method = "wilcox_test", label = "p = {p}", ref.group = "0") + #this geom is for pairwise comparison and is amazing!
    labs(x = expression("Treatment cycles"), y = expression("Log"[e] * " expression")) +
    ggtitle("EpCAM expression level changes in the nonresponder group over treatment cycles") +
    #geom_jitter() +
    #facet_wrap(vars(Marker)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
print(EpCAM_expr_over_cyc_nresp_p)

ggsave(filename = paste0("EpCAM_expression_over_cycles_Nonresponders", script_version, ".png"), EpCAM_expr_over_cyc_nresp_p,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 3500, height = 3000, units = "px")
rm(EpCAM_expr_over_cyc_nresp_p)


# PSMA responder
PSMA_expr_over_cyc_resp_p <- ggplot(data = dplyr::filter(surf_mark_df_long, Marker == "PSMA", RespGroup == "Responder"),
                                    aes(x = factor(Cycles), y = ExprLvl, fill = factor(Cycles))) +
    #geom_violin(adjust = 1, trim = TRUE, scale = "width") +
    geom_boxplot(show.legend = FALSE, position = position_dodge(width = 1, preserve = "single")) +
    geom_dotplot(show.legend = FALSE, binaxis = "y", stackdir = "center", binwidth = 0.05, alpha = 0.5, position = position_dodge(width = 1, preserve = "total"),
                 dotsize = 2) +
    scale_fill_tron() +
    #stat_compare_means(label.y = 10, label.x = 1.2, method = "kruskal.test") +
    geom_pwc(method = "wilcox_test", label = "p = {p}", ref.group = "0") + #this geom is for pairwise comparison and is amazing!
    labs(x = expression("Treatment cycles"), y = expression("Log"[e] * " expression")) +
    ggtitle("PSMA expression level changes in the responder group over treatment cycles") +
    #geom_jitter() +
    #facet_wrap(vars(Marker)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
print(PSMA_expr_over_cyc_resp_p)

ggsave(filename = paste0("PSMA_expression_over_cycles_Responders", script_version, ".png"), PSMA_expr_over_cyc_resp_p,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 3500, height = 3000, units = "px")
rm(PSMA_expr_over_cyc_resp_p)


# PSMA nonresponder
PSMA_expr_over_cyc_nresp_p <- ggplot(data = dplyr::filter(surf_mark_df_long, Marker == "PSMA", RespGroup == "Nonresponder"),
                                    aes(x = factor(Cycles), y = ExprLvl, fill = factor(Cycles))) +
    #geom_violin(adjust = 1, trim = TRUE, scale = "width") +
    geom_boxplot(show.legend = FALSE, position = position_dodge(width = 1, preserve = "single")) +
    geom_dotplot(show.legend = FALSE, binaxis = "y", stackdir = "center", binwidth = 0.05, alpha = 0.5, position = position_dodge(width = 1, preserve = "total"),
                 dotsize = 2) +
    scale_fill_tron() +
    #stat_compare_means(label.y = 10, label.x = 1.2, method = "kruskal.test") +
    geom_pwc(method = "wilcox_test", label = "p = {p}", ref.group = "0") + #this geom is for pairwise comparison and is amazing!
    labs(x = expression("Treatment cycles"), y = expression("Log"[e] * " expression")) +
    ggtitle("PSMA expression level changes in the nonresponder group over treatment cycles") +
    #geom_jitter() +
    #facet_wrap(vars(Marker)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
print(PSMA_expr_over_cyc_nresp_p)

ggsave(filename = paste0("PSMA_expression_over_cycles_Nonresponders", script_version, ".png"), PSMA_expr_over_cyc_nresp_p,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 3500, height = 3000, units = "px")
rm(PSMA_expr_over_cyc_nresp_p)




## Plot the gene expression levels of various markers in the clusters and selected groups


# Plot the following general markers (interesting to me)
# FOLH1, EPCAM ,SSTR2 ,TP53 ,GAPDH
Seurat::FeaturePlot(sct_CTC.obj,
            reduction = "umap",
            features = c("ENSG00000086205", "ENSG00000119888", "ENSG00000180616", "ENSG00000141510", "ENSG00000111640"),
            order = TRUE,
            label = FALSE,
            repel = FALSE)


Seurat::DotPlot(sct_CTC.obj,
        features = c("ENSG00000086205", "ENSG00000119888", "ENSG00000180616", "ENSG00000141510"))


general_markers_p <- Seurat::RidgePlot(sct_CTC.obj,
                                features = c("ENSG00000086205", #FOLH1
                                                "ENSG00000119888", #EPCAM
                                                "ENSG00000180616", #SSTR2
                                                "ENSG00000141510", #TP53
                                                "ENSG00000111640"),  #GAPDH
                                assay = "SCT",
                                cols = cluster_colors[c(6, 7)],
                                same.y.lims = TRUE,
                                log = TRUE,
                                combine = FALSE)


# Adjust the resulting RidgePlot
# NOTE: to be able to modify the plot it has to be split into list elements, so each change
# has to be applied to each list elemnt using a for loop
general_titles <- c("FOLH1", "EPCAM", "SSTR2", "TP53", "GAPDH")
for (i in seq_len(length(general_markers_p))) {
    general_markers_p[[i]] <- general_markers_p[[i]] + ggtitle(general_titles[i])
    general_markers_p[[i]] <- general_markers_p[[i]] + theme(legend.position = "none")
    general_markers_p[[i]] <- general_markers_p[[i]] + labs(y = expression("Clusters"))
    general_markers_p[[i]] <- general_markers_p[[i]] + theme(plot.title = element_text(hjust = 0.5))
    #general_markers_p[[i]] <- general_markers_p[[i]] + theme(axis.title.y = element_text(hjust = 0.5))
    #general_markers_p[[i]] <- general_markers_p[[i]] + theme(axis.title.x = element_text(hjust = 0.5))
}
general_markers_p <- cowplot::plot_grid(plotlist = general_markers_p, ncol = 2, byrow = FALSE)
general_markers_p <- general_markers_p + patchwork::plot_annotation(title = "Interesting marker gene expression",
    theme = theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold")))
print(general_markers_p)

cowplot::save_plot(filename = paste0("RidgePlot-General_markers", script_version, ".png"), plot = general_markers_p,
            ncol = 2, nrow = 3, base_height = 3.75, base_asp = 1.618, device = "png", path = paste0(script_version, "/", "Plots/"))
rm(general_markers_p, general_titles)



# Plot the following markers from the Kratochwil paper
# TP53, CHEK2, BRCA1, BRCA2, MSH2, MSH6, NBN, PMS1
kratochwil_markers_p <- Seurat::RidgePlot(sct_CTC.obj,
                                        features = c("ENSG00000141510",
                                                     "ENSG00000183765",
                                                     "ENSG00000012048",
                                                     "ENSG00000139618",
                                                     "ENSG00000095002",
                                                     "ENSG00000116062",
                                                     "ENSG00000104320",
                                                     "ENSG00000064933"),
                                        assay = "SCT",
                                        cols = cluster_colors[c(6, 7)],
                                        same.y.lims = TRUE,
                                        log = TRUE,
                                        combine = FALSE)


# Adjust the resulting RidgePlot
# NOTE: to be able to modify the plot it has to be split into list elements, so each change
# has to be applied to each list elemnt using a for loop
general_titles <- c("TP53", "CHEK2", "BRCA1", "BRCA2", "MSH2", "MSH6", "NBN", "PMS1")
for (i in seq_len(length(kratochwil_markers_p))) {
    kratochwil_markers_p[[i]] <- kratochwil_markers_p[[i]] + ggtitle(general_titles[i])
    kratochwil_markers_p[[i]] <- kratochwil_markers_p[[i]] + theme(legend.position = "none")
    kratochwil_markers_p[[i]] <- kratochwil_markers_p[[i]] + labs(y = expression("Clusters"))
    kratochwil_markers_p[[i]] <- kratochwil_markers_p[[i]] + theme(plot.title = element_text(hjust = 0.5))
    #kratochwil_markers_p[[i]] <- kratochwil_markers_p[[i]] + theme(axis.title.y = element_text(hjust = 0.5))
    #kratochwil_markers_p[[i]] <- kratochwil_markers_p[[i]] + theme(axis.title.x = element_text(hjust = 0.5))
}
kratochwil_markers_p <- cowplot::plot_grid(plotlist = kratochwil_markers_p, ncol = 2, byrow = FALSE)
kratochwil_markers_p <- kratochwil_markers_p + patchwork::plot_annotation(title = "Kratochwil et al. marker gene expression",
                                                                    theme = theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold")))
print(kratochwil_markers_p)

cowplot::save_plot(filename = paste0("RidgePlot-Kratochwil_markers", script_version, ".png"), plot = kratochwil_markers_p,
                   ncol = 2, nrow = 3, base_height = 3.75, base_asp = 1.618, device = "png", path = paste0(script_version, "/", "Plots/"))
rm(kratochwil_markers_p, general_titles)







# EMT related genes
# Source: Transcriptomic profiling of single circulating tumor cells provides insight into human metastatic gastric
# cancer
# ZEB2, MEF2D, NFKB1A, GATA1, SERPINE1, ZEB1, TWIST1, TWIST2, SNAI1, SNAI2, TGFB1, TGFB2
EMT_markers_p <- Seurat::RidgePlot(sct_CTC.obj,
            features = c("ENSG00000169554", #ZEB2
                        "ENSG00000116604", #MEF2D,
                        "ENSG00000109320", #NFKB1A
                        "ENSG00000179348", #GATA1
                        "ENSG00000106366", #SERPINE1
                        "ENSG00000092969"),  #TGFB2
                        assay = "SCT",
                        cols = cluster_colors[c(6, 7)],
                        same.y.lims = TRUE,
                        log = TRUE,
                        combine = FALSE)


# Adjust the resulting RidgePlot
# NOTE: to be able to modify the plot it has to be split into list elements, so each change
# has to be applied to each list element using a for loop
gene_names <- c("ZEB2", "MEF2D", "NFKB1A", "GATA1", "SERPINE1", "TGFB2")
for (i in seq_len(length(EMT_markers_p))) {
    EMT_markers_p[[i]] <- EMT_markers_p[[i]] + ggtitle(gene_names[i])
    EMT_markers_p[[i]] <- EMT_markers_p[[i]] + theme(legend.position = "none")
    EMT_markers_p[[i]] <- EMT_markers_p[[i]] + labs(y = expression("Clusters"))
    EMT_markers_p[[i]] <- EMT_markers_p[[i]] + theme(plot.title = element_text(hjust = 0.5))
    #EMT_markers_p[[i]] <- EMT_markers_p[[i]] + theme(axis.title.y = element_text(hjust = 0.5))
    #EMT_markers_p[[i]] <- EMT_markers_p[[i]] + theme(axis.title.x = element_text(hjust = 0.5))
}
EMT_markers_p <- cowplot::plot_grid(plotlist = EMT_markers_p, ncol = 3, byrow = FALSE)
EMT_markers_p <- EMT_markers_p + patchwork::plot_annotation(title = "EMT marker gene expression",
    theme = theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold")))
print(EMT_markers_p)

cowplot::save_plot(filename = paste0("RidgePlot-EMT_markers", script_version, ".png"), plot = EMT_markers_p,
            ncol = 2, nrow = 3, base_height = 3.75, base_asp = 1.618, device = "png", path = paste0(script_version, "/", "Plots/"))
rm(EMT_marker_p, gene_names)


# Epithelial markers
# EPCAM, KRT8, KRT10, KRT18, KRT19
epithelial_markers_p <- Seurat::RidgePlot(sct_CTC.obj,
        features = c("ENSG00000119888", #EPCAM
                    "ENSG00000170421", #KRT8
                    "ENSG00000186395", #KRT10
                    "ENSG00000111057", #KRT18
                    "ENSG00000171345"), #KRT19
                    assay = "SCT",
                    cols = cluster_colors[c(6, 7)],
                    same.y.lims = TRUE,
                    log = TRUE,
                    combine = FALSE)


gene_names <- c("EPCAM", "KRT8", "KRT10", "KRT18", "KRT19")
for (i in seq_len(length(epithelial_markers_p))) {
    epithelial_markers_p[[i]] <- epithelial_markers_p[[i]] + ggtitle(gene_names[i])
    epithelial_markers_p[[i]] <- epithelial_markers_p[[i]] + theme(legend.position = "none")
    epithelial_markers_p[[i]] <- epithelial_markers_p[[i]] + labs(y = expression("Clusters"))
    epithelial_markers_p[[i]] <- epithelial_markers_p[[i]] + theme(plot.title = element_text(hjust = 0.5))
    #epithelial_markers_p[[i]] <- epithelial_markers_p[[i]] + theme(axis.title.y = element_text(hjust = 0.5))
    #epithelial_markers_p[i]] <- epithelial_markers_p[[i]] + theme(axis.title.x = element_text(hjust = 0.5))
}
epithelial_markers_p <- cowplot::plot_grid(plotlist = epithelial_markers_p, ncol = 3, byrow = FALSE)
epithelial_markers_p <- epithelial_markers_p + patchwork::plot_annotation(title = "Epithelial marker gene expression",
    theme = theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold")))
print(epithelial_markers_p)

cowplot::save_plot(filename = paste0("RidgePlot-epithelial_markers", script_version, ".png"), plot = epithelial_markers_p,
            ncol = 2, nrow = 3, base_height = 3.75, base_asp = 1.618, device = "png", path = paste0(script_version, "/", "Plots/"))
rm(epithelial_markers_p, gene_names)


# Platelet related genes
# PF4, PPBP, ITGA2B, SPARC, ITGA2
platelet_markers_p <- Seurat::RidgePlot(sct_CTC.obj,
            features = c("ENSG00000163737",
                        "ENSG00000163736",
                        "ENSG00000005961",
                        "ENSG00000113140",
                        "ENSG00000164171"),
                        assay = "SCT",
                        cols = cluster_colors[c(6, 7)],
                        same.y.lims = TRUE,
                        log = TRUE,
                        combine = FALSE)


gene_names <- c("PF4", "PPBP", "ITGA2B", "SPARC", "ITGA2")
for (i in seq_len(length(platelet_markers_p))) {
    platelet_markers_p[[i]] <- platelet_markers_p[[i]] + ggtitle(gene_names[i])
    platelet_markers_p[[i]] <- platelet_markers_p[[i]] + theme(legend.position = "none")
    platelet_markers_p[[i]] <- platelet_markers_p[[i]] + labs(y = expression("Clusters"))
    platelet_markers_p[[i]] <- platelet_markers_p[[i]] + theme(plot.title = element_text(hjust = 0.5))
    #platelet_markers_p[[i]] <- platelet_markers_p[[i]] + theme(axis.title.y = element_text(hjust = 0.5))
    #platelet_markers_p[i]] <- platelet_markers_p[[i]] + theme(axis.title.x = element_text(hjust = 0.5))
}
platelet_markers_p <- cowplot::plot_grid(plotlist = platelet_markers_p, ncol = 3, byrow = FALSE)
platelet_markers_p <- platelet_markers_p + patchwork::plot_annotation(title = "Platelet binding related marker gene expression",
    theme = theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold")))
print(platelet_markers_p)

cowplot::save_plot(filename = paste0("RidgePlot-platelet_related_markers", script_version, ".png"), plot = platelet_markers_p,
            ncol = 2, nrow = 3, base_height = 3.75, base_asp = 1.618, device = "png", path = paste0(script_version, "/", "Plots/"))
rm(platelet_markers_p, gene_names)


# CSC markers
# ALDH1A1, TACSTD2, CD44, CD24, EZH2,
CSC_markers_p <- Seurat::RidgePlot(sct_CTC.obj,
            features = c("ENSG00000165092",
                        "ENSG00000184292",
                        "ENSG00000026508",
                        "ENSG00000272398",
                        "ENSG00000106462"),
                        assay = "SCT",
                        cols = cluster_colors[c(6, 7)],
                        same.y.lims = TRUE,
                        log = TRUE,
                        combine = FALSE)

gene_names <- c("ALDH1A1", "TACSTD2", "CD44", "CD24", "EZH2")
for (i in seq_len(length(CSC_markers_p))) {
    CSC_markers_p[[i]] <- CSC_markers_p[[i]] + ggtitle(gene_names[i])
    CSC_markers_p[[i]] <- CSC_markers_p[[i]] + theme(legend.position = "none")
    CSC_markers_p[[i]] <- CSC_markers_p[[i]] + labs(y = expression("Clusters"))
    CSC_markers_p[[i]] <- CSC_markers_p[[i]] + theme(plot.title = element_text(hjust = 0.5))
    #CSC_markers_p[[i]] <- CSC_markers_p[[i]] + theme(axis.title.y = element_text(hjust = 0.5))
    #CSC_markers_p[i]] <- CSC_markers_p[[i]] + theme(axis.title.x = element_text(hjust = 0.5))
}
CSC_markers_p <- cowplot::plot_grid(plotlist = CSC_markers_p, ncol = 3, byrow = FALSE)
CSC_markers_p <- CSC_markers_p + patchwork::plot_annotation(title = "CSC marker gene expression",
    theme = theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold")))
print(CSC_markers_p)

cowplot::save_plot(filename = paste0("RidgePlot-CSC_markers", script_version, ".png"), plot = CSC_markers_p,
            ncol = 2, nrow = 3, base_height = 3.75, base_asp = 1.618, device = "png", path = paste0(script_version, "/", "Plots/"))
rm(CSC_markers_p, gene_names)


# Proliferation markers
# UBE2C, CCNB1, TOP2A, PLK1, CCND1
proliferation_markers_p <- Seurat::RidgePlot(sct_CTC.obj,
            features = c("ENSG00000175063",
                        "ENSG00000134057",
                        "ENSG00000131747",
                        "ENSG00000166851",
                        "ENSG00000110092"),
                        assay = "SCT",
                        cols = cluster_colors[c(6, 7)],
                        same.y.lims = TRUE,
                        log = TRUE,
                        combine = FALSE)


gene_names <- c("UBE2C", "CCNB1", "TOP2A", "PLK1", "CCND1")
for (i in seq_len(length(proliferation_markers_p))) {
    proliferation_markers_p[[i]] <- proliferation_markers_p[[i]] + ggtitle(gene_names[i])
    proliferation_markers_p[[i]] <- proliferation_markers_p[[i]] + theme(legend.position = "none")
    proliferation_markers_p[[i]] <- proliferation_markers_p[[i]] + labs(y = expression("Clusters"))
    proliferation_markers_p[[i]] <- proliferation_markers_p[[i]] + theme(plot.title = element_text(hjust = 0.5))
    #proliferation_markers_p[[i]] <- proliferation_markers_p[[i]] + theme(axis.title.y = element_text(hjust = 0.5))
    #proliferation_markers_p[i]] <- proliferation_markers_p[[i]] + theme(axis.title.x = element_text(hjust = 0.5))
}
proliferation_markers_p <- cowplot::plot_grid(plotlist = proliferation_markers_p, ncol = 3, byrow = FALSE)
proliferation_markers_p <- proliferation_markers_p + patchwork::plot_annotation(title = "Proliferation related marker gene expression",
    theme = theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold")))
print(proliferation_markers_p)

cowplot::save_plot(filename = paste0("RidgePlot-Proliferation_markers", script_version, ".png"), plot = proliferation_markers_p,
            ncol = 2, nrow = 3, base_height = 3.75, base_asp = 1.618, device = "png", path = paste0(script_version, "/", "Plots/"))
rm(proliferation_markers_p, gene_names)


Seurat::FeaturePlot(sct_CTC.obj,
                    features = c("ENSG00000119888",
                    "ENSG00000086205"))

#################################################################   End Section  ##################################################################






#################################################################   RNA velocity analysis  ##################################################################
                                                                #                          #


# NOTE: the velocity analysis will be done in python, as the available R package cannot be installed for the time being


#Exporting data for RNA velocity analysis in python with the scVelo toolkit
#source: https://smorabit.github.io/tutorials/8_velocyto/

#export the necessary data to rebuild a seurat object like object (anndata object)
#in python for the analysis
export_Seurat_for_scVelo = function(seurat_object, gene_IDs_names_df, ENSEMBL_ID_col, gene_names_col) {
    ##Create an output folder
    if (dir.exists(paths = paste0(script_version, "/", "Velocity data")) == TRUE) {
        #do nothing
        message("The Velocity data folder already exists, no new folder will be created.")
    } else {
        dir.create(path = paste0(script_version, "/", "Velocity data"), recursive = TRUE)
        message("A new folder under the name: Velocity data was created.")
    }
    
    
    ##Save the seurat metadata table
    
    #defining the export object which will be modified
    scVelo_sct_CTC.obj <- seurat_object
    
    #saving the metadata table    
    scVelo_sct_CTC.obj$barcode <- colnames(seurat_object)
    scVelo_sct_CTC.obj$UMAP_1 <- seurat_object@reductions$umap@cell.embeddings[, 1]
    scVelo_sct_CTC.obj$UMAP_2 <- seurat_object@reductions$umap@cell.embeddings[, 2]
    
    #save the export table 
    write.csv(scVelo_sct_CTC.obj@meta.data, file = paste0(script_version, "/", "Velocity data", "/", "CTC_metadata_for_scVelo", script_version, ".csv"),
              quote = FALSE, row.names = FALSE)
    
    ##Save the RNA count matrix
    ##NOTE: the saved data will be the already processed sctransformed data, used for the clustering and DE, GSEA analysis
    
    #I will simply save the full RNA count matrix
    gene_count_matrix_IDs <- Seurat::GetAssayData(seurat_object, assay = "SCT", layer = "counts")
    Matrix::writeMM(gene_count_matrix_IDs, file = paste0(script_version, "/", "Velocity data", "/", "CTC_GeneCountMatrix_IDs_for_scVelo", script_version, ".mtx"))
    
    #as I'm unsure if the ENSEMBL IDs will be usable for the analysis I will save a trimmed matrix containing the entries with identified gene neames
    #missing_name_IDs <- unif_gene_names$ensembl_gene_id[is.na(unif_gene_names$external_gene_name[rownames(GetAssayData(seurat_object, assay = "SCT", slot = "counts")) %in% unif_gene_names$ensembl_gene_id])]
    missing_name_IDs <- gene_IDs_names_df[, ENSEMBL_ID_col][is.na(gene_IDs_names_df[, gene_names_col][rownames(gene_count_matrix_IDs) %in% gene_IDs_names_df[, ENSEMBL_ID_col]])]
    gene_count_matrix_names <- gene_count_matrix_IDs[!rownames(gene_count_matrix_IDs) %in% missing_name_IDs, ]
    Matrix::writeMM(gene_count_matrix_names, file = paste0(script_version, "/", "Velocity data", "/", "CTC_GeneCountMatrix_gNames_for_scVelo", script_version, ".mtx"))
    
    
    ##Save the  dimensionality reduction matrix
    ##NOTE: for the clustering analysis I used PCA as a dimensionality reduction method so this will be used now as well
    write.csv(seurat_object@reductions$pca@cell.embeddings, file = paste0(script_version, "/", "Velocity data", "/", "PCA_data_for_scVelo", script_version, ".csv"),
              quote = FALSE, row.names = FALSE)
    
    
    ##Save the ENSEMBL IDs and gene names
    ##NOTE: just like with the count matrices, I'm unsure if ENSEMBL IDs will work, so I will save the gene names too
    
    #saving the full ENSEMBL ID list
    write.table(data.frame("Gene" = rownames(gene_count_matrix_IDs)), file = paste0(script_version, "/", "Velocity data", "/", "ENSEMBL_IDs_for_scVelo", script_version, ".csv"),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    #saving the trimmed gene name list
    gene_names_to_save <- gene_IDs_names_df[, gene_names_col][gene_IDs_names_df[, ENSEMBL_ID_col] %in% rownames(gene_count_matrix_names)]
    write.table(data.frame("Gene" = gene_names_to_save), file = paste0(script_version, "/", "Velocity data", "/", "Gene_names_for_scVelo", script_version, ".csv"),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    ##Letting the user know that the process was successful
    message("\n",
            "Exporting the required data for scVelo was successful!",
            "\n",
            "\n", "The following new files were created in the Velocity data folder:", 
            "\n", paste0("CTC_metadata_for_scVelo", script_version, ".csv"),
            "\n", paste0("CTC_GeneCountMatrix_IDs_for_scVelo", script_version, ".mtx"),
            "\n", paste0("CTC_GeneCountMatrix_gNames_for_scVelo", script_version, ".mtx"),
            "\n", paste0("PCA_data_for_scVelo", script_version, ".csv"),
            "\n", paste0("ENSEMBL_IDs_for_scVelo", script_version, ".csv"),
            "\n", paste0("Gene_names_for_scVelo", script_version, ".csv")
    )
    
    
    
}
export_Seurat_for_scVelo(sct_CTC.obj, unif_gene_names, "ensembl_gene_id", "external_gene_name")


#copying only the .bam files for the analysis
#NOTE:very resource intensive... (should paralellize it if time allows)
copy_designated_files = function(current_file_folder, destination_folder, filename = "out.bam$") {
    #list the files you wish to copy
    files_to_copy <- list.files(current_file_folder, pattern = filename, recursive = TRUE, full.names = TRUE)
    
    #initialize the progress bar
    progressBar <- progress_bar$new(format = "Copying  files - (:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                    total = length(files_to_copy),
                                    complete = "=",
                                    incomplete = "-",
                                    current = ">",
                                    clear = FALSE,
                                    width = 100)
    
    #this loop will be responsible for copiying the designated files
    for (e in files_to_copy) {
        #initiate the progress bar
        progressBar$tick()
        
        #copy the files
        file.copy(from = e, to = destination_folder)
    }
    
}
copy_designated_files(current_file_folder = bam_files,
                      destination_folder = just_bam,
                      filename = "out.bam$")


#These are the ENSEMBL IDs corresponding to the dynamic genes from the velocity analysis.
#however, these are all unidentified genes
velocityGenes <- c("ENSG00000228658", "ENSG00000236531", "ENSG00000276256", "ENSG00000275464",
                   "ENSG00000255240", "ENSG00000235685", "ENSG00000259647", "ENSG00000283696",
                   "ENSG00000267011", "ENSG00000263427", "ENSG00000267968", "ENSG00000271254",
                   "ENSG00000233538", "ENSG00000240291")

velocityGene_names <- list()
for(i in seq_along(velocityGenes)){
    velocityGene_names[[i]] <- DE_Cluster_0v1_geneNames$geneNames[rownames(DE_Cluster_0v1_geneNames) %in% velocityGenes[i]]
}
glimpse(velocityGene_names)

##NOTE: as the resulting genes are all unidentified, I cannot in good faith take the result seriously, therefore the velocity analsis
##will not be continued

#################################################################   End Section  ##################################################################






#################################################################   Trajectory analysis  ##################################################################
                                                                #      with Monocle 3    #




## Prepare the seurat Object fot the analysis by re-running the sctransform, PCA and find neighbors


# SCTransform
sct_CTC.obj_TA <- SCTransform(filt_CTC.obj, assay = "RNA_seq", method = "glmGamPoi")


# An important thing to check is which PCs are responsible for the most variation. This is important
# for the downstream steps. One way to check this is to plot the PCs with an ElbowPlot
sct_CTC.obj_TA <- RunPCA(sct_CTC.obj_TA, assay = "SCT", approx = FALSE, verbose = FALSE) #approx = FALSE to run with normal vsd method
PC_weight_p <- ElbowPlot(sct_CTC.obj_TA, ndims = 50, reduction = "pca")
PC_weight_p2 <- PC_weight_p +
    labs(x = "PCs") +
    ggtitle("PC weights") +
    theme(plot.title = element_text(hjust = 0.5))
print(PC_weight_p2)
rm(PC_weight_p, PC_weight_p2)


# Finding neighbours and clustering
sct_CTC.obj_TA <- FindNeighbors(sct_CTC.obj_TA, reduction = "pca", dims = 1:35)


# Running dimensionality reduction (trying both UMAP and tSNE :) )
sct_CTC.obj_TA <- RunUMAP(sct_CTC.obj_TA, dims = 1:35, reduction = "pca", verbose = FALSE)




## Convert the seurat object to a monocle CDS (Cell Data Set) and do the rest with this package

# Converting the seurat object
sct_CTC.cds <- SeuratWrappers::as.cell_data_set(sct_CTC.obj_TA)


# As suggested by the monocle3 package, running the cluster_cells on the new cds object
sct_CTC.cds <- monocle3::cluster_cells(sct_CTC.cds, resolution = 0.2)

# Checking how the clustering looks like (adjust the resolution if needed)
monocle3::plot_cells(cds = sct_CTC.cds, show_trajectory_graph = FALSE, cell_size = 1)

monocle3::plot_cells(cds = sct_CTC.cds, show_trajectory_graph = FALSE, cell_size = 1,
                        color_cells_by = "partition")

# Next we will run the learn graph function



























































































































ggplot(sct_CTC_meta_df,
       aes(x = Treatment, fill = seurat_clusters,)) +
    #geom_jitter(width = 0.5, height = 0.5) +
    geom_bar(color = "black", position = "dodge") +
    #geom_pwc() +
    #facet_wrap(vars(RealCycle)) +
    theme_classic()


ggplot(sct_CTC_meta_df,
       aes(x = Treatment, y = after_stat(count), fill = seurat_clusters,)) +
    #geom_jitter(width = 0.5, height = 0.5) +
    geom_col(color = "black", position = "dodge") +
    #geom_pwc() +
    #facet_wrap(vars(RealCycle)) +
    theme_classic()






treatment_summary <- data.frame(Treatment_type = )












c25 <- c(
    "dodgerblue2", "#E31A1C", # red
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "black", "gold1",
    "skyblue2", "#FB9A99", # lt pink
    "palegreen2",
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown"
)





ggplot(sct_CTC_meta_df,
       aes(x = seurat_clusters, fill = patients,)) +
    #geom_jitter(width = 0.5, height = 0.5) +
    geom_bar(color = "black", position = "dodge") +
    scale_fill_manual(values = c25) +
    facet_wrap(vars(RealCycle)) +
    theme_classic()

ggplot(sct_CTC_meta_df,
       aes(x = RealCycle, fill = patients,)) +
    #geom_jitter(width = 0.5, height = 0.5) +
    geom_bar(color = "black", position = "dodge") +
    scale_fill_manual(values = c25) +
    #facet_wrap(vars(seurat_clusters)) +
    theme_classic()




multi_cyc_pat <- dplyr::filter(sct_CTC_meta_df, patients == "ALM76C6003" | patients == "ALM78C7005" | patients == "LM67C6011" |
                  patients == "ALM88C3002" | patients == "LM68C6005")

cycle_labs <- c("Cycle 1", "Cycle 2", "Cycle 3")

ggplot(multi_cyc_pat,
       aes(x = seurat_clusters, fill = patients,)) +
    #geom_jitter(width = 0.5, height = 0.5) +
    geom_bar(color = "black", position = "dodge") +
    scale_fill_startrek() +
    #scale_fill_manual(values = c25) +
    facet_wrap(vars(Treatment_cycle), strip.position = "bottom", labeller = "label_both") +
    labs(fill = "Patients", x = expression("Clusters"), y = expression("Cell numbers")) +
    ggtitle("CTC clusters - cell distribution based on the treatment cycles") +
    #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank())





run_GSEA(DE_gene_list = GSEA_DE_logP_ranked_lst_ord,
         verbose = TRUE,
         min_size = 5,
         max_size = 500,
         p_val_cutoff = 0.05,
         p_val_adjust_method = "BH",
         GSEA_method = "fgsea",
         ref_gene_lst = c2.cgp.v2023.1.Hs.symbols.gmt,
         GSEA_csv_file_name = "test_file",
         GSEA_result_object_name = "test_object",
         plot_title = "test_plot",
         plot_file_name = "test_plot",
         save_folder = "GSEA_results")



for (e in seq_along(reference_lst)) {
    run_GSEA(DE_gene_list = GSEA_DE_ranked_lst_ord,
             ref_gene_lst = get(reference_lst[e]),
             GSEA_csv_file_name = paste0(pathway_names[e], "_", "GSEA", "_", "result"),
             GSEA_result_object_name = paste0(pathway_names[e], "_", "GSEA", "_", "object"),
             plot_title = paste0("GSEA", pathway_names[e], "set"),
             plot_file_name = paste0(pathway_names[e], "_", "GSEA", "_", "plot"))
    
}




#Running GSEA using the clusterProfiler package
#the function uses permutation testing to adjust the p-values. as this is based on a random number, for the sake of
#reproducibility, one should set the seed
set.seed(1234)


#Run GSEA on the hallmark dataset and save the result
hallmark_gsea_res <- GSEA(geneList = GSEA_DE_ranked_lst_ord,
                          minGSSize = 5,
                          maxGSSize = 500,
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          TERM2GENE = h.all.v2023.1.Hs.symbols.gmt,
                          by = "fgsea")

hallmark_gsea_res_df <- as.data.frame(hallmark_gsea_res@result)
write.csv(hallmark_gsea_res_df, file = paste0(script_version, "/", "GSEA_results/", "hallmark_set_GSEA", script_version, ".csv"))
rm(hallmark_gsea_res_df)



#this function will visualize the GSEA result and saves it as a .png
print_and_save_gseaplot = function(gsea_result, plot_title, plot_file_name) {
    for (e in 1:nrow(gsea_result@result)) {
        #call the gseaplot function to plot the gsea result
        gsea_res_p <- gseaplot(x = gsea_result,
                               geneSetID = gsea_result@result[e,1],
                               title = paste0(plot_title, "- ", gsea_result@result[e,1]))
        #visualize the resulting plot
        print(gsea_res_p)
        
        #saving the resulting plot as a .png
        ggsave(filename = paste0(plot_file_name, "-", gsea_result@result[e,1], script_version, ".png"), gsea_res_p,
               device = "png", path = paste0(script_version, "/", "GSEA_results", "/", "Plots", "/"),
               width = 3500, height = 3000, units = "px")
        
    }
    
}
print_and_save_gseaplot(hallmark_gsea_res, "GSEA hallmark set", "GSEA hallmark set")



#Run GSEA on the hallmark dataset and save the result
hallmark_gsea_res <- GSEA(geneList = GSEA_DE_ranked_lst_ord,
                          minGSSize = 5,
                          maxGSSize = 500,
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          TERM2GENE = h.all.v2023.1.Hs.symbols.gmt,
                          by = "fgsea")

hallmark_gsea_res_df <- as.data.frame(hallmark_gsea_res@result)
write.csv(hallmark_gsea_res_df, file = paste0(script_version, "/", "GSEA_results/", "hallmark_set_GSEA", script_version, ".csv"))
rm(hallmark_gsea_res_df)


#this function will visualize the GSEA result and saves it as a .png
print_and_save_gseaplot = function(gsea_result, plot_title, plot_file_name) {
    for (e in 1:nrow(gsea_result@result)) {
        #call the gseaplot function to plot the gsea result
        gsea_res_p <- gseaplot(x = gsea_result,
                               geneSetID = gsea_result@result[e,1],
                               title = paste0(plot_title, "- ", gsea_result@result[e,1]))
        #visualize the resulting plot
        print(gsea_res_p)
        
        #saving the resulting plot as a .png
        ggsave(filename = paste0(plot_file_name, "-", gsea_result@result[e,1], script_version, ".png"), gsea_res_p,
               device = "png", path = paste0(script_version, "/", "GSEA_results", "/", "Plots", "/"),
               width = 3500, height = 3000, units = "px")
        
    }
    
}

print_and_save_gseaplot(hallmark_gsea_res, "GSEA hallmark set", "GSEA hallmark set")































for (e in 1:nrow(hallmark_gsea_res@result)) {
    h_gsea_res_p <- gseaplot(x = hallmark_gsea_res,
                             geneSetID = hallmark_gsea_res@result[e,1],
                             title = paste0("GSEA Hallmark set", "- ", hallmark_gsea_res@result[e,1]))
    print(h_gsea_res_p)
    ggsave(filename = paste0("GSEA hallmark set", "-", hallmark_gsea_res@result[e,1], script_version, ".png"), h_gsea_res_p,
           device = "png", path = paste0(script_version, "/", "GSEA_results", "/", "Plots", "/"),
           width = 3500, height = 3000, units = "px")
    
}
rm(h_gsea_res_p)





    


















    
    








test_GSEA <- fgsea(pathways = test_gmt,
                   stats = DE_ranked_list_ord[1:1000],
                   scoreType = "pos",
                   minSize = 5,
                   maxSize = 500,
                   nproc = 2)

dplyr::arrange(test_GSEA, padj)







ranked_DE_lst <- DE_Cluster_0v1_geneNames
head(ranked_DE_lst)
ranked_DE_lst <- dplyr::arrange(ranked_DE_lst, Ranking)
head(ranked_DE_lst)




ggplot(sct_CTC_meta_df,
       aes(x = RealCycle, fill = seurat_clusters,)) +
    geom_bar(color = "black") +
    geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, position = position_dodge(width = 0.9)) +
    #scale_x_discrete(labels = c("Double positive","EpCAM positive","PSMA positive")) +
    scale_fill_tron(labels = c("Cluster 0", "Cluster 1")) +
    labs(fill = "CTC cluster", x = expression("Treatment cycle"), y = expression("Cell number")) +
    facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
          strip.text = element_blank())


ggplot(cell_numbers_summary_df,
       aes(x = Cycles, y = as.integer(Cell_percent_c0), fill = cluster_colors[6])) +
    geom_col(color = "black") +
    geom_text(aes(label = as.integer(Cell_percent_c0)), vjust = -0.5, position = position_dodge(width = 0.9)) +
    scale_y_continuous(limits = c(0,100), breaks = c(0, 20, 40, 60, 80, 100)) +
    #scale_x_discrete(labels = c("Double positive","EpCAM positive","PSMA positive")) +
    scale_fill_manual(values = cluster_colors[6], label = "Cluster 0") +
    labs(fill = "CTC cluster", x = expression("Treatment cycle"), y = expression("Cell number percentages")) +
    #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
          strip.text = element_blank())

ggplot(cell_numbers_summary_df,
       aes(x = Cycles, y = as.integer(Cell_percent_c1), fill = cluster_colors[7])) +
    geom_col(color = "black") +
    geom_text(aes(label = as.integer(Cell_percent_c1)), vjust = -0.5, position = position_dodge(width = 0.9)) +
    scale_y_continuous(limits = c(0,100), breaks = c(0, 20, 40, 60, 80, 100)) +
    #scale_x_discrete(labels = c("Double positive","EpCAM positive","PSMA positive")) +
    scale_fill_manual(values = cluster_colors[7], label = "Cluster 0") +
    labs(fill = "CTC cluster", x = expression("Treatment cycle"), y = expression("Cell number percentages")) +
    #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
          strip.text = element_blank())






































##NetGSA


#NetGSA uses a count matrix for the gene set analysis, so I extracted the sct-transformed count matrix from the
#seurat object and substracetd the rows present in the DE gene list
#extracting the transformed count matrix from the seurat object for the NetGSA analysis
sct_count_matrix <- GetAssayData(sct_CTC.obj, layer = "data")

#adjusting the count matrix to contain the DE_lst genes
net_matrix <- sct_count_matrix[rownames(sct_count_matrix) %in% rownames(DEG_cluster_0v1_clean), ]
rm(sct_count_matrix)

net_matrix_2 <- net_matrix
rownames(net_matrix_2) <- DEG_cluster_0v1_clean$entrezIDs

net_matrix_2@Dimnames[[1]] <- str_c("ENTREZID:", net_matrix_2@Dimnames[[1]])
head(rownames(net_matrix_2))




#try again with this setup...
net_matrix <- matrix(ncol = ncol(sct_count_matrix), nrow = nrow(DE_Cluster_0v1_geneNames))
cut_sct_m <- sct_count_matrix[rownames(sct_count_matrix) %in% rownames(DE_Cluster_0v1_geneNames), ]

rownames(net_matrix) <- rownames(cut_sct_m)
colnames(net_matrix) <- colnames(sct_count_matrix)

for(e in 1:ncol(net_matrix)) {
    net_matrix[, e] <- rep(colnames(net_matrix)[e], nrow(net_matrix))
}



progressBar <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                total = nrow(net_matrix),
                                complete = "=",
                                incomplete = "-",
                                current = ">",
                                clear = FALSE,
                                width = 100)


net_matrix2 <- matrix(ncol = ncol(net_matrix), nrow = nrow(net_matrix))
for (e in 1:ncol(net_matrix)) {
    for (i in 1:nrow(net_matrix)) {
        progressBar$tick()
        net_matrix2[i, e] <- as.numeric(net_matrix[i, e]) 
    }
}

net_matrix2 <- GeneCountTable
net_matrix2$geneID <- NULL
net_matrix2 <- net_matrix2[rownames(net_matrix2) %in% rownames(DE_Cluster_0v1_geneNames), ]
net_matrix2 <- net_matrix2[, colnames(net_matrix2) %in% colnames(net_matrix)]

net_matrix2 <- data.matrix(net_matrix2[1:109])
net_mat <- Matrix::Matrix(net_matrix2, sparse = TRUE)
net_mat2 <- na.omit(net_mat[1:1000, ])

glmnet::glmnet(net_mat, net_group)


net_matrix2 <- glmnet::makeX(net_matrix2)
net_mat <- glmnet::makeX(net_matrix)
    




net_mat@Dimnames[[1]] <- str_c("ENSEMBL:", net_mat@Dimnames[[1]])
head(rownames(net_mat))








net_mat <- DEG_cluster_0v1_clean
net_mat$entrezIDs[742] <- "6607"
rownames(net_mat) <- net_mat$entrezIDs
rownames(net_mat) <- str_c("ENTREZID:", rownames(net_mat))
net_mat

net_mat <- glmnet::makeX(net_mat, na.impute = TRUE)


net_matrix_3 <- as.data.frame(net_matrix_2)
net_matrix_3 <- glmnet::makeX(net_matrix_3, na.impute = TRUE)

#setting up the grouping variable, in this case the cluster numbers
net_group <- as.numeric(sct_CTC.obj@meta.data$seurat_clusters)
table(net_group)


#getting the edges
edges <- netgsa::obtainEdgeList(genes = rownames(net_matrix_3), databases = c("kegg", "reactome"))
rownames(net_matrix_3)[1088] <- "ENTREZID:6607"





net_adjMat <- netgsa::prepareAdjMat(x = net_matrix_3[1:1000, ], group = net_group)











#converting the ENSEMBL IDs to ENTREZ IDs as the NetGSA algorythm has a problem working with ENSEMBL IDs
entrezID <- getBM(attributes = c("ensembl_gene_id","entrezgene_id", "external_gene_name"),
                  filters = "ensembl_gene_id",
                  values = rownames(DE_Cluster_0v1_geneNames),
                  mart = ensembl_connection)


#the conversion is not perfect as the ENSEMBL IDS contain IDs for pseudogenes, so the missing entrezIDs will be removed
entrezID_filt <- entrezID[!is.na(entrezID$entrezgene_id), ]


#subsetting the net_matrix and the DE_gene list based on the remaining IDs
net_matrix <- net_matrix[rownames(net_matrix) %in% entrezID_filt$ensembl_gene_id, ]
DE_Cluster_0v1_geneNames_filt <- DE_Cluster_0v1_geneNames[rownames(DE_Cluster_0v1_geneNames) %in% entrezID_filt$ensembl_gene_id, ]


#Re-defining the top upregulatd genes after the filtering
DE_Cluster_0v1_geneNames_filt$Significance <- ifelse(DE_Cluster_0v1_geneNames_filt$avg_log2FC > 0 & DE_Cluster_0v1_geneNames_filt$p_val_adj < 0.05, "Upregulated",
                                                ifelse(DE_Cluster_0v1_geneNames_filt$avg_log2FC < 0 & DE_Cluster_0v1_geneNames_filt$p_val_adj < 0.05, "Downregulated", "Not significant")) #labeling genes for coloring
upreg_filt <- DE_Cluster_0v1_geneNames_filt[DE_Cluster_0v1_geneNames_filt$Significance == "Upregulated", ] #subsetting the upregulated genes
topupreg_filt <- head(upreg_filt$geneName[order(upreg_filt$p_val_adj)], 10) #selecting the top upreg genes
downreg_filt <- DE_Cluster_0v1_geneNames_filt[DE_Cluster_0v1_geneNames_filt$Significance == "Downregulated", ] #subsetting the downregulated genes
topdownreg_filt <- head(downreg_filt$geneName[order(downreg_filt$p_val_adj)], 10) #selecting the top downreg genes


top50upreg_filt <- head(DE_Cluster_0v1_geneNames_filt$geneName[order(DE_Cluster_0v1_geneNames_filt$p_val_adj)], 50) #choosing the top 50 differentially expressed genes
DE_Cluster_0v1_geneNames_filt$DElabel <- ifelse(DE_Cluster_0v1_geneNames_filt$geneNames %in% topupreg_filt | DE_Cluster_0v1_geneNames_filt$geneNames %in% topdownreg_filt, DE_Cluster_0v1_geneNames_filt$geneNames, NA) #marking the top DEGs in the dataframe


#plotting the new DE_gene list
volcano_p <- ggplot(data = DE_Cluster_0v1_geneNames_filt, 
                    aes(x = avg_log2FC, y = -log10(p_val_adj), col = Significance, label = DElabel)) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_color_manual(values = cluster_colors[c(2, 5, 3)],
                       labels = c("Downregulated", "Not significant", "Upregulated")) +
    scale_x_continuous(limits = c(-10, 10), breaks = seq(-10, 10, 2)) +
    scale_y_continuous(limits = c(0, 14), breaks = seq(0, 14, 2)) +
    labs(color = "Expression status", x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
    ggtitle("Differentially expressed genes - filtered for the NetGSA analysis") +
    geom_text_repel(show.legend = FALSE ,max.overlaps = Inf) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
print(volcano_p)

ggsave(filename = paste0("NetGSA_filtered_DEGs_Volcano_plot_gb", script_version, ".png"), volcano_p,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 3500, height = 3000, units = "px")
rm(volcano_p, upreg, downreg, topupreg, topdownreg)








#adjusting the rwonames of the net_matrix to fit the requirements of the NetGSA (it has to contain what kid of ID the name is based on,
#in this case it is ENSEMBL ID so one must add ENSEMBL: in fornt of each ID)
rownames(net_matrix) <- str_c("ENSEMBL:", rownames(net_matrix))
head(rownames(net_matrix))

#Adding the sample names as colnames for the sake of easier identification
colnames(net_matrix) <- colnames(sct_CTC.obj)
head(colnames(net_matrix))


#setting up the frouping variable, in this case the cluster numbers
net_group <- as.numeric(sct_CTC.obj@meta.data$seurat_clusters)
table(net_group)


#getting the edges
edges <- netgsa::obtainEdgeList(genes = rownames(net_matrix), databases = c("kegg", "reactome"))

net_adjMat <- prepareAdjMat(x = net_matrix[1:1000, ], group = net_group)














upreg <- DE_Cluster_0v1_geneNames[-log10(DE_Cluster_0v1$p_val_adj) >= 1.3 & DE_Cluster_0v1$avg_log2FC >= 1,]
upreg <- upreg[is.na(upreg$geneNames) == FALSE, ]
downreg <- DE_Cluster_0v1_geneNames[-log10(DE_Cluster_0v1$p_val_adj) >= 1.3 & DE_Cluster_0v1$avg_log2FC <= -1,]
downreg <- downreg[is.na(downreg$geneNames) == FALSE, ]

write_csv(upreg, "Alpha_volcano_full_upreg_Cl0vs1.csv")
write_csv(downreg, "Alpha_volcano_full_downreg_Cl0vs1.csv")

























































































































