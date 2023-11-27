#!user/bin/Rscript



### SCS results - analysis




### Setting up, formatting and cleaning the relevant data


# Defining the required R packages for the full analysis
cran_packages <- c("tidyverse", "stringi", "BiocManager",
              "scales", "RCurl", "cowplot", "rebus", "ggsci",
              "progress", "metap", "doSNOW", "foreach", "scCustomize",
              "Matrix", "ggpubr", "R.utils", "devtools", "remotes")

bioconductor_packages <- c("Seurat", "glmGamPoi", "multtest", "biomaRt", "AnnotationDbi",
              "EnsDb.Hsapiens.v86", "EnhancedVolcano", "graphite", "netgsa",
              "org.Hs.eg.db", "fgsea", "clusterProfiler", "SPIA")

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


# Used library packages
packages <- c("tidyverse", "stringi", "BiocManager", "Seurat", "glmGamPoi", "biomaRt", "AnnotationDbi",
              "scales", "EnsDb.Hsapiens.v86", "RCurl", "cowplot", "rebus", "ggsci",
              "EnhancedVolcano", "progress", "metap", "doSNOW", "foreach", "scCustomize", "stringi",
              "scales", "ggsci", "graphite", "org.Hs.eg.db", "fgsea", "clusterProfiler",
              "scCustomize", "Matrix", "ggpubr", "SeuratWrappers", "presto", "multtest")
lapply(packages, library, character.only = TRUE)


# Setting wd + path and listing the read output tables
setwd("/Users/ramasz/Coding/Gabi_R")
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
read_RPG.tabs <- function (path = getwd(), filenames, test_set = FALSE) {
    RPG.tab_lst <- list.files(path = path, 
                              pattern = filenames, 
                              all.files = TRUE,
                              recursive = TRUE) #super neat function which goes into subfolders to find the pattern designated file
    
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
    progressBar <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
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
    
        
    samples <- str_remove(RPG.tab_lst, "^lane1")
    samples <- str_remove(samples, "_.*")
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


# Retrieving gene names using the ensembl IDs
listEnsembl() #lists the available biomart keynames, here I need "gene"
ensembl <- useEnsembl(biomart = "genes") #creates an ensembl object needed for the database connection
datasets <- listDatasets(ensembl) #lists the organism based datasest from which we need to chose one for the db connection
ensembl_connection <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl") #this connects the database and the dataset for the search query
attributes <- listAttributes(ensembl_connection) #needed for the query, describes the desired info like gene names
filters <- listFilters(ensembl_connection) # needed for the query, defined by upon which data we search for(?) like here the ensemblIds

gene_names <- getBM(attributes = c("ensembl_gene_id","external_gene_name"),
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


rm(gene_names, gene_names_2, clean_gene_names, missingIDs, missingNames, unif_missing)

# Re-ordering the GeneCountTable (GCT) by rows to match the re-ordered unif_gene_names df (by ensembl IDs)
# Removing the geneID column to allow ordering by sample (column) name and ordering the counts
#oGCT <- arrange(GeneCountTable, geneID)
#oGCT <- as.data.frame(oGCT)

#rownames(oGCT) <- unif_gene_names$ensembl_gene_id #optional, to track which row is which gene without messing with the ordering


# Loading the metadata.csv (only do this if you don't have the extended version yet)
metad <- read.csv("Index_table_reseq.csv", sep = ";")
glimpse(metad)
metad <- metad[, 2:3]
colnames(metad) <- c("sampleID", "RespGroup")
metad$sampleID <- str_replace_all(metad$sampleID, "\\(A\\)", "A")
patientID <- str_remove(metad$sampleID, "P" %R% one_or_more(ASCII_ALNUM))
metad <- mutate(metad, patientID = patientID)
metadata <- metad[order(metad$patientID), ]
rm(metad)
write_csv(metadata, "Reseq-sample_metadata.csv")


# Subsetting sample codes for additional entries to the metadata df
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
write_csv(metadata, "Re-seq_sample_extrended_metadata.csv")
# IMPORTANT NOTE: in this case the noted treatment cycles and the real cycles and treatments
# do not correspond, so the metadata has to be corrected accordingly. This was done separately
# and the corrected data will be loaded at the next section




















### Creating the seurat object and prepping for the data analysis (for this I need a count table and a prepped metadata dataframe)


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
    message("The metadta is loaded as the following object: meta")
}
read_metadata(filename = "Re-seq_sample_extrended_fixed_metadata_with_markers.csv", sep = ";")
glimpse(meta)
head(meta)
# NOTE:these steps seem to be necessary for the seurat object to utilize the metadata properly
#metadata <- as.data.frame(metadata)
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
save(CTC_reseq.obj, file = paste0(script_version, "/", "CTC_RNA_re-seq_full_seur_obj", script_version, ".RData"))


# Following the creation of the seurat object remove some not needed objects to save memory
rm(meta, GCT)



















### Quality control


# Load mitochondrial geneIDs and names (got it from ENSEMBL) for the QC
MT_genes <- read_csv("ENSEMBL_MT_genes.csv")
mito.genes <- rownames(CTC_reseq.obj)[rownames(CTC_reseq.obj) %in% MT_genes$`Gene stable ID`]
print(mito.genes)
rm(mito.genes)


# Alternatively one can load it from the scCustomize package
#ensembl_mito_id


# Add mitochondrial gene percentage and genes per nCount percentages (log10 tranformed for better visibility) to our
# seurat.obj
CTC_reseq.obj[["mt.percent"]] <- PercentageFeatureSet(CTC_reseq.obj, features = MT_genes$`Gene stable ID`)
CTC_reseq.obj[["nFeature.per.nCount"]] <- log10(CTC_reseq.obj$nFeature_RNA_seq) / log10(CTC_reseq.obj$nCount_RNA_seq)


# At this point I will sepearate the metadata from the seurat object, so I won't screw up the object
# while trying out different things
CTC_meta.data <- CTC_reseq.obj@meta.data
CTC_meta.data_clean <- CTC_meta.data[!is.na(CTC_meta.data$mt.percent), ]
CTC_meta.data_clean <- mutate(CTC_meta.data_clean, Cell_names = rownames(CTC_meta.data_clean))
rm(CTC_meta.data)


# Visualize the cell numbers / resistance group and save it
cn_p <- ggplot(CTC_meta.data_clean,
               aes(x = RespGroup, fill = RespGroup)) +
    geom_bar(color = "black") +
    geom_text(aes(label = after_stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) + #thx Bing
    scale_fill_tron() +
    labs(x = "Responde groups", y = "Cell numbers", fill = "Response groups", color = "Response groups") +
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
                aes(x = mt.percent, fill = RespGroup))+
    geom_density(alpha = 0.2) +
    scale_fill_tron() +
    labs(x = "MT gene expression percentage", y = "MT gene percent distribution") +
    theme_classic() +
    ggtitle("Mitochondrial gene expression") +
    theme(plot.title = element_text(hjust = 0.5))
print(mt_pd)

ggsave(filename = paste0("Mithocondrial gene expression", script_version, ".png"), plot = mt_pd,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 1500, height = 1000, units = "px")
rm(mt_pd)


# Visualizing the nCount numbers/treatment groups and save it
# The scales package helps to remove scientific notation form the axises
nCount_p <- ggplot(CTC_meta.data_clean,
                   aes(x = nCount_RNA_seq, fill = RespGroup))+
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
               aes(x = nFeature_RNA_seq, fill = RespGroup))+
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
               aes(x = nFeature.per.nCount, fill = RespGroup))+
    geom_density(alpha = 0.2) +
    scale_fill_tron() +
    labs(x = "nFeatures/nCounts", y = expression("log" [10]* " density"), fill = "Response group") +
    theme_classic() +
    ggtitle("Cell complexity distribution (nFeture/nCount)") +
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
    labs( x = "mRNA counts", y = "Gene counts", color = "MT gene expression percentage") +
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


# Clean objects
rm(CTC_meta.data_clean)




















### Data filtering


# I set the lower cutoff here to 1000 genes (anaysis tutorial, 250 genes, other sources 200 genes
# Jonathan lower 500) and the upper cutoff to 7000 (Jonathan upper 5000) there seems to be no consensus
# regarding this metric. I ultimately would not drow an upper limit as the population with
# high gen/mRNA counts seem well integrated to the corresponding population. Regardless the advice is
# that one should look the 3 parameter (nCount, nGene, MTratio) together and set the 
# thresholds accordingly (that is why I love the UvGplot)


# Filter the trimmed Seurat object according to the above decided parameters
filt_CTC.obj <- subset(CTC_reseq.obj, 
                           subset = nCount_RNA_seq >= 4000 & nFeature_RNA_seq >= 1000 &
                                    nFeature_RNA_seq <= 6500)


# Save the filtered seurat object for later use
saveRDS(object = filt_CTC.obj, file = paste0(getwd(), "/", script_version, "/", "Filtered_CTC_Seurat_object", script_version, ".rds"))


# Load the filtered surat object if not loaded yet
filt_CTC.obj <- readRDS(paste0(getwd(), "/", script_version, "/", "Filtered_CTC_Seurat_object", script_version, ".rds"))










## Another advised step is to remove genes with 0 counts across all cells


#There are two ways to do this:


# 1st

# For this first extract the counts from the seurat object using the GetAssayData function
#counts <- GetAssayData(filt_CTC.obj, slot = "counts")
#head(counts)

# Next, create a boolian matrix, describing if the count is bigger than 0 or not
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


# Checking the new cell numbers
cn_p <- ggplot(filt_CTC.obj@meta.data,
               aes(x = RespGroup, fill = RespGroup)) +
    geom_bar(color = "black") +
    geom_text(aes(label = after_stat(count)), stat = "count", position = position_dodge(width = 0.9), vjust = -0.5) +
    scale_fill_tron() +
    labs(x = "Response groups", y = "Cell numbers", fill = "Response groups") +
    ggtitle("Filtered cell numbers after sequencing") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
print(cn_p)

ggsave(filename = paste0("Filtered cell numbers", script_version, ".png"), plot = cn_p,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 1500, height = 1000, units = "px")
rm(cn_p)


# Visualizing the mt.percent across cells using density plot and save it
mt_pd <- ggplot(filt_CTC.obj@meta.data,
                aes(x = mt.percent, fill = RespGroup))+
    geom_density(alpha = 0.2) +
    scale_fill_tron() +
    theme_classic() +
    labs(x = "Cell percentage", y = "MT gene percent distribution", fill = "Response groups") +
    ggtitle("Filtered mitochondrial gene expression") +
    theme(plot.title = element_text(hjust = 0.5))
print(mt_pd)

ggsave(filename = paste0("Filtered mithocondrial gene expression", script_version, ".png"), plot = mt_pd,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 1500, height = 1000, units = "px")
rm(mt_pd)


# Visualizing the nCount numbers/treatment groups and save it
nCount_p <- ggplot(filt_CTC.obj@meta.data,
                   aes(x = nCount_RNA_seq, fill = RespGroup))+
    geom_density(alpha = 0.2) +
    scale_fill_tron() +
    labs(x = "nCounts", y = expression("log"[10]* " Count density"), fill = "Response group") +
    scale_x_log10(labels = scales::comma) + 
    theme_classic() +
    #geom_vline(xintercept = 2000, color = "red") +
    ggtitle("Filtered mRNA count distribution") +
    theme(plot.title = element_text(hjust = 0.5))
print(nCount_p)

ggsave(filename = paste0("Filtered mRNA count distribution", script_version, ".png"), plot = nCount_p,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 1500, height = 1000, units = "px")
rm(nCount_p)


# Visualizing the nFeature numbers/treatment groups and save it
nFeature_p <- ggplot(filt_CTC.obj@meta.data,
                     aes(x = nFeature_RNA_seq, fill = RespGroup))+
    geom_density(alpha = 0.2) +
    scale_fill_tron() +
    labs(x = "nFeatures", y = expression("log"[10]* " Feature density"), fill = "Response group") +
    scale_x_log10(labels = scales::comma) + 
    theme_classic() +
    ggtitle("Filtered gene count distribution") +
    #geom_vline(xintercept = 1000, color = "red") +
    #geom_vline(xintercept = 6500, color = "orange") +
    theme(plot.title = element_text(hjust = 0.5))
print(nFeature_p)

ggsave(filename = paste0("Filtered gene count distribution", script_version, ".png"), plot = nFeature_p,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 1500, height = 1000, units = "px")
rm(nFeature_p)


# Cell complexity plotting by checking the nFeature/nCount ratio (the higher the nFeature/nCount
# the more complex the cells are) and save it
nFpernC_p <-ggplot(filt_CTC.obj@meta.data,
                   aes(x = nFeature.per.nCount, fill = RespGroup))+
    geom_density(alpha = 0.2) +
    scale_fill_tron() +
    labs(x = "nFeatures/nCounts", y = expression("log"[10]* " density"), fill = "Response group") +
    theme_classic() +
    ggtitle("Cell complexity distribution (nFeture/nCount)") +
    theme(plot.title = element_text(hjust = 0.5))
print(nFpernC_p)

ggsave(filename = paste0("Filtered cell complexity distribution", script_version, ".png"), plot = nFpernC_p,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 1500, height = 1000, units = "px")
rm(nFpernC_p)


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
print(nCvnF_p)

ggsave(filename = paste0("Filtered nCount vs nFeature - QC plot", script_version, ".png"), plot = nCvnF_p,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 1500, height = 1000, units = "px")
rm(nCvnF_p)




















### Normalization and variance stabilization


## First things first, we want to look at factors commonly influencing clustering behaviour like cell cycle status or mt gene expression
## and see if they have enough effect on cell clustering or not in order to determine if want to regress them out or not

# Cell cycle
# First we check for cell cycle effects for this we do a rough nCount normalization with the
# NormalizeData() function (will divide the nCounts with the cell number and does a log10 transform)
ccPhase_CTC.obj <- NormalizeData(filt_CTC.obj)


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


# Testing the significace of PCs with the Jackstraw method and plotting the PCs
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

ggsave(filename = paste0("Cell cylce phase clustering (split)", script_version, ".png"), CellCyc_split_p2,
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

ggsave(filename = paste0("Cell cylce phase clustering", script_version, ".png"), CellCyc_p2,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 1500, height = 1000, units = "px")
rm(CellCyc_p, CellCyc_p2)


# Removing the ccPhase and ither unnecessary objects
rm(ccPhase_CTC.obj, CTC_reseq.obj)


# Mithocondrial gene expression
# First we take a look at the mitotic gene expression percentage and check the quartile values
#summary(ccPhase_CTC.obj@meta.data$mt.percent)


# Then, we turn the mt.percent into a chategorical factor vector, based on the quartile values
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


# Overlay the mt.percent using the Featureplot
#mt.per_p <- FeaturePlot(ccPhase_CTC.obj,
#                        features = "mt.percent")
#print(mt.per_p)
#rm(mt.per_p)

# In this particular case both the cell cycle and mitotic gene expression data scatters together, and there are no visible clusters
# ased on the cell cycle phase or mitotic gene expression so they do not have to be regressed out




















### Final normalization and regressing out sources of unwanted variation


## Normalization and clustering without integration

# As for the final more accurate method of normalization SCTransform, which is estimating the variance of the raw filtered data,
# and identifying the most variable genes. First i want to try and look at the data without integration.

# Try both the newer and older SCTransform version to see which one works better
sct_CTC.obj <- SCTransform(filt_CTC.obj, assay = "RNA_seq", method = "glmGamPoi")

#sct_CTC.obj <- SCTransform(filt_CTC.obj, assay = "RNA_seq")
# NOTE: The older version seems to work better, as it gives more defined clusters


# Save the generated sCTransformed seurat object
saveRDS(object = sct_CTC.obj, file = paste0(getwd(), "/", script_version, "/", "SCTransformed_CTC_Seurat_object", script_version, ".rds"))


# If not loaded, load the saved SCTransformed CTC Seurat object
sct_CTC_.obj <- readRDS(file = paste0(getwd(), "/", script_version, "/", "SCTransformed_CTC_Seurat_object", script_version, ".rds"))


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


# Another way to determine how many PCs one wnats to use is by using the PCAtools package
# NOTE: for the PCAtools to work one cannot use the seurat object itself, so some prep work is required


# For this check to work we will need to re-do the PCA with the PCATools. In order to reproduce the same PCA
# as we got with seurat we have to set the same seed as Seurat uses for the PCATools package and 
# pull the variance data from the SCTransformed object (from the SCT slot)
set.seed(42)
sct_variance <- sct_CTC.obj@assays$SCT@scale.data[VariableFeatures(sct_CTC.obj), ]


# Now that we set the seed and pulled the variance data we can run the pca
# NOTE: in order to reproduce the Seurat PCA proprly wee need to specifying the approximation method Seurat uses (IrlbaParam()) 
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

# IMPORTANT NOTE: the various methos usually do not agree on how many PCs one should keep, so you probably should try a couple and
# see if the UMAP and clustering makes sense based on what you set or not...


# Running dimensionality reduction (trying both UMAP and tSNE :) )
sct_CTC.obj <- RunUMAP(sct_CTC.obj, dims = 1:20, reduction = "pca", verbose = FALSE)
sct_CTC.obj <- RunTSNE(sct_CTC.obj, dims = 1:20, reduction = "pca", perplexity = 5)


# Finding neighbours and clustering
sct_CTC.obj <- FindNeighbors(sct_CTC.obj, reduction = "pca", dims = 1:20)
sct_CTC.obj <- FindClusters(sct_CTC.obj, resolution = c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2))


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

UMAP_p <- DimPlot(sct_CTC.obj, reduction = "umap")
print(UMAP_p)
rm(UMAP_p)

TSNE_p <- DimPlot(sct_CTC.obj, reduction = "tsne")
print(TSNE_p)
rm(TSNE_p)

PCA_p <- PCAPlot(sct_CTC.obj, reduction = "pca")
print(PCA_p)
rm(PCA_p)

# Re-run FindCluster with the prefered seeing
sct_CTC.obj <- FindClusters(sct_CTC.obj, resolution = 0.6)


# Saving the fully prepared CTC seurat object
saveRDS(object = sct_CTC.obj, file = paste0(getwd(), "/", script_version, "/", "SCTransformed_fully_prepared_CTC_Seurat_object", script_version, ".rds"))


# If not loaded, load the saved SCTransformed CTC Seurat object
sct_CTC_.obj <- readRDS(file = paste0(getwd(), "/", script_version, "/", "SCTransformed_fully_prepared_CTC_Seurat_object", script_version, ".rds"))


# Cluster visualization (PCA, UMAP and t-SNE plots)
cluster_colors <- pal_tron(palette = "legacy")(7)
cluster_colors_2 <- pal_futurama(palette = "planetexpress")(12)
cluster_colors_3 <- pal_startrek(palette = "uniform")(7)
cluster_colors_4 <- pal_rickandmorty(palette = "schwifty")(12)

PCA_p <- DimPlot(sct_CTC.obj, reduction = "pca", group.by = "seurat_clusters")
PCA_p2 <- PCA_p +
    scale_fill_manual(values = cluster_colors[c(6, 7)]) +
    scale_color_manual(values = cluster_colors[c(6, 7)]) +
    labs(x = "PC 1", y = "PC 2", color = "Clusters", fill = "Clusters") +
    ggtitle("PCA plot")
print(PCA_p2)

ggsave(filename = paste0("PCA_clustering_gb", script_version, ".png"), PCA_p2,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 1500, height = 1000, units = "px")
rm(PCA_p, PCA_p2)


UMAP_p <- DimPlot(sct_CTC.obj, reduction = "umap", group.by = "seurat_clusters")
UMAP_p2 <- UMAP_p +
    scale_fill_manual(values = cluster_colors[c(6, 7)]) +
    scale_color_manual(values = cluster_colors[c(6, 7)]) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Clusters", fill = "Clusters") +
    ggtitle("UMAP plot")
print(UMAP_p2)

ggsave(filename = paste0("UMAP_clustering_gb", script_version, ".png"), UMAP_p2,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 1500, height = 1000, units = "px")
rm(UMAP_p, UMAP_p2)


TSNE_p <- DimPlot(sct_CTC.obj, reduction = "tsne", group.by = "seurat_clusters")
TSNE_p2 <- TSNE_p +
    scale_fill_manual(values = cluster_colors[c(6, 7)]) +
    scale_color_manual(values = cluster_colors[c(6, 7)]) +
    labs(x = "tSNE 1", y = "tSNE 2", color = "Clusters", fill = "Clusters") +
    ggtitle("tSNE plot")
print(TSNE_p2)

ggsave(filename = paste0("TSNE_clustering_gb", script_version, ".png"), TSNE_p2,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 1500, height = 1000, units = "px")
rm(TSNE_p, TSNE_p2)

# Removing unnecessary objects
rm(g2m.IDs, s.IDs, MT_genes, QC_PCA, QC_PCA_par)










## Plotting different metrics over the UMAP clusters to if they correlate

# Overlaying nFeature
grad_col <- pal_locuszoom("default")(7)

UMAP_nF_p <- FeaturePlot(sct_CTC.obj,
                         reduction = "umap",
                         features = "nFeature_SCT")
UMAP_nF_p2 <- UMAP_nF_p +
    #scale_color_gradient(low = "lightgrey", high = cluster_colors[7]) +
    scale_color_gradient(low = "lightgrey", high = grad_col[5]) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Clusters", fill = "Clusters") +
    ggtitle("UMAP plot - nFeatures")
print(UMAP_nF_p2)

ggsave(filename = paste0("UMAP_nFeature_gb", script_version, ".png"), UMAP_nF_p2,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 1500, height = 1000, units = "px")
rm(UMAP_nF_p, UMAP_nF_p2)


# Overlaying nCount
UMAP_nC_p <- FeaturePlot(sct_CTC.obj,
                         reduction = "umap",
                         features = "nCount_SCT")
UMAP_nC_p2 <- UMAP_nC_p +
    scale_color_gradient(low = "lightgrey", high = grad_col[5]) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Clusters", fill = "Clusters") +
    ggtitle("UMAP plot - nCounts")
print(UMAP_nC_p2)

ggsave(filename = paste0("UMAP_nCount_gb", script_version, ".png"), UMAP_nC_p2,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 1500, height = 1000, units = "px")
rm(UMAP_nC_p, UMAP_nC_p2)


# Overlay mt.percent
UMAP_mt.percent <- FeaturePlot(sct_CTC.obj,
                               reduction = "umap",
                               features = "mt.percent")
UMAP_mt.percent2 <- UMAP_mt.percent +
    scale_color_gradient(low = "lightgrey", high = grad_col[5]) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Clusters", fill = "Clusters") +
    ggtitle("UMAP plot - MT gene expression")
print(UMAP_mt.percent2)

ggsave(filename = paste0("UMAP_mt.percent_gb", script_version, ".png"), UMAP_mt.percent2,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 1500, height = 1000, units = "px")
rm(UMAP_mt.percent, UMAP_mt.percent2)


# Overlay cell cycle
UMAP_CellCyc <- DimPlot(sct_CTC.obj,
                        reduction = "umap",
                        group.by = "Phase")
UMAP_CellCyc2 <- UMAP_CellCyc +
    scale_fill_manual(values = cluster_colors_3[c(4, 5, 7)]) +
    scale_color_manual(values = cluster_colors_3[c(4, 5 ,7)]) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Cell cycle phase", fill = "Cell cycle phase") +
    ggtitle("UMAP plot - Cell cycle phases")
print(UMAP_CellCyc2)

ggsave(filename = paste0("UMAP_Cell_cycle_gb", script_version, ".png"), UMAP_CellCyc2,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 1500, height = 1000, units = "px")
rm(UMAP_CellCyc, UMAP_CellCyc2)


# Overlay treatment cycle
UMAP_TreatCyc <- DimPlot(sct_CTC.obj,
                         reduction = "umap",
                         group.by = "RealCycle")
UMAP_TreatCyc2 <- UMAP_TreatCyc +
    scale_fill_manual(values = cluster_colors_3[c(4, 5, 7)]) +
    scale_color_manual(values = cluster_colors_3[c(4, 5, 7)]) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Treatment cycles", fill = "Treatment cyceles") +
    ggtitle("UMAP plot - Treatment cycles")
print(UMAP_TreatCyc2)

ggsave(filename = paste0("UMAP_Treatment_cycles_gb", script_version, ".png"), UMAP_TreatCyc2,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 1500, height = 1000, units = "px")
rm(UMAP_TreatCyc, UMAP_TreatCyc2)


#overlay treatment type
UMAP_TreatType <- DimPlot(sct_CTC.obj,
                         reduction = "umap",
                         group.by = "Treatment")
UMAP_TreatType2 <- UMAP_TreatType +
    scale_fill_manual(values = cluster_colors_3[c(4, 5, 7)]) +
    scale_color_manual(values = cluster_colors_3[c(4, 5, 7)]) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Treatment", fill = "Treatment") +
    ggtitle("UMAP plot - Treatment types")
print(UMAP_TreatType2)

ggsave(filename = paste0("UMAP_Treatment_type_gb", script_version, ".png"), UMAP_TreatType2,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 1500, height = 1000, units = "px")
rm(UMAP_TreatType, UMAP_TreatType2)


# Overlay RespGroup
UMAP_RespG <- DimPlot(sct_CTC.obj,
                      reduction = "umap",
                      group.by = "RespGroup")
UMAP_RespG2 <- UMAP_RespG +
    scale_fill_manual(values = cluster_colors[c(1, 2)]) +
    scale_color_manual(values = cluster_colors[c(1, 2)]) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Response group", fill = "Response group") +
    ggtitle("UMAP plot - Response groups")
print(UMAP_RespG2)

ggsave(filename = paste0("UMAP_Response_group_gb", script_version, ".png"), UMAP_RespG2,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 1500, height = 1000, units = "px")
rm(UMAP_RespG, UMAP_RespG2)


#overlay Luthetium-177 doses
UMAP_Lu <- DimPlot(sct_CTC.obj,
                      reduction = "umap",
                      group.by = "Lu177_GBq")
UMAP_Lu2 <- UMAP_Lu +
    scale_fill_manual(values = cluster_colors_4[c(4, 5, 6, 7, 8, 9)]) +
    scale_color_manual(values = cluster_colors_4[c(4, 5, 6, 7, 8, 9)]) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Lu-177 doses (GBq)", fill = "Lutetium-177 doses (GBq)") +
    ggtitle("UMAP plot - Lu-177 treatment doses")
print(UMAP_Lu2)

ggsave(filename = paste0("UMAP_Lu-177_treatment_doses_gb", script_version, ".png"), UMAP_Lu2,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 1500, height = 1000, units = "px")
rm(UMAP_Lu, UMAP_Lu2)


# Overlay Actinium-225 doses
UMAP_Ac <- DimPlot(sct_CTC.obj,
                   reduction = "umap",
                   group.by = "Ac225_MBq")
UMAP_Ac2 <- UMAP_Ac +
    scale_fill_manual(values = cluster_colors_4[c(4, 5, 6, 7)]) +
    scale_color_manual(values = cluster_colors_4[c(4, 5, 6, 7)]) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Ac-225 doses (MBq)", fill = "Ac-225 doses (MBq)") +
    ggtitle("UMAP plot - Ac-225 treatment doses")
print(UMAP_Ac2)

ggsave(filename = paste0("UMAP_Ac-225_treatment_doses_gb", script_version, ".png"), UMAP_Ac2,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 1500, height = 1000, units = "px")
rm(UMAP_Ac, UMAP_Ac2)


# Overlay marker status
UMAP_MarkerStat <- DimPlot(sct_CTC.obj,
                          reduction = "umap",
                          group.by = "Marker.status")
UMAP_MarkerStat2 <- UMAP_MarkerStat +
    scale_fill_manual(values = cluster_colors_3[c(4, 5, 7)]) +
    scale_color_manual(values = cluster_colors_3[c(4, 5, 7)]) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Marker status", fill = "Marker status") +
    ggtitle("UMAP plot - Surface marker distribution")
print(UMAP_MarkerStat2)

ggsave(filename = paste0("UMAP_Marker_Status_gb", script_version, ".png"), UMAP_MarkerStat2,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 1500, height = 1000, units = "px")
rm(UMAP_MarkerStat, UMAP_MarkerStat2)










## Visualizing the overlay data in a more understandable manner, using bar charts

# Preparing the main dataframe for the visualization
head(sct_CTC.obj@meta.data)

sct_CTC_meta_df <- sct_CTC.obj@meta.data
sct_CTC_meta_df$patients <- str_remove(sct_CTC_meta_df$patientID, "." %R% END) 

# Visualizing data from the variables side
# Visualizing the treatment cycle data
# Creatiing the plotting df
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
                                      Cell_percent_c0 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 0, seurat_clusters == 0))/
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 0)) * 100,
                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 1, seurat_clusters == 0))/
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 1)) * 100,
                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 2, seurat_clusters == 0))/
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 2)) * 100),
                                      Cell_percent_c1 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 0, seurat_clusters == 1))/
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 0)) * 100,
                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 1, seurat_clusters == 1))/
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 1)) * 100,
                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 2, seurat_clusters == 1))/
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, RealCycle == 2)) * 100))
print(Treatment_cycles_summary_df)


# Transforming the plotting df to a long format
Treatment_cycles_summary_df_long <- pivot_longer(Treatment_cycles_summary_df, cols = c(Cell_percent_c0, Cell_percent_c1), names_to = "Clusters", values_to = "Percentages")
print(Treatment_cycles_summary_df_long)


# Plotting the treatment cycle distribution
treatment_cyc_bar_p <- ggplot(Treatment_cycles_summary_df_long,
                              aes(x = Cycles, y = Percentages, fill = Clusters)) +
    geom_col(color = "black", position = "dodge") +
    geom_text(aes(label = as.integer(Percentages)), vjust = -0.5, position = position_dodge(width = 0.9)) +
    scale_y_continuous(limits = c(0,100), breaks = c(0, 20, 40, 60, 80, 100)) +
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


# Visualizing the treatment type data
# Creatiing the plotting df
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
                                      Cell_percent_c0 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "AL", seurat_clusters == 0))/
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "AL")) * 100,
                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "L", seurat_clusters == 0))/
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "L")) * 100,
                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "NT", seurat_clusters == 0))/
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "NT")) * 100),
                                      Cell_percent_c1 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "AL", seurat_clusters == 1))/
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "AL")) * 100,
                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "L", seurat_clusters == 1))/
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "L")) * 100,
                                                          nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "NT", seurat_clusters == 1))/
                                                              nrow(dplyr::filter(.data = sct_CTC_meta_df, Treatment == "NT")) * 100))
print(Treatment_summary_df)


# Transforming the plotting df to a long format
Treatment_summary_df_long <- pivot_longer(Treatment_summary_df, cols = c(Cell_percent_c0, Cell_percent_c1), names_to = "Clusters", values_to = "Percentages")
print(Treatment_summary_df_long)


# Plotting the treatment type distribution
treatment_type_bar_p <- ggplot(Treatment_summary_df_long,
                               aes(x = factor(Treatment_type, levels = c("NT", "AL", "L")), y = Percentages, fill = Clusters)) +
    geom_col(color = "black", position = "dodge") +
    geom_text(aes(label = as.integer(Percentages)), vjust = -0.5, position = position_dodge(width = 0.9)) +
    scale_y_continuous(limits = c(0,100), breaks = c(0, 20, 40, 60, 80, 100)) +
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


# Visualizing the cell cycle data
# Creatiing the plotting df
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
                                   Cell_percent_c0 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "G1", seurat_clusters == 0))/
                                                           nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "G1")) * 100,
                                                       nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "S", seurat_clusters == 0))/
                                                           nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "S")) * 100,
                                                       nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "G2M", seurat_clusters == 0))/
                                                           nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "G2M")) * 100),
                                   Cell_percent_c1 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "G1", seurat_clusters == 1))/
                                                           nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "G1")) * 100,
                                                       nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "S", seurat_clusters == 1))/
                                                           nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "S")) * 100,
                                                       nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "G2M", seurat_clusters == 1))/
                                                           nrow(dplyr::filter(.data = sct_CTC_meta_df, Phase == "G2M")) * 100))
print(Cell_cycle_summary_df)


# Transforming the plotting df to a long format
Cell_cycle_summary_df_long <- pivot_longer(Cell_cycle_summary_df, cols = c(Cell_percent_c0, Cell_percent_c1), names_to = "Clusters", values_to = "Percentages")
print(Cell_cycle_summary_df_long)


# Plotting the cell cycle distribution
cell_cyc_bar_p <- ggplot(Cell_cycle_summary_df_long,
                         aes(x = factor(CC_phase, levels = c("G1", "S", "G2M")), y = Percentages, fill = Clusters)) +
    geom_col(color = "black", position = "dodge") +
    geom_text(aes(label = as.integer(Percentages)), vjust = -0.5, position = position_dodge(width = 0.9)) +
    scale_y_continuous(limits = c(0,100), breaks = c(0, 20, 40, 60, 80, 100)) +
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


# Visualizing the response group data
# Creatiing the plotting df
Response_summary_df <- data.frame(response_group = c("Responder", "Nonresponder"),
                                    Total_cell_numbers = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Responder")),
                                                           nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Nonresponder"))),
                                    Cell_numbers_c0 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Responder", seurat_clusters == 0)),
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Nonresponder", seurat_clusters == 0))),
                                    Cell_numbers_c1 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Responder", seurat_clusters == 1)),
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Nonresponder", seurat_clusters == 1))),
                                    Cell_percent_c0 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Responder", seurat_clusters == 0))/
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Responder")) * 100,
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Nonresponder", seurat_clusters == 0))/
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Nonresponder")) * 100),
                                    Cell_percent_c1 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Responder", seurat_clusters == 1))/
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Responder")) * 100,
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Nonresponder", seurat_clusters == 1))/
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, RespGroup == "Nonresponder")) * 100))
print(Response_summary_df)


# Transforming the plotting df to a long format
Response_summary_df_long <- pivot_longer(Response_summary_df, cols = c(Cell_percent_c0, Cell_percent_c1), names_to = "Clusters", values_to = "Percentages")
print(Response_summary_df_long)


# Plotting the treatment response distribution from response view
response_group_bar_p <- ggplot(Response_summary_df_long,
                               aes(x = factor(response_group, levels = c("Responder", "Nonresponder")), y = Percentages, fill = Clusters)) +
    geom_col(color = "black", position = "dodge") +
    geom_text(aes(label = as.integer(Percentages)), vjust = -0.5, position = position_dodge(width = 0.9)) +
    scale_y_continuous(limits = c(0,100), breaks = c(0, 20, 40, 60, 80, 100)) +
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


# Visualizing the sureface marker distribution data
# Creatiing the plotting df form the markers side
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
                                    Cell_percent_c0 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "EpCAM+", seurat_clusters == 0))/
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "EpCAM+")) * 100,
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "EpCAM+PSMA+", seurat_clusters == 0))/
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "EpCAM+PSMA+")) * 100,
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "PSMA+", seurat_clusters == 0))/
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "PSMA+")) * 100),
                                    Cell_percent_c1 = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "EpCAM+", seurat_clusters == 1))/
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "EpCAM+")) * 100,
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "EpCAM+PSMA+", seurat_clusters == 1))/
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "EpCAM+PSMA+")) * 100,
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "PSMA+", seurat_clusters == 1))/
                                                            nrow(dplyr::filter(.data = sct_CTC_meta_df, Marker.status == "PSMA+")) * 100))
print(Surface_marker_summary_df)


# Transform the plotting df to the long format for plotting
Surface_marker_summary_df_long <- pivot_longer(Surface_marker_summary_df, cols = c(Cell_percent_c0, Cell_percent_c1), names_to = "Clusters", values_to = "Percentages")
print(Surface_marker_summary_df_long)


# Plotting the surface marker distribution
surface_marker_bar_p <- ggplot(Surface_marker_summary_df_long,
                               aes(x = factor(marker_status, levels = c("EpCAM+", "EpCAM+PSMA+", "PSMA+")), y = Percentages, fill = Clusters)) +
    geom_col(color = "black", position = "dodge") +
    geom_text(aes(label = as.integer(Percentages)), vjust = -0.5, position = position_dodge(width = 0.9)) +
    scale_y_continuous(limits = c(0,100), breaks = c(0, 20, 40, 60, 80, 100)) +
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





# Visualizing data from the cluster side
# Visualizing the treatment cycle data
# Creatiing the plotting df
Surface_marker_summary_df_clust <- data.frame(Clusters = c("Cluster 0", "Cluster 1"),
                                            Total_cell_numbers = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")), 
                                                                    nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1"))),
                                            Cycle_0_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", RealCycle == "0")),
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", RealCycle == "0"))),
                                            Cycle_1_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", RealCycle == "1")),
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", RealCycle == "1"))),
                                            Cycle_2_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", RealCycle == "2")),
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", RealCycle == "2"))),
                                            Cycle_0_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", RealCycle == "0"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", RealCycle == "0"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100),
                                            Cycle_1_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", RealCycle == "1"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", RealCycle == "1"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) *100),
                                            Cycle_2_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", RealCycle == "2"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", RealCycle == "2"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100))
print(Surface_marker_summary_df_clust)


# Tranform the plotting df to a long format for plotting
Treatment_cycles_summary_df_clust_long <- pivot_longer(Treatment_cycles_summary_df_clust, cols = c(Cycle_0_p, Cycle_1_p, Cycle_2_p), names_to = "Cycles", values_to = "Percentages")
print(Treatment_cycles_summary_df_clust_long)

# Plotting the treatment cycle distribution
treatment_cyc_bar_p <- ggplot(Treatment_cycles_summary_df_clust_long,
                              aes(x = Clusters, y = Percentages, fill = as.factor(Cycles))) +
    geom_col(color = "black", position = "dodge") +
    geom_text(aes(label = as.integer(Percentages)), vjust = -0.5, position = position_dodge(width = 0.9)) +
    scale_y_continuous(limits = c(0,100), breaks = c(0, 20, 40, 60, 80, 100)) +
    scale_x_discrete(labels = c("Cluster 0","Cluster 1","Cluster 2")) +
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


# Visualizing the treatment type data
# Creatiing the plotting df
Treatment_summary_df_clust <- data.frame(Clusters = c("Cluster 0", "Cluster 1"),
                                            Total_cell_numbers = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")), 
                                                                    nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1"))),
                                            AL_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Treatment == "AL")),
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Treatment == "AL"))),
                                            L_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Treatment == "L")),
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Treatment == "L"))),
                                            NT_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Treatment == "NT")),
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Treatment == "NT"))),
                                            AL_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Treatment == "AL"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Treatment == "AL"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100),
                                            L_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Treatment == "L"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Treatment == "L"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) *100),
                                            NT_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Treatment == "NT"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Treatment == "NT"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100))
print(Treatment_summary_df_clust)


# Transforming the plotting df to a long format
Treatment_summary_df_clust_long <- pivot_longer(Treatment_summary_df_clust, cols = c(AL_p, L_p, NT_p), names_to = "Treatment_type", values_to = "Percentages")
print(Treatment_summary_df_clust_long)


# Plotting the treatment type distribution
treatment_type_bar_p <- ggplot(Treatment_summary_df_clust_long,
                               aes(x = Clusters, y = Percentages, fill = factor(Treatment_type, levels = c("NT_p", "AL_p", "L_p")))) +
    geom_col(color = "black", position = "dodge") +
    geom_text(aes(label = as.integer(Percentages)), vjust = -0.5, position = position_dodge(width = 0.9)) +
    scale_y_continuous(limits = c(0,100), breaks = c(0, 20, 40, 60, 80, 100)) +
    scale_x_discrete(labels = c("Cluster 0","Cluster 1","Cluster 2")) +
    scale_fill_manual(values = cluster_colors[c(4, 6, 7)], label = c("NT", "AL", "L")) +
    labs(fill = "Treatment types", x = expression("CTC custers"), y = expression("Cell number percentages")) +
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


# Visualizing the cell cycle data
# Creatiing the plotting df
Cell_cycle_summary_df_clust <- data.frame(Clusters = c("Cluster 0", "Cluster 1"),
                                            Total_cell_numbers = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")), 
                                                                    nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1"))),
                                            S_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Phase == "S")),
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Phase == "S"))),
                                            G1_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Phase == "G1")),
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Phase == "G1"))),
                                            G2M_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Phase == "G2M")),
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Phase == "G2M"))),
                                            S_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Phase == "S"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Phase == "S"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100),
                                            G1_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Phase == "G1"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Phase == "G1"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) *100),
                                            G2M_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Phase == "G2M"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Phase == "G2M"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100))
print(Cell_cycle_summary_df_clust)


# Transforming the plotting df to a long format
Cell_cycle_summary_df_clust_long <- pivot_longer(Cell_cycle_summary_df_clust, cols = c(S_p, G1_p, G2M_p), names_to = "CC_phase", values_to = "Percentages")
print(Cell_cycle_summary_df_clust_long)


# Plotting the cell cycle distribution
cell_cyc_bar_p <- ggplot(Cell_cycle_summary_df_clust_long,
                         aes(x = Clusters, y = Percentages, fill = factor(CC_phase, levels = c("G1_p", "S_p", "G2M_p")))) +
    geom_col(color = "black", position = "dodge") +
    geom_text(aes(label = as.integer(Percentages)), vjust = -0.5, position = position_dodge(width = 0.9)) +
    scale_y_continuous(limits = c(0,100), breaks = c(0, 20, 40, 60, 80, 100)) +
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


# Visualizing the response group data
# Creatiing the plotting df
Response_summary_df_clust <- data.frame(Clusters = c("Cluster 0", "Cluster 1"),
                                            Total_cell_numbers = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")), 
                                                                    nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1"))),
                                            Resp_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", RespGroup == "Responder")),
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", RespGroup == "Responder"))),
                                            Nonresp_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", RespGroup == "Nonresponder")),
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", RespGroup == "Nonresponder"))),
                                            Resp_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", RespGroup == "Responder"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", RespGroup == "Responder"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100),
                                            Nonresp_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", RespGroup == "Nonresponder"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", RespGroup == "Nonresponder"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) *100))
print(Response_summary_df_clust)


# Transforming the plotting df to a long format
Response_summary_df_clust_long <- pivot_longer(Response_summary_df_clust, cols = c(Resp_p, Nonresp_p), names_to = "Response_group", values_to = "Percentages")
print(Response_summary_df_clust_long)


# Plotting the treatment response distribution from response view
response_group_bar_p <- ggplot(Response_summary_df_clust_long,
                               aes(x = Clusters, y = Percentages, fill = factor(Response_group, levels = c("Resp_p", "Nonresp_p")))) +
    geom_col(color = "black", position = "dodge") +
    geom_text(aes(label = as.integer(Percentages)), vjust = -0.5, position = position_dodge(width = 0.9)) +
    scale_y_continuous(limits = c(0,100), breaks = c(0, 20, 40, 60, 80, 100)) +
    scale_x_discrete(labels = c("Cluster 0", "Cluster 1", "Cluster 2")) +
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


# Visualizing the marker status data
# Creatiing the plotting df form the cluster side
Surface_marker_summary_df_clust <- data.frame(Clusters = c("Cluster 0", "Cluster 1"),
                                            Total_cell_numbers = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")), 
                                                                    nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1"))),
                                            EpCAM_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Marker.status == "EpCAM+")),
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Marker.status == "EpCAM+"))),
                                            EpCAM_PSMA_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Marker.status == "EpCAM+PSMA+")),
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Marker.status == "EpCAM+PSMA+"))),
                                            PSMA_n = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Marker.status == "PSMA+")),
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Marker.status == "PSMA+"))),
                                            EpCAM_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Marker.status == "EpCAM+"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Marker.status == "EpCAM+"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100),
                                            EpCAM_PSMA_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Marker.status == "EpCAM+PSMA+"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Marker.status == "EpCAM+PSMA+"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) *100),
                                            PSMA_p = c(nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0", Marker.status == "PSMA+"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "0")) * 100,
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1", Marker.status == "PSMA+"))/
                                                        nrow(dplyr::filter(.data = sct_CTC_meta_df, seurat_clusters == "1")) * 100))
print(Surface_marker_summary_df_clust)


# Tranform the plotting df to a long format for plotting
Surface_marker_summary_df_clust_long <- pivot_longer(Surface_marker_summary_df_clust, cols = c(EpCAM_p, EpCAM_PSMA_p, PSMA_p), names_to = "Markers", values_to = "Percentages")
print(Surface_marker_summary_df_clust_long)

surface_m_bar_p <- ggplot(Surface_marker_summary_df_clust_long,
                              aes(x = Clusters, y = Percentages, fill = as.factor(Markers))) +
    geom_col(color = "black", position = "dodge") +
    geom_text(aes(label = as.integer(Percentages)), vjust = -0.5, position = position_dodge(width = 0.9)) +
    scale_y_continuous(limits = c(0,100), breaks = c(0, 20, 40, 60, 80, 100)) +
    scale_x_discrete(labels = c("Cluster 0","Cluster 1","Cluster 2")) +
    scale_fill_manual(values = cluster_colors[c(4, 6, 7)], label = c("Cycle 0", "Cycle 1", "Cycle 2")) +
    labs(fill = "Treatment cycles", x = expression("CTC clusters"), y = expression("Cell number percentages")) +
    ggtitle("CTC clusters - cell distribution based on treatment cycle") +
    #facet_wrap(vars(seurat_clusters), strip.position = "bottom") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
          strip.text = element_blank())
print(surface_m_bar_p)

ggsave(filename = paste0("Cluster_composition-Surface_markers", script_version, ".png"), surface_m_bar_p,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 1500, height = 1000, units = "px")
rm(surface_m_bar_p)


# Remove the plotting dataframes
rm(sct_CTC_meta_df, Cell_cycle_summary_df, Cell_cycle_summary_df_long, Cell_cycle_summary_df_clust, Cell_cycle_summary_df_clust_long,
    Treatment_cycles_summary_df, Treatment_cycles_summary_df_long, Treatment_cycles_summary_df_clust, Treatment_cycles_summary_df_clust_long,
    Response_summary_df, Response_summary_df_long, Response_summary_df_clust, Response_summary_df_clust_long,
    Surface_marker_summary_df, Surface_marker_summary_df_long, Surface_marker_summary_df_clust, Surface_marker_summary_df_clust_long,
    Treatment_summary_df, Treatment_summary_df_long, Treatment_summary_df_clust, Treatment_summary_df_clust_long)














### Differential gene expression analysis

# Finding differentially expressed (DE) genes using the FindMarkers function to compare two
# clusters (if the Log2FC is positive, then it is upregulated in ident.1 vs ident.2 if it is negative then it is downregulated ident.1 vs ident.2)
DE_Cluster_0v1 <- FindMarkers(sct_CTC.obj, ident.1 = "0", ident.2 = "1")
head(DE_Cluster_0v1)


# Checking if there is differentil expression between the resp groups (not according to the clustering)
Idents(sct_CTC.obj) <- "RespGroup"
DE_RespvNonresp <- FindMarkers(sct_CTC.obj, ident.1 = "Responder", ident.2 = "Nonresponder")
head(DE_RespvNonresp)


# Looking for conserved markers between the clusters (the function needs a grouping variable, which I set to the same across all cells)
Idents(sct_CTC.obj) <- "seurat_clusters"
sct_CTC.obj[["grouping.var"]] <- "RNA"
CON_cluster_0v1 <- FindConservedMarkers(sct_CTC.obj, ident.1 = "0", ident.2 = "1", assay = "SCT", grouping.var = "grouping.var")
head(CON_cluster_0v1)


# Attaching gene names to the DE genes using the following fucntions
# This is the original, well trusted for loop based implementation
geneName_matching = function (DE_lst, gene_lst, ID_col, name_col) {
    progressBar <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                    total = nrow(DE_lst),
                                    complete = "=",
                                    incomplete = "-",
                                    current = ">",
                                    clear = FALSE,
                                    width = 100)
    position <- c()
    for(e in seq_along(rownames(DE_lst))) {
        progressBar$tick()
        
        position <-  grep(rownames(DE_lst)[e], gene_lst[, ID_col])
        DE_lst$geneName[e] <- gene_lst[, name_col][position]
        
    }
    
    original_name <- deparse(substitute(DE_lst))
    output_name <- paste0(original_name, "geneNames")
    assign(output_name, DE_lst, envir = .GlobalEnv)
    message("The gene name matched DE_cluster_0v1 list was assigned to a new object named:", "/n", output_name)
}
geneName_matching(DE_Cluster_0v1, unif_gene_names, "ensembl_gene_id", "external_gene_name")


# NOTE: this for loop is highly time and resource intensive and will take over 5 minutes to run.
# in exchange however it will insert a positionally accurate name for the DE genes in a well understood manner


# This is the parallelized version of the matching function
geneName_matching_par = function (DE_lst, gene_lst, ID_col, name_col, n_cores = 2) {
    progressBar <- progress_bar$new(format = "Progress = (:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
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

    opts <- list(progress = progress)
    
    cores <- n_cores
    clust <- parallel::makeCluster(cores)
    doSNOW::registerDoSNOW(cl = clust)
    
    position <- c()
    out <- c()
    geneNames <- foreach::foreach(i = seq_len(nrow(DE_lst)), .options.snow = opts, .combine = c) %dopar% {
        position <- grep(rownames(DE_lst)[i], gene_lst[, ID_col])
        out[i] <- gene_lst[, name_col][position]
    }
    
    parallel::stopCluster(cl = clust)
    
    out_lst <- cbind(DE_lst, geneNames)
    original_name <- deparse(substitute(DE_lst))
    output_name <- paste0(original_name, "_", "geneNames")
    assign(output_name, out_lst, envir = .GlobalEnv)
    message("The gene name matched DE_cluster_0v1 list was assigned to a new object named:", "\n", output_name)
}
geneName_matching_par(DE_Cluster_0v1, unif_gene_names, "ensembl_gene_id", "external_gene_name")


# NOTE: this version is less time intensive (runs under 5 minutes) as it runs the process on multiple cores in parallel, however
# it relies on on a foreach "loop" which is more complicated then a simple for loop




# Attaching some additional info to the DE gene list and saving the whole list
DE_Cluster_0v1_geneNames$ExprChange <- ifelse(DE_Cluster_0v1_geneNames$avg_log2FC >= 0 & DE_Cluster_0v1_geneNames$p_val_adj <= 0.05, "Upregulated",
                                    ifelse(DE_Cluster_0v1_geneNames$avg_log2FC <= -0 & DE_Cluster_0v1_geneNames$p_val_adj <= 0.05, "Downregulated", "Not significant"))

# NOTE: the ranking here will be solely based on the p_adj value (as the sign() function will turn all log2FC into a negative, 0 or positive 1)
# if you don't want that (which we want to do here) you can omit the sign() function
DE_Cluster_0v1_geneNames$Ranking <- sign(DE_Cluster_0v1_geneNames$avg_log2FC) * (-log10(DE_Cluster_0v1_geneNames$p_val_adj))


write_csv(x = DE_Cluster_0v1_geneNames, file = paste0(script_version, "/", "DE_clusters_0v1_bg", script_version, ".csv"))










## Visulazizing DE genes using a volcano plot

# Annotated visualization with ggplot2
DE_Cluster_0v1_geneNames$Significance <- ifelse(DE_Cluster_0v1$avg_log2FC > 0 & DE_Cluster_0v1$p_val_adj < 0.05, "Upregulated",
                                    ifelse(DE_Cluster_0v1$avg_log2FC < 0 & DE_Cluster_0v1$p_val_adj < 0.05, "Downregulated", "Not significant")) #labeling genes for coloring
upreg <- DE_Cluster_0v1_geneNames[DE_Cluster_0v1_geneNames$Significance == "Upregulated", ] #subsetting the upregulated genes
topupreg <- head(upreg$geneName[order(upreg$p_val_adj)], 10) #selecting the top upreg genes
downreg <- DE_Cluster_0v1_geneNames[DE_Cluster_0v1_geneNames$Significance == "Downregulated", ] #subsetting the downregulated genes
topdownreg <- head(downreg$geneName[order(downreg$p_val_adj)], 10) #selecting the top downreg genes


top50upreg <- head(DE_Cluster_0v1_geneNames$geneName[order(DE_Cluster_0v1_geneNames$p_val_adj)], 50) #choosing the top 50 differentially expressed genes
DE_Cluster_0v1_geneNames$DElabel <- ifelse(DE_Cluster_0v1_geneNames$geneNames %in% topupreg | DE_Cluster_0v1_geneNames$geneNames %in% topdownreg, DE_Cluster_0v1_geneNames$geneNames, NA) #marking the top DEGs in the dataframe

volcano_p <- ggplot(data = DE_Cluster_0v1_geneNames, 
                aes(x = avg_log2FC, y = -log10(p_val_adj), col = Significance, label = DElabel)) +
            geom_point() +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
            scale_color_manual(values = cluster_colors[c(2, 5, 3)],
                            labels = c("Downregulated", "Not significant", "Upregulated")) +
            scale_x_continuous(limits = c(-10, 10), breaks = seq(-10, 10, 2)) +
            scale_y_continuous(limits = c(0, 14), breaks = seq(0, 14, 2)) +
            labs(color = "Expression status", x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
            ggtitle("Differentially expressed genes") +
            geom_text_repel(show.legend = FALSE ,max.overlaps = Inf) +
            theme_classic() +
            theme(plot.title = element_text(hjust = 0.5))
print(volcano_p)

ggsave(filename = paste0("DEGs_Volcano_plot_gb", script_version, ".png"), volcano_p,
       device = "png", path = paste0(script_version, "/", "Plots/"),
       width = 3500, height = 3000, units = "px")
rm(volcano_p, upreg, downreg, topupreg, topdownreg)

# Alternative colors: "#EBCA5B", "grey", "#3E7CD6"


# Creating an interactive plot
library(plotly)
volcano <- plot_ly(data = DE_Cluster_0v1,
                   type = "scatter",
                   mode = "markers",
                   x = DE_Cluster_0v1$avg_log2FC,
                   y = -log10(DE_Cluster_0v1$p_val_adj),
                   text = ~paste("Gene:", geneName),
                   color = ~ExprChange,
                   colors = c("#EBCA5B", "grey", "#3E7CD6"))
volcano %>%
    add_lines(x = 1, y = range(-log10(DE_Cluster_0v1$p_val_adj)), showlegend = FALSE, inherit = FALSE,
              line = list(color = "black", dash = "dot")
    ) %>%
    add_lines(x = -1, y = range(-log10(DE_Cluster_0v1$p_val_adj)), showlegend = FALSE, inherit = FALSE,
              line = list(color = "black", dash = "dot")
    ) %>%
    add_lines(y = -log10(0.05), x = range(DE_Cluster_0v1$avg_log2FC), showlegend = FALSE, inherit = FALSE,
              line = list(color = "black", dash = "dot")
    )

save_image(p = volcano, file = paste0(script_version, "/", "Plots/", "DE_volcano_plot", script_version, ".png"), width = 1500, height = 1000, scale = 1)










## Overlaying interesting genes (DE or not)

# Surface markers (split form the main file for better customization)
Surf_mark_df <- data.frame(EpCAM = sct_CTC.obj[["SCT"]]@data[c("ENSG00000119888"), ],
                           PSMA = sct_CTC.obj[["SCT"]]@data[c("ENSG00000086205"), ],
                           Cluster = sct_CTC.obj$seurat_clusters,
                           RespGroup = sct_CTC.obj$RespGroup,
                           Cycles = sct_CTC.obj$RealCycle)

Surf_mark_df_long <- pivot_longer(data = Surf_mark_df, cols = c(EpCAM, PSMA), names_to = "Marker", values_to = "ExprLvl")

# EpCAM responder
EpCAM_expr_over_cyc_resp_p <- ggplot(data = dplyr::filter(Surf_mark_df_long, Marker == "EpCAM", RespGroup == "Responder"),
                                aes(x = factor(Cycles), y = ExprLvl, fill = factor(Cycles))) +
    #geom_violin(adjust = 1, trim = TRUE, scale = "width") +
    geom_boxplot(show.legend = FALSE, position = position_dodge(width = 1, preserve = "single")) +
    geom_dotplot(show.legend = FALSE, binaxis = "y", stackdir = "center", binwidth = 0.05, alpha = 0.5, position = position_dodge(width = 1, preserve = "total"),
                 dotsize = 2) +
    scale_fill_tron() +
    #stat_compare_means(label.y = 10, label.x = 1.2, method = "kruskal.test") +
    geom_pwc(method = "wilcox_test", label = "p = {p}", ref.group = "0") + #this geom is for pairwise comparison and is amazing!
    labs(x = expression("Treatment cycles"), y = expression("Log"[e]*" expression")) +
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
EpCAM_expr_over_cyc_nresp_p <- ggplot(data = dplyr::filter(Surf_mark_df_long, Marker == "EpCAM", RespGroup == "Nonresponder"),
                                      aes(x = factor(Cycles), y = ExprLvl, fill = factor(Cycles))) +
    #geom_violin(adjust = 1, trim = TRUE, scale = "width") +
    geom_boxplot(show.legend = FALSE, position = position_dodge(width = 1, preserve = "single")) +
    geom_dotplot(show.legend = FALSE, binaxis = "y", stackdir = "center", binwidth = 0.05, alpha = 0.5, position = position_dodge(width = 1, preserve = "total"),
                 dotsize = 2) +
    scale_fill_tron() +
    #stat_compare_means(label.y = 10, label.x = 1.2, method = "kruskal.test") +
    geom_pwc(method = "wilcox_test", label = "p = {p}", ref.group = "0") + #this geom is for pairwise comparison and is amazing!
    labs(x = expression("Treatment cycles"), y = expression("Log"[e]*" expression")) +
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
PSMA_expr_over_cyc_resp_p <- ggplot(data = dplyr::filter(Surf_mark_df_long, Marker == "PSMA", RespGroup == "Responder"),
                                    aes(x = factor(Cycles), y = ExprLvl, fill = factor(Cycles))) +
    #geom_violin(adjust = 1, trim = TRUE, scale = "width") +
    geom_boxplot(show.legend = FALSE, position = position_dodge(width = 1, preserve = "single")) +
    geom_dotplot(show.legend = FALSE, binaxis = "y", stackdir = "center", binwidth = 0.05, alpha = 0.5, position = position_dodge(width = 1, preserve = "total"),
                 dotsize = 2) +
    scale_fill_tron() +
    #stat_compare_means(label.y = 10, label.x = 1.2, method = "kruskal.test") +
    geom_pwc(method = "wilcox_test", label = "p = {p}", ref.group = "0") + #this geom is for pairwise comparison and is amazing!
    labs(x = expression("Treatment cycles"), y = expression("Log"[e]*" expression")) +
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
PSMA_expr_over_cyc_nresp_p <- ggplot(data = dplyr::filter(Surf_mark_df_long, Marker == "PSMA", RespGroup == "Nonresponder"),
                                    aes(x = factor(Cycles), y = ExprLvl, fill = factor(Cycles))) +
    #geom_violin(adjust = 1, trim = TRUE, scale = "width") +
    geom_boxplot(show.legend = FALSE, position = position_dodge(width = 1, preserve = "single")) +
    geom_dotplot(show.legend = FALSE, binaxis = "y", stackdir = "center", binwidth = 0.05, alpha = 0.5, position = position_dodge(width = 1, preserve = "total"),
                 dotsize = 2) +
    scale_fill_tron() +
    #stat_compare_means(label.y = 10, label.x = 1.2, method = "kruskal.test") +
    geom_pwc(method = "wilcox_test", label = "p = {p}", ref.group = "0") + #this geom is for pairwise comparison and is amazing!
    labs(x = expression("Treatment cycles"), y = expression("Log"[e]*" expression")) +
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


# Overlay interesting genes - PSMA(FOLH1)
# FOLH1, EPCAM ,SSTR2 ,TP53 ,GAPDH
FeaturePlot(sct_CTC.obj,
            reduction = "umap",
            features = c("ENSG00000086205", "ENSG00000119888", "ENSG00000180616", "ENSG00000141510", "ENSG00000111640"),
            order = TRUE,
            label = FALSE,
            repel = FALSE)


DotPlot(sct_CTC.obj,
        features = c("ENSG00000086205", "ENSG00000119888", "ENSG00000180616", "ENSG00000141510"))


RidgePlot(sct_CTC.obj,
            features = c("ENSG00000086205", #FOLH1
                        "ENSG00000119888", #EPCAM
                        "ENSG00000180616", #SSTR2
                        "ENSG00000141510", #TP53
                        "ENSG00000111640"),  #GAPDH
            cols = cluster_colors[c(4, 6, 7)])


# Transcriptomic profiling of single circulating tumor cells provides insight into human metastatic gastric
# cancer

# EMT related genes
# ZEB2, MEF2D, NFKB1A, GATA1, SERPINE1, ZEB1, TWIST1, TWIST2, SNAI1, SNAI2, TGFB1, TGFB2
RidgePlot(sct_CTC.obj,
            features = c("ENSG00000169554", #ZEB2
                         "ENSG00000116604", #MEF2D,
                         "ENSG00000109320", #NFKB1A
                         "ENSG00000179348", #GATA1
                         "ENSG00000106366", #SERPINE1
                         "ENSG00000092969"),  #TGFB2
                         cols = cluster_colors[c(4, 6, 7)])


#Epithelial markers
#EPCAM, KRT8, KRT10, KRT18, KRT19
RidgePlot(sct_CTC.obj,
        features = c("ENSG00000119888", #EPCAM
                     "ENSG00000170421", #KRT8
                     "ENSG00000186395", #KRT10
                     "ENSG00000111057", #KRT18
                     "ENSG00000171345"), #KRT19
                     cols =  cluster_colors[c(4, 6, 7)]) 

# Platelet related genes
# PF4, PPBP, ITGA2B, SPARC, ITGA2
RidgePlot(sct_CTC.obj,
            features = c("ENSG00000163737",
                         "ENSG00000163736",
                         "ENSG00000005961",
                         "ENSG00000113140",
                         "ENSG00000164171"),
                         cols =  cluster_colors[c(4, 6, 7)])


# CSC markers
# ALDH1A1, TACSTD2, CD44, CD24, EZH2,
RidgePlot(sct_CTC.obj,
            features = c("ENSG00000165092",
                        "ENSG00000184292",
                        "ENSG00000026508",
                        "ENSG00000272398",
                        "ENSG00000106462"),
                         cols =  cluster_colors[c(4, 6, 7)])


# Proliferation markers
# UBE2C, CCNB1, TOP2A, PLK1, CCND1
RidgePlot(sct_CTC.obj,
            features = c("ENSG00000175063",
                        "ENSG00000134057",
                        "ENSG00000131747",
                        "ENSG00000166851",
                        "ENSG00000110092"),
                         cols =  cluster_colors[c(4, 6, 7)])




Seurat::FeaturePlot(sct_CTC.obj,
                    features = c("ENSG00000119888",
                    "ENSG00000086205"))















### RNA velocity analysis
### NOTE: the velocity analysis will be done in python, as the available R package cannot be installed for the time being


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




















###Trajectory analysis with Monocle 3

##Prepare the seurat Object fot the analysis by re-running the sctransform, PCA and find neighbors

#SCTransform
sct_CTC.obj_TA <- SCTransform(filt_CTC.obj, assay = "RNA_seq", method = "glmGamPoi")


#an important thing to check is which PCs are responsible for the most variation. This is important
#for the downstream steps. One way to check this is to plot the PCs with an ElbowPlot
sct_CTC.obj_TA <- RunPCA(sct_CTC.obj_TA, assay = "SCT", approx = FALSE, verbose = FALSE) #approx = FALSE to run with normal vsd method
PC_weight_p <- ElbowPlot(sct_CTC.obj_TA, ndims = 50, reduction = "pca")
PC_weight_p2 <- PC_weight_p +
    labs(x = "PCs") +
    ggtitle("PC weights") +
    theme(plot.title = element_text(hjust = 0.5))
print(PC_weight_p2)
rm(PC_weight_p, PC_weight_p2)


#Finding neighbours and clustering
sct_CTC.obj_TA <- FindNeighbors(sct_CTC.obj_TA, reduction = "pca", dims = 1:35)


#Running dimensionality reduction (trying both UMAP and tSNE :) )
sct_CTC.obj_TA <- RunUMAP(sct_CTC.obj_TA, dims = 1:35, reduction = "pca", verbose = FALSE)




##Convert the seurat object to a monocle CDS (Cell Data Set) and do the rest with this package

#Converting the seurat object
sct_CTC.cds <- SeuratWrappers::as.cell_data_set(sct_CTC.obj_TA)


#as suggested by the monocle3 package, running the cluster_cells on the new cds object
sct_CTC.cds <- monocle3::cluster_cells(sct_CTC.cds, resolution = 0.2)

#Checking how the clustering looks like (adjust the resolution if needed)
monocle3::plot_cells(cds = sct_CTC.cds, show_trajectory_graph = FALSE, cell_size = 1)

monocle3::plot_cells(cds = sct_CTC.cds, show_trajectory_graph = FALSE, cell_size = 1,
                        color_cells_by = "partition")

#Next we will run the learn graph function



















##GSEA (Gene Set Enrichment Analysis)

#for the gene set enrichment analysis I will use the fgsea r package


#creating a ranked list of genes (a named vector basically) for the fgsa algorithm
#note: this will include all genes in the set, not just the DE genes.
Idents(sct_CTC.obj) <- "seurat_clusters"
GSEA_DE_lst <- FindMarkers(sct_CTC.obj, 
                           ident.1 = "0", 
                           ident.2 = "1", 
                           logfc.threshold = -Inf, 
                           min.pct = -Inf, 
                           min.diff.pct = -Inf, 
                           verbose = TRUE)


#Attaching gene names to the DE genes using the following fucntions
#This is the parallelized version of the matching function
#' @export 
#' @importFrom foreach %dopar%
    devtools::document()

geneName_matching_par = function (DE_lst, gene_lst, ID_col, name_col, n_cores = 2) {
    progressBar <- progress_bar$new(format = "Progress = (:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
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
    
    opts <- list(progress = progress)
    
    cores <- n_cores
    clust <- parallel::makeCluster(cores)
    doSNOW::registerDoSNOW(cl = clust)
    
    position <- c()
    out <- c()
    geneNames <- foreach::foreach(i = seq_len(nrow(DE_lst)), .options.snow = opts, .combine = c) %dopar% {
        position <- grep(rownames(DE_lst)[i], gene_lst[, ID_col])
        out[i] <- gene_lst[, name_col][position]
    }
    
    parallel::stopCluster(cl = clust)
    
    out_lst <- cbind(DE_lst, geneNames)
    original_name <- deparse(substitute(DE_lst))
    output_name <- paste0(original_name, "_", "geneNames")
    assign(output_name, out_lst, envir = .GlobalEnv)
    message("The gene name matched DE_cluster_0v1 list was assigned to a new object named:", "\n", output_name)
}
geneName_matching_par(GSEA_DE_lst, unif_gene_names, "ensembl_gene_id", "external_gene_name")

#NOTE: this version is less time intensive (runs under 5 minutes) as it runs the process on multiple cores in parallel, however
#it relies on on a foreach "loop" which is more complicated then a simple for loop


#note: the ranking here will be solely based on the p_adj value (as the sign() function will turn all log2FC into a negative, 0 or positive 1)
#if you don't want that (which we want to do here) you can omit the sign() function
GSEA_DE_lst_geneNames$Ranking <- sign(GSEA_DE_lst_geneNames$avg_log2FC) * (-log10(GSEA_DE_lst_geneNames$p_val_adj))


#note: the ranking here will take into account the magnitude of the log2 fold change, not only the the statistical significance and direction
GSEA_DE_lst_geneNames$Log_and_p_Ranking <- GSEA_DE_lst_geneNames$avg_log2FC * (-log10(GSEA_DE_lst_geneNames$p_val_adj))


#Adding an additional column with the expression change significance
GSEA_DE_lst_geneNames$Significance <- ifelse(GSEA_DE_lst_geneNames$avg_log2FC > 0 & GSEA_DE_lst_geneNames$p_val_adj < 0.05, "Upregulated",
                                                ifelse(GSEA_DE_lst_geneNames$avg_log2FC < 0 & GSEA_DE_lst_geneNames$p_val_adj < 0.05, "Downregulated", "Not significant")) #labeling genes for coloring


#to run run GSEA we should remove any NAs (not converted IDs) and duplicated entries
summary(is.na(GSEA_DE_lst_geneNames$geneNames))
GSEA_DE_lst_geneNames <- GSEA_DE_lst_geneNames[!is.na(GSEA_DE_lst_geneNames$geneNames), ]

#NOTE:checking for duplicated entries. In the case of duplicated entries I will keep the ones
#with the higher rank
summary(duplicated(GSEA_DE_lst_geneNames$geneNames))
GSEA_DE_lst_geneNames[duplicated(GSEA_DE_lst_geneNames$geneNames), ]

dup_geneName <- GSEA_DE_lst_geneNames$geneNames[duplicated(GSEA_DE_lst_geneNames$geneNames)]
dplyr::filter(GSEA_DE_lst_geneNames, geneNames %in% dup_geneName)

GSEA_DE_lst_geneNames <- GSEA_DE_lst_geneNames[-grep("ENSG00000285437", rownames(GSEA_DE_lst_geneNames)), ]
summary(duplicated(GSEA_DE_lst_geneNames$geneNames))


#creating a raked vector, which is the input format of the FGSEA and clusterProfiler package
GSEA_DE_ranked_lst <- GSEA_DE_lst_geneNames$Ranking
names(GSEA_DE_ranked_lst) <- GSEA_DE_lst_geneNames$geneNames
head(GSEA_DE_ranked_lst)

GSEA_DE_logP_ranked_lst <- GSEA_DE_lst_geneNames$Log_and_p_Ranking
names(GSEA_DE_logP_ranked_lst) <- GSEA_DE_lst_geneNames$geneNames
head(GSEA_DE_logP_ranked_lst)


#to check if there are any infinite values present or not
summary(GSEA_DE_ranked_lst)

summary(GSEA_DE_logP_ranked_lst)

#to check how the values are distributed, and how they look like
GSEA_DE_ranked_lst_ord <- sort(GSEA_DE_ranked_lst, decreasing = TRUE)
plot(GSEA_DE_ranked_lst_ord)

GSEA_DE_logP_ranked_lst_ord <- sort(GSEA_DE_logP_ranked_lst, decreasing = TRUE)
plot(GSEA_DE_logP_ranked_lst_ord)


#I will make a DE_ranked_df for ggplot
GSEA_DE_ranked_df <- GSEA_DE_lst_geneNames[, c(6, 7)]
GSEA_DE_ranked_df_ord <- dplyr::arrange(GSEA_DE_ranked_df, desc(Ranking))
head(GSEA_DE_ranked_df_ord)

GSEA_DE_logP_ranked_df <- GSEA_DE_lst_geneNames[, c(6, 7)]
GSEA_DE_logP_ranked_df_ord <- dplyr::arrange(GSEA_DE_logP_ranked_df, desc(Log_and_p_Ranking))
head(GSEA_DE_logP_ranked_df_ord)



#Cluster 0 top 50 upregulated
ggplot(GSEA_DE_ranked_df_ord[1:50, ],
       aes(x = geneNames, y = Ranking)) +
    geom_point() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(GSEA_DE_logP_ranked_df_ord[1:50, ],
       aes(x = geneNames, y = Log_and_p_Ranking)) +
    geom_point() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



#cluster 0 top 50 downregulated
ggplot(GSEA_DE_ranked_df_ord[nrow(GSEA_DE_ranked_df_ord):(nrow(GSEA_DE_ranked_df_ord) - 50), ],
       aes(x = geneNames, y = Ranking)) +
    geom_point() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(GSEA_DE_logP_ranked_df_ord[nrow(GSEA_DE_logP_ranked_df_ord):(nrow(GSEA_DE_logP_ranked_df_ord) - 50), ],
       aes(x = geneNames, y = Log_and_p_Ranking)) +
    geom_point() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#Prepare the .gmt (reference sets) and other files/data for the analysis
gmt_files <- paste0(getwd(), "/", "GSEA_ref_sets")


#List the reference pathways
pathways <- list.files(gmt_files)

#this function will read and trim the refrence pathway files (.gmt) from the proper folder,
#and assigns them into new objects in the .GlobalEnv
load_trim_gmt = function(path, DE_list) {
    #listing the files to load
    file_names <- list.files(path = path)
    
    #creates a progress bar for the file loading
    progressBar <- progress_bar$new(format = "Reading and trimming gmt - (:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                    total = length(file_names),
                                    complete = "=",
                                    incomplete = "-",
                                    current = ">",
                                    clear = FALSE,
                                    width = 100)
    
    #this loop, will ensure loading and processing of the gmt files in a sequential order
    tmp_gmt <- list()
    trimmed_gmt <- list()
    name_to_assign <- c()
    for (e in seq_along(file_names)) {
        #initiate progess bar
        progressBar$tick()
        message("\n" ,"processing .gmt file: ", e)
        
        #creating the name for the new pathway object
        name_to_assign <- file_names[e]
        
        #loading a gmt file into a temporary object
        tmp_gmt <- gmtPathways(paste0(gmt_files, "/", file_names[e]))
        
        #this lapply will 
        trimmed_gmt <- lapply(tmp_gmt, function(x)
            x[x %in% names(DE_ranked_list_ord)])
    
        #assigning the trimmed files into new objects in the .GlobalEnv
        assign(name_to_assign, trimmed_gmt, envir = .GlobalEnv)
        message("A reference pathway (.gmt) was assigned to a new object:", name_to_assign)
        
        
    }
    
    #creating a trimmed file to assign (which hopefully allow me assign all files one by one)
    gmt_to_assign <- trimmed_gmt
    
    
    
}
#load_trim_gmt(gmt_files, DE_ranked_list_ord)


#this function should load the refeence pathways (.gmt files) WItHOUT trimming
#NOTE: after talking with Jonathan and checking more references it seems that the pre trimming
#of the refence pathway list based on my gene list is not a commonly used (or good) practice
load_gmt = function(path, DE_list) {
    #listing the files to load
    file_names <- list.files(path = path)
    
    #creates a progress bar for the file loading
    progressBar <- progress_bar$new(format = "Reading and trimming gmt - (:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                    total = length(file_names),
                                    complete = "=",
                                    incomplete = "-",
                                    current = ">",
                                    clear = FALSE,
                                    width = 100)
    
    #this loop, will ensure loading and processing of the gmt files in a sequential order
    tmp_gmt <- list()
    name_to_assign <- c()
    for (e in seq_along(file_names)) {
        #initiate progess bar
        progressBar$tick()
        message("\n" ,"loading .gmt file: ", e)
        
        #creating the name for the new pathway object
        name_to_assign <- file_names[e]
        
        #loading a gmt file into a temporary object
        tmp_gmt <- gmtPathways(paste0(gmt_files, "/", file_names[e]))
        
        #assigning the trimmed files into new objects in the .GlobalEnv
        assign(name_to_assign, tmp_gmt, envir = .GlobalEnv)
        message("A reference pathway (.gmt) was assigned to a new object:", name_to_assign)
        
        
    }
    
    
    
    
    
}
#load_gmt(gmt_files, DE_ranked_list_ord)


#to bulk read in .gmt files for the clusterProfiler package
read_gmt = function(path) {
    #listing the files to load
    file_names <- list.files(path = path)
    
    #creates a progress bar for the file loading
    progressBar <- progress_bar$new(format = "Reading gmt - (:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                    total = length(file_names),
                                    complete = "=",
                                    incomplete = "-",
                                    current = ">",
                                    clear = FALSE,
                                    width = 100)
    
    #this loop, will ensure loading and processing of the gmt files in a sequential order
    tmp_gmt <- list()
    name_to_assign <- c()
    for (e in seq_along(file_names)) {
        #initiate progess bar
        progressBar$tick()
        message("\n" ,"loading .gmt file: ", e)
        
        #creating the name for the new pathway object
        name_to_assign <- file_names[e]
        
        #loading a gmt file into a temporary object
        tmp_gmt <- clusterProfiler::read.gmt(paste0(gmt_files, "/", file_names[e]))
        
        #assigning the trimmed files into new objects in the .GlobalEnv
        assign(name_to_assign, tmp_gmt, envir = .GlobalEnv)
        message("A reference pathway (.gmt) was assigned to a new object:", name_to_assign)
        
        
    }
    
    
    
    
    
}
read_gmt(gmt_files)


#This function is a wrapper function, which will do the GSEA run, assigns the result into a new object,
#saves the result as a .csv file, prints the GSEA plots and saves them as .png files
run_GSEA = function(random_seed = 1234,
                    verbose = TRUE,
                    DE_gene_list,
                    min_size = 5,
                    max_size = 500,
                    p_val_cutoff = 0.05,
                    p_val_adjust_method = "BH",
                    GSEA_method = "fgsea",
                    nPermSimple = 10000,
                    ref_gene_lst,
                    GSEA_csv_file_name,
                    GSEA_result_object_name,
                    plot_title,
                    plot_file_name,
                    save_folder = "GSEA_results") {
    
    #An if statement to check if a custom save folder was requested. If yes, it creates the folder and subfolder and saves the results there.
    if (save_folder == "GSEA_results") {
        #do nothing
    } else if (file.exists(paste0(script_version, "/", save_folder)) == TRUE) {
        message("The requested custom save folder ", save_folder, "already exists. Files will be saved there")
    } else {
        dir.create(path = paste0(script_version, "/", save_folder), recursive = TRUE)
        dir.create(path = paste0(script_version, "/", save_folder, "/", "Plots"), recursive = TRUE)
        
        message("A custom save folder under the name: ", save_folder, "was requested and created.", "\n", 
                "the GSEA results will be saved into the new folder.")
    }
    
    #Running GSEA using the clusterProfiler package
    #the function uses permutation testing to adjust the p-values. as this is based on a random number, for the sake of
    #reproducibility, one should set the seed
    set.seed(random_seed)
    
    #Run GSEA on selected dataset and save the result as a .csv and assign the result as an object in the .GlobalEnv
    GSEA_result <- clusterProfiler::GSEA(geneList = DE_gene_list,
                              minGSSize = min_size,
                              maxGSSize = max_size,
                              pvalueCutoff = p_val_cutoff,
                              pAdjustMethod = p_val_adjust_method,
                              TERM2GENE = ref_gene_lst,
                              by = GSEA_method,
                              verbose = verbose,
                              nPermSimple = nPermSimple)
    
    GSEA_result_df <- as.data.frame(GSEA_result)
    write.csv(GSEA_result_df, file = paste0(script_version, "/", save_folder, "/", GSEA_csv_file_name, script_version, ".csv"))
    message("A GSEA run was carried out using the gene list:", deparse(substitute(DE_gene_list)), "\n",
            "and the reference gene list:", deparse(substitute(ref_gene_lst)), "\n",
            "The result was saved as a .csv file:", GSEA_csv_file_name, "\n",
            "The result is assigned to the object:", GSEA_result_object_name)
    assign(GSEA_result_object_name, GSEA_result, envir = .GlobalEnv)
    
    #initiating progress bar for the plotting loop
    progressBar <- progress_bar$new(format = "Plotting GSEA - (:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                    total = nrow(GSEA_result@result),
                                    complete = "=",
                                    incomplete = "-",
                                    current = ">",
                                    clear = FALSE,
                                    width = 100)
    
    #Initiating the plotting with a messgae
    message("Plotting and saving the GSEA results:")
    
    #this loop will print and save the GSEA plot as .png files
    for (e in seq_len(nrow(GSEA_result@result))) {
        #initiating the progress bar
        progressBar$tick()
        
        #call the gseaplot function to plot the gsea result
        gsea_res_p <- clusterProfiler::gseaplot(x = GSEA_result,
                               geneSetID = GSEA_result@result[e,1],
                               title = paste0(plot_title, "- ", GSEA_result@result[e,1]))
        #visualize the resulting plot
        print(gsea_res_p)
        
        #saving the resulting plot as a .png
        ggplot2::ggsave(filename = paste0(plot_file_name, "-", GSEA_result@result[e,1], script_version, ".png"), gsea_res_p,
               device = "png", path = paste0(script_version, "/", save_folder, "/", "Plots", "/"),
               width = 3500, height = 3000, units = "px")
        
    }
        message("The resulting GSEA plots were printed and saved into the folder: ", paste0(script_version, "/", save_folder, "/", "Plots"))
    
    
    
}


#I modified the run_GSEA function to take a string containing multiple gene sets and run the analysis along the list one by one.
run_multi_GSEA = function(reference_lst,
                          random_seed_multi = 1234,
                          DE_gene_list_multi,
                          min_size_multi = 5,
                          max_size_multi = 500,
                          p_val_cutoff_multi = 0.05,
                          p_val_adjust_method_multi = "BH",
                          GSEA_method_multi = "fgsea",
                          nPermSimple_multi = 10000,
                          GSEA_csv_file_name_multi,
                          GSEA_result_object_name_multi,
                          plot_file_name_multi,
                          save_folder_multi = "GSEA_results") {
    
    
    #extracting parts of the reference pathway names to be used for result names
    pathway_names <- stringr::str_extract(reference_lst, rebus::one_or_more(ASCII_ALNUM) %R% "." %R% rebus::one_or_more(ASCII_ALNUM) %R% "." %R% rebus::one_or_more(ASCII_ALNUM))
    
    #creating the progress bar
    progressBar <- progress_bar$new(format = "Running GSEA - (:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                    total = length(reference_lst),
                                    complete = "=",
                                    incomplete = "-",
                                    current = ">",
                                    clear = FALSE,
                                    width = 100)
    
    #This for loop should run the GSEA function over the refenecence gene lists one by one
    for (e in seq_along(reference_lst)) {
        #initiating the progress bar
        progressBar$tick()
        message("\n" ,"Initiating GSEA on set ", reference_lst[e])
        
        #running the the GSEA wrapper function with the tryCatch function wrapped around it
        #NOTE: the tryCatch function will carry on executing the loop even if the GSEA function stops due to the lack of
        #enrichment. otherwise the whole loop stops moving on.
        try(run_GSEA(DE_gene_list = DE_gene_list_multi,
                 random_seed = random_seed_multi,
                 min_size = min_size_multi,
                 max_size = max_size_multi,
                 p_val_cutoff = p_val_cutoff_multi,
                 p_val_adjust_method = p_val_adjust_method_multi,
                 GSEA_method = GSEA_method_multi,
                 nPermSimple = nPermSimple_multi,
                 ref_gene_lst = get(reference_lst[e]),
                 GSEA_csv_file_name = paste0(pathway_names[e], GSEA_csv_file_name_multi),
                 GSEA_result_object_name = paste0(pathway_names[e], GSEA_result_object_name_multi),
                 plot_title = paste0("GSEA", pathway_names[e], "set"),
                 plot_file_name = paste0(pathway_names[e], plot_file_name_multi),
                 save_folder = save_folder_multi))
        
    }
    
    message("The multi GSEA run was successfully completed!")
    
    
}

run_multi_GSEA(reference_lst = pathways,
               DE_gene_list_multi = GSEA_DE_ranked_lst_ord,
               GSEA_csv_file_name_multi = "_GSEA_result",
               GSEA_result_object_name_multi = "_GSEA_object",
               plot_file_name_multi = "_GSEA_plot")

run_multi_GSEA(reference_lst = pathways,
               DE_gene_list_multi = GSEA_DE_logP_ranked_lst_ord,
               GSEA_csv_file_name_multi = "_GSEA_logP_result",
               GSEA_result_object_name_multi = "_GSEA_logP_object",
               plot_file_name_multi = "_GSEA_logP_plot",
               save_folder_multi = "GSEA_logP_results")



##SPIA (Single Pathway Impact Analysis)








































































































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
#seurat object and substracetd the rows preent in the DE gene list
#extracting the transformed count matrix from the seurat object for the NetGSA analysis
sct_count_matrix <- GetAssayData(sct_CTC.obj, slot = "data")

#adjusting the count matrix to contain the DE_lst genes
net_matrix <- sct_count_matrix[rownames(sct_count_matrix) %in% rownames(DE_Cluster_0v1_geneNames), ]
rm(sct_count_matrix)







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

#setting up the frouping variable, in this case the cluster numbers
net_group <- as.numeric(sct_CTC.obj@meta.data$seurat_clusters)
table(net_group)


#getting the edges
edges <- netgsa::obtainEdgeList(genes = rownames(net_matrix2), databases = c("kegg", "reactome"))

net_adjMat <- prepareAdjMat(x = net_mat2, group = net_group)











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

























































































































