### This is a first attempt on doing RNA velocity anylis in Python
### The script is written in pyCharm following a great tutorial on github
### source: https://smorabit.github.io/tutorials/8_velocyto/


## Import necessary packages for Seurat data reconstruction
import matplotlib
from matplotlib.pyplot import xlabel, xscale, xticks
import scanpy
import anndata
import numpy
import pandas
import scipy.stats
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import os
import loompy



## Get and set the current working directory
print(os.getcwd())
os.chdir("/home/gabi/RNA_velo_analysis_new")



## Seurat data reconstruction

# source data path
source_data_path = os.getcwd()

# load sparse matrix
sparse_matrix = io.mmread(source_data_path + "/" + "CTC_GeneCountMatrix_gNames_for_scVelo_V3.9.mtx")
print(sparse_matrix)

sparse_matrix_IDs = io.mmread(source_data_path + "/" + "CTC_GeneCountMatrix_IDs_for_scVelo_V3.9.mtx")
print(sparse_matrix_IDs)


# create an anndata object (as far as I gather this is similar to a seurat object)
adata = anndata.AnnData(X=sparse_matrix.transpose().tocsr())
print(adata)

adata_IDs = anndata.AnnData(X=sparse_matrix_IDs.transpose().tocsr())
print(adata_IDs)


# load the seurat metadata
CTC_metadata = pandas.read_csv(source_data_path + "/" + "CTC_metadata_for_scVelo_V3.9.csv")
CTC_metadata.head()


# load the gene names / ENSEMBL IDs
with open(source_data_path + "/" + "Gene_names_for_scVelo_V3.9.csv", "r") as file:
    gene_names = file.read().splitlines()
print(gene_names)

with open(source_data_path + "/" + "ENSEMBL_IDs_for_scVelo_V3.9.csv", "r") as file:
    ENSEMBL_IDs = file.read().splitlines()
print(ENSEMBL_IDs)


# set up the anndata observations and index observations by barcodes, and variables by gene names / ENSEMBL IDs
adata.obs = CTC_metadata
adata.obs.index = adata.obs["barcode"]
adata.var.index = gene_names

adata_IDs.obs = CTC_metadata
adata_IDs.obs.index = adata_IDs.obs["barcode"]
adata_IDs.var.index = ENSEMBL_IDs


## I was thinking maybe I could add the gene names together with the ENSEMBL IDs to the anndata object for the analysis, os I will try that now

# load the .csv containing the ENSEMBL IDs with missing gene names
missing_names = pandas.read_csv("/mnt/c/WSL_shared/missing_name_to_ID_V3.9.csv")
print(missing_names)

# This for loop lists out the IDs beloning to the missing gene names from the actual ID list
IDs_to_remove = []
for item in ENSEMBL_IDs:
    if missing_names["ensembl_gene_id"].isin([item]).any():
        IDs_to_remove.append(item)
print(IDs_to_remove)

# This for loop creates a new list condtaining only the IDs with gene names to them
ENSEMBL_IDs_trim = []
for item in ENSEMBL_IDs:
    if item not in IDs_to_remove:
        ENSEMBL_IDs_trim.append(item)
print(ENSEMBL_IDs_trim)

# creating a new integrated adata object and assining it the metadata to .obs and the IDs and gene names to .var
adata_integr = adata

adata.obs = CTC_metadata
adata.obs.index = adata.obs["barcode"]
adata.var.index = ENSEMBL_IDs_trim
adata.var["geneNames"] = gene_names


# load the dimesionality reduction data
pca = pandas.read_csv(source_data_path + "/" + "PCA_data_for_scVelo_V3.9.csv")
pca.head()


# set PCA and UMAP
adata_integr.obsm["X_pca"] = pca.to_numpy()
adata_integr.obsm["X_umap"] = numpy.vstack((adata.obs["UMAP_1"].to_numpy(), adata.obs["UMAP_2"].to_numpy())).T


# Plot a UMAP test, colored by the sample ID
scanpy.pl.umap(adata=adata_integr, color=["seurat_clusters"], frameon=False, save=True)


# Save the dataset as anndata format
#adata_integr.write("CTC_scSeq_velocity_anndata_obj.h5ad")


# Reload the saved anndata object
#adata = scanpy.read_h5ad("CTC_scSeq_velocity_anndata_obj.h5ad")



## Importing the necessary packages for the scVelo analysis
import scvelo
import cellrank



## scVelo analysis

# adjust a couple of settings
scvelo.settings.verbosity = 3  # this setting will show show errors(0), warnings(1), info(2), hints(3)
scvelo.settings.set_figure_params(style=None, facecolor="white", dpi=300, frameon=True, transparent=True, format="png")  # this will set the plot parameters
cellrank.settings.verbosity = 2
matplotlib.rcParams["text.color"] = "black"
matplotlib.rcParams["axes.labelcolor"] = "black"

# Load the saved anndata object if you did not do it yet
def load_anndata_object(object_name, file_name):

    # check if the data is already loaded and if not it will load it as a new object
    if object_name in globals():
        print("The Anndata object:" + file_name + "is already loaded" +
              "\n" +
              "under the name:" + object_name)
    else:
        print("The Anndata object is not loaded yet. Loading now..." +
              "\n")
        anndata_object = scanpy.read_h5ad(file_name)
    return anndata_object


#ctc_scSeq_an_obj = load_anndata_object(object_name="ctc_scSeq_an_obj", file_name="CTC_scSeq_velocity_anndata_obj.h5ad")
#ctc_scSeq_ID_an_obj = load_anndata_object(object_name="ctc_scSeq_ID_an_obj", file_name="CTC_scSeq_velocity_ENSEMBL_IDs_anndata_obj.h5ad")

ctc_scSeq_adata_integr = adata_integr

# load the loom file containing the spliced/unspliced matrices for the whole det
splice_matrix = scvelo.read("CTC_scSeq_velocity_splice_matrices.loom", cache=True, validate=False)



# The splice matrix will have to be reformatted, trimmed and cleaned up a bit so I gonna do that here
splice_matrix_copy = splice_matrix
splice_matrix_copy.obs
splice_matrix_copy.var

# remove the unncesarry parts from the sample names, so they match the exproted seurat object
splice_matrix_copy.obs.index = splice_matrix_copy.obs.index.str.replace("CTC_scRNA_seq_velocity_analysis:sorted_lane1", "").str.replace("_AAC2HT2HVAligned.sortedByCoord.out.bam", "")
splice_matrix_copy.obs

# trimming the loom file so it is equal to the exported seurat object
# Trimming the sample names (cols)
splice_obs_indexes = splice_matrix_copy.obs.index
sample_mask = splice_obs_indexes.isin(ctc_scSeq_adata_integr.obs.index)
splice_matrix_trim = splice_matrix_copy[sample_mask]
splice_matrix_trim.obs

# Trimming the gene IDs (rows)
splice_matrix_IDs = splice_matrix_trim.var["Accession"]
ensID_mask = splice_matrix_IDs.isin(ctc_scSeq_adata_integr.var.index)
splice_matrix_trim = splice_matrix_trim[: ,ensID_mask]
#NOTE: after trimming the loom file using the ENSEMBL IDs from the exported seurat object I noticed that the number of genes
#are not the same as in the seurat object. This should not be the case but I can't figure out what is causig  it.
#If I want to merge the two file I will have to remove the aforementioned genes/IDs from the seurat object as well.

# Trimming the seurat object for the merge
final_ensIDs = splice_matrix_IDs[splice_matrix_IDs.isin(ctc_scSeq_adata_integr.var.index)]
final_ensIDs_mask = ctc_scSeq_adata_integr.var.index.isin(final_ensIDs)
ctc_scSeq_adata_integr_trim = ctc_scSeq_adata_integr[: ,final_ensIDs_mask]


# Merge the splice data with the anndata object
ctc_scSeq_adata_integr_trim_merged = scvelo.utils.merge(ctc_scSeq_adata_integr_trim, splice_matrix_trim)
# Note: if you get an error stating: UnboundLocalError: cannot access local variable 'id_length' where it is not associated with a value
# run the following lines ro clean up the names in both files (no specifics were given to this fix), then re run the merge
scvelo.utils.clean_obs_names(ctc_scSeq_an_obj)
scvelo.utils.clean_obs_names(splice_matrix)


# Save the merged anndata object for later use
ctc_scSeq_adata_integr_trim_merged.write("CTC_scSeq_integrated_merged_dataframe_with_IDs_geneNames.h5ad")


# load the integrated and prepared anndata object for analysis (if not loaded yet)
ctc_scSeq_adata_integr_trim_merged = scanpy.read_h5ad("CTC_scSeq_integrated_merged_dataframe_with_IDs_geneNames.h5ad")


# make a test umap plot
scanpy.pl.umap(adata=ctc_scSeq_adata_integr_trim_merged, frameon=False, color="seurat_clusters", legend_loc="on data",
               title="RNA velocity analysis - trial plot", save="RNA_vel_trial.png")


# Computing RNA velocity
ctc_scSeq_adata_integr_trim_merged.obs["seurat_clusters"] = ctc_scSeq_adata_integr_trim_merged.obs["seurat_clusters"].astype("category") # First I will change the clusters into factors
scvelo.pl.proportions(ctc_scSeq_adata_integr_trim_merged, groupby="seurat_clusters")




## Pre-processing step
scvelo.pp.filter_and_normalize(ctc_scSeq_adata_integr_trim_merged)
scvelo.pp.moments(ctc_scSeq_adata_integr_trim_merged, n_pcs=30, n_neighbors=30)




## Comupte stohastic velocities
scvelo.tl.velocity(ctc_scSeq_adata_integr_trim_merged, mode="stochastic")
scvelo.tl.velocity_graph(ctc_scSeq_adata_integr_trim_merged)




## Visualize stohastic velocities

#visualize the velocity fields
scvelo.pl.velocity_embedding(ctc_scSeq_adata_integr_trim_merged, basis="umap", color="seurat_clusters", alpha=1, frameon=False, save="velocity_embedding_plot.png")
scvelo.pl.velocity_embedding_grid(ctc_scSeq_adata_integr_trim_merged, basis="umap", color="seurat_clusters", alpha=1, save="ctc_embedding_grid.png",
                                  title="", scale=1, size=300)
scvelo.pl.velocity_embedding_stream(ctc_scSeq_adata_integr_trim_merged, basis="umap", color="seurat_clusters", save="ctc_embedding_stream.png",
                                  title="", size=300, alpha=1)
scvelo.pl.velocity(ctc_scSeq_ID_merged, var_names=["ENSG00000086205"])
scvelo.tl.rank_velocity_genes(ctc_scSeq_ID_merged, groupby="seurat_clusters")




## Dynamic modelling

# Calculate dynamic velocities
scvelo.tl.recover_dynamics(ctc_scSeq_adata_integr_trim_merged)
scvelo.tl.velocity(ctc_scSeq_adata_integr_trim_merged, mode="dynamical")
scvelo.tl.velocity_graph(ctc_scSeq_adata_integr_trim_merged)
scvelo.tl.recover_latent_time(ctc_scSeq_adata_integr_trim_merged)

scvelo.tl.rank_dynamical_genes(ctc_scSeq_adata_integr_trim_merged, groupby="seurat_clusters")
scvelo.tl.rank_velocity_genes(ctc_scSeq_adata_integr_trim_merged, vkey="velocity", groupby="seurat_clusters")

## Visualize dynamic velocities

# Set plot parameters
scvelo.settings.set_figure_params(style="scvelo", facecolor="white", dpi=300, frameon=False, transparent=False, format="png", figsize=[10, 5])  # this will set the plot parameters
matplotlib.rcParams["text.color"] = "black"
matplotlib.rcParams["axes.labelcolor"] = "black"
matplotlib.rcParams["xtick.color"] = "black"
matplotlib.rcParams["ytick.color"] = "black"


# Plot velocity embeddings on different gorups
scvelo.pl.velocity_embedding_stream(ctc_scSeq_adata_integr_trim_merged, basis="umap", save="ctc_embedding_stream_dynamic_clusters.png",
                                    title="Velocity streams - clusters", 
                                    legend_loc="right margin",
                                    color="seurat_clusters", 
                                    size=300,
                                    alpha=1,
                                    colorbar=True,
                                    min_mass=2)

scvelo.pl.velocity_embedding_stream(ctc_scSeq_adata_integr_trim_merged, basis="umap", save="ctc_embedding_stream_dynamic_cycles.png",
                                    title="Velocity streams - treatment cycles", 
                                    legend_loc="right margin",
                                    color="RealCycle", 
                                    size=300,
                                    alpha=1,
                                    colorbar=True,
                                    min_mass=2)

scvelo.pl.velocity_embedding_stream(ctc_scSeq_adata_integr_trim_merged, basis="umap", save="ctc_embedding_stream_dynamic_latent_time.png",
                                    title="Velocity streams - latent time", 
                                    legend_loc="right margin",
                                    color="latent_time", 
                                    size=300,
                                    alpha=1,
                                    colorbar=True,
                                    min_mass=2)




## Calculate kinetic rate parameters

# Create a new dataframe contaiing the velocity genes
velocity_genes_df = ctc_scSeq_adata_integr_trim_merged.var
velocity_genes_df.head
velocity_genes_df = velocity_genes_df[(velocity_genes_df["fit_likelihood"] > 0.1) & velocity_genes_df["velocity_genes"] == True]
velocity_genes_df.head


kwargs = dict(xscale="log", fontsize=10)
with scvelo.GridSpec(ncols=3) as plot:
    plot.hist(velocity_genes_df["fit_alpha"], xlabel="Transcription rate", **kwargs)
    plot.hist(velocity_genes_df["fit_beta"] * velocity_genes_df["fit_scaling"], xlabel="Splicing rate", **kwargs)
    plot.hist(velocity_genes_df["fit_gamma"], xlabel="Degradation rate", **kwargs)

scvelo.get_df(ctc_scSeq_adata_integr_trim_merged, "fit*", dropna=True).head()


# Looking at the top genes
top_genes = ctc_scSeq_adata_integr_trim_merged.var["fit_likelihood"].sort_values(ascending=False).index.dropna()[:300]
print(top_genes)

scvelo.pl.heatmap(ctc_scSeq_adata_integr_trim_merged, var_names=top_genes, sortby="latent_time", col_colors="seurat_clusters", n_convolve=50)

top_genes = ctc_scSeq_adata_integr_trim_merged.var["fit_likelihood"].sort_values(ascending=False).dropna().index
scvelo.pl.scatter(ctc_scSeq_adata_integr_trim_merged, basis=top_genes[:5], ncols=5, frameon=False, color="seurat_clusters", )
scvelo.pl.scatter(ctc_scSeq_adata_integr_trim_merged, x="latent_time", y=top_genes[:5], frameon=False, color="seurat_clusters")

df = scvelo.get_df(ctc_scSeq_adata_integr_trim_merged, "rank_dynamical_genes/names")

for clusters in ["0", "1"]:
    scvelo.pl.scatter(ctc_scSeq_adata_integr_trim_merged, df[clusters][:5], ylabel=clusters, frameon=False)


## NOTE:
## the top gene ENSEMBL IDs are all linked to unidentified genes, indicating that his analysi cannot be interpreted
## therefore I will stop pursuing this angle :(


















































# load the dimesionality reduction data
pca = pandas.read_csv(source_data_path + "/" + "PCA_data_for_scVelo_V3.9.csv")
pca.head()


# set PCA and UMAP
adata.obsm["X_pca"] = pca.to_numpy()
adata.obsm["X_umap"] = numpy.vstack((adata.obs["UMAP_1"].to_numpy(), adata.obs["UMAP_2"].to_numpy())).T

adata_IDs.obsm["X_pca"] = pca.to_numpy()
adata_IDs.obsm["X_umap"] = numpy.vstack((adata_IDs.obs["UMAP_1"].to_numpy(), adata_IDs.obs["UMAP_2"].to_numpy())).T


# Plot a UMAP test, colored by the sample ID
scanpy.pl.umap(adata=adata, color=["seurat_clusters"], frameon=False, save=True)

scanpy.pl.umap(adata=adata_IDs, color=["seurat_clusters"], frameon=False, save=True)


# Save the dataset as anndata format
adata.write("CTC_scSeq_velocity_anndata_obj.h5ad")

adata_IDs.write("CTC_scSeq_velocity_ENSEMBL_IDs_anndata_obj.h5ad")

# Reload the saved anndata object
adata = scanpy.read_h5ad("CTC_scSeq_velocity_anndata_obj.h5ad")



## Importing the necessary packages for the scVelo analysis
import scvelo
import cellrank



## scVelo analysis

# adjust a couple of settings
scvelo.settings.verbosity = 3  # this setting will show show errors(0), warnings(1), info(2), hints(3)
scvelo.settings.set_figure_params(style="scvelo", facecolor="white", dpi=300, frameon=False, transparent=True)  # this will set the plot parameters
cellrank.settings.verbosity = 2


# Load the saved anndata object if you did not do it yet
def load_anndata_object(object_name, file_name):

    # check if the data is already loaded and if not it will load it as a new object
    if object_name in globals():
        print("The Anndata object:" + file_name + "is already loaded" +
              "\n" +
              "under the name:" + object_name)
    else:
        print("The Anndata object is not loaded yet. Loading now..." +
              "\n")
        anndata_object = scanpy.read_h5ad(file_name)
    return anndata_object


ctc_scSeq_an_obj = load_anndata_object(object_name="ctc_scSeq_an_obj", file_name="CTC_scSeq_velocity_anndata_obj.h5ad")
ctc_scSeq_ID_an_obj = load_anndata_object(object_name="ctc_scSeq_ID_an_obj", file_name="CTC_scSeq_velocity_ENSEMBL_IDs_anndata_obj.h5ad")



# load the loom file containing the spliced/unspliced matrices for the whole det
splice_matrix = scvelo.read("CTC_scSeq_velocity_splice_matrices.loom", cache=True, validate=False)



# The splice matrix will have to be reformatted, trimmed and cleaned up a bit so I gonna do that here
splice_matrix_copy = splice_matrix
splice_matrix_copy.obs
splice_matrix_copy.var

# remove the unncesarry parts from the sample names, so they match the exproted seurat object
splice_matrix_copy.obs.index = splice_matrix_copy.obs.index.str.replace("CTC_scRNA_seq_velocity_analysis:sorted_lane1", "").str.replace("_AAC2HT2HVAligned.sortedByCoord.out.bam", "")
splice_matrix_copy.obs

# trimming the loom file so it is equal to the exported seurat object
# Trimming the sample names (cols)
splice_obs_indexes = splice_matrix_copy.obs.index
sample_mask = splice_obs_indexes.isin(ctc_scSeq_ID_an_obj.obs.index)
splice_matrix_trim = splice_matrix_copy[sample_mask]
splice_matrix_trim.obs

# Trimming the gene IDs (rows)
splice_matrix_IDs = splice_matrix_trim.var["Accession"]
ensID_mask = splice_matrix_IDs.isin(ctc_scSeq_ID_an_obj.var.index)
splice_matrix_trim = splice_matrix_trim[: ,ensID_mask]
#NOTE: after trimming the loom file using the ENSEMBL IDs from the exported seurat object I noticed that the number of genes
#are not the same as in the seurat object. This should not be the case but I can't figure out what is causig  it.
#If I want to merge the two file I will have to remove the aforementioned genes/IDs from the seurat object as well.

# Trimming the seurat object for the merge
final_ensIDs = splice_matrix_IDs[splice_matrix_IDs.isin(ctc_scSeq_ID_an_obj.var.index)]
final_ensIDs_mask = ctc_scSeq_ID_an_obj.var.index.isin(final_ensIDs)
ctc_scSeq_ID_trim = ctc_scSeq_ID_an_obj[: ,final_ensIDs_mask]


# Merge the splice data with the anndata object
ctc_scSeq_ID_merged = scvelo.utils.merge(ctc_scSeq_ID_trim, splice_matrix_trim)
# Note: if you get an error stating: UnboundLocalError: cannot access local variable 'id_length' where it is not associated with a value
# run the following lines ro clean up the names in both files (no specifics were given to this fix), then re run the merge
scvelo.utils.clean_obs_names(ctc_scSeq_an_obj)
scvelo.utils.clean_obs_names(splice_matrix)


# Save the merged anndata object for later use
ctc_scSeq_ID_merged.write("CTC_scSeq_Merged_dataframe_with_IDs.h5ad")


# make a test umap plot
scanpy.pl.umap(adata=ctc_scSeq_ID_merged, frameon=False, color="seurat_clusters", legend_loc="on data",
               title="RNA velocity analysis - trial plot", save="RNA_vel_trial.pdf")

# Computing RNA velocity
ctc_scSeq_ID_merged.obs["seurat_clusters"] = ctc_scSeq_ID_merged.obs["seurat_clusters"].astype("category") # First I will change the clusters into factors
scvelo.pl.proportions(ctc_scSeq_ID_merged, groupby="seurat_clusters")


# Pre-processing step
scvelo.pp.filter_and_normalize(ctc_scSeq_ID_merged, n_top_genes=3000)
scvelo.pp.moments(ctc_scSeq_ID_merged, n_pcs=30, n_neighbors=30)


# Comupte velocity
scvelo.tl.velocity(ctc_scSeq_ID_merged, mode="stochastic")
scvelo.tl.velocity_graph(ctc_scSeq_ID_merged)




## Visualize velocities

#visualize the velocity fields
scvelo.pl.velocity_embedding(ctc_scSeq_ID_merged, basis="umap", color="seurat_clusters", alpha=1, frameon=False, save="velocity_embedding_plot.pdf")
scvelo.pl.velocity_embedding_grid(ctc_scSeq_ID_merged, basis="umap", color="seurat_clusters", alpha=1, save="ctc_embedding_grid.pdf",
                                  title="", scale=1, size=300)
scvelo.pl.velocity_embedding_stream(ctc_scSeq_ID_merged, basis="umap", color="seurat_clusters", save="ctc_embedding_stream.pdf",
                                  title="", size=300, alpha=1)
scvelo.pl.velocity(ctc_scSeq_ID_merged, var_names=["ENSG00000086205"])
scvelo.tl.rank_velocity_genes(ctc_scSeq_ID_merged, groupby="seurat_clusters")




# Dynamic modelling
scvelo.tl.recover_dynamics(ctc_scSeq_ID_merged)
scvelo.tl.velocity(ctc_scSeq_ID_merged, mode="dynamical")
scvelo.tl.velocity_graph(ctc_scSeq_ID_merged)
scvelo.pl.velocity_embedding_stream(ctc_scSeq_ID_merged, basis="umap", color="seurat_clusters", save="ctc_embedding_stream.pdf",
                                  title="", size=300, alpha=1)



















