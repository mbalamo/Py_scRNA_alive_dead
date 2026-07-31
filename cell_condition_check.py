# check the scRNA cell quality condition by various gene expression metrics

# %% 
# package import
import anndata as ad
import scanpy as sc
import scipy.sparse as sp
# for ENSG to gene name conversion
import mygene
# a curated, non-arbitrary list, the Hallmark Apoptosis gene set (MSigDB, HALLMARK_APOPTOSIS) is the 
# standard reference most pipelines use for scoring -- it contains ~161 genes and is more comprehensive 
# and curated than any manual list. access the MSigDB API. download the Hallmark collection (h.all), 
# extract the apoptosis gene list.
import decoupler as dc


# %%
# choose the file to open and downsample if desired.
# downloaded or synthetic ones created based on known gene signatures

file_name = "GSM4816047"
full_name = file_name + ".h5ad"
adata = ad.read_h5ad(full_name)

# down sample -- replace current frame
#adata = adata[:2000, :].copy()
# save the updated h5ad file
#full_name = file_name + "_2k.h5ad"
#adata.write_h5ad(full_name, compression="gzip")

# built in example data sets - unfiltered 10x PBMC dataset, which includes low-quality and dying cells
# read in from scanpy database.
#adata = sc.datasets.pbmc3k()

# %% 
# update symbols to gene names -- if needed
#
# Identify which var_names are Ensembl gene IDs (maybe a mix of ENSG and gene symbols)
is_ensg = adata.var_names.str.match(r'^ENSG\d+')
# Extract ENSG IDs from adata.var (Remove version numbers like .X if they are present)
#ensg_ids = adata.var_names.str.split('.').str[0]
ensg_ids = adata.var_names[is_ensg].str.split('.').str[0]

# only if it find something to map, do it. else skip the process or query will crash.
if len(ensg_ids) != 0:
    # check the mygene library to map ENSG IDs to gene symbols
    mg = mygene.MyGeneInfo()
    #  maps -- output the dup counts for names and no hit counts
    results = mg.querymany(ensg_ids.tolist(), scopes='ensembl.gene', fields='symbol', species='human', as_dataframe=True)
    # clean the results and assign them back to adata.var -- remove duplicates
    results = results[~results.index.duplicated(keep='first')]

    #print("Duplicates before:", adata.var_names.duplicated().sum())

    # this version does not handle NaN
    # adata.var['gene_name'] = adata.var_names.map(results['symbol'].to_dict())

    # Map using the STRIPPED ids (this matches how querymany was called)
    symbol_map = results['symbol'].to_dict()
    # Start with the original var_names — this preserves existing gene symbols untouched
    gene_names = adata.var_names.to_series()
    # Only overwrite the ENSG entries with their mygene lookup
    looked_up = ensg_ids.map(symbol_map)
    looked_up.index = adata.var_names[is_ensg]  # align back to original var_names index
    gene_names.loc[is_ensg] = looked_up
    #
    #
    # apply map to the stripped names, not the original var names.
    #gene_names = ensg_ids.map(symbol_map)
    # Fill any NaN symbols with the original ENSG id, and force everything to str
    #gene_names = gene_names.where(gene_names.notna(), adata.var_names.to_series().values)
    #
    #
    gene_names = gene_names.astype(str)
    adata.var['gene_name'] = gene_names.values
    # set gene_name as your var_names for plotting 
    adata.var.set_index('gene_name', inplace=True)
    # Make sure the new index is unique (duplicate symbols can still occur)
    adata.var_names_make_unique()

    # save the updated h5ad file
    #full_name = file_name + "_named.h5ad"
    #adata.write_h5ad(full_name)

# %% 
# check structure to confirm cleaned anndata
#
#print(adata.var.head())
#print(adata.var.keys())
print(adata.obs_names[:10])
print(adata.var_names[:10])

# print out col names for .obs and .var in adddata
#print(adata.obs.columns.tolist())
#print(adata.var.columns.tolist())

# Check if the matrix is sparse
print(isinstance(adata.X, sp.spmatrix))
#print(adata.obs["sample_id"].value_counts())

# How many genes in the data set
print(f"remaining cells: {adata.shape[0]}, remaining TOTAL genes: {adata.shape[1]}")


# %% 
# defining the gene groups for quality condition testing
#
# make all caps to match gene_names
marker_genes = {
    "apo_core": ["CASP3", "CASP7", "CASP9", "BAX", "BAK1", "BBC3", "FAS", "PMAIP1", "DDIT3"],   
    "dmg_core": [g.upper() for g in ["Fosb", "Fos", "Jun", "Junb", "Jund", "Atf3", "Egr1", "Hspa1a", "Hspa1b", "Hsp90ab1", "Hspa8", "Hspb1", "Ier3", "Ier2", "Btg1", "Btg2", "Dusp1"]],
}
# download MSigDB resource dataframe and filter for human Hallmark Apoptosis
msigdb = dc.op.resource("MSigDB")
apoptosis_df = msigdb[(msigdb['geneset'] == 'HALLMARK_APOPTOSIS')]
apoptosis_genes = apoptosis_df['genesymbol'].tolist()
# keep only those present in the data set
valid_apo_genes = [g for g in apoptosis_genes if g in adata.var_names]


# %% 
# check all cell expression against gene lists -- with visual output
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["apo_total"] = adata.var_names.isin(valid_apo_genes)
adata.var["apo_core"] = adata.var_names.isin(marker_genes["apo_core"])
adata.var["dmg_core"] = adata.var_names.isin(marker_genes["dmg_core"])

# adata.var[gene_tag] is a boolean column. we can SUM it 
print("mito total genes found: ", adata.var["mt"].sum())
print("apo total genes found: ", adata.var["apo_total"].sum())
print("apo core genes found: ", adata.var["apo_core"].sum())
print("dmg core genes found: ", adata.var["dmg_core"].sum())

# already log1p -- confirm
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "apo_total", "apo_core", "dmg_core"], inplace=True, log1p=False)
sc.pl.violin(adata, ["pct_counts_mt", "pct_counts_apo_total", "pct_counts_dmg_core"], jitter=0.4, multi_panel=True,)

# %%
# set a threshold for "alive" or "dead" or "damaged" 

dmg_threshold = 5
apo_threshold = 1
counts_by_cell_threshold = 300

adata.obs["likely_damaged"] = ((adata.obs["pct_counts_mt"] > dmg_threshold))
adata.obs["likely_apoptotic"] = ((adata.obs["pct_counts_apo_total"] > apo_threshold))

print("cells likely damaged: ", adata.obs["likely_damaged"].sum())
print("cells likely apoptotic: ", adata.obs["likely_apoptotic"].sum())

# %% 
# 2D scatter
#sc.pl.scatter(adata, x="pct_counts_mt", y="pct_counts_apo_total", color="n_genes_by_counts")

import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(6, 5))
sc_plot = ax.scatter(
    adata.obs["pct_counts_mt"],
    adata.obs["pct_counts_apo_total"],
    c=adata.obs["n_genes_by_counts"],
    cmap="viridis",
    s=3,
    alpha=0.7,
    zorder=1
)

# overlay low-count cells with a distinct marker color
low_count_mask = adata.obs["n_genes_by_counts"] < counts_by_cell_threshold
ax.scatter(
    adata.obs["pct_counts_mt"][low_count_mask],
    adata.obs["pct_counts_apo_total"][low_count_mask],
    color="red",
    s=3,
    zorder=2,
    label=f"< {counts_by_cell_threshold} genes (n={low_count_mask.sum()})",
)

ax.set_xlabel("percent counts mitochondrial")
ax.set_ylabel("percent counts apoptotic")
ax.set_title("Mito vs Apoptotic cells, colored by gene count")

#ax.legend(loc="upper right")

fig.legend(
    bbox_to_anchor=(1, 1),
    bbox_transform=fig.transFigure  
)

cbar = fig.colorbar(sc_plot, ax=ax)
cbar.set_label("")
# damage threshold
ax.axvline(dmg_threshold, color="red", linestyle="--", linewidth=1)
# apoptotic threshold
ax.axhline(apo_threshold, color="orange", linestyle="--", linewidth=1)
plt.show()


# %%
# RUN