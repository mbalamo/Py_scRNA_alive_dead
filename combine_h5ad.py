# combine two (or more) h5ad files and downsample if needed

# %% 
# package import
import anndata as ad


# %% 
# open two files and combine
adata1 = ad.read_h5ad("cardiomyocytes.h5ad")
adata2 = ad.read_h5ad("neural_crest.h5ad")
# subselect rows
sub1 = adata1[:5000, :].copy()
sub2 = adata2[:5000, :].copy()

print(sub1.obs.columns.tolist())
print(sub2.obs.columns.tolist())

# combine along the row axis. 
# join="outer": keps all genes from both files. If a gene is missing in one file, it fills it with NaN values.
# join="inner": keeps only overlap.
# appends "-0" to file1 names and "-1" to file2 names
combined = ad.concat([sub1, sub2], axis=0, join="inner", label="batch", keys=["file1", "file2"])

# write
# clear empty tracking references to prevent save errors
combined.uns.clear()
combined.write_h5ad("combined_output.h5ad")
