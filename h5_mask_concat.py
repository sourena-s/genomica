import numpy as np
import h5py
import os
from scipy.stats import norm
import pandas as pd
import gc
import nibabel as nib

list_path =  "all_h5_rho_point5.txt"
#list_path = "curated_phenotype_list_491.proc"

h5_base_path = "/projects/0/einf2700/sourena/genica-2025/gwas-databases/panukbb/sumstats/EUR/"
#var_catalogue_file = "/projects/0/einf2700/sourena/genica-2025/gwas-databases/panukbb/full_variant_qc_metrics.sorted.sst.rsid.txt"

#var_catalogue_file = '/home/ssoheili/genomica-code/clumps/psy_gen_corr_r2_point5/clumps.sorted.mask.in.var.catalogue'
#var_catalogue_file = '/home/ssoheili/genomica-code/clumps/manual-r2point3/clumps.sorted.mask.in.var.catalogue'
#var_catalogue_file = '/home/ssoheili/genomica-code/clumps/manual-r2point3-maf-1percent/clumps.sorted.mask.in.var.catalogue'
#var_catalogue_file = '/home/ssoheili/genomica-code/clumps/merged_manual-r2point3_psy_gen_corr_r2_point5/clumps.sorted.mask.in.var.catalogue'

var_catalogue_file = '/home/ssoheili/genomica-code/clumps/hapmap3/clumps.sorted.mask.in.var.catalogue'

df = pd.read_csv(var_catalogue_file, sep=' ', header=None)

catalogue = df.iloc[:, 0].values              # rsid (col1)
assert set(df.iloc[:,1].unique()) <= {0,1}
mask = df.iloc[:, 1].values.astype(bool)  #clump mask (col 2)

with open(list_path) as f:
    files = [os.path.join(h5_base_path, line.strip()) for line in f if line.strip()]

vectors = []
batch_id=0
max_abs = -1

for i, path in enumerate(files, start=1):
    with h5py.File(path, "r") as h5f:
        v = np.array(h5f["data"])
        if v.shape[0] != mask.shape[0]:
            raise ValueError(f"Dimension mismatch! maskdim: {mask.shape[0]} and last dataset dim: {v.shape[0]}")

        v = v[mask]          # shape (M,)
        vectors.append(v)    # list of (M,)
        print(f'v shape { v.shape}')
        n_nan = np.isnan(v).sum()
        mean = np.nanmean(v)
        std  = np.nanstd(v)
        print(
            f"[{i}/{len(files)}] "
            f"shape={v.shape} "
            f"mean={mean:.4g} std={std:.4g} n_nan={n_nan}"
        )

X = np.stack(vectors, axis=1).astype(np.float32)   # shape (N, M)

print(f'stacked X dim: {X.shape}')

#assert X.shape[0] == 63 * 47 * 191

X4d = X[:,None,None,:]

n_nan = np.isnan(X4d).sum()
n_inf = np.isinf(X4d).sum()
n_fin = np.isfinite(X4d).sum()


print(f"X4d: nan={n_nan}, inf={n_inf}, finite={n_fin}, total={X4d.size}")
X4d = np.nan_to_num(X4d, nan=0.0)

print(f"X4d after nan removal: nan={np.isnan(X4d).sum()}")

print(f'Final Matrix dim: {X4d.shape}')

affine = np.eye(4, dtype=np.float32)
img = nib.Nifti2Image(X4d, affine)
#img = nib.Nifti1Image(X4d, affine)
img.set_qform(affine, code=1)
img.set_sform(affine, code=1)

#nib.save(img, "matrix_4d_r2psy.nii.gz")
#nib.save(img, "matrix_4d_manual_rpoint3_maf_1percent.nii.gz")
#nib.save(img, "matrix_4d_manual_rpoint3.nii.gz")

#nib.save(img, "matrix_4d_r2psy_mergedmask.nii.gz")
nib.save(img, "/projects/0/einf2700/sourena/genica-2025/melodic_results/matrix_4d_hapmap_r2psy_traits.nii.gz")
