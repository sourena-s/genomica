import numpy as np
import h5py
import os
from scipy.stats import norm
import pandas as pd
import gc
import nibabel as nib

list_path = "curated_phenotype_list_491.proc"
list_path =  "all_h5_rho_point5.txt"

h5_base_path = "/projects/0/einf2700/sourena/genica-2025/gwas-databases/panukbb/sumstats/EUR/"
#var_catalogue_file = "/projects/0/einf2700/sourena/genica-2025/gwas-databases/panukbb/full_variant_qc_metrics.sorted.sst.rsid.txt"

var_catalogue_file = '/home/ssoheili/genomica-code/clumps/manual/clumps.sorted.mask.in.var.catalogue'
var_catalogue_file = '/home/ssoheili/genomica-code/clumps/psy_gen_corr_r2_point5/clumps.sorted.mask.in.var.catalogue'

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
        print(f'[{i}/{len(files)}] Loaded, array {v.shape}')
        if v.shape[0] != mask.shape[0]:
            raise ValueError(f'Dimension mismatch! maskdim: { mask.shape[0]} and last dataset dim: {v.shape[0]}. Aborting.')
        
        v = v[mask][np.newaxis,:]
        vectors.append(v)


X = np.concatenate(vectors, axis=0).astype(np.float32)


X4d = X.T[:, None, None, :]

print(f'Final Matrix dim: {X4d.shape}')

affine = np.eye(4, dtype=np.float32)
img = nib.Nifti2Image(X4d, affine)
img.set_qform(affine, code=1)
img.set_sform(affine, code=1)
#nib.save(img, "matrix_4d_manual.nii.gz")
nib.save(img, "matrix_4d_r2psy.nii.gz")


