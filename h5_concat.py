import numpy as np
import h5py
import os
from scipy.stats import norm
import pandas as pd
import gc

list_path = "all_h5.txt"
h5_base_path = "/home/ssoheili/genetic-data/genica/gwas-databases/panukbb/sumstats/EUR/" 
var_catalogue_file = '/home/ssoheili/genetic-data/genica/gwas-databases/panukbb/full_variant_qc_metrics.sorted.sst.rsid.txt'

with open(list_path) as f:
    files = [os.path.join(h5_base_path, line.strip()) for line in f if line.strip()]

vectors = []
global_dim = -1
running_batch_size=25
batch_id=0
max_abs = -1

for i, path in enumerate(files, start=1):
    with h5py.File(path, "r") as h5f:
        v = np.array(h5f["data"])
        print(f'[{i}/{len(files)}] Loaded, array {v.shape}')
        if global_dim != v.shape[0] and global_dim != -1:
            raise ValueError(f'Dimension mismatch! global dim: {global_dim} and last dataset dim: {v.shape[0]}. Aborting.')
        else:
            global_dim = v.shape[0]
        
        v = v[np.newaxis,:]
        batch_id= batch_id+1
        vectors.append(v)
        if batch_id >= running_batch_size:
            print(f'Calculating max z across {batch_id} vectors..')

            if np.isscalar(max_abs):
                max_abs = np.full((global_dim,), np.nan)

            result = np.concatenate(vectors, axis=0)
            max_abs = np.nanmax(np.concatenate([np.abs(result), max_abs[np.newaxis, :]], axis=0), axis=0)
            #release the memory
            vectors = []
            del result
            gc.collect()
            batch_id=0

#leftover?
if batch_id>0:
    result = np.concatenate(vectors, axis=0)
    max_abs = np.nanmax(np.concatenate([np.abs(result), max_abs[np.newaxis, :]], axis=0), axis=0)

print(f"Result shape: {max_abs.shape}")  # (n_phenotypes, n_variants)

neg_log10_p= -(np.log(2) + norm.logsf(abs(max_abs, 0, 1)) / np.log(10)

var_catalogue = pd.read_csv(var_catalogue_file, sep=' ', header=None, dtype=str)

if var_catalogue.shape[0] != p_values.shape[0]:
    raise ValueError(f"Row count mismatch: var_catalogue={var_catalogue.shape[0]}, p_values={p_values.shape[0]}")

var_catalogue['z'] = max_abs
var_catalogue['neg_log10_p'] = neg_log10_p

var_catalogue.to_csv('variants_min_pval.txt', sep=' ', index=False)


