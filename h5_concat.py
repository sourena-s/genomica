import numpy as np
import h5py
import os
from scipy.stats import norm


list_path = "test_all_h5.txt"
h5_base_path = "/home/ssoheili/genetic-data/genica/gwas-databases/panukbb/sumstats/EUR/" 

with open(list_path) as f:
    files = [os.path.join(h5_base_path, line.strip()) for line in f if line.strip()]

vectors = []
global_dim = -1
for i, path in enumerate(files, start=1):
    with h5py.File(path, "r") as h5f:
        v = np.array(h5f["data"])
        print(f'[{i}/{len(files)}] Loaded, array {v.shape}')
        if global_dim != v.shape[0] and global_dim != -1:
            print(f'Dimension mismatch! global dim: {global_dim} and last dataset dim: {v.shape[0]}. Aborting.')
            sys.exit(1)
        else:
            global_dim = v.shape[0]
        
        v = v[np.newaxis,:]
        vectors.append(v)

result = np.concatenate(vectors, axis=0)

print(f"Result shape: {result.shape}")  # (n_phenotypes, n_variants)

max_abs = np.max(np.abs(result), axis=0, keepdims=False)

p_values = 2 * (1 - norm.cdf(np.abs(max_abs)))

np.savetxt("test_min_pvalues.txt", p_values, fmt="%.6e")

