#!/usr/bin/env python3
import numpy as np
import os
import sys
import h5py

if len(sys.argv) < 2:
    print(f"Usage: {sys.argv[0]} <input_file>")
    sys.exit(1)

file_path = sys.argv[1]

data = np.genfromtxt(file_path, dtype=float, missing_values="nan", filling_values=np.nan, comments=None)

num_nan = np.isnan(data).sum()
mean_val = np.nanmean(data)
std_val = np.nanstd(data)

output_file = file_path + ".h5"
log_file = file_path + ".h5log"

with h5py.File(output_file, "w") as f:
    f.create_dataset("data", data=data, compression="gzip", compression_opts=9)

with open(log_file, "w") as log:
    log.write(f"file {output_file} shape {data.shape} , nan count  {num_nan} {(num_nan/data.shape[0])*100:.2f}%, mean {mean_val:.4f} sd {std_val:.4f} \n")

