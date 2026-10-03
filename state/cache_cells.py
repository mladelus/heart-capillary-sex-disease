"""One-time conversion of the per-donor R files (dl/*.rds) to fast .npz files with the endothelial count matrix."""
import os, sys, glob
import numpy as np
from multiprocessing import Pool
from common import DIRS, OUT, read_rds
CACHE = os.path.join(OUT, "cells"); os.makedirs(CACHE, exist_ok=True)
def conv(a):
    co, f = a; o = os.path.join(CACHE, co + "__" + os.path.basename(f)[:-4] + ".npz")
    if os.path.exists(o): return 0
    x = read_rds(f); m = x["ec_counts"]
    np.savez_compressed(o, dataset_id=str(x["dataset_id"][0]), donor_id=str(x["donor_id"][0]), genes=np.asarray(m["Dimnames"][0], dtype=str),
                        cells=np.asarray(m["Dimnames"][1], dtype=str), x=m["x"].astype(np.float32), i=m["i"].astype(np.int32), p=m["p"].astype(np.int64), dim=np.asarray(m["Dim"], dtype=np.int64))
    return 1
if __name__ == "__main__":
    jobs = [(co, f) for co, d in DIRS.items() for f in sorted(glob.glob(os.path.join(d, "dl", "*.rds")))]
    with Pool(int(os.environ.get("NCPU", "4"))) as p: n = sum(p.imap_unordered(conv, jobs, chunksize=2))
    print("converted", n, "of", len(jobs))
