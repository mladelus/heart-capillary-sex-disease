"""Read the saved R result files (.rds) of both projects into one Python cache (D.pkl) for the figure scripts.
Folders: healthy-heart project (HCAP_DATA, default ~/Documents/heart_capillary_project) and
disease project (HCAP2_DATA, default ~/Documents/heart_capillary_disease). Run after scripts 00-20."""
import os, glob, pickle, sys
import pandas as pd
from rds import read_rds
HERE = os.path.dirname(os.path.abspath(__file__))
P1 = os.path.expanduser(os.environ.get("HCAP_DATA", "~/Documents/heart_capillary_project"))
P2 = os.path.expanduser(os.environ.get("HCAP2_DATA", "~/Documents/heart_capillary_disease"))
if len(sys.argv) == 3: P1, P2 = sys.argv[1:3]
SKIP = {"pb", "ec_subtypes", "inv_cells", "gene_map", "stage3_umap"}          # large files not needed for figures
D = {}
for folder, prefix in ((P1, "p1_"), (P2, "")):
    for f in sorted(glob.glob(os.path.join(folder, "*.rds"))):
        k = os.path.basename(f)[:-4]
        if k in SKIP: continue
        D[(prefix + k) if k.startswith("inv") else k] = read_rds(f)
u = os.path.join(P2, "tables", "S3L_ec_umap_cells.csv.gz")
if os.path.exists(u): D["umap"] = pd.read_csv(u)
pickle.dump(D, open(os.path.join(HERE, "D.pkl"), "wb"))
print("cached:", ", ".join(sorted(D)))
