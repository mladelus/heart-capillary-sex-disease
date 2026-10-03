"""Post hoc check (plan section M, addendum): does the flow-response program (KLF2, KLF4, NOS3, THBD) itself differ by sex
in capillary cells? Uses the donor-level score saved by R script 18 and the per-cell scores of m1_signature.py."""
from common import *
S = samples(); E = read_rds(os.path.join(P2, "stage3_endothelin.rds"))["donors"].copy()
for c in ("stratum", "cohort", "sex", "assay"): E[c] = E[c].astype(str)
r = pooled(E, "Shear_KLF2"); b = r["both"]
print(f"Flow-response score, capillary pseudobulk, female - male: {b['est']:+.3f} ({b['lo']:+.3f} to {b['hi']:+.3f}), p = {b['p']:.3f}, same direction {b['same_dir']}/{b['k']}")
C = pd.read_pickle(os.path.join(OUT, "cells.pkl")); cap = C[C.subtype == "capillary"]; g = cap.groupby("sample")
d = pd.DataFrame(dict(n=g.size(), mean_flow=g.flow_score.mean())); d = d[d.n >= MIN_CAP].join(S).reset_index()
b2 = pooled(d, "mean_flow")["both"]; print(f"Mean per-cell flow-response score: {b2['est']:+.3f} ({b2['lo']:+.3f} to {b2['hi']:+.3f}), p = {b2['p']:.3f}")
pickle.dump(dict(b, per_cell=b2), open(os.path.join(OUT, "flow_sex.pkl"), "wb"))
