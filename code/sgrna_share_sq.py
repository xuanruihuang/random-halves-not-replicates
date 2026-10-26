"""sgRNA cell shares per target for the three-part decomposition: mean over targets of w0^2 + w1^2, where w_g is the
share of the target's cells carrying sgRNA g (n8 codes gp*4 + r1*2 + r2), over the targets of clone_tests_cd4.py
(>= 8 cells in every quarter side); expressed targets (control mean >= 0.5) where sum8 is available, else all targets."""
import json, os, sys
import numpy as np
out = {}
for d in ("D1", "D2", "D3", "D4"):
    for st, lab in (("Rest", "rest"), ("Stim48hr", "stim48")):
        p = os.environ["STORE_ROOT"] + "/cd4_sc/cells_%s_%s" % (d, st)
        if d == "D1" and os.path.exists(p + "_v2"): p = p + "_v2"
        meta = json.load(open(p + "/meta.json")); groups = meta["groups"]; n8 = np.load(p + "/n8.npy").astype(float)
        tids = [i for i in range(len(groups) - 1) if n8[i, [0, 4, 2, 6]].min() >= 8]
        w0 = n8[tids, :4].sum(1) / n8[tids].sum(1); f_all = float(np.mean(w0 ** 2 + (1 - w0) ** 2))
        f_exp = None
        if os.path.exists(p + "/sum8.npy"):
            S = np.load(p + "/sum8.npy", mmap_mode="r"); genes = meta["genes"]; ci = len(groups) - 1
            cm = np.asarray(S[ci].sum(0)) / n8[ci].sum(); gpos = {g: j for j, g in enumerate(genes)}
            e = np.array([groups[i] in gpos and cm[gpos[groups[i]]] >= 0.5 for i in tids])
            f_exp = float(np.mean((w0 ** 2 + (1 - w0) ** 2)[e])); n_exp = int(e.sum())
        out["%s_%s" % (d, lab)] = dict(f_all=f_all, f_expressed=f_exp, n_targets=len(tids), dir=p)
        print(d, lab, round(f_all, 4), f_exp and round(f_exp, 4), len(tids), flush=True)
json.dump(dict(rules=__doc__, shares=out), open(os.environ["STORE_ROOT"] + "/review_tests/sgrna_share_sq.json", "w"), indent=1)
