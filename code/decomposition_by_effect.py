"""Three-part decomposition of the random-half signal of expressed CD4 targets, stratified by effect size (A100).
Per donor and state, the per-target random-half (r1 = 0 vs 1) and sgRNA-half (gp = 0 vs 1) cross-half signals are
computed exactly as in clone_tests_cd4.py (same genes, targets, control referencing, target gene zeroed, referenced to the
mean of the unexpressed targets). Targets are stratified by an effect size measured without shared clones: the mean
sgRNA-half cross-half signal of the same target in the other three donors (same state), in tertiles over the expressed
targets present in all four donors. In each tertile: reproduced = sum G / sum R; sgRNA-specific = f * ratio * reproduced
(f from sgrna_share_sq.json, ratio from guide_donor.json); within one transduction = the rest.
Usage: decomposition_by_effect.py <out.json>"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from clone_tests_cd4 import DESIGNS

GD = {"rest": 0.80, "stim48": 0.48}   # guide_specific_over_diff_guide (guide_donor.json, released pseudobulk of all four donors)


def per_target(d):
    meta = json.load(open(d + "/meta.json")); genes, groups = meta["genes"], meta["groups"]
    S = np.load(d + "/sum8.npy", mmap_mode="r"); n8 = np.load(d + "/n8.npy")
    ci = len(groups) - 1
    cm = np.asarray(S[ci].sum(0)) / n8[ci].sum()
    trg = {j for j, g in enumerate(genes) if g.startswith(("TRBC", "TRAC", "TRBV", "TRAV", "TRGC", "TRDC"))}
    gk = np.array([j for j in np.where(cm >= 0.1)[0] if j not in trg])
    gpos = {g: j for j, g in enumerate(genes)}; gkpos = {j: k for k, j in enumerate(gk)}
    texp = {t: cm[gpos[t]] for t in groups[:-1] if t in gpos}
    tids = [i for i, t in enumerate(groups[:-1]) if t in texp and n8[i, [0, 4, 2, 6]].min() >= 8]
    tn = [groups[i] for i in tids]
    inert = np.array([texp[t] < 0.01 for t in tn]); expd = np.array([texp[t] >= 0.5 for t in tn])
    Sg = np.asarray(S[tids][:, :, gk], dtype=np.float64); Sc = np.asarray(S[ci][:, gk], dtype=np.float64)
    nt = n8[tids].astype(float); nc = n8[ci].astype(float)
    out = {}
    for nm in ("random_half", "guide_half"):
        a, b = DESIGNS[nm]; sides = []
        for codes in (a, b):
            m = Sg[:, codes].sum(1) / nt[:, codes].sum(1)[:, None] - Sc[codes].sum(0) / nc[codes].sum()
            for k, i in enumerate(tids):
                j = gkpos.get(gpos.get(groups[i], -1))
                if j is not None:
                    m[k, j] = 0
            sides.append(m)
        A, B = sides
        out[nm] = ((A - A[inert].mean(0)) * (B - B[inert].mean(0))).sum(1)
    return {t: (out["random_half"][k], out["guide_half"][k]) for k, t in enumerate(tn) if expd[k]}


def main():
    outp = sys.argv[1]
    f = json.load(open(os.environ["STORE_ROOT"] + "/review_tests/sgrna_share_sq.json"))["shares"]
    res = {}
    for st, lab in (("Rest", "rest"), ("Stim48hr", "stim48")):
        pt = {}
        for d in ("D1", "D2", "D3", "D4"):
            p = os.environ["STORE_ROOT"] + "/cd4_sc/cells_%s_%s" % (d, st)
            if d == "D1" and os.path.exists(p + "_v2"):
                p += "_v2"
            pt[d] = per_target(p); print(d, lab, len(pt[d]), flush=True)
        common = sorted(set.intersection(*[set(v) for v in pt.values()]))
        for d in pt:
            others = [o for o in pt if o != d]
            eff = np.array([np.mean([pt[o][t][1] for o in others]) for t in common])
            R = np.array([pt[d][t][0] for t in common]); G = np.array([pt[d][t][1] for t in common])
            cuts = np.quantile(eff, [1 / 3, 2 / 3]); tert = np.digitize(eff, cuts)
            rows = []
            for k in range(3):
                w = tert == k
                rep = G[w].sum() / R[w].sum(); gsp = f["%s_%s" % (d, lab)]["f_expressed"] * GD[lab] * rep
                rows.append(dict(tertile=k + 1, n=int(w.sum()), effect_other_donors_mean=float(eff[w].mean()),
                                 reproduced=float(rep), specific=float(gsp), within_transduction=float(1 - rep - gsp)))
            repa = G.sum() / R.sum()
            res["%s_%s" % (d, lab)] = dict(n_common=len(common), all=dict(reproduced=float(repa)), tertiles=rows)
            print(d, lab, json.dumps(rows), flush=True)
    json.dump(dict(rules=__doc__, results=res), open(outp, "w"), indent=1)


if __name__ == "__main__":
    main()
