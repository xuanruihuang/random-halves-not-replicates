"""Robustness of the CD4 null result to the gene set (A100; output of cd4_sc_extract_cells.py). For unexpressed and
expressed targets, random halves (r1) and sgRNA halves (sgRNA parity), on the same targets as clone_tests_cd4.py
(>= 8 cells in every quarter side): share of targets with a positive cross-half signal and mean cross-half signal,
for four gene sets: all genes with control mean >= 0.1 (the main analyses), the 2,000 of these with the highest control
mean, those without mitochondrial (MT-) and ribosomal-protein (RPL, RPS) genes, and genes with control mean >= 0.5. TCR
genes are excluded throughout; effects and the null-baseline reference as in clone_tests_cd4.py.
Usage: geneset_robustness_cd4.py <cells_dir> <label> [<cells_dir> <label> ...] <out.json>"""
import sys, json
import numpy as np

MINQ = 8
DES = {"random": ([0, 1, 4, 5], [2, 3, 6, 7]), "sgRNA": ([0, 1, 2, 3], [4, 5, 6, 7])}
TCR = ("TRBC", "TRAC", "TRBV", "TRAV", "TRGC", "TRDC")


def run(d, label):
    meta = json.load(open(d + "/meta.json")); genes, groups = meta["genes"], meta["groups"]
    S = np.load(d + "/sum8.npy", mmap_mode="r"); n8 = np.load(d + "/n8.npy"); ci = len(groups) - 1
    cm = np.asarray(S[ci].sum(0), dtype=np.float64) / n8[ci].sum()
    base = np.array([j for j in np.where(cm >= 0.1)[0] if not genes[j].startswith(TCR)])
    sets = {"all kept (control mean >= 0.1)": base,
            "top 2,000 by control mean": np.sort(base[np.argsort(-cm[base])[:2000]]),
            "without MT-, RPL, RPS genes": np.array([j for j in base if not genes[j].startswith(("MT-", "RPL", "RPS"))]),
            "control mean >= 0.5": np.array([j for j in base if cm[j] >= 0.5])}
    gpos = {g: j for j, g in enumerate(genes)}
    texp = {t: cm[gpos[t]] for t in groups[:-1] if t in gpos}
    tids = [i for i, t in enumerate(groups[:-1]) if t in texp and n8[i, [0, 4, 2, 6]].min() >= MINQ]
    tn = [groups[i] for i in tids]
    inert = np.array([texp[t] < 0.01 for t in tn]); expd = np.array([texp[t] >= 0.5 for t in tn])
    nt = n8[tids].astype(float); nc = n8[ci].astype(float)
    out = {"n_unexpressed": int(inert.sum()), "n_expressed": int(expd.sum()), "sets": {}}
    for nm, gk in sets.items():
        Sg = np.asarray(S[tids][:, :, gk], dtype=np.float64); Sc = np.asarray(S[ci][:, gk], dtype=np.float64)
        kpos = {j: k for k, j in enumerate(gk)}
        r = {"n_genes": int(len(gk))}
        for des, (a, b) in DES.items():
            sides = []
            for codes in (a, b):
                m = Sg[:, codes].sum(1) / nt[:, codes].sum(1)[:, None] - Sc[codes].sum(0) / nc[codes].sum()
                for k, t in enumerate(tn):
                    j = kpos.get(gpos.get(t, -1))
                    if j is not None:
                        m[k, j] = 0
                sides.append(m)
            A, B = sides; A = A - A[inert].mean(0); B = B - B[inert].mean(0); x = (A * B).sum(1)
            for s, mm in (("unexpressed", inert), ("expressed", expd)):
                r["%s_%s_pos" % (des, s)] = float((x[mm] > 0).mean()); r["%s_%s_mean" % (des, s)] = float(x[mm].mean())
        out["sets"][nm] = r
        print(label, nm, json.dumps(r), flush=True)
    return out


def main():
    a = sys.argv[1:]; res = {"rules": __doc__}
    for i in range(0, len(a) - 1, 2):
        res[a[i + 1]] = run(a[i], a[i + 1])
    json.dump(res, open(a[-1], "w"), indent=1)


if __name__ == "__main__":
    main()
