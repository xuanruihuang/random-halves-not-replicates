"""Check: do the T-cell receptor genes (TRBC1, TRBC2, TRAC, TRGC2; clone markers by construction) drive the transfer
scores? The baselines of pilot_baselines.py (copy, additive, linear) are scored with all pilot genes and with the TCR
genes left out (pilot_common.score gene_mask); for the linear model the TCR genes are also left out of the fit.
Usage: PILOT_DONOR=<D> pilot_tcr_check.py <data_dir> <out.json>"""
import sys, json
import numpy as np
from pilot_common import effects, score, test_pairs, STATES

TCR = ("TRBC", "TRAC", "TRBV", "TRAV", "TRGC", "TRDC")


def ridge_fit(Z, Y, lam):
    return np.linalg.solve(Z.T @ Z + lam * np.eye(Z.shape[1]), Z.T @ Y)


def preds_for(R, S, train, tests, ri, si, keep):
    Xtr = np.stack([R["eff"][ri[l]][keep] for l in train]); Ytr = np.stack([S["eff"][si[l]][keep] for l in train])
    out = {"copy": {}, "additive": {}, "linear": {}}
    shift = (Ytr - Xtr).mean(0); mu_x, mu_y = Xtr.mean(0), Ytr.mean(0)
    U, s, Vt = np.linalg.svd(Xtr - mu_x, full_matrices=False); V = Vt[:50].T; Z = (Xtr - mu_x) @ V
    rng = np.random.default_rng(0); fold = rng.integers(0, 5, len(train)); cv = {}
    for lam in (1e-2, 1e-1, 1, 10, 100, 1000):
        err = 0.0
        for k in range(5):
            tr, te = fold != k, fold == k
            W = ridge_fit(Z[tr], Ytr[tr] - mu_y, lam); err += ((Z[te] @ W + mu_y - Ytr[te]) ** 2).sum()
        cv[lam] = err
    W = ridge_fit(Z, Ytr - mu_y, min(cv, key=cv.get))
    for g in tests:
        x = R["eff"][ri[g]][keep]; full = np.zeros(len(keep))
        for m, v in (("copy", x), ("additive", x + shift), ("linear", ((x - mu_x) @ V) @ W + mu_y)):
            p = full.copy(); p[keep] = v; out[m][g] = p
    return out


def main():
    d, outp = sys.argv[1:3]
    des, res = effects(d)
    R, S = res[STATES[0]], res[STATES[1]]; genes = R["genes"]
    ri, si = {l: i for i, l in enumerate(R["labels"])}, {l: i for i, l in enumerate(S["labels"])}
    train = [l for l in S["labels"] if l.split("|")[0] in set(des["train"]) and l in ri]
    tests = [g for _, _, g, _ in test_pairs(des) if g in ri]
    tcr = np.array([g.startswith(TCR) for g in genes]); keep_all = np.ones(len(genes), bool)
    out = {"tcr_genes": [g for g in genes if g.startswith(TCR)], "with_tcr": {}, "without_tcr": {}}
    for lab, keep in (("with_tcr", keep_all), ("without_tcr", ~tcr)):
        P = preds_for(R, S, train, tests, ri, si, keep)
        for m, p in P.items():
            out[lab][m] = score(p, des, res, genes, gene_mask=keep)
    json.dump(dict(rules=__doc__, **out), open(outp, "w"), indent=1)
    for m in ("copy", "additive", "linear"):
        for k in ("inert", "strong"):
            a, b = out["with_tcr"][m][k], out["without_tcr"][m][k]
            print(m, k, "A-B with %.3f without %.3f | A %.3f/%.3f B %.3f/%.3f" % (a["leak"], b["leak"], a["pearson_A"], b["pearson_A"], a["pearson_B"], b["pearson_B"]), flush=True)


if __name__ == "__main__":
    main()
