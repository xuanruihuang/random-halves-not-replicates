"""Supplementary Table S4c: cell-eval 0.8.2 ceilings of the seven metrics of Fig. 5b in each CD4 donor (resting cells),
random halves (mean of seeds 0-2) / sgRNA halves (mean of the two directions); n.d., not defined (split-half
correlation at most 0). Donors ordered by clone concentration (TRBC readout, Supplementary Table S2b).
Inputs: v1 celleval.json (donor 1), results/celleval/celleval_d{2,3,4}.json (run_celleval_ceiling.py).
Usage: table_s4c_celleval.py <v1 data dir> <results dir> <out.md>"""
import os, sys, json
import numpy as np

METS = [("overlap_at_N", "DE overlap"), ("precision_at_N", "DE precision"), ("de_sig_genes_recall", "DE recall"),
        ("de_direction_match", "DE direction match*"), ("pearson_delta", "Pearson Δ"),
        ("discrimination_score_l1", "Discrimination (L1)*"), ("discrimination_score_cosine", "Discrimination (cosine)*")]
SUBS = [("all", "All perturbations"), ("inert", "Unexpressed targets"), ("strong", "Strength-matched knockdowns")]


def main():
    v1, res, out = sys.argv[1:4]
    src = [("4", os.path.join(res, "celleval", "celleval_d4.json")), ("1", os.path.join(v1, "celleval.json")),
           ("2", os.path.join(res, "celleval", "celleval_d2.json")), ("3", os.path.join(res, "celleval", "celleval_d3.json"))]
    data = {d: json.load(open(f)) for d, f in src}
    rows = ["| Metric | " + " | ".join("Donor %s" % d for d, _ in src) + " |", "|---" * (len(src) + 1) + "|"]
    n_cmp = n_higher = 0
    for sub, lab in SUBS:
        rows.append("| **%s** |" % lab + " |" * len(src))
        for m, name in METS:
            cells = []
            for d, _ in src:
                c = data[d]["shared"]
                r = np.nanmean([c["random_%d" % s][sub]["ceiling"].get(m) or np.nan for s in (0, 1, 2)])
                gv = [x for x in (c[k][sub]["ceiling"].get(m) for k in ("guide_0real_1pred_0", "guide_1real_0pred_0")) if x is not None and not np.isnan(x)]
                g = float(np.mean(gv)) if gv else None
                n_cmp += 1; n_higher += int(g is None or r > g)
                cells.append("%.2f / %s" % (r, "n.d." if g is None else "%.2f" % g))
            rows.append("| %s | %s |" % (name, " | ".join(cells)))
    rows.append("| Perturbations (n) | " + " | ".join(str(data[d]["n_shared_perturbations"]) for d, _ in src) + " |")
    open(out, "w", encoding="utf-8").write("\n".join(rows) + "\n")
    print("random-half ceiling higher in %d of %d comparisons" % (n_higher, n_cmp))


if __name__ == "__main__":
    main()
