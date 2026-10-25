"""Transfer-model scores in the four CD4 donors, collected from the per-donor score files (results/models4/<D>_<name>.json,
copied from the 5090: pilot_baselines.py, pilot_state_score.py, pilot_pie_score.py, pilot_scgen.py) for Fig. 7a and
Supplementary Table S6a. Score = mean Pearson r over test sgRNAs between predicted and answer effects; A, answer from the
same sgRNA (shares clones with the input); B, from the other sgRNA of the target (another transduction).
Usage as a module: collect(res) -> {model: {donor: {kind: (A, B, n)}}}; as a script: table_s6a(res, out_md)."""
import os, sys, json

DONORS = ("D1", "D2", "D3", "D4")
MODELS = [("mean", "Context mean (baseline)"), ("copy", "Copy resting effect"), ("additive", "Additive"), ("linear", "Linear (ridge)"),
          ("arc_perturb_mean", "Arc perturb_mean"), ("arc_pseudobulk", "Arc pseudobulk"), ("scgen", "scGen"),
          ("state_st", "STATE"), ("state_se", "State (SE-600M + STATE)"), ("pie", "PIE (expression shift)")]
FILES = {"arc_perturb_mean": "arc_perturb_mean_scores", "arc_pseudobulk": "arc_pseudobulk_scores", "scgen": "scgen_scores",
         "state_st": "state_st_scores", "state_se": "state_se_scores", "pie": "pie_scores"}


def collect(res):
    d4 = os.path.join(res, "models4"); out = {m: {} for m, _ in MODELS}
    for d in DONORS:
        f = os.path.join(d4, "%s_baselines.json" % d)
        if os.path.exists(f):
            b = json.load(open(f))["models"]
            for m in ("mean", "copy", "additive", "linear"):
                out[m][d] = {k: (v["pearson_A"], v["pearson_B"], v["n"]) for k, v in b[m].items() if k in ("inert", "strong")}
        for m, name in FILES.items():
            f = os.path.join(d4, "%s_%s.json" % (d, name))
            if not os.path.exists(f):
                continue
            s = json.load(open(f))["scores"]
            if m == "pie":
                s = s["delta_p_pred"]
            out[m][d] = {k: (v["pearson_A"], v["pearson_B"], v["n"]) for k, v in s.items() if k in ("inert", "strong")}
    return out


def table_s6a(res, outp):
    c = collect(res)
    rows = ["| Model | " + " | ".join("Donor %s" % d[1] for d in DONORS) + " |", "|---" * (len(DONORS) + 1) + "|"]
    for kind, lab in (("inert", "Unexpressed targets: A / B"), ("strong", "Strength-matched knockdowns: A / B")):
        rows.append("| **%s** |" % lab + " |" * len(DONORS))
        for m, name in MODELS:
            cells = []
            for d in DONORS:
                v = c[m].get(d, {}).get(kind)
                cells.append("–" if v is None else ("%.3f / %.3f" % (v[0], v[1])).replace("-", "−"))
            rows.append("| %s | %s |" % (name, " | ".join(cells)))
    n = [c["copy"][d]["inert"][2] if d in c["copy"] else "–" for d in DONORS]
    k = [c["copy"][d]["strong"][2] if d in c["copy"] else "–" for d in DONORS]
    rows.append("| Test sgRNAs, unexpressed / knockdowns (n) | " + " | ".join("%s / %s" % (a, b) for a, b in zip(n, k)) + " |")
    open(outp, "w", encoding="utf-8").write("\n".join(rows) + "\n")


if __name__ == "__main__":
    table_s6a(sys.argv[1], sys.argv[2])
