"""Supplementary Table S2b: the clone tests of the main text in each CD4 donor (rows = quantities, columns = donor and
state). Inputs (all computed elsewhere; this script only collects them):
  clone_tests_D<d>.json      clone_tests_cd4.py: TRBC var(z), TRBC deviation r between random and sgRNA halves
  trbc_state_keff.json (v1, donors 1-2), trbc_state_keff_D34.json: upper bound on K_eff (median n / z^2, resting
                             and stimulated), r of the TRBC deviation between states
  ci_cd4.json, ci_cd4_D2.json (v1), ci_D3.json, ci_D4.json: ci_bootstrap.py cd4 - share of unexpressed targets with a
                             positive cross-half signal (random halves r1 / sgRNA halves), reliability labels (resting)
  trbc_split/<D>_<state>.json  trbc_split_cd4.py: constant-region halves
  clone_D1.json, clone_D2.json (v1), clone_tests_D3/D4.json, guide_donor.json (v1): the three parts of Fig. 5a
Usage: table_per_donor_cd4.py <v1 data dir> <results dir> <out.md>"""
import os, sys, json

DON = ("D1", "D2", "D3", "D4"); ST = (("rest", "Rest", "rest"), ("stim48", "Stim48", "stim"))


def J(*p):
    return json.load(open(os.path.join(*p)))


def main():
    v1, res, out = sys.argv[1:4]
    ct = {d: J(res, "clone_tests_%s.json" % d) for d in DON}
    ke = J(v1, "trbc_state_keff.json"); ke.update(J(res, "trbc_state_keff_D34.json"))
    ci = {"D1": J(v1, "ci_cd4.json"), "D2": J(v1, "ci_cd4_D2.json"), "D3": J(res, "ci_D3.json"), "D4": J(res, "ci_D4.json")}
    gd = J(v1, "guide_donor.json")
    cl = {"D1": J(v1, "clone_D1.json"), "D2": J(v1, "clone_D2.json"), "D3": ct["D3"], "D4": ct["D4"]}
    cols = [(d, s) for d in DON for s in ST]
    rows = []

    def row(name, f, fmt):
        vals = []
        for d, s in cols:
            v = f(d, s)
            vals.append("–" if v is None else fmt(v))
        rows.append("| %s | %s |" % (name, " | ".join(vals)))

    pc = lambda v: "%d" % round(100 * v)
    rows.append("| **TRBC clone marker** |" + " |" * len(cols))
    row("var(z) of the TRBC1 share (1 if no clones)", lambda d, s: ct[d][s[0]]["trbc"]["overdispersion"]["var_z"], lambda v: "%.0f" % v)
    row("r of the TRBC1 deviation, random halves", lambda d, s: ct[d][s[0]]["trbc"]["reproducibility_random_half"]["pearson"], lambda v: "%.2f" % v)
    row("r of the TRBC1 deviation, sgRNA halves", lambda d, s: ct[d][s[0]]["trbc"]["reproducibility_guide_half"]["pearson"], lambda v: ("%.2f" % v).replace("-0.00", "0.00").replace("-", "−"))
    row("r of the TRBC1 deviation, resting vs stimulated", lambda d, s: ke[d]["pearson_rest_vs_stim"] if s[0] == "rest" else None, lambda v: "%.2f" % v)
    row("Upper bound on K~eff~ (median n/z^2^)", lambda d, s: ke[d]["keff_%s" % ("rest" if s[0] == "rest" else "stim")]["keff_upper_median"], lambda v: "%.1f" % v)
    rows.append("| **Unexpressed targets, positive cross-half signal (%)** |" + " |" * len(cols))
    row("Random halves", lambda d, s: ci[d]["random_r1"]["inert_frac_pos_%s" % s[2]][0], pc)
    row("sgRNA halves", lambda d, s: ci[d]["guide"]["inert_frac_pos_%s" % s[2]][0], pc)
    row("Targets (n)", lambda d, s: ci[d]["n_inert"], lambda v: "{:,}".format(v))
    rows.append("| **Constant-region halves, unexpressed targets** |" + " |" * len(cols))
    tr = lambda d, s: J(res, "trbc_split", "%s_%s.json" % (d, s[1]))["groups"]
    row("Positive, random halves (%)", lambda d, s: tr(d, s)["inert"]["share_pos_random_r1"], pc)
    row("Positive, constant-region halves (%)", lambda d, s: tr(d, s)["inert"]["share_pos_trbc"], pc)
    row("Cross-half signal kept (%)", lambda d, s: tr(d, s)["inert"]["mean_trbc"] / tr(d, s)["inert"]["mean_random_r1"], pc)
    row("Targets higher under random halves (%)", lambda d, s: tr(d, s)["inert"]["share_random_gt_trbc"], pc)
    row("Targets (n)", lambda d, s: tr(d, s)["inert"]["n_targets"], lambda v: "{:,}".format(v))
    row("Expressed targets: signal kept (%)", lambda d, s: tr(d, s)["expressed"]["mean_trbc"] / tr(d, s)["expressed"]["mean_random_r1"], pc)

    def parts(d, s):
        x = cl[d][s[0]]["depth"]; R = x["random_half"]["expressed"]["cross_per_target"]; G = x["guide_half"]["expressed"]["cross_per_target"]
        rep = G / R; gsp = gd["guide_specific_over_diff_guide_" + s[2]] * rep
        return rep, gsp, 1 - rep - gsp
    rows.append("| **Expressed targets, random-half signal (%)** |" + " |" * len(cols))
    row("Reproduced by the other sgRNA", lambda d, s: parts(d, s)[0], pc)
    row("Specific to the sgRNA", lambda d, s: parts(d, s)[1], pc)
    row("Only within one transduction", lambda d, s: parts(d, s)[2], pc)
    row("Targets (n)", lambda d, s: cl[d][s[0]]["n_expressed"], lambda v: "{:,}".format(v))
    rows.append("| **Labelled reliable, resting cells (%)** |" + " |" * len(cols))
    rl = lambda d, s, k, g: ci[d]["wang_reliable_share_rest"][k][g] if s[0] == "rest" else None
    p1 = lambda v: "%.1f" % (100 * v) if v < 0.1 else "%d" % round(100 * v)
    row("Unexpressed, random halves", lambda d, s: rl(d, s, "random", "inert"), p1)
    row("Unexpressed, sgRNA halves", lambda d, s: rl(d, s, "guide", "inert"), p1)
    row("Expressed, random halves", lambda d, s: rl(d, s, "random", "expressed"), p1)
    row("Expressed, sgRNA halves", lambda d, s: rl(d, s, "guide", "expressed"), p1)
    row("Targets, unexpressed / expressed (n)", lambda d, s: "%s / %s" % ("{:,}".format(ci[d]["n_inert"]), "{:,}".format(ci[d]["n_expressed"])) if s[0] == "rest" else None, str)
    head = "| Quantity | " + " | ".join("%s, %s" % (d.replace("D", "Donor "), s[2]) for d, s in cols) + " |"
    md = "\n".join([head, "|" + "---|" * (len(cols) + 1)] + rows) + "\n"
    open(out, "w", encoding="utf-8").write(md)


if __name__ == "__main__":
    main()
