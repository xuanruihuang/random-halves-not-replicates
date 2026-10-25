"""Null subtraction as an approximate correction without independent replicates (CD4 donors 1-4, resting and
stimulated): mean random-half cross-half signal of expressed targets minus that of unexpressed targets, compared with
the sgRNA-half signal plus the sgRNA-specific part (sgRNA-half signal times the across-donor sgRNA-specific ratio of
guide_donor.json), as in Fig. 5a. Inputs: the v1 result files clone_D1.json, clone_D2.json, guide_donor.json;
for donors 3 and 4, the same quantities from clone_tests_cd4.py (model_eval/results/clone_tests_D3.json, _D4.json).
Usage: null_subtraction_cd4.py <v1 data dir> <out.json> [<donor>=<clone_tests json> ...]"""
import sys, os, json


def main():
    d, out = sys.argv[1:3]
    files = [(os.path.join(d, "clone_D1.json"), "D1"), (os.path.join(d, "clone_D2.json"), "D2")]
    files += [(a.split("=", 1)[1], a.split("=", 1)[0]) for a in sys.argv[3:]]
    gd = json.load(open(os.path.join(d, "guide_donor.json")))
    res = {"rules": __doc__, "rows": []}
    for f, donor in files:
        c = json.load(open(f))
        for st, rk in (("rest", "rest"), ("stim48", "stim")):
            x = c[st]["depth"]
            R = x["random_half"]["expressed"]["cross_per_target"]; G = x["guide_half"]["expressed"]["cross_per_target"]
            N = x["random_half"]["inert"]["cross_per_target"]
            sg = G * (1 + gd["guide_specific_over_diff_guide_" + rk])
            res["rows"].append(dict(donor=donor, state=st, random=R, sgRNA_half=G, null_random=N, random_minus_null=R - N,
                                    sgRNA_related=sg, random_over_sgRNA_half=R / G, corrected_over_sgRNA_related=(R - N) / sg))
    for r in res["rows"]:
        print(r["donor"], r["state"], "random/sgRNA %.1f, corrected/sgRNA-related %.2f" % (r["random_over_sgRNA_half"], r["corrected_over_sgRNA_related"]))
    json.dump(res, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
