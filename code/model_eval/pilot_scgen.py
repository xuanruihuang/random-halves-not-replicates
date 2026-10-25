"""Model evaluation: scGen (theislab/scgen 2.1.1 from its GitHub repository, scvi-tools 1.3.3), the classic
context-transfer model (latent-space arithmetic), on the same cells, genes and few-shot split as the other models.
Data: pilot_common.build_adata (log1p CP10K over the pilot genes). scGen's cell type = the sgRNA label (controls
'ctrl'), its condition = the state (rest / stim). Training cells: all cells except the stimulated cells of the test
sgRNAs (as in scGen's own hold-out of one cell type in the stimulated condition); the resting cells of the test sgRNAs
are in training. Training: the settings of scGen's perturbation-prediction tutorial (SCGEN defaults, max_epochs 100,
batch_size 32, early stopping with patience 25 on scvi-tools' own validation split).
Prediction as SCGEN.predict(restrict_arithmetic_to='all'): delta = mean latent of the stimulated cells minus mean
latent of the resting cells after scGen's balancer (every sgRNA resampled to the size of the largest); the balanced
means are computed exactly as the equal-weight means over sgRNAs of the per-sgRNA latent means, instead of by
resampling (the resampled set would hold about 2,500 x 6,000 cells). Each test sgRNA's resting cells are encoded,
shifted by delta and decoded (module.generative, as in predict); its predicted effect is the mean decoded expression
minus that of the control cells' resting cells shifted and decoded in the same way. Scored with pilot_common.score.
Usage: CUDA_VISIBLE_DEVICES=<gpu> PILOT_DONOR=<D> pilot_scgen.py <data_dir> <out_json> [max_epochs]"""
import sys, json
import numpy as np, torch
from pilot_common import build_adata, effects, score, test_targets


def group_mean(Z, lab):
    u, inv = np.unique(lab, return_inverse=True)
    M = np.zeros((len(u), Z.shape[1])); np.add.at(M, inv, Z)
    return u, M / np.bincount(inv)[:, None]


def main():
    d, outp = sys.argv[1:3]; me = int(sys.argv[3]) if len(sys.argv) > 3 else 100
    import scgen, scvi
    scvi.settings.seed = 0
    des, a, genes = build_adata(d)
    test = sorted({"%s|g%d" % (t, p) for t in test_targets(des) for p in (0, 1)})
    a.obs["sgrna"] = a.obs["condition"].astype(str)
    hold = (a.obs["state"] == "stim") & a.obs["sgrna"].isin(test)
    tr = a[~hold.to_numpy()].copy()
    tr.obs = tr.obs[["sgrna", "state"]].astype("category")
    scgen.SCGEN.setup_anndata(tr, batch_key="state", labels_key="sgrna")
    m = scgen.SCGEN(tr)
    m.train(max_epochs=me, batch_size=32, early_stopping=True, early_stopping_patience=25)
    hist = m.history
    Z = m.get_latent_representation(tr)
    st = tr.obs["state"].to_numpy(); lab = tr.obs["sgrna"].astype(str).to_numpy()
    _, Mr = group_mean(Z[st == "rest"], lab[st == "rest"]); _, Ms = group_mean(Z[st == "stim"], lab[st == "stim"])
    delta = Ms.mean(0) - Mr.mean(0)
    m.module.eval()

    def decode(zs):
        with torch.no_grad():
            return m.module.generative(torch.Tensor(zs + delta).to(m.device))["px"].cpu().numpy()
    rest = st == "rest"
    base = decode(Z[rest & (lab == "ctrl")]).mean(0)
    pred = {}
    for g in test:
        w = rest & (lab == g)
        if w.sum():
            pred[g] = decode(Z[w]).mean(0) - base
    _, res = effects(d)
    sc = score(pred, des, res, genes)
    ep = len(hist["elbo_train"]) if "elbo_train" in hist else None
    json.dump(dict(model="scGen", n_pred_guides=len(pred), epochs_run=ep, n_train_cells=int(tr.n_obs), scores=sc, rules=__doc__),
              open(outp, "w"), indent=1)
    print("epochs", ep, len(pred), json.dumps(sc))


if __name__ == "__main__":
    main()
