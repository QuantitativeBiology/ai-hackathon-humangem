"""W4 re-check: does the A4 genomic ('loss') variant depend on having run the 'none' variant first in the same
process? F0 and F2 (GLPK), R1, ACH-002024 and ACH-003383: job genomic=('none','loss') vs job genomic=('loss',)."""
import hashlib, json
import cobra
cobra.Configuration().solver = "glpk"
import pandas as pd
from nextgem import context_inputs as ci, context_run as cr, paths
W4 = paths.results_dir() / "w4"


def main():
    m = cobra.io.read_sbml_model(str(paths.data_root() / "Human-GEM-main/model/Human-GEM.xml"))
    genes = ci.gene_table(m)
    lines = ["ACH-002024", "ACH-003383"]
    rna = ci.read_expression(genes, lines)
    acc = json.loads((W4 / "tmp" / "accepted_losses_panel.json").read_text())
    st = ci.ThresholdStats("rna", 0, 0, pd.Series(dtype=float), 0, 0)
    items = []
    for line in lines:
        for fam in ("F0", "F2"):
            for gen in (("none", "loss"), ("loss",)):
                j = cr.Job(build_id=f"{line}|{fam}|{'+'.join(gen)}", line=line, family=fam, rule="R1", modality="rna",
                           variant="real", scores=ci.signed_scores(line, "R1", "rna", genes, rna, st), genomic=gen,
                           accepted=tuple(acc[line]), tasks=False, meta={"solver": "glpk"})
                items.append((j.build_id, "glpk", j, None, (line, fam, gen)))
    out = {}

    def done(key, status, res, payload, wall):
        assert status == "ok", res
        row = [r for r in res["rows"] if r["genomic"] == "loss"][0]
        ko = res["ko"]["loss"]
        out[key] = (row["open_mask"], hashlib.sha256(pd.util.hash_pandas_object(ko[["gene", "ratio_r"]], index=False).values.tobytes()).hexdigest()[:16],
                    row["n_reopened"], ko.set_index("gene").ratio_r)
    cr.Runner({"glpk": 4}, 1800, None).run_all(items, done)
    res = {}
    for line in lines:
        for fam in ("F0", "F2"):
            a, b = out[f"{line}|{fam}|none+loss"], out[f"{line}|{fam}|loss"]
            ma = int.from_bytes(bytes.fromhex(a[0]), "big"); mb = int.from_bytes(bytes.fromhex(b[0]), "big")
            res[f"{line}|{fam}"] = {"loss_mask_equal": a[0] == b[0], "n_open_reactions_differ": bin(ma ^ mb).count("1"),
                                    "loss_ko_equal": a[1] == b[1], "n_reopened": [a[2], b[2]],
                                    "ko_n_differ": int((a[3] - b[3]).abs().gt(1e-9).sum()),
                                    "lethal_differ": int(((a[3] < 0.01) != (b[3] < 0.01)).sum())}
    (W4 / "review" / "r12_loss_variant_path.json").write_text(json.dumps(res, indent=1)); print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
