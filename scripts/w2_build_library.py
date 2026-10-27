"""W2 step 1: machine-readable media library, exchange map, per-line medium assignment (label-free).

Reads data/media/* (four Gibco/ATCC sheets, HPLM mmc1.xlsx, Plasmax aau7314_sm.pdf), DepMap Model.csv,
ModelCondition.csv, ScreenSequenceMap.csv, CRISPRScreenMap.csv and the three SLC1A* expression columns. Opens no
label-bearing file.

Outputs (results/w2/media_library/):
  compositions.tsv           every formulation row on disk (raw component, concentration, source row)
  canonical_composition.tsv  medium x canonical metabolite (uM), mapping status
  exchange_map.tsv           raw component -> canonical (stoichiometry) -> exchange id per parent, status, note
  unmapped_by_medium.tsv     every component of every medium that is not mapped to an exchange (and why)
  formulation_strings.tsv    every distinct DepMap medium string: base, class, parts, serum, supplements, flags
  line_media.tsv             per screened line: tier, text, stratum, attribution flag, medium used, medium_key, gate
  distinct_media.tsv         medium_key -> representative text, open canonical set, flags, n_lines
  mf_crosscheck.tsv          MF-ID family x base medium (text) agreement
  shuffled_media.tsv         negative control: medium_key permuted within lineage x serum-free, SEEDS[0:5]
  serum_model.tsv            serum components, rule, scaling, representability
  library_summary.json       counts (tiers, strata, attribution set by base, unmapped per medium)

Usage: .venv/bin/python scripts/w2_build_library.py
"""

from __future__ import annotations

import json

import pandas as pd

from nextgem import media
from nextgem.config import SEEDS
from nextgem.w2_runtime import PARENTS, Manifest, atomic_write_tsv, load_parent, w2_dir

ON_DISK = set(media.SHEET_FILES)  # base media whose sheet is on disk


def main() -> None:
    out = w2_dir() / "media_library"
    out.mkdir(parents=True, exist_ok=True)
    dm = media.paths.depmap_dir()
    inputs = [media.media_dir() / f for f, _ in media.SHEET_FILES.values()] + [
        media.media_dir() / media.HPLM_FILE, media.media_dir() / media.PLASMAX_FILE] + [
        dm / f for f in ("Model.csv", "ModelCondition.csv", "ScreenSequenceMap.csv", "CRISPRScreenMap.csv",
                         "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")]
    with Manifest("scripts/w2_build_library.py", solver="glpk", inputs=inputs) as man:
        rows = media.read_compositions()
        comp = media.canonical_composition(rows)
        comp["status"] = comp.canonical.map(media.mapping_status)
        comp["exchange"] = comp.canonical.map(lambda c: media.EXCHANGES.get(c, ("",))[0])
        sheets = media.sheet_media(comp)

        # ---- exchange map, verified per parent
        checks = {}
        for p in PARENTS:
            m = load_parent(p, solver="glpk")  # verification only; no Gurobi session needed
            checks[p] = media.verify_exchange_map(m).set_index("canonical")
            del m
        em = []
        for raw in sorted(set(rows.raw_component)):
            mapped = media._canon_rows(raw)
            if isinstance(mapped, str):
                em.append(dict(raw_component=raw, canonical=f"<{mapped}>", stoich=None, status="non_nutrient",
                               exchange_human1="", exchange_human2="", note=f"{mapped}: not a nutrient"))
                continue
            for canon, n in mapped:
                st = media.mapping_status(canon)
                row = dict(raw_component=raw, canonical=canon, stoich=n, status=st,
                           note=media.MAPPING_NOTES.get(canon, ""))
                for p in PARENTS:
                    if st == "mapped":
                        c = checks[p].loc[canon]
                        assert c.is_exchange and c.name_ok, (p, canon, c.to_dict())
                        row[f"exchange_{p}"] = c.exchange
                        row[f"met_name_{p}"] = c.met_name
                    else:
                        row[f"exchange_{p}"] = ""
                        row[f"met_name_{p}"] = ""
                em.append(row)
        em = pd.DataFrame(em)
        em.insert(0, "mapping_version", media.MAPPING_VERSION)
        # supplement-only canonical metabolites (from formulation strings), verified the same way
        sup_canon = sorted({c for c in media.EXCHANGES} - set(em.canonical))
        sup = pd.DataFrame([dict(mapping_version=media.MAPPING_VERSION, raw_component="(supplement/serum/floor/"
                                 "background)", canonical=c, stoich=1, status="mapped",
                                 note=media.MAPPING_NOTES.get(c, ""),
                                 **{f"exchange_{p}": checks[p].loc[c].exchange for p in PARENTS},
                                 **{f"met_name_{p}": checks[p].loc[c].met_name for p in PARENTS})
                            for c in sup_canon if all(checks[p].loc[c].name_ok for p in PARENTS)])
        em = pd.concat([em, sup], ignore_index=True)
        bad = [c for c in media.EXCHANGES if not all(checks[p].loc[c].name_ok and checks[p].loc[c].is_exchange
                                                      for p in PARENTS)]
        assert not bad, f"exchange map fails verification: {bad}"

        unm = comp[comp.status != "mapped"][["medium", "canonical", "status", "uM", "raw_components"]].copy()
        unm["note"] = unm.canonical.map(lambda c: media.MAPPING_NOTES.get(c, "non-nutrient (buffer/dye)"
                                                                           if c.startswith("<") else ""))

        # ---- formulation strings and line assignment
        lines = media.assign_line_media(ON_DISK)
        texts = sorted(set(lines.text.dropna()) | set(lines.onboarded_text.dropna()) | {media.IMPUTED_TEXT})
        fs, meds, keys = [], {}, {}
        for t in texts:
            p = media.parse_formulation(t)
            buildable = bool(p.parts) and all(b in ON_DISK for b, _ in p.parts)
            k = None
            if buildable:
                med = media.build_medium(p, sheets)
                k = media.medium_key(p, med, sheets)
                med.key = k
                meds.setdefault(k, (t, med))
                keys[t] = k
            fs.append(dict(text=t, base=p.base, base_class=p.base_class,
                           parts=";".join(f"{b}:{f:.4g}" for b, f in p.parts), serum_pct=p.serum_pct,
                           serum_type=p.serum_type, serum_free=p.serum_free,
                           supplements=" | ".join(f"{s.token.strip()} -> {s.kind}:" + ",".join(
                               f"{c}={'' if v is None else round(v, 4)}" for c, v in s.components.items())
                               + (f" [{s.note}]" if s.note else "") for s in p.supplements),
                           flags=";".join(p.flags), buildable_from_disk=buildable, medium_key=k,
                           n_lines_screen_or_onboarded=int((lines.text == t).sum())))
        fs = pd.DataFrame(fs)
        lines["medium_key"] = lines.medium_text_used.map(keys)
        assert lines.medium_key.notna().all()
        gate = media.gate_states(lines.ModelID.tolist())
        lines = lines.merge(gate, on="ModelID", how="left")
        lines["gate_open"] = lines.gate_open.astype("boolean")
        lines["gate_effective"] = lines.gate_open.fillna(True).astype(bool)
        lines["gate_flag"] = lines.gate_open.isna().map({True: "no_expression_gate_open_by_default", False: ""})

        dist = []
        for k, (t, med) in sorted(meds.items()):
            dist.append(dict(medium_key=k, representative_text=t, serum=med.serum, allow_cysteine=med.allow_cysteine,
                             n_open_canonical=len([c for c in med.present() if c in media.EXCHANGES]),
                             open_canonical=",".join(sorted(c for c in med.present() if c in media.EXCHANGES)),
                             unmapped_canonical=",".join(sorted(c for c in med.present() if c not in media.EXCHANGES)),
                             flags=";".join(med.flags),
                             n_lines=int((lines.medium_key == k).sum()),
                             n_attribution=int(((lines.medium_key == k) & lines.attribution_set).sum())))
        dist = pd.DataFrame(dist)

        mf = media.mf_crosscheck(lines)
        shuf = media.shuffled_assignment(lines, SEEDS[0:5])

        serum = pd.DataFrame(
            [dict(component=c, exchange=media.EXCHANGES[c][0], rule="pFBA uptake at >=99% WT growth with c open "
                  "at 1000 on the reference medium (RPMI sheet + floor at 1000)", bound="f x req_c",
                  f_levels=",".join(map(str, media.SERUM_F_LEVELS)), representable=True) for c in media.SERUM_F_SCALED]
            + [dict(component=c, exchange=media.EXCHANGES[c][0], rule="availability (trace metal)", bound="1000",
                    f_levels="", representable=True) for c in media.SERUM_TRACE_METALS]
            + [dict(component="selenite", exchange="", rule="availability (trace element)", bound="none",
                    f_levels="", representable=False)]
            + [dict(component=c, exchange=media.EXCHANGES[c][0], rule="floor component (not serum, not f-scaled)",
                    bound="10 x req_c (FVA-min)", f_levels="", representable=True)
               for c in ("linoleate", "linolenate")]
        )
        serum["excluded_by_design"] = ""
        serum = pd.concat([serum, pd.DataFrame([dict(component=c, rule="excluded (DESIGN W2: no nucleosides and no "
                                                     "hypoxanthine in serum)", bound="0", representable=True,
                                                     excluded_by_design="yes")
                                                for c in ("hypoxanthine", "uridine", "thymidine")])],
                          ignore_index=True)
        serum["fbs_reference_conc"] = "PENDING domain agent (DESIGN W2: concentrations + citations, then frozen)"

        files = {
            "compositions.tsv": rows, "canonical_composition.tsv": comp, "exchange_map.tsv": em,
            "unmapped_by_medium.tsv": unm, "formulation_strings.tsv": fs, "line_media.tsv": lines,
            "distinct_media.tsv": dist, "mf_crosscheck.tsv": mf, "shuffled_media.tsv": shuf,
            "serum_model.tsv": serum,
        }
        for name, df in files.items():
            man.add(atomic_write_tsv(df, out / name))

        summary = dict(
            mapping_version=media.MAPPING_VERSION, cysteine_rule=media.CYSTEINE_RULE,
            sheet_bases_on_disk=sorted(ON_DISK),
            n_components_by_medium=rows.groupby("medium").size().to_dict(),
            unmapped_by_medium={m: sorted(g.canonical) for m, g in unm.groupby("medium")},
            tier_counts=lines.tier.value_counts().to_dict(), stratum_counts=lines.stratum.value_counts().to_dict(),
            attribution_set_n=int(lines.attribution_set.sum()),
            attribution_set_by_base=lines[lines.attribution_set].base_class.value_counts().to_dict(),
            n_distinct_strings=len(fs), n_distinct_medium_keys=len(dist),
            n_serum_free_screened=int(lines.serum_free.astype("boolean").fillna(False).sum()),
            gate_open_counts=lines.gate_open.astype("string").fillna("NA").value_counts().to_dict(),
            mf_disagreements=mf[~mf.agrees].to_dict("records"),
        )
        js = out / "library_summary.json"
        js.write_text(json.dumps(summary, indent=1, default=str))
        man.add(js)
        print(json.dumps({k: summary[k] for k in ("unmapped_by_medium", "stratum_counts", "attribution_set_by_base",
                                                  "n_distinct_strings", "n_distinct_medium_keys",
                                                  "gate_open_counts")}, indent=1))


if __name__ == "__main__":
    main()
