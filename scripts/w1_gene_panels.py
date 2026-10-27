"""M1 leftover: generate ``data/controls/gene_panels.tsv`` from ``notes/apriori_lists_final.md`` §1 (DESIGN.md
§5.4, Appendix A, amendment A2). Label-free: reads the list file, the baseline gene maps, the models' GPRs
and the CRISPR *header* through ``labels.ids()`` (identifiers only).

The list file is the authority; this script only resolves its text into gene rows. Rules:
- Each subsection's *membership text* is taken as the file states it (the gene list of a §1.1 item is
  its first sentence after the dash; a §1.2 row's first table cell; a §1.4/§1.5/§1.6 item's paragraph,
  minus citation brackets). Genes named only as exceptions ("X moves to §1.2") or in the "Deliberately
  excluded" line are not members.
- Shorthand is expanded mechanically: ``A/B/C`` suffix lists (ATP5F1A/B/C, COQ2/3/6/7, CYP4F2/3/11/12),
  ranges (SDHA-D, MT-ND1-6, AK1-AK4, SPCS1-3, ATP6V0A1-4), and family wildcards (ATP6V*, PIG*, ALG*,
  NDUFAF*, COA*, POLR1/2/3*, UGT1A*, COQ*, CIAO*, xenobiotic families) against the gene vocabulary.
- A named Human-GEM reaction (``MARxxxxx``) that the list uses as a gene set (§1.4d pseudo-reaction GPRs,
  §1.4g MAR06435) is resolved to the genes of that reaction's GPR in each parent.
- Rules the list states as computations are computed: §1.5a (gene closes no reaction: from
  results/w1/gene_features.tsv), §1.6b (Ensembl 116 paralog >= 95% identity, protein-coding part only:
  the local paralog table has no pseudogene identities, stated as a caveat), §1.6e multi-Entrez genes.
- Descriptive phrases with no gene names ("the core complex I and complex III subunits", "OXPHOS
  partials", "carnitine synthesis") are NOT resolved; they are listed in gene_panels_unresolved.tsv.
Every token that is not a gene symbol in the vocabulary (DepMap CRISPR header, both models, Ensembl 116
gene map) is dropped and logged.

Usage: .venv/bin/python scripts/w1_gene_panels.py
Outputs: data/controls/gene_panels.tsv (DESIGN path), results/w1/gene_panels_unresolved.tsv,
results/w1/gene_panels_summary.tsv
"""

from __future__ import annotations

import re
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from w1_common import BASELINE, OUT_DIR, PROJECT, atomic_write_tsv, current, run_main, sha256  # noqa: E402

SRC = PROJECT / "notes" / "apriori_lists_final.md"
OUT = PROJECT / "data" / "controls" / "gene_panels.tsv"
TOKEN = re.compile(r"[A-Z][A-Z0-9]*(?:-[A-Z0-9]+)*(?:/[A-Z0-9]+)*\*?(?:-[A-Z0-9]+)?|[A-Z]+\d+[A-Z]*\*")
TAGS = ("strong", "moderate", "model-masked", "weaker", "uncertain")

FAMILY_RE = {  # wildcard -> regex over the vocabulary (anchored)
    "ATP6V*": r"ATP6V\w+", "PIG*": r"PIG[A-Z]+", "ALG*": r"ALG\d+[A-Z]?", "NDUFAF*": r"NDUFAF\d+",
    "COA*": r"COA\d+", "POLR1/2/3*": r"POLR[123][A-Z0-9]+", "UGT1A*": r"UGT1A\d+", "COQ*": r"COQ\d+[A-Z]?",
    "CIAO*": r"CIAO\d+[A-Z]?", "UGT2B*": r"UGT2B\d+", "CYP4F*": r"CYP4F\d+", "PPIAL4*": r"PPIAL4[A-Z]", "USP17L*": r"USP17L\d+",
    "*ARS": r"[A-Z]ARS1", "*ARS2": r"[A-Z]ARS2",
    "CYP": r"CYP\d+[A-Z]\d+", "UGT": r"UGT\d+[A-Z]\d+", "GST": r"GST[AKMOPTZ]\d+[A-Z]?",
    "SULT": r"SULT\d+[A-Z]\d+", "NAT": r"NAT[12]",
}


def vocabulary():
    from nextgem import labels

    rows = []
    for c in labels.ids()["gene_columns"]:
        m = labels.GENE_COL.match(c)
        rows.append((m["symbol"], int(m["entrez"]), "crispr"))
    for name in ("human1", "human2"):
        mp = pd.read_csv(PROJECT / "results" / "baseline" / f"{name}_gene_mapping.tsv", sep="\t", dtype=str)
        for s, e in zip(mp.symbol, mp.entrez):
            if isinstance(s, str):
                rows.append((s, int(e) if isinstance(e, str) else None, name))
    em = pd.read_csv(PROJECT / "data" / "paralogs" / "ensembl_gene_map.tsv", sep="\t", dtype=str)
    for s, e in zip(em["Gene name"], em["NCBI gene (formerly Entrezgene) ID"]):
        if isinstance(s, str):
            rows.append((s, int(e) if isinstance(e, str) else None, "ensembl"))
    v = pd.DataFrame(rows, columns=["symbol", "entrez", "source"])
    prio = {"crispr": 0, "human2": 1, "human1": 2, "ensembl": 3}
    v = v.assign(p=v.source.map(prio)).sort_values(["symbol", "p"])
    v = v.dropna(subset=["entrez"]).drop_duplicates("symbol").set_index("symbol").entrez.astype(int)
    return v


def expand(tok: str) -> list[str]:
    tok = tok.strip().rstrip(".,;:")
    tok = re.sub(r"\((?:H1|H2)\)", "", tok)
    out = []
    for piece in [tok] if "/" not in tok or tok.endswith("*") else _slash(tok):
        m = re.fullmatch(r"(.*?)(\d+)-(?:\1)?(\d+)", piece)  # MT-ND1-6, SPCS1-3, AK1-AK4, ATP6V0A1-4
        if m and int(m[3]) > int(m[2]) and int(m[3]) - int(m[2]) < 20:
            out += [f"{m[1]}{i}" for i in range(int(m[2]), int(m[3]) + 1)]
            continue
        m = re.fullmatch(r"(.*?)([A-Z])-(?:\1)?([A-Z])", piece)  # SDHA-D, AMY1A-C
        if m and m[1] and ord(m[3]) > ord(m[2]):
            out += [f"{m[1]}{chr(c)}" for c in range(ord(m[2]), ord(m[3]) + 1)]
            continue
        out.append(piece)
    return out


def _slash(tok: str) -> list[str]:
    """``A/B/C`` lists: a digit-only part replaces the trailing digits (COQ2/3/6/7), a single letter the
    trailing letter (ATP5F1A/B/C), a part starting with a digit the same-length tail (SULT1A1/1A2);
    anything else is a full symbol (GYG1/GYG2, HK1/HK2/HK3/HKDC1)."""
    parts = tok.split("/")
    first = parts[0]
    res = [first]
    for p in parts[1:]:
        if p.isdigit():
            res.append(re.sub(r"\d+$", "", first) + p)
        elif len(p) == 1 and p.isalpha():
            res.append(first[:-1] + p)
        elif p[0].isdigit():
            res.append(first[: len(first) - len(p)] + p)
        else:
            res.append(p)
    return res


CLUSTER = {"UGT1A": "UGT1A*", "UGT2B": "UGT2B*", "CYP4F": "CYP4F*"}  # "the UGT1A cluster", "UGT2B", "the CYP4F cluster"


def tokens(text: str) -> list[str]:
    text = text.replace("**", "")
    text = re.sub(r"\[(?:LIT|MS|PRED|OMICS)[^\]]*\]", " ", text)
    text = re.sub(r"`MAR\d+[^`]*`", " ", text)
    text = re.sub(r"\(`[a-z-]+`[^)]*\)", " ", text)
    return [CLUSTER.get(t, t) for t in TOKEN.findall(text)]


def section_text(md: str, start: str, end: str | None) -> str:
    i = md.index(start)
    j = md.index(end, i + len(start)) if end else len(md)
    return md[i:j]


def main() -> int:
    t0, c0 = time.time(), time.process_time()
    md = SRC.read_text()
    vocab = vocabulary()
    import cobra

    from nextgem.model import load_model

    cobra.Configuration().solver = "glpk"
    models = {n: load_model(n) for n in ("human1", "human2")}
    sym = {n: pd.read_csv(PROJECT / "results" / "baseline" / f"{n}_gene_mapping.tsv", sep="\t", dtype=str)
           .set_index("gene").symbol for n in models}
    gf = pd.read_csv(OUT_DIR / "gene_features.tsv", sep="\t")
    U = set(gf.entrez)
    in_h = {n: set(pd.read_csv(PROJECT / "results" / "baseline" / f"{n}_gene_mapping.tsv", sep="\t", dtype=str)
                   .entrez.dropna().astype(int)) for n in models}
    rows, unresolved, dropped = [], [], []

    def add(panel, sub, syms, how, tag="", role="", section=""):
        for s in syms:
            if s not in vocab.index:
                dropped.append({"section": section, "panel": panel, "subpanel": sub, "token": s})
                continue
            rows.append({"section": section, "panel": panel, "subpanel": sub, "symbol": s, "entrez": int(vocab[s]),
                         "role": role, "tag": tag, "how": how})

    def add_tokens(panel, sub, text, tag="", role="", section=""):
        for t in tokens(text):
            if t.endswith("*") or t in FAMILY_RE:
                fam = FAMILY_RE.get(t)
                if fam is None:
                    unresolved.append({"section": section, "subpanel": sub, "phrase": t, "reason": "unknown wildcard"})
                    continue
                add(panel, sub, [s for s in vocab.index if re.fullmatch(fam, s)], f"wildcard:{t}", tag, role, section)
            else:
                add(panel, sub, expand(t), "token", tag, role, section)

    def gpr_of(panel, sub, rid, section, tag=""):
        for n, m in models.items():
            if rid in m.reactions:
                add(panel, sub, sorted({sym[n].get(g.id) for g in m.reactions.get_by_id(rid).genes} - {None}),
                    f"gpr_of:{rid}:{n}", tag, "", section)

    def phrase(section, sub, p):
        unresolved.append({"section": section, "subpanel": sub, "phrase": p, "reason": "descriptive phrase, no gene names"})

    # ---- 1.1 pan-medium essentials: first sentence after the dash of each lettered item
    s11 = section_text(md, "## 1.1 ", "## 1.2 ")
    for m in re.finditer(r"\*\*([A-K])\. ([^*]+?)(?: \(([^)]*)\))?\*\* — (.+)", s11):
        letter, title, head, body = m[1], m[2], m[3] or "", m[4]
        first = re.split(r"\.\s+(?=[A-Z*])", body, maxsplit=1)[0]
        tag = next((t for t in TAGS if t in head), "")
        for part in re.split(r";|, plus ", first):
            ptag = "moderate" if "`moderate`" in part else tag
            if re.search(r"\bsubunits\b", part) and not re.search(r"[A-Z]{2,}\d", part):
                phrase("1.1", letter, part.strip())
                continue
            if "subunits" in part and "ATP5F1A/B/C" in part:
                phrase("1.1", letter, "the ATP-synthase subunits")
            add_tokens("pan_medium", f"{letter}:{title.strip()}", part, ptag or "strong", section="1.1")

    # ---- 1.2 medium-dependent: first cell of each table row; rescuer = second cell
    s12 = section_text(md, "## 1.2 ", "## 1.3 ")
    for line in s12.splitlines():
        if not line.startswith("| ") or line.startswith("| Gene") or line.startswith("|---"):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        genes, rescuer = cells[0], re.sub(r"\s+", " ", cells[1])
        tag = "weaker" if "weaker" in genes else ""
        if "ETC subunits" in genes:
            phrase("1.2", rescuer, "ETC subunits")
        if "selenoproteins" in genes:
            phrase("1.2", rescuer, "selenoproteins (GPX4 named)")
        if "sterol branch" in genes:
            genes = genes.replace("sterol branch", "")
        # role: (+ SLC1A3 as the gate) and (modifier NAPRT) are modifiers
        mod = re.findall(r"\(\+ ([A-Z0-9]+)|modifier ([A-Z0-9]+)", genes)
        mods = {a or b for a, b in mod}
        for t in tokens(genes):
            for s in expand(t):
                add("medium_dependent", rescuer[:60], [s], "token", tag if s == "DHCR24" else "",
                    "modifier" if s in mods else "member", "1.2")

    # ---- 1.3 collateral lethality (trigger -> target) and paralog-buffering sets
    s13 = section_text(md, "## 1.3 ", "## 1.4 ")
    coll = s13.split("**Collateral lethality.**")[1].split("**Paralog-buffering pairs**")[0]
    for m in re.finditer(r"\((\d)\) ([^;\[]+?)(?:\[|;|$)", coll):
        k, txt = m[1], m[2]
        if "→" not in txt:
            continue
        trig, targ = txt.split("→", 1)
        add_tokens("collateral_lethality", f"C{k}", trig, role="trigger", section="1.3")
        add_tokens("collateral_lethality", f"C{k}", targ, role="target", section="1.3")
    buf = s13.split("**Paralog-buffering pairs**")[1].split("**Model GPR error")[0]
    buf = buf.split("]:", 1)[1]
    for i, grp in enumerate(re.split(r",\s+(?![^()]*\))", buf)):
        grp = re.sub(r"\(.*?\)", "", grp).strip(" .\n")
        if grp:
            add_tokens("paralog_buffering", f"B{i + 1:02d}:{grp[:40]}", grp, section="1.3")
    add("gpr_error", "RNR_AND (RRM2B false positive)", ["RRM2B"], "token", section="1.3")

    # ---- 1.4 false-positive classes (a)-(k), 1.5 false-negative classes (a)-(i), 1.6 artefact strata (a)-(g)
    def items(sec_start, sec_end):
        txt = section_text(md, sec_start, sec_end)
        return re.findall(r"^\(([a-k])\) \*\*([^*]+)\*\*(.*)$", txt, flags=re.M)

    for letter, title, body in items("## 1.4 ", "## 1.5 "):
        sub = f"{letter}:{title.strip()}"
        if letter == "d":
            for rid in re.findall(r"MAR\d{5}", body):
                gpr_of("fp_class", sub, rid, "1.4")
            continue
        if letter == "g":
            gpr_of("fp_class", sub, "MAR06435", "1.4")
            add_tokens("fp_class", sub, body.split("the 31-gene")[0], section="1.4")
            phrase("1.4", sub, "the 31-gene nucleoporin AND; the RNR AND (genes as in 1.4d / 1.3)")
            continue
        if letter == "i":
            for fam in ("CYP", "UGT", "GST", "SULT", "NAT"):
                add("fp_class", sub, [s for s in vocab.index if re.fullmatch(FAMILY_RE[fam], s)], f"wildcard:{fam}", section="1.4")
            continue
        if letter in ("e", "k"):
            phrase("1.4", sub, "descriptive pathway names (sterol branch, heme pathway, carnitine synthesis / bile-acid and steroid-hormone synthesis)")
        if letter == "h":
            phrase("1.4", sub, "the OXPHOS partials")
        tag = "uncertain" if letter == "a" else ""
        for t in tokens(body):
            if t in FAMILY_RE or t.endswith("*"):
                add_tokens("fp_class", sub, t, section="1.4")
            else:
                for s in expand(t):
                    add("fp_class", sub, [s], "token", "uncertain" if (tag and s == "ECHS1") else "", section="1.4")

    for letter, title, body in items("## 1.5 ", "## 1.6 "):
        sub = f"{letter}:{title.strip()}"
        if letter == "a":
            nonaddr = gf[(gf.in_human1 & ~gf.addressable_human1) | (gf.in_human2 & ~gf.addressable_human2)]
            for _, r in nonaddr.iterrows():
                rows.append({"section": "1.5", "panel": "fn_class", "subpanel": sub, "symbol": "", "entrez": int(r.entrez),
                             "role": "h1" if (r.in_human1 and not r.addressable_human1) else "h2", "tag": "",
                             "how": "computed:closes_no_reaction"})
            continue
        if letter == "i":
            phrase("1.5", sub, "OXPHOS and glycolysis genes (no list)")
            continue
        if letter == "g":
            add("fn_class", sub, [s for s in vocab.index if re.fullmatch(FAMILY_RE["ATP6V*"], s)], "wildcard:ATP6V*", section="1.5")
        if letter == "e":
            phrase("1.5", sub, "nucleotide-pool turnover, mitochondrial translation and OXPHOS assembly factors, NAD+/GSH turnover")
        add_tokens("fn_class", sub, body, section="1.5")

    for letter, title, body in items("## 1.6 ", "---"):
        sub = f"{letter}:{title.strip()}"
        if letter == "b":
            seeds = body.split("Seed loci always included:")[1]
            add_tokens("crispr_artefact", sub, seeds, role="seed", section="1.6")
            par = _paralog95(vocab)
            add("crispr_artefact", sub, par, "computed:ensembl_paralog_id>=95", role="rule", section="1.6")
            phrase("1.6", sub, "pseudogene >= 95% identity part of the rule (no local pseudogene identity table)")
            continue
        if letter == "e":
            add_tokens("crispr_artefact", sub, body.split("11 multi-Entrez")[0], section="1.6")
            for n in models:
                mp = pd.read_csv(PROJECT / "results" / "baseline" / f"{n}_gene_mapping.tsv", sep="\t", dtype=str)
                add("crispr_artefact", sub, sorted(mp.loc[mp.map_status == "one_to_many", "symbol"]), f"computed:one_to_many:{n}",
                    section="1.6")
            continue
        if letter == "g":
            for fam in ("*ARS", "*ARS2", "ATP6V*", "PIG*", "ALG*", "NDUFAF*", "COA*"):
                add("crispr_artefact", sub, [s for s in vocab.index if re.fullmatch(FAMILY_RE[fam], s)], f"wildcard:{fam}", section="1.6")
            body2 = re.sub(r"`[^`]*`", " ", body.split("plus the §1.4d")[0])
            add_tokens("crispr_artefact", sub, body2, section="1.6")
            for rid in ("MAR13086", "MAR07161", "MAR06978", "MAR10062"):
                gpr_of("crispr_artefact", sub, rid, "1.6")
            continue
        if letter == "c":
            body = body.split("[OMICS]:")[1]
        add_tokens("crispr_artefact", sub, body, section="1.6")

    df = pd.DataFrame(rows)
    df = df.drop_duplicates(["panel", "subpanel", "entrez", "role"]).reset_index(drop=True)
    df["symbol"] = [s if s else next((k for k, v in vocab.items() if v == e), "") for s, e in zip(df.symbol, df.entrez)] \
        if (df.symbol == "").any() else df.symbol
    df["in_human1"] = df.entrez.isin(in_h["human1"])
    df["in_human2"] = df.entrez.isin(in_h["human2"])
    df["in_U"] = df.entrez.isin(U)
    df["source_sha256"] = sha256(SRC)[:16]
    df = df[["section", "panel", "subpanel", "symbol", "entrez", "in_human1", "in_human2", "in_U", "role", "tag", "how",
             "source_sha256"]]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_tsv(df, OUT, float_format=None)
    un = pd.concat([pd.DataFrame(unresolved).assign(kind="unresolved_phrase"),
                    pd.DataFrame(dropped).drop_duplicates().rename(columns={"token": "phrase"}).assign(kind="dropped_token")],
                   ignore_index=True)
    atomic_write_tsv(un, OUT_DIR / "gene_panels_unresolved.tsv", float_format=None)
    summ = df.groupby(["section", "panel", "subpanel"]).agg(n=("entrez", "nunique"), n_in_U=("in_U", "sum")).reset_index()
    atomic_write_tsv(summ, OUT_DIR / "gene_panels_summary.tsv", float_format=None)
    current().output(OUT, OUT_DIR / "gene_panels_unresolved.tsv", OUT_DIR / "gene_panels_summary.tsv")
    print(f"{len(df)} rows, {df.entrez.nunique()} genes, {len(summ)} subpanels; "
          f"{len(unresolved)} unresolved phrases, {len(set(d['token'] for d in dropped))} dropped tokens; sha256 {sha256(OUT)}")
    return 0


def _paralog95(vocab) -> list[str]:
    em = pd.read_csv(PROJECT / "data" / "paralogs" / "ensembl_gene_map.tsv", sep="\t", dtype=str)
    ens2sym = em.dropna(subset=["Gene name"]).drop_duplicates("Gene stable ID").set_index("Gene stable ID")["Gene name"]
    model_ens = set()
    for n in ("human1", "human2"):
        model_ens |= set(pd.read_csv(PROJECT / "results" / "baseline" / f"{n}_gene_mapping.tsv", sep="\t", dtype=str).gene)
    par = pd.read_csv(PROJECT / "data" / "paralogs" / "ensembl_human_paralogs.tsv", sep="\t", usecols=[0, 1, 3, 4],
                      names=["a", "b", "id_ab", "id_ba"], header=0)
    par = par[par.a.isin(model_ens) & ((par.id_ab >= 95) | (par.id_ba >= 95)) & (par.a != par.b)]
    return sorted({ens2sym[a] for a in par.a if a in ens2sym.index} & set(vocab.index))


if __name__ == "__main__":
    raise SystemExit(run_main(main, inputs=[SRC, BASELINE / "human1_gene_mapping.tsv", BASELINE / "human2_gene_mapping.tsv",
        PROJECT / "data" / "paralogs" / "ensembl_gene_map.tsv", PROJECT / "data" / "paralogs" / "ensembl_human_paralogs.tsv",
        OUT_DIR / "gene_features.tsv", PROJECT / "data" / "labels" / "label_ids.json",
        PROJECT / "data" / "Human1" / "model" / "Human-GEM.xml", PROJECT / "data" / "Human-GEM-main" / "model" / "Human-GEM.xml"],
        solver="glpk (structure only)"))
