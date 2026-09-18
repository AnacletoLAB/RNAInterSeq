"""Every number quoted in the manuscript, computed from the current release.

Reads the newest data/interactions_with_sequences_RNA_RNA_union_*.pkl, the
by-key entity tables of 00b and the curated (pre-filter) interaction table,
and writes data/paper_numbers_<release date>.json. The evidence-method
categories reproduce the hand-curated normalisation of mRNAInter_stats.ipynb
(cells 24-25) verbatim, quirks included, so that the tables stay comparable
with the earlier release.
"""
import glob, json, os, re, sys, time
from collections import Counter
from itertools import combinations
import numpy as np, pandas as pd

t0 = time.time()
def log(*a): print(f"[{time.time()-t0:7.1f}s]", *a, flush=True)

PATH = sorted(glob.glob("data/interactions_with_sequences_RNA_RNA_union_*.pkl"))[-1]
DATE = re.search(r"union_(.*)\.pkl", PATH).group(1)
OUT = f"data/paper_numbers_{DATE}.json"
R = {"release": os.path.basename(PATH)}
def nv(x): return [s for s in x if isinstance(s, str) and len(s) > 0] if isinstance(x, list) else []

df = pd.read_pickle(PATH); log("release rows:", len(df), "cols:", len(df.columns))
R["n_columns"] = len(df.columns)
is_m = df.Category1.eq("mRNA") | df.Category2.eq("mRNA")
R["interactions"] = dict(total=len(df), mrna=int(is_m.sum()), ncrna=int((~is_m).sum()),
                         mrna_pct=round(100 * is_m.mean(), 1), ncrna_pct=round(100 * (~is_m).mean(), 1))
keys = pd.concat([df.key1, df.key2]).drop_duplicates()
R["entities"] = dict(total=len(keys), mrna=int(sum(k[1] == "mRNA" for k in keys)), ncrna=int(sum(k[1] != "mRNA" for k in keys)))
R["entities_by_biotype"] = dict(Counter(k[1] for k in keys).most_common())
tax = set(df.Species1_TaxID) | set(df.Species2_TaxID); names = set(df.Species1) | set(df.Species2)
R["species"] = dict(taxids=len(tax), names=len(names))
R["biotypes"] = len(set(df.Category1) | set(df.Category2))
# distinct sequences and their lengths
seqs = set();
for col in ("ent1_sequences", "ent2_sequences"):
    for L in df[col]: seqs.update(nv(L))
ln = np.array([len(s) for s in seqs])
R["distinct_sequences"] = dict(n=len(seqs), min=int(ln.min()), median=float(np.median(ln)), mean=round(float(ln.mean()), 1), max=int(ln.max()))
del seqs, ln
# transcript-level pairs
n1 = df.ent1_sequences.apply(lambda x: len(nv(x))); n2 = df.ent2_sequences.apply(lambda x: len(nv(x)))
tp = n1 * n2
nc = df[~is_m]
dp, dpo = set(), set()                                    # unordered / ordered distinct sequence pairs
for a, b in zip(nc.ent1_sequences, nc.ent2_sequences):
    for s in nv(a):
        for t in nv(b): dp.add((s, t) if s <= t else (t, s)); dpo.add((s, t))
R["transcript_pairs"] = dict(total=int(tp.sum()), mrna=int(tp[is_m].sum()), ncrna=int(tp[~is_m].sum()),
                             ncrna_distinct_unordered=len(dp), ncrna_distinct_ordered=len(dpo))
del dp, dpo
# evidence
sb, wb, pb = df.strong.notna(), df.weak.notna(), df.predict.notna()
R["evidence_all"] = dict(strong=int(sb.sum()), weak=int(wb.sum()), predict=int(pb.sum()), validated=int((sb | wb).sum()))

def normalize_method(name):                              # mRNAInter_stats.ipynb cell 24, verbatim
    name = name.lower().strip()
    name = name.replace('real-time qpc', 'qrt-pcr')
    if 'qrt' in name or 'real-time' in name or 'qPCR' in name: return 'qRT-PCR'
    elif 'western blot' in name: return 'Western blot'
    elif 'luciferase' in name: return 'Luciferase reporter assay'
    elif 'rip' in name: return 'RIP'
    elif 'chip' in name: return 'ChIP'
    elif 'fish' in name: return 'FISH'
    elif 'pull-down' in name or 'pull down' in name: return 'RNA pull-down assay'
    elif 'clap' in name or 'clip' in name: return 'CLIP'
    elif 'ihc' in name: return 'IHC'
    else: return name.title()
normalize_map = {'qPCR': 'qRT-PCR', 'Qpcr': 'qRT-PCR', 'Rt-Pcr': 'qRT-PCR', 'Luciferase reporter assay': 'Reporter assay'}
category_map = {
    'qRT-PCR': 'PCR-based', 'Rna-Seq': 'High-throughput sequencing', 'Microarray': 'High-throughput sequencing',
    'Mario': 'High-throughput sequencing', 'Pare': 'High-throughput sequencing', 'Gmuct': 'High-throughput sequencing',
    'Ligr-Seq': 'High-throughput sequencing', 'Clash': 'High-throughput sequencing', 'RIP': 'RNA-protein binding',
    'Risc-Ip': 'RNA-protein binding', 'Hur-Ip': 'RNA-protein binding', 'Ago2-Ip': 'RNA-protein binding',
    'CLIP': 'RNA-protein binding', 'ChIP': 'DNA/RNA-protein binding', 'Chirp': 'DNA/RNA-protein binding',
    'Reporter assay': 'Reporter assay', 'Egfp Reporter Assay': 'Reporter assay',
    'Dual Fluorescent Reporter Assay': 'Reporter assay', 'Luciferase reporter assay': 'Reporter assay',
    'Western blot': 'Protein-based', 'Immunoblot': 'Protein-based', 'ELISA': 'Protein-based',
    'Northern Blot': 'RNA detection', 'RNA pull-down assay': 'RNA-protein binding', 'RNA interference': 'Functional assay',
    'FISH': 'Imaging', 'IHC': 'Imaging', 'In Situ Hybridization': 'Imaging', 'Immunofluorescence Assay': 'Imaging',
    'MTT Assay': 'Functional assay', 'Transwell Assay': 'Functional assay', 'Brdu Incorporation Assay': 'Functional assay',
    'Chromatography Technology': 'Other', 'Sequencing': 'Other', 'Rpa': 'Other', 'CRISPR-Cas9': 'Other',
    'Rnascope': 'Other', 'Race': 'Other', 'Rlm-Race': 'Other', 'Psilac': 'Other', 'Sybr Green Qpcr': 'PCR-based',
    'Co-Ip': 'RNA-protein binding', 'Chirp-Qpcr': 'PCR-based', 'Wst-1 Assay': 'Functional assay',
    'Rescue Assay': 'Functional assay', 'RNA Hybrid': 'Other', 'EMSA': 'Other', 'IFA': 'Imaging', 'IFS': 'Imaging'}

def evidence_block(sub):
    s, w, p = sub.strong.notna(), sub.weak.notna(), sub.predict.notna(); n = len(sub); v = s | w
    out = dict(n=n, strong=int(s.sum()), weak=int(w.sum()), predict=int(p.sum()), validated=int(v.sum()),
               pred_and_val=int((p & v).sum()), pred_only=int((p & ~v).sum()), val_only=int((v & ~p).sum()),
               strong_and_weak=int((s & w).sum()))
    for k in ("strong", "weak", "predict", "validated"): out[k + "_pct"] = round(100 * out[k] / n, 2)
    meth = Counter()
    for entry in list(sub.strong.dropna()) + list(sub.weak.dropna()):
        for m in entry.split("//"): meth[normalize_method(m)] += 1
    ms = pd.DataFrame(meth.items(), columns=["Method", "Count"])
    ms["Method_normalized"] = ms["Method"].replace(normalize_map)
    summ = ms.groupby("Method_normalized")["Count"].sum()
    cat = summ.groupby(summ.index.map(category_map).fillna("Other")).sum().sort_values(ascending=False)
    out["method_categories"] = [(c, int(v), round(100 * v / n, 2)) for c, v in cat.items()]
    return out
R["evidence_mrna"] = evidence_block(df[is_m]); R["evidence_ncrna"] = evidence_block(nc)
log("evidence done")

# GO / Rfam coverage (01_interaction_figures.py logic)
NS = {"biological_process": "BP", "molecular_function": "MF", "cellular_component": "CC"}
manual = {"GO:0007569": "BP", "GO:0016458": "BP", "GO:0071241": "BP"}
go_ns, cur, ns = {}, [], None
for line in open("go-basic.obo"):
    line = line.rstrip("\n")
    if line in ("[Term]", "[Typedef]"):
        for i in cur:
            if ns: go_ns[i] = NS.get(ns, "Unknown")
        cur, ns = [], None
    elif line.startswith("id: GO:"): cur.append(line[4:])
    elif line.startswith("alt_id: GO:"): cur.append(line[8:])
    elif line.startswith("namespace: "): ns = line[11:]
for i in cur:
    if ns: go_ns[i] = NS.get(ns, "Unknown")
def cat_of(g): return go_ns.get(g, manual.get(g, "Unknown"))
def go_mg(entry, et):
    if not isinstance(entry, dict): return []
    gt = (entry.get("type_of_gene") or "").lower(); et = str(et).lower()
    is_mrna = ("mrna" in et) or ("protein" in et)
    if gt == "protein-coding" and not is_mrna: return []
    if gt != "protein-coding" and is_mrna: return []
    g = entry.get("go", {}) or {}; out = []
    for c in ["BP", "MF", "CC"]:
        v = g.get(c)
        if not v: continue
        if isinstance(v, list): out += [x["id"] for x in v if "id" in x]
        elif isinstance(v, dict) and "id" in v: out.append(v["id"])
    return out
def has_rfam(x):
    if x is None: return False
    if isinstance(x, (np.ndarray, pd.Series)): return x.size > 0 and pd.Series(x).notna().any()
    if isinstance(x, (list, tuple, set)): return any(pd.notna(i) for i in x)
    return bool(pd.notna(x))
def go_block(sub):
    agg = Counter(); st = Counter(); rf = Counter(); cache = {}
    for rc1, mg1, c1, rc2, mg2, c2, r1, r2 in zip(sub.ent1_GO_RNAcentral, sub.ent1_mygene_raw, sub.Category1,
                                                  sub.ent2_GO_RNAcentral, sub.ent2_mygene_raw, sub.Category2,
                                                  sub.ent1_Rfam, sub.ent2_Rfam):
        pres = []
        for rc, mg, cc in ((rc1, mg1, c1), (rc2, mg2, c2)):
            s = {g for g in (rc if isinstance(rc, list) else []) + go_mg(mg, cc) if pd.notna(g)}
            for g in s:
                if g not in cache: cache[g] = cat_of(g)
                agg[cache[g]] += 1
            pres.append(bool(s))
        st["both" if all(pres) else ("one" if any(pres) else "none")] += 1
        h = [has_rfam(r1), has_rfam(r2)]
        rf["both" if all(h) else ("one" if any(h) else "none")] += 1
    tot = sum(agg[c] for c in ("BP", "MF", "CC", "Unknown") if agg.get(c))
    return dict(go_status=dict(st), rfam_status=dict(rf), go_domain_fraction={c: round(agg[c] / tot, 3) for c in ("BP", "MF", "CC") if agg.get(c)},
                go_domain_counts={c: agg[c] for c in ("BP", "MF", "CC", "Unknown") if agg.get(c)})
R["annotation_mrna"] = go_block(df[is_m]); R["annotation_ncrna"] = go_block(nc); log("GO/Rfam done")

# biotype / species pairs and interaction-weighted transcript lengths (01 logic; dog merged by taxid)
NCBI_NAME = {"9615": "Canis lupus familiaris"}
for c in ("1", "2"):
    df[f"Species{c}"] = [NCBI_NAME.get(str(t), n) for t, n in zip(df[f"Species{c}_TaxID"], df[f"Species{c}"])]
def landscape(sub):
    tp_ = pd.Series([tuple(sorted(t)) for t in zip(sub.Category1, sub.Category2)], index=sub.index)
    sp_ = pd.Series([tuple(sorted(t)) for t in zip(sub.Species1, sub.Species2)], index=sub.index)
    vc = tp_.value_counts(); keep = tp_.isin(vc[vc >= 100].index)
    L = []
    for col in ("ent1_sequences", "ent2_sequences"):
        for l in sub.loc[keep, col]: L.extend(len(s) for s in nv(l))
    L = np.array(L)
    return dict(biotype_pairs_all=int(tp_.nunique()), species_pairs_all=int(sp_.nunique()),
                n_filtered=int(keep.sum()), biotype_pairs_filtered=int(tp_[keep].nunique()), species_pairs_filtered=int(sp_[keep].nunique()),
                top_biotype_pair=vc.index[0], top_species_pair=sp_.value_counts().index[0],
                weighted_length_mean=round(float(L.mean())), weighted_length_median=round(float(np.median(L))))
R["landscape_mrna"] = landscape(df[is_m]); R["landscape_ncrna"] = landscape(nc); log("landscape done")
del df, nc

# entity-level statistics from the by-key tables of 00b
def q(a): return dict(n=len(a), mean=round(float(np.mean(a))), median=round(float(np.median(a))), sd=round(float(np.std(a, ddof=1))),
                      min=int(np.min(a)), max=int(np.max(a)), q1=round(float(np.percentile(a, 25))), q3=round(float(np.percentile(a, 75))))
def bins(n): return "1" if n == 1 else "2" if n == 2 else "3-5" if 3 <= n <= 5 else "6-10" if 6 <= n <= 10 else ">10"
E = {}
for tag in ("mrna", "ncrna"):
    e = pd.read_pickle(f"data/entities_bykey_{tag}.pkl")
    e["n"] = e.seq_list.apply(lambda x: len(nv(x))); e["ml"] = e.seq_list.apply(lambda x: float(np.mean([len(s) for s in nv(x)])))
    d = dict(n_entities=len(e), length=q(e.ml.values),
             multiplicity=dict(mean=round(float(e.n.mean()), 2), median=float(e.n.median()), sd=round(float(e.n.std()), 2), max=int(e.n.max()),
                               bins={b: int(v) for b, v in e.n.apply(bins).value_counts().items()}))
    if tag == "mrna":
        m = e[e.category == "mRNA"]; ncp = e[e.category != "mRNA"]
        d["length_mrna_entities"] = q(m.ml.values); d["n_mrna_entities"] = len(m)
        d["length_noncoding_partners"] = q(ncp.ml.values); d["n_noncoding_partners"] = len(ncp)
        d["overlap_with_ncrna_subset"] = int(e.entity_key.isin(pd.read_pickle("data/entities_bykey_ncrna.pkl").entity_key).sum())
    st = pd.read_pickle(f"data/df_step3_bykey_{tag}_with_similarity.pkl")
    st["bin"] = st.n_transcripts_filtered.apply(bins)
    d["similarity"] = dict(n=len(st), mean=round(float(st.mean_sim.mean()), 3), median=round(float(st.mean_sim.median()), 3),
                           min=round(float(st.mean_sim.min()), 3), max=round(float(st.mean_sim.max()), 3),
                           var_mean=round(float(st.var_sim.mean()), 3), var_median=round(float(st.var_sim.median()), 3),
                           by_bin={b: dict(n=int(len(g)), mean=round(float(g.mean_sim.mean()), 3), median=round(float(g.mean_sim.median()), 3),
                                           var_median=round(float(g.var_sim.median()), 3), pairs_median=float(g.n_pairs.median()))
                                   for b, g in st.groupby("bin")})
    E[tag] = d; log(tag, "entity stats done")
R["entities_mrna_subset"] = E["mrna"]; R["entities_ncrna_subset"] = E["ncrna"]

# curated table: species and the interactor-field losses
cur = pd.read_pickle("data/interactions_processed_noSeq.pkl")
R["curated"] = dict(n=len(cur), taxids=len(set(cur.Species1_TaxID) | set(cur.Species2_TaxID)), names=len(set(cur.Species1) | set(cur.Species2)))
json.dump(R, open(OUT, "w"), indent=1, default=str)
log("saved", OUT)
print(json.dumps(R, indent=1, default=str))
