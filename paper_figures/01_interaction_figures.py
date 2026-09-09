"""Regenerate the paper figures flagged for regeneration.

Changes with respect to the notebook versions:
  * histogram/bar panels report RELATIVE FREQUENCIES instead of absolute counts
    (EC and GV): Fig1 type distribution, Fig2 species landscape, Fig4 transcript
    length, GO_distribution, transcript count panels, transcript_pairwise_id;
  * predicted_evidence_type: tight bounding box so the method labels are no
    longer truncated on the left margin;
  * GO panels keep the ordering/colour and computed-label fixes already applied.
Panels that are not frequencies (medians, boxplots, entropy, scatter) keep their
original scale; they are re-exported so the whole set comes from one run.
"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib as mpl
import seaborn as sns, pandas as pd, numpy as np, glob, time, os
from collections import Counter
from scipy.stats import entropy

t0=time.time()
def log(*a): print(f"[{time.time()-t0:7.1f}s]", *a, flush=True)

MD="data/mrna_stats_union/"; ND="data/ncrna_stats_union/"
os.makedirs(MD, exist_ok=True); os.makedirs(ND, exist_ok=True)
mpl.rcParams.update({"figure.figsize":(7,5),"figure.dpi":300,"axes.labelsize":18,
    "axes.titlesize":20,"xtick.labelsize":14,"ytick.labelsize":14,
    "axes.spines.top":False,"axes.spines.right":False,"font.family":"sans-serif"})
sns.set_style("white")
COLOR_TYPE="#4C72B0"; COLOR_SPECIES="#DD8452"; COLOR_LENGTH="#6E6E6E"
COLOR_ENTROPY="#937860"; PASTEL_COLOR="#A6CEE3"; PASTEL_SIM="#B2DF8A"; PASTEL_BOX="#FDBF6F"
GO_ORDER=["BP","MF","CC","Unknown"]
GO_COLORS={"BP":"skyblue","MF":"salmon","CC":"lightgreen","Unknown":"lightgrey"}

PATH=sorted(glob.glob("data/interactions_with_sequences_RNA_RNA_union_*.pkl"))[-1]
log("loading", os.path.basename(PATH))
df_all=pd.read_pickle(PATH)
is_m=df_all.Category1.eq("mRNA")|df_all.Category2.eq("mRNA")
log("rows:", len(df_all))

# GO namespace map from go-basic.obo (avoids a goatools dependency)
NS={"biological_process":"BP","molecular_function":"MF","cellular_component":"CC"}
manual={"GO:0007569":"BP","GO:0016458":"BP","GO:0071241":"BP"}
go_ns={}; cur=[]; ns=None
for line in open("go-basic.obo"):
    line=line.rstrip("\n")
    if line in ("[Term]","[Typedef]"):
        for i in cur:
            if ns: go_ns[i]=NS.get(ns,"Unknown")
        cur=[]; ns=None
    elif line.startswith("id: GO:"): cur.append(line[4:])
    elif line.startswith("alt_id: GO:"): cur.append(line[8:])
    elif line.startswith("namespace: "): ns=line[11:]
for i in cur:
    if ns: go_ns[i]=NS.get(ns,"Unknown")
log("go-basic.obo:", len(go_ns), "ids")
def cat_of(g): return go_ns.get(g, manual.get(g,"Unknown"))

def go_mg(entry, et):
    if not isinstance(entry,dict): return []
    gt=(entry.get("type_of_gene") or "").lower(); et=str(et).lower()
    is_mrna=("mrna" in et) or ("protein" in et)
    if gt=="protein-coding" and not is_mrna: return []
    if gt!="protein-coding" and is_mrna: return []
    g=entry.get("go",{}) or {}; out=[]
    for c in ["BP","MF","CC"]:
        v=g.get(c)
        if not v: continue
        if isinstance(v,list): out+=[x["id"] for x in v if "id" in x]
        elif isinstance(v,dict) and "id" in v: out.append(v["id"])
    return out

def nvalid(x): return len([s for s in x if isinstance(s,str) and len(s)>0]) if isinstance(x,list) else 0

for tag, OUT, sub in [("mRNA", MD, df_all[is_m]), ("ncRNA", ND, df_all[~is_m])]:
    log(f"===== {tag}: {len(sub)} interactions")
    sub=sub.copy()
    sub["type_pair"]=[tuple(sorted(t)) for t in zip(sub.Category1, sub.Category2)]
    sub["sp_pair"]=[tuple(sorted(t)) for t in zip(sub.Species1, sub.Species2)]
    sub["n_tx_total"]=sub.ent1_sequences.apply(nvalid)+sub.ent2_sequences.apply(nvalid)
    vc=sub["type_pair"].value_counts()
    df_f=sub[sub["type_pair"].isin(vc[vc>=100].index)]
    N=len(df_f); log(" filtered:", N)

    # --- Fig1 interaction-biotype distribution -> PROPORTION ---
    tc=(df_f["type_pair"].value_counts()/N).sort_values()
    plt.figure()
    lab=[f"{a}–{b}" for a,b in tc.index]
    sns.barplot(x=tc.values, y=lab, hue=lab, palette=[COLOR_TYPE]*len(lab), legend=False)
    plt.xscale("log"); plt.xlim(tc.values.min()*0.5, 1.0)
    plt.xlabel("Fraction of interactions (log scale)")
    plt.tight_layout(); plt.savefig(OUT+"Fig1_type_distribution.pdf", bbox_inches="tight"); plt.close()
    log("  Fig1 done (proportion)")

    # --- Fig2 species landscape -> PROPORTION ---
    sc=(df_f["sp_pair"].value_counts()/N).head(15).sort_values()
    fig,ax=plt.subplots(figsize=(7,5))
    lab=[f"{a}–{b}" for a,b in sc.index]
    sns.barplot(x=sc.values, y=lab, hue=lab, palette=[COLOR_SPECIES]*len(lab), legend=False, ax=ax)
    ax.set_xscale("log"); ax.set_xlim(sc.values.min()*0.5, 1.0)
    ax.set_xlabel("Fraction of interactions (log scale)")
    fig.subplots_adjust(left=0.35, right=0.95, top=0.95, bottom=0.15)
    plt.savefig(OUT+"Fig2_species_landscape.pdf", bbox_inches="tight"); plt.close()
    log("  Fig2 done (proportion)")

    # --- Fig3b median transcripts per biotype pair (a median: unchanged) ---
    tm=df_f.groupby("type_pair")["n_tx_total"].median().sort_values()
    plt.figure()
    lab=[f"{a}–{b}" for a,b in tm.index]
    sns.barplot(x=tm.values, y=lab, hue=lab, palette=["#C44E52"]*len(lab), legend=False)
    plt.xlabel("Median transcripts per interaction")
    plt.tight_layout(); plt.savefig(OUT+"Fig3b_transcript_median_per_type.pdf", bbox_inches="tight"); plt.close()
    log("  Fig3b done")

    # --- Fig4 transcript length -> PROPORTION ---
    L=[]
    for col in ["ent1_sequences","ent2_sequences"]:
        for l in df_f[col]:
            if isinstance(l,list): L.extend(len(s) for s in l if isinstance(s,str))
    L=np.array(L); L=L[L>0]
    log(f"  lengths n={len(L)} mean={L.mean():.1f} median={np.median(L):.1f}")
    plt.figure()
    sns.histplot(np.log10(L), bins=30, color=COLOR_LENGTH, stat="probability")
    for v,c in [(4096,"#C44E52"),(L.mean(),"#4C72B0"),(np.median(L),"#55A868")]:
        plt.axvline(np.log10(v), color=c, linestyle="--", linewidth=2)
    plt.xlabel("log10 transcript length (nt)"); plt.ylabel("Fraction of transcripts")
    plt.tight_layout(); plt.savefig(OUT+"Fig4_transcript_length_log_distribution.pdf", bbox_inches="tight"); plt.close()
    log("  Fig4 done (proportion)")

    # --- Fig7 species entropy (entropy: unchanged) ---
    sp15=sub["sp_pair"].value_counts().head(15).index
    ev=[entropy(sub[sub["sp_pair"]==s]["type_pair"].value_counts()) for s in sp15]
    es=pd.Series(ev, index=[f"{a}–{b}" for a,b in sp15]).sort_values()
    plt.figure()
    sns.barplot(x=es.values, y=list(es.index), hue=list(es.index),
                palette=[COLOR_ENTROPY]*len(es), legend=False)
    plt.xlabel("Interaction-biotype entropy")
    plt.tight_layout(); plt.savefig(OUT+"Fig7_species_entropy.pdf", bbox_inches="tight"); plt.close()
    log("  Fig7 done")

    # --- GO: aggregate and per-interaction status ---
    agg=Counter(); n_go=np.zeros(len(sub),dtype=int); status=[]; cache={}
    for i,(rc1,mg1,c1,rc2,mg2,c2) in enumerate(zip(
            sub.ent1_GO_RNAcentral, sub.ent1_mygene_raw, sub.Category1,
            sub.ent2_GO_RNAcentral, sub.ent2_mygene_raw, sub.Category2)):
        tot=0; pres=[]
        for rc,mg,cc in ((rc1,mg1,c1),(rc2,mg2,c2)):
            rcl=rc if isinstance(rc,list) else []
            s={g for g in rcl+go_mg(mg,cc) if pd.notna(g)}
            for g in s:
                if g not in cache: cache[g]=cat_of(g)
                agg[cache[g]]+=1
            tot+=len(s); pres.append(bool(s))
        n_go[i]=tot; status.append("both" if all(pres) else ("one" if any(pres) else "none"))
    sub["total_GO"]=n_go; sub["GO_status"]=status
    log("  GO status:", dict(Counter(status)))

    # GO_distribution -> PROPORTION of GO annotations
    cats=[c for c in GO_ORDER if agg.get(c)]
    tot_go=sum(agg[c] for c in cats)
    plt.figure()
    plt.bar(cats,[agg[c]/tot_go for c in cats], color=[GO_COLORS[c] for c in cats])
    plt.ylabel("Fraction of GO annotations"); plt.title(f"{tag} GO annotations by domain")
    plt.tight_layout(); plt.savefig(OUT+"GO_distribution.pdf", bbox_inches="tight"); plt.close()
    log("  GO_distribution done (proportion):", {c:round(agg[c]/tot_go,4) for c in cats})

    # GO boxplot (distribution, labels computed from data)
    sc2=pd.Series(status).value_counts()
    plt.figure(figsize=(8,6))
    ax=sns.boxplot(x="GO_status", y="total_GO", data=sub, order=["one","both"],
                   hue="GO_status", hue_order=["one","both"],
                   palette=["#66C2A5","#FC8D62"], legend=False)
    ax.set_xticks([0,1])
    ax.set_xticklabels([f"one ({sc2.get('one',0):,} interactions)",
                        f"both ({sc2.get('both',0):,} interactions)"])
    plt.xlabel(""); plt.ylabel("GO annotations per interaction")
    plt.tight_layout(); plt.savefig(OUT+"GO_annotations_boxplot.pdf", bbox_inches="tight"); plt.close()
    log("  GO_annotations_boxplot done")

    # --- evidence-method panels (already percentages; fix truncation) ---
    for field, fname, pal in [("predict","predicted_evidence_type.pdf","coolwarm")]:
        meth=sub[field].fillna("").str.split("//")
        ex=meth[sub[field].notna()].explode()
        cnt=ex.value_counts().rename_axis("Method").reset_index(name="Count")
        cnt["Percentage"]=100*cnt["Count"]/len(sub)
        cnt=cnt.sort_values("Count", ascending=True)
        plt.figure(figsize=(10,6))
        sns.barplot(x="Percentage", y="Method", data=cnt, hue="Method",
                    palette=pal, legend=False)
        for i,c in enumerate(cnt["Count"]):
            plt.text(cnt["Percentage"].iloc[i]+0.5, i, f"{c}", va="center", fontsize=9)
        plt.title("Predicted interactions by computational method")
        plt.xlabel("Percentage of interactions"); plt.ylabel("")
        plt.xlim(0, cnt["Percentage"].max()+5)
        plt.tight_layout()
        plt.savefig(OUT+fname, bbox_inches="tight")   # <- fixes the truncated labels
        plt.close()
        log(f"  {fname} done (tight bbox, {len(cnt)} methods)")
    del sub, df_f

del df_all
log("interaction-level figures done")
