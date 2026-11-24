"""Generate reproducible figures, tables, perturbation controls and a run audit."""
import argparse
import csv
import json
import math
import os
from collections import Counter
from pathlib import Path
os.environ.setdefault("MPLCONFIGDIR", "/tmp/erdos99-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle
import numpy as np
from .geometry import normalize, pairs, squared_distances, diagnostics
from .starts import known_small
from .verify import verify_points

# Historical reference values, NOT claimed to be current records or proven minima.
# Audet, Fournier, Hansen & Messine (2010), Table 1, DOI 10.1007/s11590-010-0185-y.
HISTORICAL = {4:math.sqrt(2),5:(1+math.sqrt(5))/2,6:2*math.sin(2*math.pi/5),7:2.,
    8:1/(2*math.sin(math.pi/14)),9:2.569236979951405,10:2.777310779181221,
    11:2.867454500369003,12:math.sqrt(5+2*math.sqrt(3)),13:3.151195929924547,
    14:3.315828389056499,15:3.469645266872202,16:3.590157365310607,
    17:3.753733968657966,18:3.837352745620864,19:2*math.sqrt(2+math.sqrt(3)),
    20:4.090785438142876,21:4.216039902557467,22:4.365086999476002,
    23:4.476583810221013,24:math.sqrt(14+4*math.sqrt(3)),25:4.703357120488387,
    26:4.776039181929329,27:4.863703305156273,28:4.998143544423982,
    29:5.092001689651576,30:5.184522011522479}
# n=27 is calculated from the largest root of the paper's cubic, avoiding OCR.
HISTORICAL[27] = float(max(np.roots([1,-3,-12,16]).real))


def savefig(fig,path):
    # Append extensions: with_suffix would discard '.001' in margin_0.001.
    fig.savefig(Path(str(path)+".png"),dpi=150,bbox_inches="tight")
    fig.savefig(Path(str(path)+".svg"),bbox_inches="tight")
    plt.close(fig)


def draw(ax,points,title,disks=False,labels=False):
    p=normalize(points,safety=0)
    p-=p.mean(axis=0)
    i,j=pairs(len(p)); d=np.sqrt(squared_distances(p))
    contacts=d <= 1+1e-6
    diam=d >= d.max()-1e-6
    if disks:
        for x,y in p:
            ax.add_patch(Circle((x,y),.5,facecolor="#e8eff6",edgecolor="#b9c9da",lw=.7))
    ax.add_collection(LineCollection(np.stack([p[i[diam]],p[j[diam]]],axis=1),
                      colors="#e3a587",linestyles="dashed",linewidths=.8,zorder=2))
    ax.add_collection(LineCollection(np.stack([p[i[contacts]],p[j[contacts]]],axis=1),
                      colors="#176b91",linewidths=1.3,zorder=3))
    ax.scatter(*p.T,s=18,c="#112b42",zorder=4)
    if labels:
        for k,(x,y) in enumerate(p):
            ax.annotate(str(k),(x,y),xytext=(4,4),textcoords="offset points",fontsize=8)
    ax.set_aspect("equal"); ax.autoscale(); ax.margins(.15 if not disks else .2)
    ax.set_title(title,fontsize=10); ax.axis("off")


def perturbation_control(out):
    p=known_small(7)
    rng=np.random.default_rng(99)
    direction=rng.normal(size=p.shape)
    rows=[]
    for amplitude in (0.,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10):
        q=normalize(p+amplitude*direction)
        v=verify_points(q.tolist())
        rows.append(dict(amplitude=amplitude,diameter=v["normalized_diameter"],
                         diameter_excess=v["normalized_diameter"]-2,
                         triangle_margin=v["triangle_margin"],
                         triangles_1e8=v["near_unit_triangles"]["1e-8"],
                         triangles_1e6=v["near_unit_triangles"]["1e-6"],
                         triangles_1e4=v["near_unit_triangles"]["1e-4"]))
    with (out/"perturbation.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=rows[0]); w.writeheader(); w.writerows(rows)
    fig,ax=plt.subplots(figsize=(7,4))
    ax.loglog([r["amplitude"] for r in rows[1:]],[r["diameter_excess"] for r in rows[1:]],"o-",label="diameter minus 2")
    ax.loglog([r["amplitude"] for r in rows[1:]],[r["triangle_margin"] for r in rows[1:]],"s-",label="minimum triangle margin")
    ax.set(xlabel="size of random displacement",ylabel="change in length",title="Tiny perturbations remove exact triangles at tiny cost")
    ax.legend(); ax.grid(alpha=.2); fig.tight_layout()
    savefig(fig,out/"figures"/"perturbation")


def rhombic_control(out):
    rows=[]
    for n in (7,19,30):
        for degrees in (60.,60.01,60.1,61.,65.,75.,90.):
            theta=np.deg2rad(degrees)
            r=int(np.ceil(np.sqrt(n)))+2
            ij=np.array([(a,b) for a in range(-r,r+1) for b in range(-r,r+1)])
            p=np.column_stack([ij[:,0]+ij[:,1]*np.cos(theta),ij[:,1]*np.sin(theta)])
            p=normalize(p[np.argsort(np.sum(p*p,axis=1),kind="stable")[:n]])
            v=verify_points(p.tolist())
            rows.append(dict(n=n,angle_degrees=degrees,diameter=v["normalized_diameter"],
                             triangle_margin=v["triangle_margin"],triangles_1e6=v["near_unit_triangles"]["1e-6"],
                             leading_diameter_ratio=math.sqrt(math.sin(theta)/math.sin(math.pi/3))))
    with (out/"rhombic.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=rows[0]); w.writeheader(); w.writerows(rows)


def analyze(out,skip_individual=False):
    figures=out/"figures"; figures.mkdir(exist_ok=True)
    configs=[json.loads(p.read_text()) for p in sorted((out/"configurations").glob("*.json"))]
    trials=[json.loads(line) for line in (out/"trials.jsonl").read_text().splitlines()]
    best={c["n"]:c for c in configs if c["margin"] == 0}
    ns=sorted(best)
    # Refresh richer diagnostics without rewriting the actual search record.
    for c in configs:
        p=np.asarray(c["points"])
        c["diagnostics"]=diagnostics(p)
        if skip_individual:
            continue
        name=f"n{c['n']:02d}_"+("unrestricted" if c["margin"] == 0 else f"margin_{c['margin']:g}")
        fig,ax=plt.subplots(figsize=(5,5))
        draw(ax,p,f"n={c['n']}; D={c['diameter']:.9f}; margin={c['margin']:g}",True,True)
        fig.text(.5,.03,"blue: distance 1 within 1e-6   ·   dashed: longest distances",ha="center",fontsize=8)
        savefig(fig,figures/name)
    fig,axes=plt.subplots(5,6,figsize=(17,14))
    for ax,n in zip(axes.flat,ns):
        c=best[n]; d=c["diagnostics"]
        draw(ax,c["points"],f"n={n}    D={c['diameter']:.6f}\n{d['near_unit_triangles']} near-unit triangles")
    for ax in list(axes.flat)[len(ns):]: ax.axis("off")
    fig.suptitle("Best configurations found · n=4–30\nBlue: near-unit contacts (tolerance 1e-6); dashed orange: diameters",fontsize=15)
    fig.tight_layout(rect=(0,0,1,.95)); savefig(fig,figures/"gallery")
    fig,axes=plt.subplots(2,3,figsize=(12,8))
    for ax,n in zip(axes.flat,[4,6,7,8,9,19]):
        draw(ax,best[n]["points"],f"n={n}    D={best[n]['diameter']:.9f}",True)
    fig.suptitle("Six configurations to inspect\nBlue: near-unit contacts; dashed orange: longest distances")
    fig.tight_layout(rect=(0,0,1,.93)); savefig(fig,figures/"highlights")
    margins=sorted({c["margin"] for c in configs})
    fig,(ax,gap)=plt.subplots(2,1,figsize=(10,8),sharex=True)
    for m in margins:
        cs=sorted([c for c in configs if c["margin"] == m],key=lambda c:c["n"])
        label="unrestricted best found" if m == 0 else f"triangle margin ≥ {m:g} (within 1e-8)"
        ax.plot([c["n"] for c in cs],[c["diameter"] for c in cs],"o-",ms=3,label=label)
        if m:
            gap.plot([c["n"] for c in cs],[c["diameter"]-best[c["n"]]["diameter"] for c in cs],"o-",ms=3,label=label)
    ax.plot(ns,[HISTORICAL[n] for n in ns],"k.",label="published upper bounds (2010; not all proven optimal)")
    # Stronger planar lower bound of Bezdek–Fodor, as recorded by BHS (2014).
    ax.plot(ns,[max(math.sqrt(n)-1,math.sqrt(n/(math.pi/math.sqrt(12)))-1) for n in ns],
            "--",color="gray",label="rigorous general lower bound")
    ax.set(ylabel="diameter",title="Numerical upper bounds and the remaining proof gap")
    gap.set(xlabel="number of points n",ylabel="extra diameter found")
    ax.legend(fontsize=8); gap.legend(fontsize=8)
    for a in (ax,gap): a.grid(alpha=.2)
    fig.tight_layout(); savefig(fig,figures/"comparison")
    # Starting-family distributions expose different local minima.
    fig,ax=plt.subplots(figsize=(10,4))
    families=sorted({t["family"] for t in trials if t["margin"]==0 and t["optimizer"].get("method")!="construction"})
    for family in families:
        ts=[t for t in trials if t["margin"]==0 and t["family"]==family]
        ax.scatter([t["n"] for t in ts],[t["diameter"]-best[t["n"]]["diameter"] for t in ts],s=12,label=family,alpha=.7)
    ax.set(xlabel="n",ylabel="diameter above best found",title="Different starts reach different answers")
    ax.legend(fontsize=7,ncol=4); ax.grid(alpha=.2); fig.tight_layout(); savefig(fig,figures/"starts")
    perturbation_control(out); rhombic_control(out)
    # Verify exported JSON again, reading decimal tokens directly, not binary floats.
    from decimal import Decimal
    verified=[]
    for file in sorted((out/"configurations").glob("*.json")):
        c=json.loads(file.read_text(),parse_float=Decimal)
        v=verify_points(c["points"],float(c["margin"]))
        verified.append(dict(file=file.name,**v))
    solved=[t for t in trials if "success" in t["optimizer"]]
    bad=[t for t in solved if not t["optimizer"]["success"]]
    audit=dict(trial_count=len(trials),solver_calls=len(solved),solver_failures=len(bad),
               failure_messages=dict(Counter(t["optimizer"]["message"] for t in bad)),
               fallback_count=sum(t["optimizer"].get("fallback",False) for t in solved),
               all_trials_feasible=all(t["verification"]["exact_saved_feasible"] for t in trials),
               selected_configurations=len(verified),
               all_selected_feasible=all(v["exact_saved_feasible"] for v in verified),
               requested_margin_exact_failures=sum(not v["exact_requested_margin"] for v in verified),
               all_selected_margin_within_1e8=all(v["margin_within_acceptance"] for v in verified),
               selected_verification=verified)
    (out/"audit.json").write_text(json.dumps(audit,indent=2)+"\n")
    (out/"structure.json").write_text(json.dumps([
        dict(n=c["n"],margin=c["margin"],**c["diagnostics"]) for c in configs],indent=2)+"\n")
    lines=["# Initial experiment table", "", "All values are best found in this run, not proofs of global optimality.",
           "Margin ε means every triple has a side ≥ 1+ε after minimum-distance normalization (acceptance tolerance 1e-8).", "",
           "Historical column: upper bounds in Audet et al. (2010), Table 1; these are not asserted to be current records.",
           "Source: https://doi.org/10.1007/s11590-010-0185-y", "",
           "| n | Best found D | ε=0.001 | ε=0.01 | Near-unit triangles | 2010 upper bound |", "|---:|---:|---:|---:|---:|---:|"]
    comparison=[]
    for n in ns:
        c=best[n]
        a=next(x for x in configs if x["n"]==n and x["margin"]==.001)
        b=next(x for x in configs if x["n"]==n and x["margin"]==.01)
        lines.append(f"| {n} | {c['diameter']:.9f} | {a['diameter']:.9f} | {b['diameter']:.9f} | {c['diagnostics']['near_unit_triangles']} | {HISTORICAL[n]:.9f} |")
        comparison.append(dict(n=n,best_found=c["diameter"],margin_0001=a["diameter"],margin_001=b["diameter"],
                               historical_2010=HISTORICAL[n],gap_to_historical=c["diameter"]-HISTORICAL[n]))
    (out/"TABLE.md").write_text("\n".join(lines)+"\n")
    with (out/"historical_comparison.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=comparison[0]); w.writeheader(); w.writerows(comparison)
    print(json.dumps({k:v for k,v in audit.items() if k!="selected_verification"},indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("output",type=Path,nargs="?",default=Path("results/initial"))
    p.add_argument("--skip-individual",action="store_true",help="Reuse already generated individual figures")
    args=p.parse_args()
    analyze(args.output,args.skip_individual)


if __name__=="__main__": main()
