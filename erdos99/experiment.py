"""Run: python -m erdos99.experiment --help. Every trial is persisted as JSONL."""
import argparse
import csv
import json
import os
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import scipy
from .geometry import normalize, squared_distances, diagnostics
from .optimize import optimize
from .starts import generate
from .verify import verify_points


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, allow_nan=False)+"\n")


def run(args):
    if args.n_min < 4 or args.n_max < args.n_min:
        raise ValueError("Use 4 <= n-min <= n-max")
    if any(m <= 0 or m > 1 for m in args.margins):
        raise ValueError("Margins must lie in (0,1]")
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    if (out/"trials.jsonl").exists():
        raise ValueError("Output already contains trials; choose a fresh output directory")
    metadata = dict(command=sys.argv, started_utc=datetime.now(timezone.utc).isoformat(),
                    python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
                    platform=platform.platform(), seeds=args.seeds, margins=args.margins,
                    maxiter=args.maxiter, ftol=1e-10, contact_tolerance=1e-6,
                    acceptance_margin_tolerance=1e-8, safety_inflation=2e-12,
                    thread_env={k:os.environ.get(k) for k in
                                ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS")})
    write_json(out/"metadata.json", metadata)
    all_rows = []
    begin = time.monotonic()
    with (out/"trials.jsonl").open("w") as log:
        for n in range(args.n_min, args.n_max+1):
            starts = list(generate(n,args.seeds))
            if args.custom:
                data = json.loads(args.custom.read_text())
                cp = np.asarray(data["points"] if isinstance(data,dict) else data)
                if len(cp) == n:
                    starts.append(("custom", 0, normalize(cp)))
            trials = []

            def record(points, family, seed, margin, info):
                p = normalize(points)
                v = verify_points(p.tolist(), margin)
                if not v["exact_saved_feasible"]:
                    raise RuntimeError("Independent exact-rational feasibility check failed")
                if not v["margin_within_acceptance"]:
                    raise RuntimeError("Independent margin check failed")
                item = dict(n=n,family=family,seed=seed,margin=margin,
                            diameter=v["normalized_diameter"],points=p.tolist(),
                            verification=v,optimizer=info)
                log.write(json.dumps(item,allow_nan=False)+"\n")
                log.flush()
                trials.append(item)
                return item

            for family, seed, p in starts:
                record(p, family+"_baseline", seed, 0, dict(method="construction"))
                q, info = optimize(p,maxiter=args.maxiter)
                record(q,family,seed,0,info)
            # Include the best few distinct outcomes as well as geometric branches.
            unrestricted = min(trials,key=lambda t:t["diameter"])
            restricted_starts = [("best_unrestricted",0,np.asarray(unrestricted["points"]))]
            restricted_starts += starts
            for margin in args.margins:
                for family,seed,p in restricted_starts:
                    q,info = optimize(p,maxiter=args.maxiter,margin=margin)
                    record(q,family,seed,margin,info)
            # A restricted trial is ALSO a valid unrestricted competitor.
            best = min(trials,key=lambda t:t["diameter"])
            modes = [("unrestricted",best)]
            for margin in args.margins:
                eligible = [t for t in trials if t["verification"]["triangle_margin"] >= margin-1e-8]
                modes.append((f"margin_{margin:g}",min(eligible,key=lambda t:t["diameter"])))
            for label,t in modes:
                item = dict(t)
                # Record requested comparison margin even if winner came from another mode.
                item["source_margin"] = item["margin"]
                item["margin"] = 0. if label == "unrestricted" else float(label.split("_")[1])
                name = f"n{n:02d}_{label}"
                item["verification"] = verify_points(item["points"],item["margin"],out/"pairs"/(name+".csv"))
                item["diagnostics"] = diagnostics(item["points"])
                write_json(out/"configurations"/(name+".json"),item)
                d = item["diagnostics"]
                all_rows.append(dict(n=n,mode=label,diameter=item["diameter"],
                    certified_upper=item["verification"]["certified_diameter_upper"],
                    gap_to_best=item["diameter"]-best["diameter"],
                    contacts=d["contacts"],triangles_tol_1e6=d["near_unit_triangles"],
                    triangle_margin=d["triangle_margin"],contact_rank=d["contact_rank"],
                    bar_flex_dimension=d["bar_flex_dimension"],rotational_orders=str(d["rotational_orders"]),
                    family=item["family"],seed=item["seed"],
                    exact_requested_margin=item["verification"]["exact_requested_margin"],
                    exact_saved_feasible=item["verification"]["exact_saved_feasible"]))
            with (out/"summary.csv").open("w",newline="") as f:
                w=csv.DictWriter(f,fieldnames=list(all_rows[0]))
                w.writeheader(); w.writerows(all_rows)
            failures = sum(t["optimizer"].get("success") is False for t in trials)
            print(f"n={n:2d} best={best['diameter']:.9f} "
                  +" ".join(f"{label}={t['diameter']:.9f}" for label,t in modes[1:])
                  +f" trials={len(trials)} solver_failures={failures}",flush=True)
    metadata.update(finished_utc=datetime.now(timezone.utc).isoformat(),
                    elapsed_seconds=time.monotonic()-begin)
    write_json(out/"metadata.json",metadata)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-min",type=int,default=4)
    parser.add_argument("--n-max",type=int,default=30)
    parser.add_argument("--seeds",type=int,nargs="+",default=[99,2026,31415])
    parser.add_argument("--margins",type=float,nargs="+",default=[.001,.01])
    parser.add_argument("--maxiter",type=int,default=600)
    parser.add_argument("--output",type=Path,default=Path("results/initial"))
    parser.add_argument("--custom",type=Path,help="JSON array of [x,y] or object with points")
    run(parser.parse_args())


if __name__=="__main__":
    main()
