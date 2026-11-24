"""Independent stdlib-only verifier. No optimizer geometry routines or NumPy.

Exact rational arithmetic certifies feasibility of the *saved decimal coordinates*.
It cannot establish optimality or interpret near-unit contacts as exact contacts.
"""
import argparse
import csv
import json
import math
from decimal import Decimal, localcontext, ROUND_CEILING
from fractions import Fraction
from itertools import combinations
from pathlib import Path


def verify_points(points, requested_margin=0., pair_csv=None):
    if len(points) < 2 or any(len(p) != 2 for p in points):
        raise ValueError("Expected at least two planar points")
    p = [[Fraction(str(x)) for x in row] for row in points]
    q = {}
    for i, j in combinations(range(len(p)), 2):
        q[i,j] = sum((a-b)**2 for a,b in zip(p[i],p[j]))
    lo, hi = min(q.values()), max(q.values())
    if lo <= 0:
        raise ValueError("Coincident points")
    triple_q = [max(q[a,b],q[a,c],q[b,c]) for a,b,c in combinations(range(len(p)),3)]
    margin_q = min(triple_q)/lo if triple_q else None
    with localcontext() as ctx:
        ctx.prec = 60
        def dec(f):
            return Decimal(f.numerator)/Decimal(f.denominator)
        ds = {ij: dec(v).sqrt() for ij,v in q.items()}
        mind, maxd = dec(lo).sqrt(), dec(hi).sqrt()
        # Round diameter up to 10 decimal places, then check that bound exactly.
        upper = maxd.quantize(Decimal("1e-10"), rounding=ROUND_CEILING)
        if Fraction(upper)**2 < hi:
            upper += Decimal("1e-10")
        thresholds = ("1e-8", "1e-6", "1e-4")
        triangles = {}
        for tol in thresholds:
            cap = (1+Fraction(tol))**2*lo
            triangles[tol] = sum(t <= cap for t in triple_q)
        minimum = str(mind)
        maximum = str(maxd)
        ratio = float(dec(hi/lo).sqrt())
        gap = float(dec(margin_q).sqrt()-1) if margin_q is not None else None
        margin_lower = (dec(margin_q).sqrt()-1).quantize(Decimal("1e-10"), rounding="ROUND_FLOOR") if margin_q is not None else None
        if margin_lower is not None:
            while (1+Fraction(margin_lower))**2 > margin_q:
                margin_lower -= Decimal("1e-10")
        if pair_csv:
            path = Path(pair_csv)
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["i", "j", "distance_decimal60"])
                writer.writerows((i,j,str(d)) for (i,j),d in ds.items())
    # Allow a documented 1e-8 numerical margin residual for optimizer acceptance;
    # separately expose the exact requested-margin check, with NO such tolerance.
    exact_margin = margin_q is None or margin_q >= (1+Fraction(str(requested_margin)))**2
    return dict(n=len(p), pair_count=len(q), min_distance_decimal60=minimum,
                max_distance_decimal60=maximum, normalized_diameter=ratio,
                exact_saved_feasible=bool(lo >= 1), certified_diameter_upper=str(upper),
                exact_saved_unit_triangles=sum(t == 1 and q[a,b] == 1 and q[a,c] == 1 and q[b,c] == 1
                    for (a,b,c),t in zip(combinations(range(len(p)),3),triple_q)),
                near_unit_triangles=triangles, triangle_margin=gap,
                certified_triangle_margin_lower=str(margin_lower) if margin_lower is not None else None,
                requested_margin=requested_margin, exact_requested_margin=bool(exact_margin),
                margin_within_acceptance=bool(gap is None or gap >= requested_margin-1e-8))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path)
    parser.add_argument("--pair-csv", type=Path)
    args = parser.parse_args()
    data = json.loads(args.path.read_text(), parse_float=Decimal)
    points = data["points"] if isinstance(data, dict) else data
    margin = float(data.get("margin", 0)) if isinstance(data, dict) else 0
    result = verify_points(points, margin, args.pair_csv)
    print(json.dumps(result, indent=2))
    if not result["exact_saved_feasible"] or not result["margin_within_acceptance"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
