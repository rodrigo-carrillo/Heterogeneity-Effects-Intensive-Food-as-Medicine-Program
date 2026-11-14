"""
Agent A (rerun) comparison: committed outputs vs outputs regenerated in a clean copy.

    python compare_outputs.py <committed_root> <clean_root> <out_csv> [subdir ...]

Compares every file under <root>/<subdir> (default: outputs/v3 outputs/replication) plus
manuscript/figure1.png.
  CSV  : cell by cell. Numeric cells -> abs/rel difference. Text cells -> every number embedded
         in the text is compared (e.g. "-0.05 (-0.44 to +0.33)"); non-numeric text must match.
  JSON : flattened leaf by leaf (lists indexed).
  MD/TXT: line by line; numbers on each line compared token by token, the rest of the line must
         match after numbers are masked. console_log.txt: if the two logs contain exactly the
         same lines in a different order (parallel worker completion messages), this is noted
         and not counted as a difference.
  PNG  : byte hash, then pixel difference (max abs channel difference, % pixels differing).
Writes one row per differing item (abs diff > 1e-6 or any non-numeric mismatch) to <out_csv>
and a per-file summary to <out_csv minus .csv>_files.csv.
"""
import hashlib
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

TOL = 1e-6
NUM = re.compile(r"[-+−]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")


def to_f(s):
    return float(s.replace("−", "-"))


def nums(s):
    return [to_f(m) for m in NUM.findall(str(s))]


def mask(s):
    return NUM.sub("#", str(s))


class Cmp:
    def __init__(self):
        self.rows, self.files = [], []

    def start(self, f, kind):
        self.cur = dict(file=f, kind=kind, n_compared=0, n_diff=0, max_abs=0.0, max_rel=0.0,
                        text_mismatch=0, note="")

    def num(self, loc, a, b):
        self.cur["n_compared"] += 1
        if (isinstance(a, float) and math.isnan(a)) and (isinstance(b, float) and math.isnan(b)):
            return
        if a == b:
            return
        if any(isinstance(x, float) and math.isnan(x) for x in (a, b)):
            d, r = float("inf"), float("inf")
        else:
            d = abs(a - b)
            r = d / max(abs(a), abs(b)) if max(abs(a), abs(b)) > 0 else 0.0
        self.cur["max_abs"] = max(self.cur["max_abs"], d)
        self.cur["max_rel"] = max(self.cur["max_rel"], r)
        if d > TOL:
            self.cur["n_diff"] += 1
            self.rows.append(dict(file=self.cur["file"], location=loc, committed=a, clean=b,
                                  abs_diff=d, rel_diff=r))

    def text(self, loc, a, b):
        self.cur["n_compared"] += 1
        if a != b:
            self.cur["text_mismatch"] += 1
            self.rows.append(dict(file=self.cur["file"], location=loc, committed=str(a)[:300],
                                  clean=str(b)[:300], abs_diff=np.nan, rel_diff=np.nan))

    def mixed(self, loc, a, b):
        """Compare two strings: numbers numerically, the rest textually."""
        if mask(a) != mask(b):
            self.text(loc, a, b)
            return
        na, nb = nums(a), nums(b)
        if not na:
            self.cur["n_compared"] += 1
            return
        for i, (x, y) in enumerate(zip(na, nb)):
            self.num(f"{loc} [num {i}]", x, y)

    def end(self, note=""):
        self.cur["note"] = note
        self.files.append(self.cur)


def cmp_csv(c, a, b):
    da, db = pd.read_csv(a, keep_default_na=True), pd.read_csv(b, keep_default_na=True)
    if list(da.columns) != list(db.columns) or da.shape != db.shape:
        c.text("shape/columns", f"{da.shape} {list(da.columns)}", f"{db.shape} {list(db.columns)}")
        return "SHAPE DIFFERS"
    for col in da.columns:
        for i in range(len(da)):
            x, y = da[col].iloc[i], db[col].iloc[i]
            loc = f"row {i} ({da.iloc[i, 0]} | {da.iloc[i, 1] if da.shape[1] > 1 else ''}) col {col}"
            if pd.isna(x) and pd.isna(y):
                c.cur["n_compared"] += 1
            elif isinstance(x, (int, float, np.integer, np.floating)) and \
                    isinstance(y, (int, float, np.integer, np.floating)):
                c.num(loc, float(x), float(y))
            else:
                c.mixed(loc, "" if pd.isna(x) else str(x), "" if pd.isna(y) else str(y))
    return ""


def flatten(o, p=""):
    if isinstance(o, dict):
        for k, v in o.items():
            yield from flatten(v, f"{p}.{k}")
    elif isinstance(o, list):
        for i, v in enumerate(o):
            yield from flatten(v, f"{p}[{i}]")
    else:
        yield p, o


def cmp_json(c, a, b):
    fa, fb = dict(flatten(json.loads(a.read_text()))), dict(flatten(json.loads(b.read_text())))
    if set(fa) != set(fb):
        c.text("keys", sorted(set(fa) - set(fb)), sorted(set(fb) - set(fa)))
    for k in fa:
        if k not in fb:
            continue
        x, y = fa[k], fb[k]
        if isinstance(x, bool) or isinstance(y, bool) or x is None or y is None:
            c.text(k, x, y)
        elif isinstance(x, (int, float)) and isinstance(y, (int, float)):
            c.num(k, float(x), float(y))
        else:
            c.mixed(k, str(x), str(y))
    return ""


def cmp_text(c, a, b):
    la = a.read_text(encoding="utf-8", errors="replace").splitlines()
    lb = b.read_text(encoding="utf-8", errors="replace").splitlines()
    note = "" if len(la) == len(lb) else f"line count {len(la)} vs {len(lb)}"
    if a.name == "console_log.txt" and la != lb and sorted(la) == sorted(lb):
        # Parallel cross-fitting workers print "repeat k/n done" in completion order, which
        # varies run to run; compare the log as a multiset of lines instead.
        c.cur["n_compared"] += len(la)
        n_moved = sum(x != y for x, y in zip(la, lb))
        return f"identical set of lines; {n_moved} lines in different order (worker completion order)"
    for i, (x, y) in enumerate(zip(la, lb)):
        c.mixed(f"line {i + 1}", x, y)
    if len(la) != len(lb):
        c.text("line count", len(la), len(lb))
    return note


def cmp_png(c, a, b):
    ha = hashlib.sha256(a.read_bytes()).hexdigest()
    hb = hashlib.sha256(b.read_bytes()).hexdigest()
    c.cur["n_compared"] += 1
    if ha == hb:
        return "byte-identical"
    from PIL import Image
    ia, ib = np.asarray(Image.open(a).convert("RGB")).astype(int), \
        np.asarray(Image.open(b).convert("RGB")).astype(int)
    if ia.shape != ib.shape:
        c.text("image size", ia.shape, ib.shape)
        return "size differs"
    d = np.abs(ia - ib)
    frac = float((d.max(axis=2) > 0).mean())
    c.cur["max_abs"] = float(d.max())
    if d.max() > 0:
        c.cur["n_diff"] += 1
        c.rows.append(dict(file=c.cur["file"], location="pixels", committed=ha[:12], clean=hb[:12],
                           abs_diff=float(d.max()), rel_diff=frac))
    return f"bytes differ; pixels differing {100 * frac:.3f}%, max channel diff {d.max()}"


def main():
    ra, rb, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    subs = sys.argv[4:] or ["outputs/v3", "outputs/replication"]
    c = Cmp()
    files = []
    for s in subs:
        files += sorted(p.relative_to(ra) for p in (ra / s).rglob("*") if p.is_file())
    files.append(Path("manuscript/figure1.png"))
    for rel in files:
        a, b = ra / rel, rb / rel
        kind = rel.suffix.lower().lstrip(".")
        c.start(rel.as_posix(), kind)
        if not b.exists():
            c.text("file", "present", "MISSING in clean run")
            c.end("MISSING")
            continue
        fn = {"csv": cmp_csv, "json": cmp_json, "md": cmp_text, "txt": cmp_text,
              "png": cmp_png}.get(kind)
        c.end(fn(c, a, b) if fn else "not compared")
    # files produced only in the clean run
    for s in subs:
        for p in (rb / s).rglob("*"):
            if p.is_file() and not (ra / p.relative_to(rb)).exists():
                c.start(p.relative_to(rb).as_posix(), p.suffix)
                c.text("file", "MISSING in committed", "present")
                c.end("EXTRA in clean run")
    pd.DataFrame(c.rows, columns=["file", "location", "committed", "clean", "abs_diff",
                                  "rel_diff"]).to_csv(out, index=False)
    fs = pd.DataFrame(c.files)
    fs.to_csv(out.with_name(out.stem + "_files.csv"), index=False)
    with pd.option_context("display.width", 250, "display.max_rows", 500,
                           "display.max_colwidth", 60):
        print(fs.to_string(index=False))
    print(f"\n{len(c.rows)} differing items written to {out}")


if __name__ == "__main__":
    main()
