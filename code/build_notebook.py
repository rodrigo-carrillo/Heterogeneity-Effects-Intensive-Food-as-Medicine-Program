"""
Build notebooks/hte_food_as_medicine.ipynb from code/hte_food_as_medicine.py.

The .py script is the single source of truth. Cells are delimited by
"# %%" lines; "# %% [markdown]" cells become markdown (leading "# " removed).
Run:  python code/build_notebook.py
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "code" / "hte_food_as_medicine.py"
DST = ROOT / "notebooks" / "hte_food_as_medicine.ipynb"


def split_cells(text):
    cells, kind, buf = [], "code", []
    for line in text.splitlines():
        if line.startswith("# %%"):
            cells.append((kind, buf))
            kind, buf = ("markdown" if "[markdown]" in line else "code"), []
        else:
            buf.append(line)
    cells.append((kind, buf))
    return cells


def to_nb_cell(kind, lines):
    while lines and not lines[-1].strip():
        lines = lines[:-1]
    while lines and not lines[0].strip():
        lines = lines[1:]
    if not lines:
        return None
    if kind == "markdown":
        lines = [l[2:] if l.startswith("# ") else l.lstrip("#") for l in lines]
        return {"cell_type": "markdown", "metadata": {}, "source": "\n".join(lines)}
    return {"cell_type": "code", "metadata": {}, "execution_count": None,
            "outputs": [], "source": "\n".join(lines)}


cells = []
for i, (kind, lines) in enumerate(split_cells(SRC.read_text(encoding="utf-8"))):
    if i == 0:
        # module docstring -> markdown intro (strip the triple quotes)
        body = "\n".join(lines).strip().strip('"').strip()
        cells.append({"cell_type": "markdown", "metadata": {},
                      "source": "```\n" + body + "\n```\n\n*Generated from "
                                "`code/hte_food_as_medicine.py` - edit the script, "
                                "then run `python code/build_notebook.py`.*"})
        continue
    c = to_nb_cell(kind, lines)
    if c:
        cells.append(c)

nb = {"cells": cells, "metadata": {"kernelspec": {"display_name": "Python 3",
      "language": "python", "name": "python3"}, "language_info": {"name": "python"}},
      "nbformat": 4, "nbformat_minor": 5}
DST.parent.mkdir(exist_ok=True)
DST.write_text(json.dumps(nb, indent=1), encoding="utf-8")
print(f"Wrote {DST.relative_to(ROOT)} ({len(cells)} cells)")
