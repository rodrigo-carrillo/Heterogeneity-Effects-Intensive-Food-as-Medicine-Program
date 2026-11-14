#!/bin/bash
# Runs the authors' Stata code and the diagnostics in docs/07_stata_check.md, in a throw-away
# directory, and writes the summaries to verification/stata/results/.
#   bash verification/stata/run_stata_checks.sh [path/to/stata-se]
# Runs on macOS and on Windows (Git Bash). Environment overrides: JOBS (Stata processes at a
# time), PYTHON (interpreter with pandas). On Windows the default is a network Stata,
# one process at a time: it is a 10-seat shared licence.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
case "$(uname -s)" in
    MINGW*|MSYS*|CYGWIN*)
        ST="${1:-//<server>/<share>/StataSE-64.exe}"
        BATCH=/e; JOBS="${JOBS:-1}"
        PY="${PYTHON:-$HOME/.venvs/ete-fam/Scripts/python.exe}"
        export MSYS2_ARG_CONV_EXCL='*' ;;   # keep Git Bash from turning /e into a path
    *)
        ST="${1:-/Applications/StataNow/StataSE.app/Contents/MacOS/stata-se}"
        BATCH=-b; JOBS="${JOBS:-8}"; PY="${PYTHON:-python3}" ;;
esac
stata() { "$ST" "$BATCH" do "$@"; }
# stata_watch <dir> <do-file> [args]: run Stata in <dir>. The network Stata on Windows sometimes
# hangs at random (0 CPU, log frozen), so a run whose log has not changed for STALL seconds is
# killed and retried, up to 3 attempts. Returns 1 if every attempt stalled.
STALL="${STALL:-600}"
kill_proc() {  # on Windows, kill -9 only ends Git Bash's stub, not Stata itself
    local w
    if w=$(cat "/proc/$1/winpid" 2>/dev/null); then taskkill /F /T /PID "$w" >/dev/null 2>&1 || true; fi
    kill -9 "$1" 2>/dev/null || true
}
stata_watch() {
    local dir="$1"; shift
    local log="$dir/$(basename "$1" .do).log" try pid m prev last
    for try in 1 2 3; do
        rm -f "$log"
        ( cd "$dir" && exec "$ST" "$BATCH" do "$@" ) & pid=$!
        prev=""; last=$(date +%s)
        while kill -0 "$pid" 2>/dev/null; do
            sleep 15
            m=$(date -r "$log" +%s 2>/dev/null || echo none)
            if [ "$m" != "$prev" ]; then prev=$m; last=$(date +%s)
            elif (( $(date +%s) - last > STALL )); then
                echo "stalled: $dir (attempt $try), killing" >&2
                kill_proc "$pid"; break
            fi
        done
        wait "$pid" 2>/dev/null && return 0
        (( $(date +%s) - last > STALL )) || return 0   # exited by itself (Stata errors exit too)
    done
    echo "gave up after 3 stalled attempts: $dir" >&2; return 1
}
throttle() { while (( $(jobs -rp | wc -l) >= JOBS )); do wait -n || true; done; }
SRC="$ROOT/data/raw/dataverse_files"
HERE="$ROOT/verification/stata"
OUT="$HERE/results"
# STATA_WORK=<dir> keeps the work directory, so an interrupted run resumes: finished runs are
# skipped. Without it, a temporary directory is used and deleted at the end.
if [ -n "${STATA_WORK:-}" ]; then WORK="$STATA_WORK"; mkdir -p "$WORK"
else WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT; fi
mkdir -p "$OUT"
done_() { [ -f "$WORK/$1/output/T4_regressions.csv" ]; }   # the do-file's last output

# The only edit needed to run the posted code on current Stata (see docs/07 §2): the patient id.
patch() { sed 's/egen n_subject_id = max(_n), by(subject_id)/egen n_subject_id = group(subject_id)/' "$1"; }
setup() {  # setup <dir> <extra lines after insheet>
    done_ "$1" && return 0
    mkdir -p "$WORK/$1/output"; cp "$SRC/analysis_sample.csv" "$WORK/$1/"
    printf '%s' "$2" > "$WORK/$1/extra.txt"     # lines inserted right after the data are read
    patch "$SRC/01_published_tables.do" \
        | awk -v f="$WORK/$1/extra.txt" '{print} /^insheet using/ {while ((getline l < f) > 0) print l}' \
        > "$WORK/$1/01_published_tables.do"
}
run() { done_ "$1" || stata_watch "$WORK/$1" 01_published_tables.do || true; }

# 0. Unmodified code (expected to stop with r(451) on Stata >= 16.1)
mkdir -p "$WORK/orig/output"; cp "$SRC/analysis_sample.csv" "$SRC/01_published_tables.do" "$WORK/orig/"
stata_watch "$WORK/orig" 01_published_tables.do || true
grep -m1 -B1 "^r(" "$WORK/orig/01_published_tables.log" > "$OUT/unmodified_run_error.txt" || echo "no error" > "$OUT/unmodified_run_error.txt"

# 1. Patched code, default state, twice (determinism)
for r in a b; do setup default_$r ""; run default_$r; done
cp "$WORK/default_a/output/"*.csv "$OUT/" 2>/dev/null
for f in "$WORK"/default_a/output/*.csv; do cmp -s "$f" "$WORK/default_b/output/$(basename "$f")" || echo "differs: $(basename "$f")"; done > "$OUT/determinism.txt"
[ -s "$OUT/determinism.txt" ] || echo "two runs identical" > "$OUT/determinism.txt"

# 2. Authors' code under 40 sort seeds (imputation seed unchanged)
for k in $(seq 1 40); do setup ss$k "set sortseed $k"; done
for k in $(seq 1 40); do throttle; run ss$k & done; wait

# 3. Same, with weight in pounds (to test whether the adjusted weight was run in lb)
for k in $(seq 1 40); do setup lb$k "set sortseed $k
foreach m in 0 6 12 {
    replace weight_val_\`m'mo = weight_lb_val_\`m'mo
}"; done
for k in $(seq 1 40); do throttle; run lb$k & done; wait

# 4. Proper per-imputation baseline, under 10 sort seeds
for k in $(seq 1 10); do
    done_ proper$k && continue
    setup proper$k "set sortseed $k"
    "$PY" - "$WORK/proper$k/01_published_tables.do" <<'PY'
import sys; p = sys.argv[1]; s = open(p).read()
old = """            sort n_subject_id post_period
            bysort n_subject_id: replace `var'_baseline = `var'_panel[1] if ///
                mi(`var'_baseline)

            gen post_baseline = post_period * `var'_baseline
"""
new = """            mi passive: by n_subject_id: egen bfill = total(`var'_panel * (post_period == 0))
            mi passive: gen post_baseline = post_period * bfill
"""
assert s.count(old) == 1; open(p, "w").write(s.replace(old, new))
PY
done
for k in $(seq 1 10); do throttle; run proper$k & done; wait

# 5. Diagnostics: what the fill step does (4 sort states) and observed-baseline model
mkdir -p "$WORK/diag"; cp "$SRC/analysis_sample.csv" "$HERE/weight_baseline_fill.do" "$HERE/weight_observed_baselines.do" "$WORK/diag/"
: > "$OUT/weight_baseline_fill.txt"
for k in "" 1 2 3; do
    stata_watch "$WORK/diag" weight_baseline_fill.do $k
    { echo "== sortseed: ${k:-Stata default}"; sed -n '/n_subject_id   treat_now/,/^$/p' "$WORK/diag/weight_baseline_fill.log"; } >> "$OUT/weight_baseline_fill.txt"
done
stata_watch "$WORK/diag" weight_observed_baselines.do
grep "^RESULT" "$WORK/diag/weight_observed_baselines.log" > "$OUT/weight_observed_baselines.txt"

# 6. Summaries
winpath() { command -v cygpath >/dev/null && cygpath -m "$1" || echo "$1"; }  # Windows Python needs C:/ paths
"$PY" "$(winpath "$HERE/summarize_stata.py")" "$(winpath "$WORK")" "$(winpath "$OUT")" "$(winpath "$ROOT")"
echo "done: $OUT"
