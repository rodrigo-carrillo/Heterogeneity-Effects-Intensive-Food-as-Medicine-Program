"""Build prompts/prompt_log.md: the author's prompts to Claude Code, in sequence, from session transcripts.

Journals may ask authors who used a large language model in the research to report the prompts and
their sequence. This script keeps only what the author typed (and the author's answers to the model's
multiple-choice questions), and of that only the prompts directly related to the analysis and the
derivation of the manuscript: each prompt is marked include or exclude in prompts/prompt_review.txt (kept
in the internal repository only); prompts not yet reviewed are left out. It drops model output, tool output (which can contain rows of the trial
data), system messages, background-task notifications and automatic context summaries. Directory
paths, session IDs and e-mail addresses are replaced with [path], [id] and [e-mail]; the terms listed
in prompts/redact_terms.txt (one per line; kept in the internal repository only) with [redacted].

Usage (one or more transcripts; sessions are merged and sorted by time, duplicates removed):
    python prompts/extract_prompts.py ~/.claude/projects/<project-folder>/*.jsonl
Claude Code keeps transcripts in ~/.claude/projects/<folder named after the project path>/<session>.jsonl.
"""
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "prompt_log.md"
TERMS_FILE = HERE / "redact_terms.txt"
REVIEW_FILE = HERE / "prompt_review.txt"
SKIP_PREFIXES = ("<task-notification>", "<local-command", "<command-name>", "<command-message>",
                 "[Request interrupted", "Caveat:", "<bash-stdout>", "<bash-stderr>")
ROOT = r"(?<![\w.~/\\])(?:[A-Za-z]:[\\/]|~[\\/]|/(?:[A-Za-z]|Users|home|mnt|tmp|var)/)"  # start of an absolute path
PATH_RES = [
    re.compile(r'"%s[^"\n]*"' % ROOT),                       # quoted path (may contain spaces)
    re.compile(r"'%s[^'\n]*'" % ROOT),
    re.compile(r"%s[^\"'<>|\n]*" % ROOT),                    # unquoted: up to a quote, >, <, | or line end
]


def load_terms():
    if not TERMS_FILE.exists():
        print(f"warning: {TERMS_FILE.name} not found; no names redacted", file=sys.stderr)
        return []
    terms = [t.strip() for t in TERMS_FILE.read_text(encoding="utf-8").splitlines()]
    return sorted((t for t in terms if t and not t.startswith("#")), key=len, reverse=True)


TERMS = load_terms()


def clean(text):
    text = re.sub(r"<system-reminder>.*?</system-reminder>", "", text, flags=re.S)
    text = re.sub(r"</?pasted_content[^>]*>", "", text)                          # paste wrappers
    for rx in PATH_RES:
        text = rx.sub(lambda m: m.group(0)[0] + "[path]" + m.group(0)[0] if m.group(0)[0] in "\"'"
                      else "[path]", text)
    text = re.sub(r"\b[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}\b", "[id]", text)
    text = re.sub(r"[\w.+-]+@[\w-]+\.[\w.-]+", "[e-mail]", text)
    text = re.sub(r"[\w.-]+\.docx\b", "[file]", text)                          # document file names
    for t in TERMS:
        text = re.sub(re.escape(t), "[redacted]", text, flags=re.I)
    return text.strip()


def prompts(path):
    for line in open(path, encoding="utf-8"):
        try:
            o = json.loads(line)
        except json.JSONDecodeError:
            continue
        if o.get("type") != "user" or o.get("isMeta") or o.get("isCompactSummary") or o.get("isSidechain"):
            continue
        ts = o.get("timestamp", "")
        c = o.get("message", {}).get("content")
        res = o.get("toolUseResult")
        if isinstance(res, dict) and isinstance(res.get("answers"), dict):       # answers to the model's questions
            answers = "\n".join(f"- {q} -> **{a}**" for q, a in res["answers"].items())
            yield ts, "Answer to the model's questions", clean(answers)
            continue
        if isinstance(c, list):
            if any(b.get("type") == "tool_result" for b in c):
                continue
            c = "\n".join(b.get("text", "") for b in c if b.get("type") == "text")
        if not isinstance(c, str):
            continue
        t = clean(c)
        if not t or t.startswith(SKIP_PREFIXES):
            continue
        m = re.fullmatch(r"<bash-input>(.*)</bash-input>", t, flags=re.S)       # a shell command typed by the author
        if m:
            yield ts, "Shell command typed by the author", m.group(1).strip()
            continue
        yield ts, "Prompt", t


def load_review():
    """{transcript timestamp: 'include' | 'exclude'} from prompt_review.txt."""
    review = {}
    if REVIEW_FILE.exists():
        for line in REVIEW_FILE.read_text(encoding="utf-8").splitlines():
            parts = line.split("#", 1)[0].split()
            if len(parts) == 2 and parts[0] in ("include", "exclude"):
                review[parts[1]] = parts[0]
    else:
        print(f"warning: {REVIEW_FILE.name} not found; no prompts listed", file=sys.stderr)
    return review


def main(paths):
    review = load_review()
    every = sorted({r for p in paths for r in prompts(p)})
    unreviewed = [r for r in every if r[0] not in review]
    for ts, kind, text in unreviewed:
        print(f"not reviewed (left out): {ts}  # {' '.join(text.split())[:70]}", file=sys.stderr)
    rows = [r for r in every if review.get(r[0]) == "include"]
    lines = ["# Prompt log", "",
             "The author's prompts to Claude Code (Anthropic; Claude Opus 5.5), in sequence. Generated by",
             "`prompts/extract_prompts.py` from the session transcripts; model and tool output are not",
             "included. Times are UTC.", "",
             "Only the sessions on the author's Windows computer (from 22 October 2026) are included; the",
             "prompts of the earlier sessions (June 2026 and 21 October 2026, on another computer) were",
             "not retained.", "",
             "Only the prompts directly related to the analysis and the derivation of the manuscript are",
             "listed; other prompts (for example, status checks, troubleshooting, and file, repository or",
             "administrative tasks) are left out. Document file names are shown as [file].", ""]
    sessions = sorted({Path(p).stem for p in paths})
    lines += [f"Transcripts included: {len(sessions)} session(s), {len(rows)} entries, "
              f"{rows[0][0][:10] if rows else '-'} to {rows[-1][0][:10] if rows else '-'}.", ""]
    day = None
    for ts, kind, text in rows:
        if ts[:10] != day:
            day = ts[:10]
            lines += [f"## {day}", ""]
        body = text.replace("\n", "\n> ")
        lines += [f"**{ts[11:16]} UTC, {kind}**", "", f"> {body}", ""]
    OUT.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {OUT} ({len(rows)} of {len(every)} entries listed, {len(unreviewed)} not reviewed; "
          f"{len(sessions)} session(s))")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
