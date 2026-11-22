"""P84 real-repository downstream scan; classify active and historical consumers.

The uploaded scanner's /tmp/kit_p78 listing is NOT a repository audit.
This script must run from the repository root in CI, with the REAL
working tree. Does not mutate caches, past sources or evidence.
"""
from pathlib import Path
from collections import Counter
import argparse
import hashlib
import json
import re

SCHEMA="NUMERIA.P84.S3.REAL_REPOSITORY.DOWNSTREAM.1"
OLD=("17.21","19.31","21.77","24.83","29.01","35.35","46.06")
NEW=("17.49","19.67","22.25","25.48","29.89","36.59","51.65")
PATTERN=re.compile(r"\b(?:TABLE_V2|TABLE_S2_PUBLIC_ROBUST|TABLE_V2_S2_PUBLIC_ROBUST|S2_ROBUST_B1|row_audit|B1_replacement|B1_upper_dyadic)\b|(?:\.b1\b)")
EXCLUDES={".git",".venv","venv","node_modules","__pycache__",".mypy_cache",".pytest_cache"}
ROOT=Path(__file__).resolve().parents[1]
UPLOADED=Path(__file__).resolve().parent/"input_unverified"/"p78_s3_downstream_scan.json"

def role(name,content):
    if name in ("NUMERIA_P76_COEFFICIENTS/table.py","NUMERIA_P75_CS2/cs_v2.py"):
        return "VERSIONED_TABLE_OR_CONSUMER"
    if name.startswith(("NUMERIA_P77_B1/","NUMERIA_P78_RAMARE_PUBLIC/",
                        "NUMERIA_P79_B1_PUBLIC/","NUMERIA_P75_CS2/")):
        return "HISTORICAL_SOURCE_OR_OLD_CERTIFICATE"
    if name.startswith("NUMERIA_P84_S3_REPO/"):
        return "NEW_VERSIONED_ADOPTION"
    if "TABLE_V2" in content or "TABLE" in content and "B1" in content:
        return "POTENTIAL_RUNTIME_CONSUMER"
    return "REFERENCE_ONLY"

def scan(root=ROOT):
    root=Path(root).resolve()
    if not (root/"NUMERIA_P76_COEFFICIENTS"/"table.py").is_file() or not (root/"NUMERIA_P75_CS2"/"cs_v2.py").is_file():
        raise ValueError("expected a real scicalc repository root")
    examined=0
    findings=[]
    for p in sorted(root.rglob("*")):
        if not p.is_file() or p.suffix not in (".py",".json",".md",".yml",".yaml"):
            continue
        if any(x in EXCLUDES for x in p.relative_to(root).parts):
            continue
        if p.stat().st_size>2_000_000:
            continue
        name=p.relative_to(root).as_posix()
        # Original legacy S3 inputs are data: do not mistake their old B1
        # literals for active runtime consumers.
        if "/input_unverified/" in name:
            continue
        examined+=1
        try:
            s=p.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        matches=list(PATTERN.finditer(s))
        old=[v for v in OLD if re.search(r'(?<![0-9.])'+re.escape(v)+r'(?![0-9.])',s)]
        new=[v for v in NEW if re.search(r'(?<![0-9.])'+re.escape(v)+r'(?![0-9.])',s)]
        if not(matches or old or new):
            continue
        first_lines=sorted(set(s.count("\n",0,m.start())+1 for m in matches[:100]))
        findings.append({
            "file":name,
            "classification":role(name,s),
            "B1_related_reference_lines":first_lines[:25],
            "old_B1_literals":old,
            "new_S2_B1_literals":new,
            "sha256":hashlib.sha256(p.read_bytes()).hexdigest(),
        })
    must_have=("NUMERIA_P76_COEFFICIENTS/table.py","NUMERIA_P75_CS2/cs_v2.py")
    seen={f["file"] for f in findings}
    if not set(must_have).issubset(seen):
        raise ArithmeticError("missing versioned table or consumer")
    old=json.loads(UPLOADED.read_text()) if UPLOADED.exists() else {}
    oldroot=old.get("summary",{}).get("root","")
    if oldroot==str(root):
        raise ArithmeticError("legacy kit scanner not independent of this repository scan")
    summary=Counter(x["classification"] for x in findings)
    return {
        "schema":SCHEMA,
        "scan_root":str(root),
        "files_examined":examined,
        "files_with_B1_references":len(findings),
        "classifications":dict(summary),
        "historical_kit_scan_root":oldroot,
        "original_kit_scan_does_not_certify_repo":True,
        "reproduction_policy":"Never rewrite P75/P77/P79 certificates; opt in to versioned S2 tables and produce distinct artifacts",
        "unreviewed_consumer_count":summary["POTENTIAL_RUNTIME_CONSUMER"],
        "requires_follow_up_if_potential_consumers":bool(summary["POTENTIAL_RUNTIME_CONSUMER"]),
        "findings":findings
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--repo-root",type=Path,default=ROOT)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    data=scan(a.repo_root)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(data,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({k:data[k] for k in ("files_examined","files_with_B1_references","classifications","unreviewed_consumer_count","historical_kit_scan_root")},indent=2))

if __name__=="__main__":
    main()
