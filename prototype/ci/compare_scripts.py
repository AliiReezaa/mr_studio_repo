"""Cross-OS check of the generated scripts and their results (contracts v1.1, C6.2 and C6.5).
For every operating system the pipeline saves the pair folder that the script tests produced.
This compares, against the first operating system (alphabetically):
  1. the generated R scripts: must be byte-identical (they contain no date, user or path);
  2. the analysis data tables: numbers within 1e-12 relative;
  3. the result tables written by running the scripts: numbers within 1e-8 relative.
Usage: python3 compare_scripts.py <reports_dir>"""
import csv, glob, os, sys

TOL_DATA, TOL_RESULT = 1e-12, 1e-8
root = sys.argv[1]
pairs = {}
for d in sorted(glob.glob(os.path.join(root, "report-*"))):
    os_id = os.path.basename(d).replace("report-", "")
    hits = glob.glob(os.path.join(d, "**", "script-pair"), recursive=True)
    if hits:
        pairs[os_id] = hits[0]
if len(pairs) < 2:
    print("fewer than two operating systems saved a script pair; nothing to compare"); sys.exit(1)

def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.reader(f))

def is_num(x):
    try:
        float(x); return True
    except ValueError:
        return False

def compare_tables(ref_path, other_path, tol, label):
    a, b = read_csv(ref_path), read_csv(other_path)
    if len(a) != len(b) or a[0] != b[0]:
        print(f"FAIL {label}: shape or header differs"); return 1
    worst, bad = 0.0, 0
    for ra, rb in zip(a[1:], b[1:]):
        for x, y in zip(ra, rb):
            if is_num(x) and is_num(y):
                fx, fy = float(x), float(y)
                rel = abs(fx - fy) / max(1.0, abs(fx)); worst = max(worst, rel)
            elif x != y:
                bad += 1
    status = "PASS" if (worst <= tol and bad == 0) else "FAIL"
    print(f"{status} {label}: max relative difference {worst:.2e}, text mismatches {bad}")
    return 0 if status == "PASS" else 1

ref_os = sorted(pairs)[0]
ref = pairs[ref_os]
failures = 0
for other, path in pairs.items():
    if other == ref_os:
        continue
    for sub, kind in (("02_scripts", "script"), ("01_data", "data"), ("04_tables", "result")):
        names = sorted(os.path.basename(p) for p in glob.glob(os.path.join(ref, sub, "*")))
        names_o = sorted(os.path.basename(p) for p in glob.glob(os.path.join(path, sub, "*")))
        if names != names_o:
            print(f"FAIL {sub}: file lists differ between {ref_os} and {other}"); failures += 1; continue
        for n in names:
            pa, pb = os.path.join(ref, sub, n), os.path.join(path, sub, n)
            label = f"{sub}/{n} ({ref_os} vs {other})"
            if kind == "script":
                same = open(pa, "rb").read() == open(pb, "rb").read()
                print(f"{'PASS' if same else 'FAIL'} {label}: byte-identical" if same else f"FAIL {label}: scripts differ")
                failures += not same
            else:
                failures += compare_tables(pa, pb, TOL_DATA if kind == "data" else TOL_RESULT, label)
print(f"\n{failures} failure(s) across {len(pairs)} operating systems")
sys.exit(1 if failures else 0)
